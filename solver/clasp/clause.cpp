// 
// Copyright (c) 2006-2011, Benjamin Kaufmann
// 
// This file is part of Clasp. See http://www.cs.uni-potsdam.de/clasp/ 
// 
// Clasp is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// Clasp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Clasp; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
#include <clasp/clause.h>
#include <clasp/solver.h>
#include <clasp/util/misc_types.h>
#include <algorithm>

namespace Clasp { namespace Detail {

struct GreaterLevel {
	GreaterLevel(const Solver& s) : solver_(s) {}
	bool operator()(const Literal& p1, const Literal& p2) const {
		assert(solver_.value(p1.var()) != value_free && solver_.value(p2.var()) != value_free);
		return solver_.level(p1.var()) > solver_.level(p2.var());
	}
private:
	GreaterLevel& operator=(const GreaterLevel&);
	const Solver& solver_;
};

struct Sink { 
	explicit Sink(SharedLiterals* c) : clause(c) {}
	~Sink() { if (clause) clause->release(); }
	SharedLiterals* clause;   
};

// returns an abstraction of p's decision level, s.th
// abstraction(any true literal) > abstraction(any free literal) > abstraction(any false literal).
static uint32 levelAbstraction(const Solver& s, Literal p) {
	ValueRep value_p = s.value(p.var());
	// DL+1,  if isFree(p)
	// DL(p), if isFalse(p)
	// ~DL(p),if isTrue(p)
	uint32   abstr_p = value_p == value_free ? s.decisionLevel()+1 : s.level(p.var()) ^ -(value_p==trueValue(p));
	assert(abstr_p > 0 || (s.isFalse(p) && s.level(p.var()) == 0));
	return abstr_p;
}


} // namespace Detail

/////////////////////////////////////////////////////////////////////////////////////////
// SharedLiterals
/////////////////////////////////////////////////////////////////////////////////////////
SharedLiterals* SharedLiterals::newShareable(const Literal* lits, uint32 size, ConstraintType t, uint32 numRefs) {
	#pragma message TODO("replace with CACHE_LINE_ALIGNED alloc")
	void* m = ::operator new(sizeof(SharedLiterals)+(size*sizeof(Literal)));
	return new (m) SharedLiterals(lits, size, t, numRefs);
}

SharedLiterals::SharedLiterals(const Literal* a_lits, uint32 size, ConstraintType t, uint32 refs) 
	: size_type_( (size << 2) + t ) {
	refCount_ = refs;
	std::memcpy(lits_, a_lits, size*sizeof(Literal));
}

uint32 SharedLiterals::simplify(Solver& s) {
	bool   removeFalse = unique();
	uint32   newSize   = 0;
	Literal* r         = lits_;
	Literal* e         = lits_+size();
	ValueRep v;
	for (Literal* c = r; r != e; ++r) {
		if ( (v = s.value(r->var())) == value_free ) {
			if (c != r) *c = *r;
			++c; ++newSize;
		}
		else if (v == trueValue(*r)) {
			newSize = 0; break;
		}
		else if (!removeFalse) ++c;
	}
	if (removeFalse && newSize != size()) {
		size_type_ = (newSize << 2) | (size_type_ & uint32(3));
	}
	return newSize;
}

void SharedLiterals::release() {
	if (--refCount_ == 0) {
		this->~SharedLiterals();
		::operator delete(this);
	}
}
SharedLiterals* SharedLiterals::share() {
	++refCount_;
	return this;
}
/////////////////////////////////////////////////////////////////////////////////////////
// ClauseCreator
/////////////////////////////////////////////////////////////////////////////////////////
ClauseCreator::ClauseCreator(Solver* s) 
: solver_(s) { }

void ClauseCreator::init(ConstraintType t) {
	assert(!solver_ || solver_->decisionLevel() == 0 || t != Constraint_t::static_constraint);
	literals_.clear();
	extra_     = ClauseInfo().setType(t);
	fwLev_     = 0;
	swLev_     = 0;
	asserting_ = false;
}

ClauseCreator& ClauseCreator::start(ConstraintType t) {
	init(t);
	return *this;
}

ClauseCreator& ClauseCreator::startAsserting(ConstraintType t, const Literal& p) {
	init(t);
	literals_.push_back(p);
	fwLev_     = UINT32_MAX;
	asserting_ = true;
	return *this;
}

ClauseCreator& ClauseCreator::add(const Literal& p) {
	assert(solver_); const Solver& s = *solver_;
	uint32 abstr_p = Detail::levelAbstraction(s, p);
	if (abstr_p != 0) {
		literals_.push_back(p);
		// watch the two literals with highest abstraction, i.e.
		// prefer true over free over false literals
		if (abstr_p > fwLev_) {
			std::swap(abstr_p, fwLev_);
			std::swap(literals_[0], literals_.back());
		}
		if (abstr_p > swLev_) {
			std::swap(literals_[1], literals_.back());
			swLev_= abstr_p;
		}
	}
	// else: drop uncoditionally false literals 
	return *this;
}

void ClauseCreator::simplify() {
	if (literals_.empty()) return;
	swLev_    = 0;
	uint32 sw = 0;
	solver_->markSeen(literals_[0]);
	LitVec::size_type i, j = 1;
	for (i = 1; i != literals_.size(); ++i) {
		Literal x = literals_[i];
		if (solver_->seen(~x)) { fwLev_ = UINT32_MAX; asserting_ = false; }
		if (!solver_->seen(x)) {
			solver_->markSeen(x);
			if (i != j) { literals_[j] = literals_[i]; }
			if (swLev_ != uint32(-1)) {
				uint32 xLev = solver_->value(x.var()) == value_free ? uint32(-1) : solver_->level(x.var());
				if (swLev_ < xLev) {
					swLev_ = xLev;
					sw     = (uint32)j;
				}
			}
			++j;
		}
	}
	literals_.erase(literals_.begin()+j, literals_.end());
	for (LitVec::iterator it = literals_.begin(), end = literals_.end(); it != end; ++it) {
		solver_->clearSeen(it->var());
	}
	if (literals_.size() >= 2 && sw != 1) {
		std::swap(literals_[1], literals_[sw]);
	}
}

ClauseCreator::Status ClauseCreator::status() const {
	if (fwLev_ > varMax) {
		if (asserting_)           return status_unit;
		if (fwLev_ == UINT32_MAX) return status_subsumed;
		return status_sat;
	}
	else if (literals_.empty()) {
		return status_empty;
	}
	else if (solver_->isFalse(literals_[0])) {
		return status_conflicting;
	}
	else if (literals_.size()==1 || solver_->isFalse(literals_[1])) {
		return status_unit;
	}
	return status_open;
}

ClauseCreator::Result ClauseCreator::end() {
	assert(solver_);
	Result ret(0, true);
	Solver& s     = *solver_;
	Status  stat  = status();
	if (stat == status_subsumed || stat == status_sat) {
		return ret;
	}
	if (stat == status_empty) { 
		ret.ok = s.addUnary(negLit(0), extra_.type());
		return ret; 
	}
	if (stat == status_conflicting) {
		ret.ok = false;
		return ret;
	}
	if (extra_.type() != Constraint_t::static_constraint && !isSentinel(s.sharedContext()->tagLiteral())) {
		if (std::find(literals_.begin(), literals_.end(), ~s.sharedContext()->tagLiteral()) != literals_.end()) {
			extra_.tag();
		}
	}
	return ClauseCreator::create(s, literals_, 0, extra_);
}

void ClauseCreator::initWatches(Solver& s, LitVec& lits, bool allFree) {
	uint32 size = static_cast<uint32>(lits.size());
	if (size < 2) return;
	uint32 fw, sw;
	if (allFree) {
		if (s.strategies().randomWatches) {
			RNG& r = s.strategies().rng;
			fw = r.irand(size);
			while ( (sw = r.irand(size)) == fw) {/*intentionally empty*/}
		}
		else {
			uint32 watch[2] = {0, size-1};
			uint32 count[2] = {s.numWatches(~lits[0]), s.numWatches(~lits[size-1])};
			uint32 maxCount = count[0] < count[1];
			for (uint32 x = 1, end = size-1; count[maxCount] > 0u && x != end; ++x) {
				uint32 cxw = s.numWatches(~lits[x]);
				if (cxw < count[maxCount]) {
					if (cxw < count[1-maxCount]) {
						watch[maxCount]   = watch[1-maxCount];
						count[maxCount]   = count[1-maxCount];
						watch[1-maxCount] = x;
						count[1-maxCount] = cxw;
					}
					else {
						watch[maxCount]   = x;
						count[maxCount]   = cxw;
					}
				}
			}
			fw = watch[0];
			sw = watch[1];
		}
		std::swap(lits[0], lits[fw]);
		std::swap(lits[1], lits[sw]);
	}
	else {
		// first watch: free/sat or highest dl
		// second watch: free/sat or second highest dl
		fw = 0; sw = 1;
		uint32 abstr_fw = Detail::levelAbstraction(s, lits[0]);
		uint32 abstr_sw = 0;
		for (LitVec::size_type i = 1; i != lits.size(); ++i) {
			uint32 abstr_p = Detail::levelAbstraction(s, lits[i]);
			if (abstr_p > abstr_fw) {
				assert(abstr_fw >= abstr_sw);
				sw        = fw;
				abstr_sw  = abstr_fw;
				fw        = i;
				abstr_fw  = abstr_p;
			}
			else if (abstr_p > abstr_sw) {
				sw       = i;
				abstr_sw = abstr_p;
			}
		}
		if (sw == 0) { std::swap(lits[0], lits[1]); }
		std::swap(lits[0], lits[fw]);
		std::swap(lits[1], lits[sw]);
	}
}

ClauseHead* ClauseCreator::newProblemClause(Solver& s, LitVec& lits, const ClauseInfo& e, uint32 flags) {
	ClauseHead* ret;
	initWatches(s, lits, true);
	if (lits.size() <= Clause::MAX_SHORT_LEN || !s.sharedContext()->shareConstraints()) {
		ret = Clause::newClause(s, &lits[0], (uint32)lits.size(), e);
	}
	else {
		ret = Clause::newShared(s, SharedLiterals::newShareable(lits, e.type(), 1), e, &lits[0], false);
	}
	if ( (flags & clause_no_add) == 0 ) {
		s.add(ret);
	}
	return ret;
}

ClauseHead* ClauseCreator::newLearntClause(Solver& s, LitVec& lits, const ClauseInfo& e, uint32 flags) {
	ClauseHead* ret;
	Detail::Sink sharedPtr(0);
	sharedPtr.clause = s.sharedContext()->distribute(s, &lits[0], static_cast<uint32>(lits.size()), e);
	if (lits.size() <= Clause::MAX_SHORT_LEN || sharedPtr.clause == 0) {
		if (!s.isFalse(lits[1]) || lits.size() < s.strategies().compress()) {
			ret = Clause::newLearntClause(s, lits, e);
		}
		else {
			std::stable_sort(lits.begin()+2, lits.end(), Detail::GreaterLevel(s));
			ret = Clause::newContractedClause(s, lits, e, 2, true);
		}
	}
	else {
		ret              = Clause::newShared(s, sharedPtr.clause, e, &lits[0], false);
		sharedPtr.clause = 0;
	}
	if ( (flags & clause_no_add) == 0 ) {
		s.addLearnt(ret, (uint32)lits.size(), e.tagged());
	}
	return ret;
}

ClauseCreator::Status ClauseCreator::status(const Solver& s, Literal* lits, uint32 flags) {
	if (s.isFalse(lits[1]) && (!s.isTrue(lits[0]) || s.level(lits[1].var()) < s.level(lits[0].var()))) {
		return status_unit; // or conflicting!
	}
	else if (s.isTrue(lits[0])) {
		uint32 satLevel = s.level(lits[0].var());
		return satLevel == 0 || ((flags & clause_not_sat) != 0) || ((flags & clause_not_root_sat) != 0 && satLevel <= s.rootLevel())
			? status_subsumed
			: status_sat;
	}
	return status_open;
}

ClauseCreator::Result ClauseCreator::create(Solver& s, LitVec& lits, uint32 flags, const ClauseInfo& extra) {
	assert(s.decisionLevel() == 0 || extra.type() != Constraint_t::static_constraint);
	Result ret;
	if (extra.type() == Constraint_t::static_constraint && s.satPrepro()) {
		ret.ok = s.satPrepro()->addClause(lits);
		return ret;
	}	
	if (lits.size() > 1) {
		Status x    = status(s, &lits[0], flags);
		if (x != status_subsumed) {
			s.strategies().heuristic->newConstraint(s, &lits[0], lits.size(), extra.type());
			uint32 impL = (x == status_unit ? s.level(lits[1].var()) : UINT32_MAX);
			if (lits.size() > 3 || extra.tagged() || (flags&clause_explicit) != 0) {
				ret.local = (extra.type() == Constraint_t::static_constraint)
						? newProblemClause(s, lits, extra, flags)
						: newLearntClause(s, lits, extra, flags);
				ret.ok = x != status_unit || s.force(lits[0], impL, Antecedent(ret.local));
				return ret;
			}
			else if (lits.size() == 3) {
				s.addTernary(lits[0], lits[1], lits[2], extra.type());
				ret.ok = x != status_unit || s.force(lits[0], impL, Antecedent(~lits[1], ~lits[2]));
				return ret;
			}
			else {
				s.addBinary(lits[0], lits[1], extra.type());
				ret.ok = x != status_unit || s.force(lits[0], impL, Antecedent(~lits[1]));
				return ret;
			}
		}
		return ret;
	}
	else {
		Literal imp = !lits.empty()	? lits[0] : negLit(0);
		ret.ok = s.addUnary(imp, extra.type());
		return ret;
	}	
}
struct ClauseCreator::Watches {
	Watches() : fwLev(0), swLev(0) {
		lits[0] = lits[1] = lits[2] = negLit(0);
	}
	void init(Solver& s, const Literal* it, const Literal* end) {
		for (uint32 pos = 0; it != end; ++it) {
			Literal p      = *it;
			uint32 abstr_p = Detail::levelAbstraction(s, p);
			lits[pos]      = p;
			// watch the two literals with highest abstraction, i.e.
			// prefer true over free over false literals
			if (abstr_p > fwLev) {
				std::swap(abstr_p, fwLev);
				std::swap(lits[0], lits[pos]);
			}
			if (abstr_p > swLev) {
				std::swap(lits[1], lits[pos]);
				swLev = abstr_p;
			}
			if (++pos == Clause::MAX_SHORT_LEN) {
				pos = 2;
			}
		}
	}
	uint32  fwLev;
	uint32  swLev;
	Literal lits[Clause::MAX_SHORT_LEN];
};

ClauseCreator::Result ClauseCreator::integrate(Solver& s, SharedLiterals* clause, uint32 modeFlags, const ClauseInfo& e) {
	assert(!s.hasConflict() && "ClauseCreator::integrate() - precondition violated!");
	Detail::Sink shared( (modeFlags & clause_no_release) == 0 ? clause : 0);
	// determine state of clause
	if (clause->size() > 1) {
		Watches w;
		if ( (modeFlags & clause_known_order) == 0 ) {
			w.init(s, clause->begin(), clause->end());
		}
		else {
			std::memcpy(w.lits, clause->begin(), std::min(clause->size(), uint32(Clause::MAX_SHORT_LEN))*sizeof(Literal));
		}
		Status stat = status(s, w.lits, modeFlags);
		if (stat != status_subsumed) {
			Literal imp = stat == status_unit ? w.lits[0] : posLit(0);
			ClauseHead* x = 0;
			if (clause->size() > Clause::MAX_SHORT_LEN) {
				x = Clause::newShared(s, clause, e, w.lits, shared.clause == 0);
				shared.clause = 0;
			}
			else if (clause->size() > 3 || (modeFlags & clause_explicit) != 0) {
				x = Clause::newClause(s, w.lits, clause->size(), e);
			}
			else {
				// implicitly shared via binary/ternary implication graph;
				// only check for implication/conflict but do not create
				// a local representation for the clause
				Antecedent ante(~w.lits[1], ~w.lits[2]);
				return Result(x, imp == posLit(0) || s.force(imp, s.level(w.lits[1].var()), ante));
			}
			if ( (modeFlags & clause_no_add) == 0 ) { s.addLearnt(x, clause->size()); }
			return Result(x, imp == posLit(0) || s.force(imp, s.level(w.lits[1].var()), x));
		}
		return Result(0, true);
	}
	else {
		Literal imp = clause->size() > 0 ? *clause->begin() : negLit(0);
		return Result(0, s.addUnary(imp, Constraint_t::static_constraint));
	}
}
ClauseCreator::Result ClauseCreator::integrate(Solver& s, SharedLiterals* clause, uint32 modeFlags) { 
	return integrate(s, clause, modeFlags, ClauseInfo().setType(clause->type())); 
}

/////////////////////////////////////////////////////////////////////////////////////////
// Clause
/////////////////////////////////////////////////////////////////////////////////////////
void* Clause::alloc(Solver& s, uint32 lits) {
	if (lits <= Clause::MAX_SHORT_LEN) {
		return s.allocSmall();
	}
	#pragma message TODO("replace with CACHE_LINE_ALIGNED alloc") 
	uint32 extra = std::max((uint32)ClauseHead::HEAD_LITS, lits) - ClauseHead::HEAD_LITS; 
	return ::operator new(sizeof(Clause) + (extra)*sizeof(Literal));
}

ClauseHead* Clause::newClause(Solver& s, const Literal* lits, uint32 size, const ClauseInfo& e) {
	assert(size >= 2);
	return new (alloc(s, size)) Clause(s, e, lits, size, size, false);
}

ClauseHead* Clause::newShared(Solver& s, SharedLiterals* shared_lits, const ClauseInfo& e, const Literal* lits, bool addRef) {
	return mt::SharedLitsClause::newClause(s, shared_lits, e, lits, addRef);
}

ClauseHead* Clause::newContractedClause(Solver& s, const LitVec& lits, const ClauseInfo& e, LitVec::size_type tailStart, bool extend) {
	assert(lits.size() >= 2);
	return new (alloc(s, (uint32)lits.size())) Clause(s, e, &lits[0], static_cast<uint32>(lits.size()), static_cast<uint32>(tailStart), extend);
}

ClauseHead* Clause::newShortClause(Solver& s, void* mem, const Literal* lits, uint32 size, const ClauseInfo& extra) {
	assert(size >= 2);
	assert(size <= static_cast<uint32>(Clause::MAX_SHORT_LEN));
	return new (mem ? mem : s.allocSmall()) Clause(s, extra, lits, size, size, false);
}

Clause::Clause(Solver& s, const ClauseInfo& e, const Literal* a_lits, uint32 size, uint32 tail, bool extend) 
	: ClauseHead(e) {
	assert(tail == size || s.isFalse(a_lits[tail]));
	if (size > MAX_SHORT_LEN) {
		tail_.data.init(false);
		// copy literals
		std::memcpy(head_, a_lits, size*sizeof(Literal));
		tail           = std::max(tail, (uint32)ClauseHead::HEAD_LITS);
		tail_.data.setSize(tail);
		tail_.data.setIdx(3);         // next idx to look at propagation
		if (tail < size) {            // contracted clause
			head_[size-1].watch();      // mark last literal of clause
			Literal t = head_[tail];
			if (s.level(t.var()) > 0) {
				tail_.data.markContracted();
				if (extend) {
					s.addUndoWatch(s.level(t.var()), this);
				}
			}
		}
	}
	else {
		tail_.data.init(true);
		std::memcpy(head_, a_lits, std::min(size, (uint32)ClauseHead::HEAD_LITS)*sizeof(Literal));
		Literal* t = smallLits();
		t[0] = size > ClauseHead::HEAD_LITS   ? a_lits[ClauseHead::HEAD_LITS]   : negLit(0);
		t[1] = size > ClauseHead::HEAD_LITS+1 ? a_lits[ClauseHead::HEAD_LITS+1] : negLit(0);
		assert(tail_.data.isSmall() && Clause::size() == size);
	}
	attach(s);
}

Clause::Clause(Solver& s, const Clause& other) : ClauseHead(ClauseInfo()) {
	info_.rep      = info_.rep;
	uint32 oSize   = other.size();
	if (oSize <= MAX_SHORT_LEN) {
		tail_.data.init(true);
		std::memcpy(head_, other.head_, ClauseHead::HEAD_LITS*sizeof(Literal));
		LitRange tail  = const_cast<Clause&>(other).tail();
		smallLits()[0] = tail.first[0];
		smallLits()[1] = tail.first[1];
	}
	else {
		tail_.data.init(false);
		std::memcpy(head_, other.head_, oSize*sizeof(Literal));
		tail_.data.setSize(oSize);
		tail_.data.setIdx(3);
	}
	attach(s);
}

Constraint* Clause::cloneAttach(Solver& other) {
	assert(type() == Constraint_t::static_constraint);
	return new (alloc(other, Clause::size())) Clause(other, *this);
}

void Clause::destroy(Solver* s, bool detachFirst) {
	if (detachFirst && s) { 
		Clause::detach(*s);
	}
	void* mem   = static_cast<Constraint*>(this);
	bool  small = isSmall();
	this->~Clause();
	if (!small) { ::operator delete(mem); }
	else if (s) { s->freeSmall(mem);      }
}

void Clause::detach(Solver& s) {
	if (contracted()) {
		Literal* eoc = longEnd();
		if (s.isFalse(*eoc) && s.level(eoc->var()) != 0) {
			s.removeUndoWatch(s.level(eoc->var()), this);
		}
	}
	ClauseHead::detach(s);
}


bool Clause::updateWatch(Solver& s, uint32 pos) {
	int32 idx   = -3, bound = -1;
	if (!isSmall()) {
		bound = (int32)(tail_.data.size());
		idx   = (int32)(tail_.data.idx_);
	}
	for (;;) {
		for (int32 i = idx; i < bound; ++i) {
			if (!s.isFalse(head_[i])) {
				std::swap(head_[pos], head_[i]);
				if (i > 0) { tail_.data.setIdx(i+1); }
				return true;
			}
		}
		if (idx <= 3) { break; }
		bound = idx;
		idx   = 3;
	}
	return false;
}

void Clause::reason(Solver& s, Literal p, LitVec& out) {
	ClauseHead::bumpActivity();
	LitVec::size_type i = out.size();
	out.push_back(~head_[p == head_[0]]);
	if (!isSentinel(head_[2])) {
		out.push_back(~head_[2]);
		LitRange t = tail();
		for (const Literal* r = t.first; r != t.second; ++r) {
			out.push_back(~*r);
		}
		if (contracted()) {
			const Literal* r = t.second;
			do { out.push_back(~*r); } while (!r++->watched());
		}
	}
	if (s.strategies().updateLbd && ClauseHead::type() != Constraint_t::static_constraint) {
		setLbd( s.computeLbd(p, &out[0]+i, &out[0]+out.size()) );
	}
}

bool Clause::minimize(Solver& s, Literal p, CCMinRecursive* rec) {
	ClauseHead::bumpActivity();
	uint32 other = p == head_[0];
	bool   ok    = s.ccMinimize(~head_[other], rec) && s.ccMinimize(~head_[2], rec);
	if (ok) {
		LitRange t = tail();
		for (const Literal* r = t.first; r != t.second; ++r) {
			if (!s.ccMinimize(~*r, rec)) return false;
		}
		if (contracted()) {
			for (const Literal* r = t.second; (ok = s.ccMinimize(~*r, rec)) == true && !r->watched(); ) {
				++r;
			}
		}
	}
	return ok;
}

bool Clause::isReverseReason(const Solver& s, Literal p, uint32 maxL, uint32 maxN) {
	uint32 other = p == head_[0];
	if (!isRevLit(s, head_[other], maxL) || !isRevLit(s, head_[2], maxL)) {
		return false;
	}
	uint32 notSeen = !s.seen(head_[other].var()) + !s.seen(head_[2].var());
	LitRange t     = tail();
	for (const Literal* r = t.first; r != t.second && notSeen <= maxN; ++r) {
		if (!isRevLit(s, *r, maxL)) { return false; }
		notSeen += !s.seen(r->var());
	}
	if (contracted()) {
		const Literal* r = t.second;
		do { notSeen += !s.seen(r->var()); } while (notSeen <= maxN && !r++->watched());
	}
	return notSeen <= maxN;
}

bool Clause::simplify(Solver& s, bool reinit) {
	assert(s.decisionLevel() == 0 && s.queueSize() == 0 && !contracted());
	if (ClauseHead::satisfied(s)) {
		detach(s);
		return true;
	}
	LitRange t  = tail();
	Literal* it = !isSmall() ? head_+2 : t.first;
	Literal* j  = t.second;
	for (; it != t.second && s.value(it->var()) == value_free; ++it) { }
	if (it != t.second) {
		for (j = it; it != t.second; ++it) {
			if (s.isTrue(*it))  { Clause::detach(s); return true; }
			if (!s.isFalse(*it)){ *j++ = *it; }
		}
		for (Literal* r = j; r != t.second; ++r) { *r = negLit(0); }
	}
	if (!isSmall()) {
		uint32 size= std::max((uint32)ClauseHead::HEAD_LITS, static_cast<uint32>(j-head_));
		tail_.data.setSize(size);
		tail_.data.setIdx(3);
		if (reinit && size > 3) {
			detach(s);
			std::random_shuffle(head_, j, s.strategies().rng);
			attach(s);	
		}
	}
	else if (s.isFalse(head_[2])) {
		head_[2]   = t.first[0];
		t.first[0] = t.first[1];
		t.first[1] = negLit(0);
		--j;
	}
	return j <= t.first && ClauseHead::toImplication(s);
}

bool Clause::isSatisfied(const Solver& s, LitVec& freeLits) {
	if (ClauseHead::satisfied(s)) {
		return true;
	}
	freeLits.push_back(head_[0]);
	freeLits.push_back(head_[1]);
	if (!s.isFalse(head_[2])) freeLits.push_back(head_[2]);
	ValueRep v;
	LitRange t = tail();
	for (Literal* r = t.first; r != t.second; ++r) {
		if ( (v = s.value(r->var())) == value_free) {
			freeLits.push_back(*r);
		}
		else if (v == trueValue(*r)) {
			std::swap(head_[2], *r);
			return true;
		}
	}
	return false;
}

void Clause::undoLevel(Solver& s) {
	assert(!isSmall());
	uint32   t = tail_.data.size();
	Literal* r = head_+t;
	while (!r->watched() && s.value(r->var()) == value_free) {
		++t;
		++r;
	}
	if (r->watched() || s.level(r->var()) == 0) {
		r->clearWatch();
		t += !isSentinel(*r);
		tail_.data.clearContracted();
	}
	else {
		s.addUndoWatch(s.level(r->var()), this);
	}
	tail_.data.setSize(t);
}

void Clause::toLits(LitVec& out) const {
	out.insert(out.end(), head_, (head_+ClauseHead::HEAD_LITS)-isSentinel(head_[2]));
	LitRange t = const_cast<Clause&>(*this).tail();
	if (contracted()) {
		while (!t.second++->watched()) {;}
	}
	out.insert(out.end(), t.first, t.second);
}

ClauseHead::BoolPair Clause::strengthen(Solver& s, Literal p, bool toShort) {
	LitRange t  = tail();
	Literal* eoh= head_+ClauseHead::HEAD_LITS;
	Literal* eot= t.second;
	Literal* it = std::find(head_, eoh, p);
	BoolPair ret(false, false);
	if (it != eoh) {
		if (it == head_ || it == head_+1) {
			// update watch
			uint32 w = it != head_;
			head_[w] = head_[2];
			s.removeWatch(~p, this);
			Literal* best = head_+w;
			for (it = t.first; it != eot && s.isFalse(*best); ++it) {
				if (!s.isFalse(*it) || s.level(it->var()) > s.level(best->var())) { 
					best = it; 
				}
			}
			std::swap(head_[w], *best);
			s.addWatch(~head_[w], ClauseWatch(this));
		}	
		// replace cache literal with literal from tail
		head_[2]  = *t.first;
		eot       = removeFromTail(s, t.first, eot);
		ret.first = true;
	}
	else if ((it = std::find(t.first, eot, p)) != eot || (contracted() && *it == p)) {
		eot = removeFromTail(s, it, eot);
		ret.first = true;
	}
	else if (contracted()) {
		assert(*it != p);
		for (; !it++->watched();) {
			if (*it == p) {
				ret.first = true;
				Literal* c = it-1;
				while (!it->watched()) { *++c = *++it; }
				c->watch();
				break;
			}
		}
		eot = it;
	}
	if (ret.first && ~p == s.sharedContext()->tagLiteral()) {
		clearTagged();
	}
	ret.second = toShort && eot == t.first && toImplication(s);
	return ret;
}

Literal* Clause::removeFromTail(Solver& s, Literal* it, Literal* end) {
	if (!contracted()) {
		if (it != end) {
			*it = *--end;
			if (!isSmall()) { 
				tail_.data.setSize(tail_.data.size()-1);
				tail_.data.setIdx(3);
			}
		}
		*end = negLit(0);
		return end;
	}
	else if (!end->watched()) {
		Literal wOld = *end;
		Literal* j   = it;
		do { *j = *++it; } while (!j++->watched()); 
		*j = negLit(0);
		if (s.level(wOld.var()) != s.level(end->var())) {
			s.removeUndoWatch(s.level(wOld.var()), this);
			s.addUndoWatch(s.level(end->var()), this);
		}
		return j;
	}
	else {
		assert(it == end);
		s.removeUndoWatch(s.level(end->var()), this);
		tail_.data.clearContracted();
		*end = negLit(0);
		return end;
	}
}
uint32 Clause::size() const {
	return !isSmall() 
		? tail_.data.size()
		: 2 + !isSentinel(head_[2]) + !isSentinel(*smallLits()) + !isSentinel(smallLits()[1]);
}

/////////////////////////////////////////////////////////////////////////////////////////
// mt::SharedLitsClause
/////////////////////////////////////////////////////////////////////////////////////////
namespace mt {
ClauseHead* SharedLitsClause::newClause(Solver& s, SharedLiterals* shared_lits, const ClauseInfo& e, const Literal* lits, bool addRef) {
	return new (s.allocSmall()) SharedLitsClause(s, shared_lits, lits, e, addRef);
}

SharedLitsClause::SharedLitsClause(Solver& s, SharedLiterals* shared_lits, const Literal* w, const ClauseInfo& e, bool addRef) 
	: ClauseHead(e) {
	static_assert(sizeof(SharedLitsClause) <= 32);
	tail_.shared = addRef ? shared_lits->share() : shared_lits;
	std::memcpy(head_, w, std::min((uint32)ClauseHead::HEAD_LITS, shared_lits->size())*sizeof(Literal));
	attach(s);
}

Constraint* SharedLitsClause::cloneAttach(Solver& other) {
	return SharedLitsClause::newClause(other, tail_.shared, ClauseInfo().setType(this->type()), head_);
}

bool SharedLitsClause::updateWatch(Solver& s, uint32 pos) {
	Literal  other = head_[1^pos];
	for (const Literal* r = tail_.shared->begin(), *end = tail_.shared->end(); r != end; ++r) {
		// at this point we know that head_[2] is false so we only need to check 
		// that we do not watch the other watched literal twice!
		if (!s.isFalse(*r) && *r != other) {
			head_[pos] = *r; // replace watch
			// try to replace cache literal
			switch( std::min(static_cast<uint32>(8), static_cast<uint32>(end-r)) ) {
				case 8: if (!s.isFalse(*++r) && *r != other) { head_[2] = *r; return true; }
				case 7: if (!s.isFalse(*++r) && *r != other) { head_[2] = *r; return true; }
				case 6: if (!s.isFalse(*++r) && *r != other) { head_[2] = *r; return true; }
				case 5: if (!s.isFalse(*++r) && *r != other) { head_[2] = *r; return true; }
				case 4: if (!s.isFalse(*++r) && *r != other) { head_[2] = *r; return true; }
				case 3: if (!s.isFalse(*++r) && *r != other) { head_[2] = *r; return true; }
				case 2: if (!s.isFalse(*++r) && *r != other) { head_[2] = *r; return true; }
				default: return true;
			}
		}
	}
	return false;
}

void SharedLitsClause::reason(Solver& s, Literal p, LitVec& out) {
	ClauseHead::bumpActivity();
	LitVec::size_type i = out.size();
	for (const Literal* r = tail_.shared->begin(), *end = tail_.shared->end(); r != end; ++r) {
		assert(s.isFalse(*r) || *r == p);
		if (*r != p) { out.push_back(~*r); }
	}
	if (s.strategies().updateLbd && ClauseHead::type() != Constraint_t::static_constraint) {
		setLbd( s.computeLbd(p, &out[0]+i, &out[0]+out.size()) );
	}
}

bool SharedLitsClause::minimize(Solver& s, Literal p, CCMinRecursive* rec) {
	ClauseHead::bumpActivity();
	for (const Literal* r = tail_.shared->begin(), *end = tail_.shared->end(); r != end; ++r) {
		if (*r != p && !s.ccMinimize(~*r, rec)) { return false; }
	}
	return true;
}

bool SharedLitsClause::isReverseReason(const Solver& s, Literal p, uint32 maxL, uint32 maxN) {
	uint32 notSeen = 0;
	for (const Literal* r = tail_.shared->begin(), *end = tail_.shared->end(); r != end; ++r) {
		if (*r == p) continue;
		if (!isRevLit(s, *r, maxL)) return false;
		if (!s.seen(r->var()) && ++notSeen > maxN) return false;
	}
	return true;
}

bool SharedLitsClause::simplify(Solver& s, bool reinit) {
	if (ClauseHead::satisfied(s)) {
		detach(s);
		return true;
	}
	uint32 optSize = tail_.shared->simplify(s);
	if (optSize == 0) {
		detach(s);
		return true;
	}
	else if (optSize <= Clause::MAX_SHORT_LEN) {
		Literal lits[Clause::MAX_SHORT_LEN];
		Literal* j = lits;
		for (const Literal* r = tail_.shared->begin(), *e = tail_.shared->end(); r != e; ++r) {
			if (!s.isFalse(*r)) *j++ = *r;
		}
		uint32 rep= info_.rep;
		// detach & destroy but do not release memory
		detach(s);
		SharedLitsClause::destroy(0, false);
		// construct short clause in "this"
		ClauseHead* h = Clause::newShortClause(s, this, lits, static_cast<uint32>(j-lits), ClauseInfo());
		// restore extra data - safe because memory layout has not changed!
		info_.rep = rep;
		return h->simplify(s, reinit);
	}
	else if (s.isFalse(head_[2])) {
		// try to replace cache lit with non-false literal
		for (const Literal* r = tail_.shared->begin(), *e = tail_.shared->end(); r != e; ++r) {
			if (!s.isFalse(*r) && std::find(head_, head_+2, *r) == head_+2) {
				head_[2] = *r;
				break;
			}
		}
	}
	return false;
}

void SharedLitsClause::destroy(Solver* s, bool detachFirst) {
	if (s && detachFirst) {
		ClauseHead::detach(*s);
	}
	tail_.shared->release();
	void* mem = this;
	this->~SharedLitsClause();
	if (s) {
		// return memory to solver
		s->freeSmall(mem);
	}
}

bool SharedLitsClause::isSatisfied(const Solver& s, LitVec& freeLits) {
	if (ClauseHead::satisfied(s)) {
		return true;
	}
	Literal* head = head_;
	ValueRep v;
	for (const Literal* r = tail_.shared->begin(), *end = tail_.shared->end(); r != end; ++r) {
		if ( (v = s.value(r->var())) == value_free ) {
			freeLits.push_back(*r);
		}
		else if (v == trueValue(*r)) {
			head[2] = *r; // remember as cache literal
			return true;
		}
	}
	return false;
}

void SharedLitsClause::toLits(LitVec& out) const {
	out.insert(out.end(), tail_.shared->begin(), tail_.shared->end());
}

ClauseHead::BoolPair SharedLitsClause::strengthen(Solver&, Literal, bool) {
	return BoolPair(false, false);
}

uint32 SharedLitsClause::size() const { return tail_.shared->size(); }

} // end namespace mt

/////////////////////////////////////////////////////////////////////////////////////////
// LoopFormula
/////////////////////////////////////////////////////////////////////////////////////////
LoopFormula::LoopFormula(Solver& s, uint32 size, Literal* bodyLits, uint32 numBodies, uint32 bodyToWatch) {
	activity_       = (uint32)s.stats.restarts + (size-numBodies);  
	end_            = numBodies + 2;
	size_           = end_+1;
	other_          = end_-1;
	lits_[0]        = Literal();  // Starting sentinel
	lits_[end_-1]   = Literal();  // Position of active atom
	lits_[end_]     = Literal();  // Ending sentinel - active part
	for (uint32 i = size_; i != size+3; ++i) {
		lits_[i] = Literal();
	}

	// copy bodies: S B1...Bn, watch one
	std::memcpy(lits_+1, bodyLits, numBodies * sizeof(Literal));
	s.addWatch(~lits_[1+bodyToWatch], this, ((1+bodyToWatch)<<1)+1);
	lits_[1+bodyToWatch].watch();
}

void LoopFormula::destroy(Solver* s, bool detach) {
	if (s && detach) {
		for (uint32 x = 1; x != end_-1; ++x) {
			if (lits_[x].watched()) {
				s->removeWatch(~lits_[x], this);
				lits_[x].clearWatch();
			}
		}
		if (lits_[end_-1].watched()) {
			lits_[end_-1].clearWatch();
			for (uint32 x = end_+1; x != size_; ++x) {
				s->removeWatch(~lits_[x], this);
				lits_[x].clearWatch();
			}
		}
	}
	void* mem = static_cast<Constraint*>(this);
	this->~LoopFormula();
	::operator delete(mem);
}


void LoopFormula::addAtom(Literal atom, Solver& s) {
	uint32 pos = size_++;
	assert(isSentinel(lits_[pos]));
	lits_[pos] = atom;
	lits_[pos].watch();
	s.addWatch( ~lits_[pos], this, (pos<<1)+0 );
	if (isSentinel(lits_[end_-1])) {
		lits_[end_-1] = lits_[pos];
	}
}

void LoopFormula::updateHeuristic(Solver& s) {
	Literal saved = lits_[end_-1];
	for (uint32 x = end_+1; x != size_; ++x) {
		lits_[end_-1] = lits_[x];
		s.strategies().heuristic->newConstraint(s, lits_+1, end_-1, Constraint_t::learnt_loop);
	}
	lits_[end_-1] = saved;
}

bool LoopFormula::watchable(const Solver& s, uint32 idx) {
	assert(!lits_[idx].watched());
	if (idx == end_-1) {
		for (uint32 x = end_+1; x != size_; ++x) {
			if (s.isFalse(lits_[x])) {
				lits_[idx] = lits_[x];
				return false;
			}
		}
	}
	return true;
}

bool LoopFormula::isTrue(const Solver& s, uint32 idx) {
	if (idx != end_-1) return s.isTrue(lits_[idx]);
	for (uint32 x = end_+1; x != size_; ++x) {
		if (!s.isTrue(lits_[x])) {
			lits_[end_-1] = lits_[x];
			return false;
		}
	}
	return true;
}

Constraint::PropResult LoopFormula::propagate(Solver& s, Literal, uint32& data) {
	if (isTrue(s, other_)) {          // ignore clause, as it is 
		return PropResult(true, true);  // already satisfied
	}
	uint32  pos   = data >> 1;
	uint32  idx   = pos;
	if (pos > end_) {
		// p is one of the atoms - move to active part
		lits_[end_-1] = lits_[pos];
		idx           = end_-1;
	}
	int     dir   = ((data & 1) << 1) - 1;
	int     bounds= 0;
	for (;;) {
		for (idx+=dir;s.isFalse(lits_[idx]);idx+=dir) {;} // search non-false literal - sentinels guarantee termination
		if (isSentinel(lits_[idx])) {             // Hit a bound,
			if (++bounds == 2) {                    // both ends seen, clause is unit, false, or sat
				if (other_ == end_-1) {
					uint32 x = end_+1;
					for (; x != size_ && s.force(lits_[x], this);  ++x) { ; }
					return Constraint::PropResult(x == size_, true);  
				}
				else {
					return Constraint::PropResult(s.force(lits_[other_], this), true);  
				}
			}
			idx   = std::min(pos, end_-1);          // halfway through, restart search, but
			dir   *= -1;                            // this time walk in the opposite direction.
			data  ^= 1;                             // Save new direction of watch
		}
		else if (!lits_[idx].watched() && watchable(s, idx)) { // found a new watchable literal
			if (pos > end_) {     // stop watching atoms
				lits_[end_-1].clearWatch();
				for (uint32 x = end_+1; x != size_; ++x) {
					if (x != pos) {
						s.removeWatch(~lits_[x], this);
						lits_[x].clearWatch();
					}
				}
			}
			lits_[pos].clearWatch();
			lits_[idx].watch();
			if (idx == end_-1) {  // start watching atoms
				for (uint32 x = end_+1; x != size_; ++x) {
					s.addWatch(~lits_[x], this, static_cast<uint32>(x << 1) + 0);
					lits_[x].watch();
				}
			}
			else {
				s.addWatch(~lits_[idx], this, static_cast<uint32>(idx << 1) + (dir==1));
			}
			return Constraint::PropResult(true, false);
		} 
		else if (lits_[idx].watched()) {          // Hit the other watched literal
			other_  = idx;                          // Store it in other_
		}
	}
}

// Body: all other bodies + active atom
// Atom: all bodies
void LoopFormula::reason(Solver&, Literal p, LitVec& lits) {
	// all relevant bodies
	for (uint32 x = 1; x != end_-1; ++x) {
		if (lits_[x] != p) {
			lits.push_back(~lits_[x]);
		}
	}
	// if p is a body, add active atom
	if (other_ != end_-1) {
		lits.push_back(~lits_[end_-1]);
	}
	++activity_;
}

bool LoopFormula::minimize(Solver& s, Literal p, CCMinRecursive* rec) {
	++activity_;
	for (uint32 x = 1; x != end_-1; ++x) {
		if (lits_[x] != p && !s.ccMinimize(~lits_[x], rec)) {
			return false;
		}
	}
	return other_ == end_-1
		|| s.ccMinimize(~lits_[end_-1], rec);
}

uint32 LoopFormula::size() const {
	return size_ - 3;
}

bool LoopFormula::locked(const Solver& s) const {
	if (other_ != end_-1) {
		return s.isTrue(lits_[other_]) && s.reason(lits_[other_]) == this;
	}
	for (uint32 x = end_+1; x != size_; ++x) {
		if (s.isTrue(lits_[x]) && s.reason(lits_[x]) == this) {
			return true;
		}
	}
	return false;
}

bool LoopFormula::isSatisfied(const Solver& s, LitVec& freeLits) {
	if (other_ != end_-1 && s.isTrue(lits_[other_])) return true;
	for (uint32 x = 1; x != end_-1; ++x) {
		if (s.isTrue(lits_[x])) {
			other_ = x;
			return true;
		}
		else if (!s.isFalse(lits_[x])) { freeLits.push_back(lits_[x]); }
	}
	bool sat = true;
	for (uint32 x = end_+1; x != size_; ++x) {
		if (s.value(lits_[x].var()) == value_free) {
			freeLits.push_back(lits_[x]);
			sat = false;
		}
		else sat &= s.isTrue(lits_[x]);
	}
	return sat;
}

bool LoopFormula::simplify(Solver& s, bool) {
	assert(s.decisionLevel() == 0);
	typedef std::pair<uint32, uint32> WatchPos;
	bool      sat = false;          // is the constraint SAT?
	WatchPos  bodyWatches[2];       // old/new position of watched bodies
	uint32    bw  = 0;              // how many bodies are watched?
	uint32    j   = 1, i;
	// 1. simplify the set of bodies:
	// - search for a body that is true -> constraint is SAT
	// - remove all false bodies
	// - find the watched bodies
	for (i = 1; i != end_-1; ++i) {
		assert( !s.isFalse(lits_[i]) || !lits_[i].watched() ); // otherwise should be asserting 
		if (!s.isFalse(lits_[i])) {
			sat |= s.isTrue(lits_[i]);
			if (i != j) { lits_[j] = lits_[i]; }
			if (lits_[j].watched()) { bodyWatches[bw++] = WatchPos(i, j); }
			++j;
		}
	}
	uint32  newEnd    = j + 1;
	uint32  numBodies = j - 1;
	j += 2;
	// 2. simplify the set of atoms:
	// - remove all determined atoms
	// - remove/update watches if necessary
	for (i = end_ + 1; i != size_; ++i) {
		if (s.value(lits_[i].var()) == value_free) {
			if (i != j) { lits_[j] = lits_[i]; }
			if (lits_[j].watched()) {
				if (sat || numBodies <= 2) {
					s.removeWatch(~lits_[j], this);
					lits_[j].clearWatch();
				}
				else if (i != j) {
					GenericWatch* w = s.getWatch(~lits_[j], this);
					assert(w);
					w->data = (j << 1) + 0;
				}
			}
			++j;
		}
		else if (lits_[i].watched()) {
			s.removeWatch(~lits_[i], this);
			lits_[i].clearWatch();
		}
	}
	size_         = j;
	end_          = newEnd;
	lits_[end_]   = Literal();
	lits_[end_-1] = lits_[end_+1];
	if (sat || numBodies < 3 || size_ == end_ + 1) {
		for (i = 0; i != bw; ++i) {
			s.removeWatch(~lits_[bodyWatches[i].second], this);
			lits_[bodyWatches[i].second].clearWatch();
		}
		if (sat || size_ == end_+1) { return true; }
		// replace constraint with short clauses
		ClauseCreator creator(&s);
		creator.start();
		for (i = 1; i != end_; ++i) { creator.add(lits_[i]); }
		for (i = end_+1; i != size_; ++i) {
			creator[creator.size()-1] = lits_[i];
			creator.end();
		}
		return true;
	}
	other_ = 1;
	for (i = 0; i != bw; ++i) {
		if (bodyWatches[i].first != bodyWatches[i].second) {
			GenericWatch* w  = s.getWatch(~lits_[bodyWatches[i].second], this);
			assert(w);
			w->data = (bodyWatches[i].second << 1) + (w->data&1);
		}
	}
	return false;
}
}
