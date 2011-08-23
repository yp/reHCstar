// 
// Copyright (c) 2006-2010, Benjamin Kaufmann
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
#include <clasp/solver.h>
#include <clasp/clause.h>

#ifdef _MSC_VER
#pragma warning (disable : 4996) // 'std::copy': Function call with parameters that may be unsafe
#endif

namespace Clasp { 

DecisionHeuristic::~DecisionHeuristic() {}

/////////////////////////////////////////////////////////////////////////////////////////
// SelectFirst selection strategy
/////////////////////////////////////////////////////////////////////////////////////////
// selects the first free literal
Literal SelectFirst::doSelect(Solver& s) {
	for (Var i = 1; i <= s.numVars(); ++i) {
		if (s.value(i) == value_free) {
			return s.sharedContext()->preferredLiteralByType(i);
		}
	}
	assert(!"SelectFirst::doSelect() - precondition violated!\n");
	return Literal();
}
/////////////////////////////////////////////////////////////////////////////////////////
// SelectRandom selection strategy
/////////////////////////////////////////////////////////////////////////////////////////
// Selects a random literal from all free literals.
class SelectRandom : public DecisionHeuristic {
public:
	SelectRandom() : randFreq_(1.0), pos_(0) {}
	void shuffle(RNG& rng) {
		std::random_shuffle(vars_.begin(), vars_.end(), rng);
		pos_ = 0;
	}
	void   randFreq(double d) { randFreq_ = d; }
	double randFreq() const   { return randFreq_; }
	void   endInit(Solver& s) {
		vars_.clear();
		for (Var i = 1; i <= s.numVars(); ++i) {
			if (s.value( i ) == value_free) {
				vars_.push_back(i);
			}
		}
		pos_ = 0;
	}
private:
	Literal doSelect(Solver& s) {
		LitVec::size_type old = pos_;
		do {
			if (s.value(vars_[pos_]) == value_free) {
				Literal l = savedLiteral(s, vars_[pos_]);
				return l != posLit(0)
					? l
					: s.sharedContext()->preferredLiteralByType(vars_[pos_]);
			}
			if (++pos_ == vars_.size()) pos_ = 0;
		} while (old != pos_);
		assert(!"SelectRandom::doSelect() - precondition violated!\n");
		return Literal();
	}
	VarVec            vars_;
	double            randFreq_;
	VarVec::size_type pos_;
};
/////////////////////////////////////////////////////////////////////////////////////////
// Post propagator list
/////////////////////////////////////////////////////////////////////////////////////////
Solver::PPList::PPList() : head(0), look(0), saved(0) { }
Solver::PPList::~PPList() {
	for (PostPropagator* r = head; r;) {
		PostPropagator* t = r;
		r = r->next;
		t->destroy();
	}
}
void Solver::PPList::add(PostPropagator* p) {
	assert(p && p->next == 0);
	uint32 prio = p->priority();
	if (!head || prio == PostPropagator::priority_highest || prio < head->priority()) {
		p->next = head;
		head    = p;
	}
	else {
		for (PostPropagator* r = head; ; r = r->next) {
			if (r->next == 0) {
				r->next = p; break;
			}
			else if (prio < r->next->priority()) {
				p->next = r->next;
				r->next = p;
				break;
			}				
		}
	}
	if (prio >= PostPropagator::priority_lookahead && (!look || prio < look->priority())) {
		look = p;
	}
	saved = head;
}
void Solver::PPList::remove(PostPropagator* p) {
	assert(p);
	if (!head) return;
	if (p == head) {
		head = head->next;
		p->next = 0;
	}
	else {
		for (PostPropagator* r = head; ; r = r->next) {
			if (r->next == p) {
				r->next = r->next->next;
				p->next = 0;
				break;
			}
		}
	}
	if (p == look) { look = look->next; }
	saved = head;
}
bool Solver::PPList::propagate(Solver& s, PostPropagator* x) {
	for (PostPropagator* p = head, *t; p != x;) {
		// just in case t removes itself from the list
		// during propagateFixpoint
		t = p;
		p = p->next;
		if (!t->propagateFixpoint(s)) { return false; }
	}
	return true;
}
void Solver::PPList::simplify(Solver& s, bool shuf) {
	for (PostPropagator* r = head; r;) {
		PostPropagator* t = r;
		r = r->next;
		if (t->simplify(s, shuf)) {
			remove(t);
		}
	}
}
bool Solver::PPList::init(Solver& s) {
	for (PostPropagator* r = head; r; r = r->next) {
		if (!r->init(s)) return false;
	}
	return true;
}
bool PostPropagator::propagateFixpoint(Solver& s) {
	bool ok = propagate(s);
	while (ok && s.queueSize() > 0) {
		ok = s.propagateUntil(this) && propagate(s);
	}
	return ok;
}
void Solver::PPList::reset()            { for (PostPropagator* r = head; r; r = r->next) { r->reset(); } }
bool Solver::PPList::isModel(Solver& s) {
	saved = head;
	for (PostPropagator* r = head; r; r = r->next) {
		if (!r->isModel(s)) return false;
	}
	return true;
}
bool Solver::PPList::nextSymModel(Solver& s, bool expand) {
	for (; saved; saved = saved->next) {
		if (saved->nextSymModel(s, expand)) return true;
	}
	saved = head;
	return false;
}
/////////////////////////////////////////////////////////////////////////////////////////
// SolverStrategies
/////////////////////////////////////////////////////////////////////////////////////////
SolverStrategies::SolverStrategies()
	: heuristic(new SelectFirst) 
	, search(use_learning)
	, saveProgress(0) 
	, reverseArcs(0)
	, otfs(0)
	, cflMinAntes(all_antes)
	, reduceAlgo(reduce_linear)
	, reduceScore(score_act)
	, reduceGlue(0)
	, strengthenRecursive(false)
	, randomWatches(false)
	, updateLbd(false)
	, compress_(250) {
}
/////////////////////////////////////////////////////////////////////////////////////////
// Solver: Construction/Destruction/Setup
////////////////////////////////////////////////////////////////////////////////////////
Solver::Solver() 
	: strategy_()
	, ccMin_(0)
	, smallAlloc_(new SmallClauseAlloc)
	, shared_(0)
	, randHeuristic_(0)
	, levConflicts_(0)
	, undoHead_(0)
	, enum_(0)
	, id_(0)
	, units_(0)
	, lastSimplify_(0)
	, rootLevel_(0)
	, btLevel_(0)
	, tagged_(0)
	, lbdTime_(0)
	, shuffle_(false) {
	Var sentVar = assign_.addVar();
	assign_.setValue(sentVar, value_true);
	markSeen(sentVar);
}

Solver::~Solver() {
	freeMem();
}

void Solver::freeMem() {
	std::for_each( constraints_.begin(), constraints_.end(), DestroyObject());
	std::for_each( learnts_.begin(), learnts_.end(), DestroyObject() );
	constraints_.clear();
	learnts_.clear();
	setEnumerationConstraint(0);
	PodVector<WatchList>::destruct(watches_);
	// free undo lists
	// first those still in use
	for (DecisionLevels::size_type i = 0; i != levels_.size(); ++i) {
		delete levels_[i].undo;
	}
	// then those in the free list
	for (ConstraintDB* x = undoHead_; x; ) {
		ConstraintDB* t = x;
		x = (ConstraintDB*)x->front();
		delete t;
	}
	delete smallAlloc_;
	delete ccMin_;
	delete levConflicts_;
	delete randHeuristic_;
	smallAlloc_   = 0;
	ccMin_        = 0;
	levConflicts_ = 0;
	randHeuristic_= 0;
}

SatPreprocessor* Solver::satPrepro() const { return shared_->satPrepro.get(); }

void Solver::reset() {
	// hopefully, no one derived from this class...
	this->~Solver();
	new (this) Solver();
}

void Solver::startInit(SharedContext* ctx, uint32 numConsGuess) {
	shared_ = ctx;
	assert(numVars() <= shared_->numVars());
	assign_.resize(shared_->numVars() + 1);
	watches_.resize(assign_.numVars()<<1);
	// pre-allocate some memory
	assign_.trail.reserve(numVars());
	constraints_.reserve(numConsGuess/2);
	levels_.reserve(25);
	if (smallAlloc_ == 0) {
		smallAlloc_ = new SmallClauseAlloc();
	}
	if (undoHead_ == 0) {
		for (uint32 i = 0; i != 25; ++i) { 
			undoFree(new ConstraintDB(10)); 
		}
	}
	strategy_.heuristic->startInit(*this);
}

bool Solver::endInit() {
	if (!propagate()) { return false; }
	strategy_.heuristic->endInit(*this);
	// Force propagation again - just in case
	// heuristic->endInit() added new information
	if (!propagate() || !simplify()) {
		return false;
	}
	if (randHeuristic_) randHeuristic_->endInit(*this);
	return true;
}

bool Solver::addUnary(Literal p, ConstraintType t) {
	undoUntil(0);
	const Antecedent ante(posLit(0));
	if (decisionLevel() != 0) {
		// remove any existing records
		for (ImpliedLits::size_type i = 0; i < impliedLits_.size(); ++i) {
			if (impliedLits_[i].lit == p) {
				impliedLits_[i] = impliedLits_.back();
				impliedLits_.pop_back();
				--i;
			}
		}
		// remember implication
		impliedLits_.push_back(ImpliedLiteral(p, 0, ante));
		if (isTrue(p)) {
			// replace antecedent with stronger constraint
			setReason(p, ante);
		}
	}
	bool ret = force(p, ante);
	if (t != Constraint_t::static_constraint && shared_->learnBinary(p, negLit(0))) {
		stats.addLearnt(1, t);
		shared_->distribute(*this, p, negLit(0), negLit(0), t);	
	}
	return ret;
}

bool Solver::addBinary(Literal p, Literal q, ConstraintType t) {
	assert(validWatch(~p) && validWatch(~q) && "ERROR: startAddConstraints not called!");
	if (t == Constraint_t::static_constraint) {
		shared_->addBinary(p, q);
	}
	else if (t != Constraint_t::static_constraint && shared_->learnBinary(p, q)) {
		stats.addLearnt(2, t);
		shared_->distribute(*this, p, q, negLit(0), t);
	}
	return true;
}

bool Solver::addTernary(Literal p, Literal q, Literal r, ConstraintType t) {
	assert(validWatch(~p) && validWatch(~q) && validWatch(~r) && "ERROR: startAddConstraints not called!");
	assert(p != q && q != r && "ERROR: ternary clause contains duplicate literal");
	if (t == Constraint_t::static_constraint) {
		shared_->addTernary(p, q, r);
	}
	else if (t != Constraint_t::static_constraint && shared_->learnTernary(p, q, r)) {
		stats.addLearnt(3, t);
		shared_->distribute(*this, p, q, r, t);
	}
	return true;
}

void Solver::add(Constraint* c) {
	assert(shared_->master() == this);
	shared_->add(c);
}

uint32 Solver::numConstraints() const {
	return static_cast<uint32>(constraints_.size())
		+ (shared_ ? shared_->numBinary()+shared_->numTernary() : 0);
}

bool Solver::popRootLevel(uint32 i, bool resolve)  {
	clearStopConflict();
	if (i > rootLevel_) i = rootLevel_;
	rootLevel_ = btLevel_ = rootLevel_ - i;
	// first, resolve any open conflicts
	if (hasConflict() && resolve && !resolveConflict()) {
		return false;
	}
	// second, go back to new root level and try to re-assert still implied literals
	undoUntil(rootLevel_);
	if (!isTrue(sharedContext()->tagLiteral())) {
		removeConditional();
	}
	return forceImplied();
}

bool Solver::clearAssumptions()  {
	return popRootLevel(rootLevel(), false)
		&& simplify();
}

void Solver::clearStopConflict() {
	if (hasConflict() && hasStopConflict()) {
		rootLevel_ = conflict_[1].asUint();
		conflict_.clear();
	}
}

bool Solver::hasStopConflict() const { return conflict_[0] == negLit(0); }
void Solver::setStopConflict() {
	if (!hasConflict()) {
		// we use the nogood {FALSE} to represent the unrecoverable conflict -
		// note that {FALSE} can otherwise never be a violated nogood because
		// TRUE is always true in every solver
		conflict_.push_back(negLit(0));
		// remember the current root-level
		conflict_.push_back(Literal::fromRep(rootLevel_));
	}
	// artificially increase root level -
	// this way, the solver is prevented from resolving the conflict
	pushRootLevel(decisionLevel());
}

void Solver::updateGuidingPath(LitVec& gpOut, LitVec::size_type& start, uint32& implied) {
	if (start < shared_->topLevelSize()) {
		start = shared_->topLevelSize();
	}
	const LitVec& tr = assign_.trail;
	LitVec::size_type end = rootLevel_ == decisionLevel() ? tr.size() : levels_[rootLevel_].trailPos;
	for (LitVec::size_type i = start; i < end; ++i) {
		if (reason(tr[i]).isNull()) {
			gpOut.push_back(tr[i]);
		}
	}
	start = end;
	uint32 implLits = 0;
	for (LitVec::size_type i = 0; i != impliedLits_.size(); ++i) {
		if (impliedLits_[i].level <= rootLevel_ && ++implLits > implied) {
			gpOut.push_back(impliedLits_[i].lit);
		}
	}
	implied = implLits;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Solver: Watch management
////////////////////////////////////////////////////////////////////////////////////////
uint32 Solver::numWatches(Literal p) const {
	assert( validVar(p.var()) );
	if (!validWatch(p)) return 0;
	return static_cast<uint32>(watches_[p.index()].size()) 
		+ shared_->numShortImplications(p);
}
	
bool Solver::hasWatch(Literal p, Constraint* c) const {
	if (!validWatch(p)) return false;
	const WatchList& pList = watches_[p.index()];
	return std::find_if(pList.right_begin(), pList.right_end(), GenericWatch::EqConstraint(c)) != pList.right_end();
}

bool Solver::hasWatch(Literal p, ClauseHead* h) const {
	if (!validWatch(p)) return false;
	const WatchList& pList = watches_[p.index()];
	return std::find_if(pList.left_begin(), pList.left_end(), ClauseWatch::EqHead(h)) != pList.left_end();
}

GenericWatch* Solver::getWatch(Literal p, Constraint* c) const {
	if (!validWatch(p)) return 0;
	const WatchList& pList = watches_[p.index()];
	WatchList::const_right_iterator it = std::find_if(pList.right_begin(), pList.right_end(), GenericWatch::EqConstraint(c));
	return it != pList.right_end()
		? &const_cast<GenericWatch&>(*it)
		: 0;
}

void Solver::removeWatch(const Literal& p, Constraint* c) {
	assert(validWatch(p));
	WatchList& pList = watches_[p.index()];
	pList.erase_right(std::find_if(pList.right_begin(), pList.right_end(), GenericWatch::EqConstraint(c)));
}

void Solver::removeWatch(const Literal& p, ClauseHead* h) {
	assert(validWatch(p));
	WatchList& pList = watches_[p.index()];
	pList.erase_left(std::find_if(pList.left_begin(), pList.left_end(), ClauseWatch::EqHead(h)));
}

void Solver::removeUndoWatch(uint32 dl, Constraint* c) {
	assert(dl != 0 && dl <= decisionLevel() );
	if (levels_[dl-1].undo) {
		ConstraintDB& uList = *levels_[dl-1].undo;
		ConstraintDB::iterator it = std::find(uList.begin(), uList.end(), c);
		if (it != uList.end()) {
			*it = uList.back();
			uList.pop_back();
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// Solver: Basic DPLL-functions
////////////////////////////////////////////////////////////////////////////////////////
void Solver::initRandomHeuristic(double randFreq) {
	randFreq = std::min(1.0, std::max(0.0, randFreq));
	if (randFreq == 0.0) {
		delete randHeuristic_;
		randHeuristic_ = 0;
		return;
	}
	if (!randHeuristic_) {
		randHeuristic_ = new SelectRandom();
		randHeuristic_->endInit(*this);
	}
	static_cast<SelectRandom*>(randHeuristic_)->shuffle(strategies().rng);
	static_cast<SelectRandom*>(randHeuristic_)->randFreq(randFreq);
}


// removes all satisfied binary and ternary clauses as well
// as all constraints for which Constraint::simplify returned true.
bool Solver::simplify() {
	if (decisionLevel() != 0) return true;
	if (hasConflict())        return false;
	if (lastSimplify_ != assign_.trail.size()) {
		LitVec::size_type old = lastSimplify_;
		if (!simplifySAT()) { return false; }
		assert(lastSimplify_ == assign_.trail.size());
		strategy_.heuristic->simplify(*this, old);
	}
	if (shuffle_) { simplifySAT(); }
	return true;
}

void Solver::removeConditional() { 
	Literal p = ~sharedContext()->tagLiteral();
	if (!isSentinel(p) && tagged_) {
		ConstraintDB::size_type i, j, end = learnts_.size();
		for (i = j = 0; i != end; ++i) {
			ClauseHead* c = static_cast<LearntConstraint*>(learnts_[i])->clause();
			if (!c || !c->tagged()) {
				learnts_[j++] = c;
			}
			else {
				assert((decisionLevel() == rootLevel() || !c->locked(*this)) && "Solver::removeConditional(): must not remove locked constraint!");
				c->destroy(this, true);
			}
		}
		learnts_.erase(learnts_.begin()+j, learnts_.end());
		tagged_ = 0;
	}
}

void Solver::strengthenConditional() { 
	Literal p = ~sharedContext()->tagLiteral();
	if (!isSentinel(p) && tagged_) {
		ConstraintDB::size_type i, j, end = learnts_.size();
		for (i = j = 0; i != end; ++i) {
			ClauseHead* c = static_cast<LearntConstraint*>(learnts_[i])->clause();
			if (!c || !c->tagged() || !c->strengthen(*this, p, true).second) {
				learnts_[j++] = c;
			}
			else {
				assert((decisionLevel() == rootLevel() || !c->locked(*this)) && "Solver::strengthenConditional(): must not remove locked constraint!");
				c->destroy(this, false);
			}
		}
		learnts_.erase(learnts_.begin()+j, learnts_.end());
		tagged_ = 0;
	}
}

bool Solver::simplifySAT() {
	if (queueSize() > 0 && !propagate()) {
		return false;
	}
	assert(assign_.qEmpty());
	assign_.front = lastSimplify_;
	while (!assign_.qEmpty()) {
		simplifyShort(assign_.qPop()); // remove satisfied binary- and ternary clauses
	}
	lastSimplify_ = assign_.front;
	if (shuffle_) {
		std::random_shuffle(constraints_.begin(), constraints_.end(), strategies().rng);
		std::random_shuffle(learnts_.begin(), learnts_.end(), strategies().rng);
	}
	simplifyDB(constraints_);
	simplifyDB(learnts_);
	post_.simplify(*this, shuffle_);
	if (enum_ && enum_->simplify(*this, shuffle_)) {
		enum_->destroy(this, false);
		enum_ = 0;
	}
	shuffle_ = false;
	return true;
}

void Solver::simplifyDB(ConstraintDB& db) {
	ConstraintDB::size_type i, j, end = db.size();
	for (i = j = 0; i != end; ++i) {
		Constraint* c = db[i];
		if (c->simplify(*this, shuffle_)) { c->destroy(this, false); }
		else                              { db[j++] = c;  }
	}
	db.erase(db.begin()+j, db.end());
}

// removes all binary clauses containing p - those are now SAT
// binary clauses containing ~p are unit and therefore likewise SAT. Those
// are removed when their second literal is processed.
// Note: Binary clauses containing p are those that watch ~p.
//
// Simplifies ternary clauses.
// Ternary clauses containing p are SAT and therefore removed.
// Ternary clauses containing ~p are now either binary or SAT. Those that
// are SAT are removed when the satisfied literal is processed. 
// All conditional binary-clauses are replaced with a real binary clause.
// Note: Ternary clauses containing p watch ~p. Those containing ~p watch p.
// Note: Those clauses are now either binary or satisfied.
void Solver::simplifyShort(Literal p) {
	releaseVec(watches_[p.index()]);
	releaseVec(watches_[(~p).index()]);
	if (shared_->unique()) {
		shared_->simplifyBtig(p);
	}
}

void Solver::setConflict(Literal p, const Antecedent& a, uint32 data) {
	conflict_.push_back(~p);
	if (strategy_.search != SolverStrategies::no_learning && !a.isNull()) {
		if (data == UINT32_MAX) {
			a.reason(*this, p, conflict_);
		}
		else {
			// temporarily replace old data with new data
			uint32 saved = assign_.data(p.var());
			assign_.setData(p.var(), data);
			// extract conflict using new data
			a.reason(*this, p, conflict_);
			// restore old data
			assign_.setData(p.var(), saved);
		}
	}
}

bool Solver::assume(const Literal& p) {
	if (value(p.var()) == value_free) {
		assert(decisionLevel() != assign_.maxLevel());
		++stats.choices;
		levels_.push_back(DLevel(numAssignedVars(), 0));
		if (levConflicts_) levConflicts_->push_back( (uint32)stats.conflicts );
		return assign_.assign(p, decisionLevel(), Antecedent());
	}
	return isTrue(p);
}

bool Solver::propagate() {
	if (unitPropagate() && post_.propagate(*this, 0)) {
		assert(queueSize() == 0);
		return true;
	}
	assign_.qReset();
	post_.reset();
	return false;
}

uint32 Solver::mark(uint32 s, uint32 e) {
	while (s != e) { markSeen(assign_.trail[s++].var()); }
	return e;
}

Constraint::PropResult ClauseHead::propagate(Solver& s, Literal p, uint32&) {
	Literal* head = head_;
	uint32 wLit   = (head[1] == ~p); // pos of false watched literal
	if (s.isTrue(head[1-wLit])) {
		return Constraint::PropResult(true, true);
	}
	else if (!s.isFalse(head[2])) {
		assert(!isSentinel(head[2]) && "Invalid ClauseHead!");
		head[wLit] = head[2];
		head[2]    = ~p;
		s.addWatch(~head[wLit], ClauseWatch(this));
		return Constraint::PropResult(true, false);
	}
	else if (updateWatch(s, wLit)) {
		assert(!s.isFalse(head_[wLit]));
		s.addWatch(~head[wLit], ClauseWatch(this));
		return Constraint::PropResult(true, false);
	}	
	return PropResult(s.force(head_[1^wLit], this), true);
}

bool Solver::unitPropagate() {
	assert(!hasConflict());
	Literal p, q, r;
	uint32 idx, ignore;
	const ShortImplicationsGraph& btig = shared_->shortImplications();
	while ( !assign_.qEmpty() ) {
		p             = assign_.qPop();
		idx           = p.index();
		WatchList& wl = watches_[idx];
		// first: short clause BCP
		if (!btig.propagate(*this, p)) {
			return false;
		}
		// second: clause BCP
		if (wl.left_size() != 0) {
			WatchList::left_iterator it, end, j = wl.left_begin(); 
			Constraint::PropResult res;
			for (it = wl.left_begin(), end = wl.left_end();  it != end;  ) {
				ClauseWatch& w = *it++;
				res = w.head->ClauseHead::propagate(*this, p, ignore);
				if (res.keepWatch) {
					*j++ = w;
				}
				if (!res.ok) {
					wl.shrink_left(std::copy(it, end, j));
					return false;
				}
			}
			wl.shrink_left(j);
		}
		// third: general constraint BCP
		if (wl.right_size() != 0) {
			WatchList::right_iterator it, end, j = wl.right_begin(); 
			Constraint::PropResult res;
			for (it = wl.right_begin(), end = wl.right_end(); it != end; ) {
				GenericWatch& w = *it++;
				res = w.propagate(*this, p);
				if (res.keepWatch) {
					*j++ = w;
				}
				if (!res.ok) {
					wl.shrink_right(std::copy(it, end, j));
					return false;
				}
			}
			wl.shrink_right(j);
		}
	}
	return decisionLevel() > 0 || (units_=mark(units_, uint32(assign_.front))) == assign_.front;
}

bool Solver::test(Literal p, Constraint* c) {
	assert(value(p.var()) == value_free && !hasConflict());
	assume(p); --stats.choices;
	freezeLevel(decisionLevel()); // can't split-off this level
	if (unitPropagate() && post_.propagate(*this, post_.look)) {
		assert(decision(decisionLevel()) == p);
		unfreezeLevel(decisionLevel());
		if (c) c->undoLevel(*this);
		undoUntil(decisionLevel()-1);
		return true;
	}
	unfreezeLevel(decisionLevel());
	assert(decision(decisionLevel()) == p);
	assign_.qReset();
	post_.reset();
	return false;
}

bool Solver::resolveConflict() {
	assert(hasConflict());
	++stats.conflicts;
	if (decisionLevel() > rootLevel_) {
		if (decisionLevel() != btLevel_ && strategy_.search != SolverStrategies::no_learning) {
			uint32 uipLevel = analyzeConflict();
			stats.updateJumps(decisionLevel(), uipLevel, btLevel_);
			undoUntil( uipLevel );
			return ClauseCreator::create(*this, cc_, 0, ccInfo_);
		}
		else {
			return backtrack();
		}
	}
	// do not count artificial conflicts that are only used
	// to stop the current search
	stats.conflicts -= hasStopConflict();
	return false;
}

bool Solver::backtrack() {
	Literal lastChoiceInverted;
	do {
		if (decisionLevel() == rootLevel_) return false;
		lastChoiceInverted = ~decision(decisionLevel());
		btLevel_ = decisionLevel() - 1;
		undoUntil(btLevel_);
	} while (!forceImplied() || !force(lastChoiceInverted, 0));
	return true;
}

bool Solver::forceImplied() {
	bool r = !hasConflict();
	const uint32 DL = decisionLevel();
	ImpliedLits::size_type j = 0;
	for (ImpliedLits::size_type i = 0; i < impliedLits_.size(); ++i) {
		if (impliedLits_[i].level <= DL) {
			r = r && force(impliedLits_[i].lit, impliedLits_[i].ante.ante(), impliedLits_[i].ante.data());
			impliedLits_[j++] = impliedLits_[i];
		}
	}
	if (DL == 0) j = 0;
	impliedLits_.erase(impliedLits_.begin()+j, impliedLits_.end());
	return r;
}

void Solver::undoUntil(uint32 level) {
	assert(btLevel_ >= rootLevel_);
	level = std::max( level, btLevel_ );
	if (level >= decisionLevel()) return;
	conflict_.clear();
	strategy_.heuristic->undoUntil( *this, levels_[level].trailPos);
	bool sp = strategy_.saveProgress > 0 && (decisionLevel() - level) > (uint32)strategy_.saveProgress;
	undoLevel(false);
	while (decisionLevel() != level) {
		undoLevel(sp);
	}
}

uint32 Solver::undoUntil(uint32 level, bool popBt) {
	bool checkImplied = false;
	if (popBt && backtrackLevel() > level && !shared_->project(decision(backtrackLevel()).var())) {
		setBacktrackLevel(level);
		checkImplied = true;
	}
	undoUntil(level);
	if (checkImplied) {
		forceImplied();
	}
	return decisionLevel();
}

uint32 Solver::estimateBCP(const Literal& p, int rd) const {
	if (value(p.var()) != value_free) return 0;
	LitVec::size_type first = assign_.assigned();
	LitVec::size_type i     = first;
	Solver& self            = const_cast<Solver&>(*this);
	self.assign_.setValue(p.var(), trueValue(p));
	self.assign_.trail.push_back(p);
	const ShortImplicationsGraph& btig = shared_->shortImplications();
	do {
		Literal x = assign_.trail[i++];  
		if (!btig.propagateBin(self.assign_, x, 0)) {
			break;
		}
	} while (i < assign_.assigned() && rd-- != 0);
	i = assign_.assigned()-first;
	while (self.assign_.assigned() != first) {
		self.assign_.undoLast();
	}
	return (uint32)i;
}
/////////////////////////////////////////////////////////////////////////////////////////
// Solver: Private helper functions
////////////////////////////////////////////////////////////////////////////////////////
Solver::ConstraintDB* Solver::allocUndo(Constraint* c) {
	if (undoHead_ == 0) {
		return new ConstraintDB(1, c);
	}
	assert(undoHead_->size() == 1);
	ConstraintDB* r = undoHead_;
	undoHead_ = (ConstraintDB*)undoHead_->front();
	r->clear();
	r->push_back(c);
	return r;
}
void Solver::undoFree(ConstraintDB* x) {
	// maintain a single-linked list of undo lists
	x->clear();
	x->push_back((Constraint*)undoHead_);
	undoHead_ = x;
}
// removes the current decision level
void Solver::undoLevel(bool sp) {
	assert(decisionLevel() != 0 && levels_.back().trailPos != assign_.trail.size() && "Decision Level must not be empty");
	assign_.undoTrail(levels_.back().trailPos, sp);
	if (levels_.back().undo) {
		const ConstraintDB& undoList = *levels_.back().undo;
		for (ConstraintDB::size_type i = 0, end = undoList.size(); i != end; ++i) {
			undoList[i]->undoLevel(*this);
		}
		undoFree(levels_.back().undo);
	}
	levels_.pop_back();
	if (levConflicts_) levConflicts_->pop_back();
}

// computes the First-UIP clause and stores it in cc_, where cc_[0] is the asserting literal (inverted UIP)
// and cc_[1] is a literal from the asserting level (if > 0)
// RETURN: asserting level of the derived conflict clause
uint32 Solver::analyzeConflict() {
	// must be called here, because we unassign vars during analyzeConflict
	strategy_.heuristic->undoUntil( *this, levels_.back().trailPos );
	uint32 onLevel  = 0;        // number of literals from the current DL in resolvent
	uint32 resSize  = 0;        // size of current resolvent
	Literal p;                  // literal to be resolved out next
	cc_.assign(1, p);           // will later be replaced with asserting literal
	Antecedent lhs, last;       // resolve operands
	const bool doOtfs = strategy_.otfs > 0;
	strategy_.heuristic->updateReason(*this, conflict_, p);
	for (;;) {
		uint32 lhsSize = resSize;
		Antecedent rhs = last;
		for (LitVec::size_type i = 0; i != conflict_.size(); ++i) {
			Literal& q = conflict_[i];
			if (!seen(q.var())) {
				++resSize;
				assert(isTrue(q) && "Invalid literal in reason set!");
				uint32 cl = level(q.var());
				assert(cl > 0 && "Top-Level implication not marked!");
				markSeen(q.var());
				if (cl == decisionLevel()) {
					++onLevel;
				}
				else {
					cc_.push_back(~q);
					markLevel(cl);
				}
			}
		}
		if (resSize != lhsSize) {
			lhs = 0;
		}
		if (conflict_.size() != resSize || lastSimplify_ != units_) {
			// If lastSimplify_ != units_ conflict_ may contain literals assigned
			// on level 0 which are not counted in resSize
			rhs = 0;
		}
		if (doOtfs && (!rhs.isNull() || !lhs.isNull())) {
			// resolvent subsumes rhs and possibly also lhs
			otfs(lhs, rhs, p, onLevel == 1);
		}
		// search for the last assigned literal that needs to be analyzed...
		while (!seen(assign_.last().var())) {
			assign_.undoLast();
		}
		if (--onLevel == 0) {
			conflict_.push_back(~p);
			break;
		}
		p = assign_.last();
		clearSeen(p.var());
		--resSize;              // p was resolved out
		last = reason(p);       // next resolve operand
		reason(p, conflict_);
		strategy_.heuristic->updateReason(*this, conflict_, p);
	}
	cc_[0] = ~assign_.last(); // store the 1-UIP
	clearSeen(cc_[0].var());
	assert(decisionLevel() == level(cc_[0].var()));
	ClauseHead* cRhs = 0;
	if (doOtfs) {
		if (!lhs.isNull()) {
			cRhs = otfsRemove(lhs.constraint()->clause(), &cc_);
		}
		else if (!last.isNull() && last.type() == Antecedent::generic_constraint) {
			cRhs = strategy_.otfs > 1 ? last.constraint()->clause() : 0;
		}
	}
	ccInfo_ = ClauseInfo()
		.setType(Constraint_t::learnt_conflict)
		.setActivity(static_cast<uint32>(1+stats.restarts));
	return finalizeConflictClause(cRhs);
}

void Solver::otfs(Antecedent& lhs, const Antecedent& rhs, Literal p, bool final) {
	ClauseHead* cLhs = 0, *cRhs = 0;
	ClauseHead::BoolPair x;
	if (!lhs.isNull() && (cLhs = lhs.constraint()->clause()) != 0) {
		x = cLhs->strengthen(*this, ~p, !final);
		if (!x.first || x.second) {
			cLhs = !x.first ? 0 : otfsRemove(cLhs, 0);
		}
	}
	lhs = cLhs;
	if (!rhs.isNull() && rhs.type() == Antecedent::generic_constraint && (cRhs = rhs.constraint()->clause()) != 0) {
		x = cRhs->strengthen(*this, p, !final);
		if (!x.first || x.second) {
			cRhs = !x.first ? 0 : otfsRemove(cRhs, 0);
		}
		if (cLhs && cRhs) {
			// lhs and rhs are now equal - only one of them is needed
			if (cLhs->type() == Constraint_t::static_constraint) {
				std::swap(cLhs, cRhs);
			}
			otfsRemove(cLhs, 0);
		}
		lhs = cRhs;
	}
}

ClauseHead* Solver::otfsRemove(ClauseHead* c, const LitVec* newC) {
	ConstraintType t = c->type();
	bool remStatic   = !newC || (newC->size() <= 3 && shared_->unique());
	if ((t == Constraint_t::static_constraint && remStatic) || 
	    (t == Constraint_t::learnt_conflict || t == Constraint_t::learnt_loop)) {
		ConstraintDB& db = (t == Constraint_t::static_constraint ? constraints_ : learnts_);
		ConstraintDB::iterator it;
		if ((it = std::find(db.begin(), db.end(), c)) != db.end()) {
			db.erase(it);
			c->destroy(this, true);
			c = 0;
		}
	}
	return c;
}

// minimizes the conflict clause in cc_ w.r.t selected strategy and
// computes asserting level and lbd of minimized clause.
// PRE:
//  - cc_ is a valid conflict clause and cc_[0] is the UIP-literal
//  - all literals in cc_ except cc_[0] are marked
//  - all decision levels of literals in cc_ are marked
// POST:
//  - literals and decision levels in cc_ are no longer marked
//  - if cc_.size() > 1: cc_[1] is a literal from the asserting level
// RETURN: asserting level of derived conflict clause
uint32 Solver::finalizeConflictClause(ClauseHead* rhs) {
	// 1. remove redundant literals from conflict clause
	LitVec::size_type stop = assign_.trail.size();
	if (strategy_.strengthenRecursive && !ccMin_) {
		ccMin_ = new CCMinRecursive;
	}
	uint32 onAssert = ccMinimize(cc_, assign_.trail, strategy_.cflMinAntes, ccMin_);
	uint32 jl       = cc_.size() > 1 ? level(cc_[1].var()) : 0;
	// clear seen flags of removed literals - keep levels marked
	for (LitVec::size_type x = assign_.trail.size(); x != stop; ) {
		clearSeen(assign_.trail[--x].var());
	}
	rhs             = rhs && cc_.size() <= conflict_.size() ? rhs : 0;
	if (rhs) {
		markSeen(cc_[0].var());
		// check if rhs is subsumed by cc_
		// NOTE: at this point rhs == ~conflict_
		uint32 lhsOpen = (uint32)cc_.size();
		for (LitVec::size_type i = 0; i != conflict_.size() && lhsOpen; ++i) {
			// NOTE: at this point the DB might not be fully simplified,
			//       e.g. because of mt or lookahead, hence we must explicitly
			//       check for literals assigned on DL 0
			lhsOpen -= level(conflict_[i].var()) > 0 && seen(conflict_[i].var());
  	}
		if ( (rhs = lhsOpen == 0 ? otfsRemove(rhs, &cc_) : 0) != 0 ) {
			// rhs is subsumed by cc_ but could not be deleted -
			// try to strengthen.
			ClauseHead::BoolPair x(true, false);
			for (LitVec::size_type i = 0; i != conflict_.size() && x.first; ++i) {
				if (!seen(conflict_[i].var()) || level(conflict_[i].var()) == 0) {
					x = rhs->strengthen(*this, ~conflict_[i], false);
				}
			}
		}
		clearSeen(cc_[0].var());
	}
	if (onAssert == 1 && strategy_.reverseArcs > 0) {
		uint32 maxN = (uint32)strategy_.reverseArcs;
		if      (maxN > 2) maxN = UINT32_MAX;
		else if (maxN > 1) maxN = static_cast<uint32>(cc_.size() / 2);
		markSeen(cc_[0].var());
		Antecedent ante;
		if (ccReverseArc(cc_[1], jl, maxN, ante)) {
			conflict_.clear();
			ante.reason(*this, ~cc_[1], conflict_);
			strategy_.heuristic->updateReason(*this, conflict_, cc_[1]);
			Literal x;
			for (LitVec::size_type i = 0; i != conflict_.size(); ++i) {
				x = conflict_[i];
				assert(isTrue(x));
				if (!seen(x.var())) {
					assert(level(x.var()) < jl);
					markLevel(level(x.var()));
					cc_.push_back(~x);
				}
			}
			clearSeen(cc_[1].var());
			unmarkLevel(level(cc_[1].var()));
			cc_[1] = cc_.back();
			cc_.pop_back();
		}
		clearSeen(cc_[0].var());
	}
	// 2. clear flags and compute lbd
	uint32  lbd         = 1;
	uint32  aLbd        = 1;
	uint32  varLevel    = 0;
	uint32  assertLevel = 0;
	uint32  assertPos   = 1;
	Literal tagLit      = ~sharedContext()->tagLiteral();
	for (LitVec::size_type i = 1; i != cc_.size(); ++i) {
		clearSeen(cc_[i].var());
		if (cc_[i] == tagLit) { ccInfo_.tag(); }
		if ( (varLevel = level(cc_[i].var())) > assertLevel ) {
			assertLevel = varLevel;
			assertPos   = static_cast<uint32>(i);
		}
		if (hasLevel(varLevel)) {
			unmarkLevel(varLevel);
			++lbd;
			aLbd += (varLevel > rootLevel());
		}
	}
	if (assertPos != 1)         { std::swap(cc_[1], cc_[assertPos]); }
	ccInfo_.setLbd(lbd, aLbd);
	// 3. clear level flags of redundant literals
	while (assign_.trail.size() != stop) {
		Var v = assign_.trail.back().var();
		assign_.trail.pop_back();
		unmarkLevel(level(v));
	}
	return assertLevel;
}

// conflict clause minimization
// PRE: 
//  - cc is an asserting clause and cc[0] is the asserting literal
//  - all literals in cc are marked as seen
//  -  if ccMin != 0, all decision levels of literals in cc are marked
// POST:
//  - redundant literals were added to removed
//  - if (cc.size() > 1): cc[1] is a literal from the asserting level
// RETURN
//  - the number of literals from the asserting level
uint32 Solver::ccMinimize(LitVec& cc, LitVec& removed, uint32 anteMask, CCMinRecursive* ccMin) {
	if (ccMin) { ccMin->init(numVars()+1); }
	// skip the asserting literal
	LitVec::size_type j = 1;
	uint32 assertLevel  = 0;
	uint32 assertPos    = 1;
	uint32 onAssert     = 0;
	uint32 varLevel     = 0;
	for (LitVec::size_type i = 1; i != cc.size(); ++i) { 
		if (!ccRemovable(~cc[i], anteMask, ccMin)) {
			if ( (varLevel = level(cc[i].var())) > assertLevel ) {
				assertLevel = varLevel;
				assertPos   = static_cast<uint32>(j);
				onAssert    = 0;
			}
			onAssert += (varLevel == assertLevel);
			cc[j++] = cc[i];
		}
		else { 
			removed.push_back(cc[i]);
		}
	}
	cc.erase(cc.begin()+j, cc.end());
	if (assertPos != 1) {
		std::swap(cc[1], cc[assertPos]);
	}
	if (ccMin) { ccMin->clear(); }
	return onAssert;
}

// returns true if p is redundant in current conflict clause
bool Solver::ccRemovable(Literal p, uint32 m, CCMinRecursive* ccMin) {
	const Antecedent& ante = reason(p);
	if (ante.isNull() || ((ante.type()+1) & m) == 0) {
		return false;
	}
	if (!ccMin) { return ante.minimize(*this, p, 0); }
	// recursive minimization
	LitVec& dfsStack = ccMin->dfsStack;
	assert(dfsStack.empty());
	CCMinRecursive::State dfsState = CCMinRecursive::state_removable;
	p.clearWatch();
	dfsStack.push_back(p);
	for (Literal x;; ) {
		x = dfsStack.back();
		dfsStack.pop_back();
		assert(!seen(x.var()) || x == p);
		if (x.watched()) {
			if (x == p) return dfsState == CCMinRecursive::state_removable;
			ccMin->markVisited(x, dfsState);
		}
		else if (dfsState != CCMinRecursive::state_poison) {
			CCMinRecursive::State temp = ccMin->state(x);
			if (temp == CCMinRecursive::state_open) {
				assert(value(x.var()) != value_free && hasLevel(level(x.var())));
				x.watch();
				dfsStack.push_back(x);
				const Antecedent& ante = reason(x);
				if (ante.isNull() || ((ante.type()+1)&m) == 0 || !ante.minimize(*this, x, ccMin)) {
					dfsState = CCMinRecursive::state_poison;
				}
			}
			else if (temp == CCMinRecursive::state_poison) {
				dfsState = temp;
			}
		}
	}
}

// checks whether there is a valid "reverse arc" for the given literal p that can be used
// to resolve p out of the current conflict clause
// PRE: 
//  - all literals in the current conflict clause are marked
//  - p is a literal of the current conflict clause and level(p) == maxLevel
// RETURN
//  - false if no reverse arc was found. Otherwise true and the antecdent is returned in out
bool Solver::ccReverseArc(const Literal p, uint32 maxLevel, uint32 maxNew, Antecedent& out) {
	assert(seen(p.var()) && isFalse(p) && level(p.var()) == maxLevel);
	const ShortImplicationsGraph& btig = shared_->shortImplications();
	WatchList& wl   = watches_[p.index()];
	if (!btig.reverseArc(*this, p, maxLevel, out)) {
		WatchList::left_iterator it, end; 
		for (it = wl.left_begin(), end = wl.left_end();  it != end;  ++it) {
			if (it->head->isReverseReason(*this, ~p, maxLevel, maxNew)) {
				out = it->head;
				return true;
			}
		}
		return false;
	}
	return true;
}

// (inefficient) default implementation 
bool Constraint::minimize(Solver& s, Literal p, CCMinRecursive* rec) {
	LitVec temp;
	reason(s, p, temp);
	for (LitVec::size_type i = 0; i != temp.size(); ++i) {
		if (!s.ccMinimize(temp[i], rec)) {
			return false;
		}
	}
	return true;
}


// Selects next branching literal. Use user-supplied heuristic if rand() < randProp.
// Otherwise makes a random choice.
// Returns false if assignment is total.
bool Solver::decideNextBranch() {
	DecisionHeuristic* heu = strategy_.heuristic.get();
	if (randHeuristic_ && strategy_.rng.drand() < static_cast<SelectRandom*>(randHeuristic_)->randFreq()) {
		heu = randHeuristic_;
	}
	return heu->select(*this);
}

// Removes up to maxRem% of the learnt nogoods but
// keeps those that are locked or are highly active.
void Solver::reduceLearnts(float maxRem) {
	uint32 oldS = numLearntConstraints();
	uint32 j    = 0;
	if (maxRem < 1.0f) {    
		uint32 remMax = static_cast<uint32>(numLearntConstraints() * std::min(1.0f, std::max(0.15f, maxRem)));
		uint32 remMin = static_cast<uint32>(numLearntConstraints() * 0.15f);
		switch (strategies().reduceAlgo) {
			case SolverStrategies::reduce_linear:   j = reduceLinear(remMin, remMax); break;
			case SolverStrategies::reduce_in_place: j = reduceSortedInPlace(remMin, remMax); break;
			case SolverStrategies::reduce_stable:   
			default: j = reduceSorted(remMin, remMax); break;
		}
	}
	else {
		// remove all nogoods that are not locked
		for (LitVec::size_type i = 0; i != learnts_.size(); ++i) {
			LearntConstraint* c = static_cast<LearntConstraint*>(learnts_[i]);
			if (c->locked(*this)) {
				c->decreaseActivity();
				learnts_[j++] = c;
			}
			else {
				c->destroy(this, true);
			}
		}
	}
	shrinkVecTo(learnts_, j);
	stats.deleted += oldS - numLearntConstraints();
}

// Removes up to maxR of the learnt nogoods.
// Keeps those that are locked or have a high activity.
uint32 Solver::reduceLinear(uint32, uint32 maxR) {
	// compute average activity
	uint64 scoreSum  = 0;
	for (LitVec::size_type i = 0; i != learnts_.size(); ++i) {
		scoreSum += strategies().score(static_cast<LearntConstraint*>(learnts_[i])->activity());
	}
	double avgAct    = (scoreSum / (double) numLearntConstraints());
	// constraints with socre > 1.5 times the average are "active"
	double scoreThresh = avgAct * 1.5;
	double scoreMax    = (double)strategies().score(maxActivity());
	if (scoreThresh > scoreMax) {
		scoreThresh = (scoreMax + (scoreSum / (double) numLearntConstraints())) / 2.0;
	}
	// remove up to maxR constraints but keep "active" and locked once
	const uint32 glue = strategies().reduceGlue;
	uint32 j = 0;
	for (LitVec::size_type i = 0; i != learnts_.size(); ++i) {
		LearntConstraint* c          = static_cast<LearntConstraint*>(learnts_[i]);
		LearntConstraint::Activity a = c->activity();
		if (maxR == 0 || c->locked(*this) || strategies().score(a) > scoreThresh || a.lbd() <= glue) {
			c->decreaseActivity();
			learnts_[j++] = c;
		}
		else {
			--maxR;
			c->destroy(this, true);
		}
	}
	return j;
}

// Removes up to maxR of the most inactive constraints
// while reordering learnts_
uint32 Solver::reduceSortedInPlace(uint32 minR, uint32 maxR) {
	assert(minR <= maxR && maxR <= learnts_.size());
	CmpScore cmp(learnts_, strategies());
	ConstraintDB::iterator hBeg = learnts_.begin();
	ConstraintDB::iterator hEnd = learnts_.begin();
	const uint32 glue           = strategies().reduceGlue;
	// move maxR least active constraints to front of learnts
	for (LitVec::size_type i = 0, end = learnts_.size(); i != end; ++i) {
		LearntConstraint* c = static_cast<LearntConstraint*>(learnts_[i]);
		if (!c->locked(*this)) {
			if (minR || (maxR && c->activity().lbd() > glue)) {
				learnts_[i] = *hEnd;
				*hEnd       = c;
				++hEnd;
				--maxR;
				if      (!minR)       { std::push_heap(hBeg, hEnd, cmp); }
				else if (--minR == 0) { std::make_heap(hBeg, hEnd, cmp); }
			}
			else if (cmp(c, learnts_[0])) {
				std::pop_heap(hBeg, hEnd, cmp);
				learnts_[i] = *(hEnd-1);
				*(hEnd-1)   = c;
				std::push_heap(hBeg, hEnd, cmp);
			}
		}
	}
	// remove all constraints in heap
	for (ConstraintDB::iterator it = hBeg; it != hEnd; ++it) {
		static_cast<LearntConstraint*>(*it)->destroy(this, true);
	}
	// copy remaining constraints down
	uint32 j = 0;
	for (ConstraintDB::iterator it = hEnd, end = learnts_.end(); it != end; ++it) {
		LearntConstraint* c = static_cast<LearntConstraint*>(*it);
		c->decreaseActivity();
		learnts_[j++] = c;
	}
	return j;
}

// Removes up to maxR of the most inactive constraints
// without reordering learnts_
uint32 Solver::reduceSorted(uint32 minR, uint32 maxR) {
	assert(minR <= maxR && maxR <= learnts_.size());
	VarVec temp; temp.reserve(maxR);
	CmpScore cmp(learnts_, strategies());
	const uint32 glue = strategies().reduceGlue;
	// add maxR most inactive constraints to heap
	for (LitVec::size_type i = 0; i != learnts_.size(); ++i) {
		LearntConstraint* c = static_cast<LearntConstraint*>(learnts_[i]);
		if (!c->locked(*this)) {
			if (minR || (maxR && c->activity().lbd() > glue)) {
				temp.push_back(i);
				--maxR;
				if      (!minR)       { std::push_heap(temp.begin(), temp.end(), cmp); }
				else if (--minR == 0) { std::make_heap(temp.begin(), temp.end(), cmp); }
			}
			else if (cmp(i, temp[0])) {
				std::pop_heap(temp.begin(), temp.end(), cmp);
				temp.back() = i;
				std::push_heap(temp.begin(), temp.end(), cmp);
			}
		}
	}
	while (!temp.empty()) {
		learnts_[temp.back()]->destroy(this, true);
		learnts_[temp.back()] = 0;
		temp.pop_back();
	}
	uint32 j = 0;
	for (LitVec::size_type i = 0; i != learnts_.size(); ++i) {
		if (learnts_[i]) {
			static_cast<LearntConstraint*>(learnts_[i])->decreaseActivity();
			learnts_[j++] = learnts_[i];
		}
	}
	return j;
}

uint32 Solver::computeLbd(Literal p, const Literal* first, const Literal* last) {
	if (++lbdTime_ == 0) { lbdStamp_.clear(); lbdTime_ = 1; }
	lbdStamp_.resize(levels_.size()+1, lbdTime_);
	uint32 lbd = 1;
	lbdStamp_[value(p.var()) != value_free ? level(p.var()) : decisionLevel()] = lbdTime_;
	for (; first != last; ++first) {
		if (lbdStamp_[level(first->var())] != lbdTime_) {
			lbdStamp_[level(first->var())] = lbdTime_;
			++lbd;
		}
	}
	return lbd;
}
/////////////////////////////////////////////////////////////////////////////////////////
// The basic DPLL-like search-function
/////////////////////////////////////////////////////////////////////////////////////////
ValueRep Solver::search(RestartLimit& rlimit, LearntLimit& dlimit, double randProp) {
	initRandomHeuristic(randProp);
	uint64 cfl     = 0; // number of conflicts
	rlimit.maxConf = std::max(uint64(1), rlimit.maxConf);
	if (rlimit.countLocal) {
		if (!levConflicts_) levConflicts_ = new VarVec();
		if (decisionLevel() == rootLevel()) { levConflicts_->clear(); }
		levConflicts_->resize(decisionLevel()+1, (uint32)stats.conflicts);
	}
	do {
		if ((hasConflict()&&!resolveConflict()) || !simplify()) { return value_false; }
		do {
			while (!propagate()) {
				++cfl;
				if (!resolveConflict() || (decisionLevel() == 0 && !simplify())) {
					return value_false;
				}
				if ( (cfl == rlimit.maxConf && !rlimit.countLocal) ||
					(rlimit.countLocal && (stats.conflicts - (*levConflicts_)[decisionLevel()]) >= rlimit.maxConf)) {
					rlimit.maxConf = 0;
					dlimit.maxConf-= std::min(cfl, dlimit.maxConf);
					return value_free;  
				}
			}
			if (numLearntConstraints()>dlimit.maxLearnt || cfl >= dlimit.maxConf) { 
				dlimit.maxConf -= std::min(cfl, dlimit.maxConf);
				rlimit.maxConf -= !rlimit.countLocal ? cfl : 0;
				return value_free; // dlimit reached
			}
		} while (decideNextBranch());
		// found a model candidate
		assert(numFreeVars() == 0 || hasConflict());
	} while (hasConflict() || !post_.isModel(*this));
	stats.addModel(decisionLevel());
	if (satPrepro()) {
		satPrepro()->extendModel(assign_, unconstr_);
	}
	dlimit.maxConf -= cfl;
	rlimit.maxConf -= !rlimit.countLocal ? cfl : 0;
	return value_true;
}

bool Solver::nextSymModel(bool expand) {
	assert(numFreeVars() == 0);
	if (expand) {
		if (satPrepro() != 0 && !unconstr_.empty()) {
			stats.addModel(decisionLevel());
			satPrepro()->extendModel(assign_, unconstr_);
			return true;
		}
		else if (post_.nextSymModel(*this, true)) {
			stats.addModel(decisionLevel());
			return true;
		}
	}
	else {
		post_.nextSymModel(*this, false);
		unconstr_.clear();
	}
	return false;
}

}
