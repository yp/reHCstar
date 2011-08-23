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

#include <clasp/weight_constraint.h>
#include <clasp/clause.h>
#include <clasp/solver.h>
#include <clasp/util/misc_types.h>
#include <algorithm>

namespace Clasp { namespace mt {

struct SharedWeightConstraintLits {
	SharedWeightConstraintLits(uint32 size, bool hasWeight) : shared_lits(size, true, hasWeight) {
		shared_counter = 1; 
	}
	std::atomic<int>     shared_counter;
	WeightConstraintLits shared_lits;
};

SharedWeightConstraintLits* toSharedLits(WeightConstraintLits* x) {
	assert(x->sharable() && "WeightConstraintLits are not shareable!");
#define OFFSETOF(S, x)  (size_t)&reinterpret_cast<const volatile unsigned char&>((((S*)0)->*&S::x))
	return reinterpret_cast<SharedWeightConstraintLits*>(
		reinterpret_cast<unsigned char*>(x) - OFFSETOF(SharedWeightConstraintLits, shared_lits)
	);
#undef OFFSETOF 
}

} // Clasp::mt

WeightConstraintLits* WeightConstraintLits::share() {
	if (sharable()) {
		++mt::toSharedLits(this)->shared_counter;
		return this;
	}
	else {
		std::size_t bytes = sizeof(WeightConstraintLits) + (size() << uint32(hasWeights()))*sizeof(Literal);
		WeightConstraintLits* x = new (::operator new(bytes)) WeightConstraintLits();
		std::memcpy(x, this, bytes);
		return x;
	}
}

void WeightConstraintLits::release() {
	if (!sharable()) {
		::operator delete (this);
	}
	else if (--mt::toSharedLits(this)->shared_counter == 0) {
		void* m = mt::toSharedLits(this);
		::operator delete (m);
	}
}

bool WeightConstraintLits::unique() const {
	return !sharable() || mt::toSharedLits(const_cast<WeightConstraintLits*>(this))->shared_counter == 1;
}
/////////////////////////////////////////////////////////////////////////////////////////
// WeightConstraint
/////////////////////////////////////////////////////////////////////////////////////////
bool WeightConstraint::newWeightConstraint(SharedContext& ctx, Literal con, WeightLitVec& lits, weight_t bound, Constraint** out) {
	assert(ctx.master()->decisionLevel() == 0 && ctx.unique());
	if (out) *out = 0;
	weight_t sumW = canonicalize(*ctx.master(), lits, bound);
	if      (bound <= 0)  { return ctx.addUnary(con); } // trivially SAT
	else if (bound > sumW){ return ctx.addUnary(~con);} // trivially UNSAT
	Solver& s = *ctx.master();
	if (s.value(con.var()) != value_free) {
		// backward propagate
		weight_t b  = s.isTrue(con) ? (sumW-bound)+1 : bound;
		bool bpTrue = s.isTrue(con);
		uint32 i;
		for (i = 0; i != lits.size() && (b - lits[i].second) <= 0; ++i) {
			sumW -= lits[i].second;
			if (bpTrue) bound -= lits[i].second;
			if (!s.force(bpTrue ? lits[i].first : ~lits[i].first, 0)) {
				return false;
			}
		}
		if (i == lits.size()) return true;
		lits.erase(lits.begin(), lits.begin()+i);
		if (lits.front().second == lits.back().second && lits.front().second != 1) {
			// disguised cardinality constraint
			weight_t w = lits.front().second;
			bound      = (bound+(w-1))/w;
			sumW       = (sumW+(w-1))/w;
		}
		if (bound == 1) { // clause
			assert(s.isTrue(con));
			ClauseCreator clause(&s); clause.reserve(lits.size());
			clause.start();
			for (i = 0; i != lits.size(); ++i) { clause.add(lits[i].first); }
			return clause.end();
		}
	}
	bool   hasW = lits.front().second != lits.back().second;
	uint32 size = (uint32)lits.size()+1;
	void* m     = ::operator new(sizeof(WeightConstraint) + (size+hasW)*sizeof(UndoInfo));
	Constraint* c = new (m) WeightConstraint(ctx, con, lits, bound, sumW);
	if (out) *out = c;
	ctx.add(c);
	return !s.hasConflict();
}
	
// removes assigned and merges duplicate/complementary literals
// return: achievable weight
// post  : lits is sorted by decreasing weights
weight_t WeightConstraint::canonicalize(Solver& s, WeightLitVec& lits, weight_t& bound) {
	// Step 1: remove assigned/superfluous literals and merge duplicate/complementary ones
	LitVec::size_type j = 0, other;
	for (LitVec::size_type i = 0; i != lits.size(); ++i) {
		Literal x = lits[i].first;
		if (s.value(x.var()) == value_free && lits[i].second > 0) {
			if (!s.seen(x.var())) { // first time we see x, keep and mark x
				if (i != j) lits[j] = lits[i];
				s.markSeen(x);
				++j;
			}
			else if (!s.seen(~x)) { // multi-occurrences of x, merge
				for (other = 0; other != j && lits[other].first != x; ++other) { ; }
				lits[other].second += lits[i].second;
			}
			else {                  // complementary literals ~x x
				for (other = 0; other != j && lits[other].first != ~x; ++other) { ; }
				bound              -= lits[i].second; // decrease by min(w(~x), w(x)) ; assume w(~x) > w(x)
				lits[other].second -= lits[i].second; // keep ~x, 
				if (lits[other].second < 0) {         // actually, w(x) > w(~x), 
					lits[other].first  = x;             // replace ~x with x
					lits[other].second = -lits[other].second;
					s.clearSeen(x.var());
					s.markSeen(x);
					bound += lits[other].second;        // and correct the bound
				}
				else if (lits[other].second == 0) {   // w(~x) == w(x) - drop both lits
					s.clearSeen(x.var());
					std::memmove(&lits[0]+other, &lits[0]+other+1, (j-other-1)*sizeof(lits[other]));
					--j;
				}
			}
		}
		else if (s.isTrue(x)) { bound -= lits[i].second; }
	}
	lits.erase(lits.begin()+j, lits.end());
	// Step 2: compute min,max, achievable weight and clear flags set in step 1
	weight_t sumW = 0;
	uint32 minW   = std::numeric_limits<uint32>::max(), maxW = 0;
	for (LitVec::size_type i = 0; i != lits.size(); ++i) {
		assert(lits[i].second > 0);
		s.clearSeen(lits[i].first.var());
		if (uint32(lits[i].second) > maxW) { maxW = lits[i].second; }
		if (uint32(lits[i].second) < minW) { minW = lits[i].second; }
		sumW += lits[i].second;
	}
	if (bound < 0) {
		bound = 0;
	}
	// Step 3: sort by decreasing weight
	if (maxW != minW) {
		std::stable_sort(lits.begin(), lits.end(), compose22(
			std::greater<weight_t>(),
			select2nd<WeightLiteral>(),
			select2nd<WeightLiteral>()));
	}
	else if (minW != 1) {
		// disguised cardinality constraint
		bound = (bound+(minW-1))/minW;
		sumW  = (sumW+(minW-1))/minW;
		for (LitVec::size_type i = 0; i != lits.size(); ++i) { lits[i].second = 1; }
	}
	return sumW;
}

WeightConstraint::WeightConstraint(SharedContext& ctx, Literal W, const WeightLitVec& lits, uint32 bound, uint32 sumWeights) {
	const bool hasWeights = lits.front().second != lits.back().second;
	if (!ctx.shareConstraints()) {
		void* m = ::operator new(sizeof(WeightConstraintLits) + ((1+lits.size()) << uint32(hasWeights))*sizeof(Literal));
		lits_   = new (m) WeightConstraintLits(uint32(1+lits.size()), false, hasWeights);
	}
	else {
		void* m = ::operator new(sizeof(mt::SharedWeightConstraintLits) + ((1+lits.size()) << uint32(hasWeights))*sizeof(Literal));
		mt::SharedWeightConstraintLits* x = new (m) mt::SharedWeightConstraintLits(uint32(1+lits.size()), hasWeights);
		lits_ = &x->shared_lits;
	}
	Literal* heu  = new Literal[1+lits.size()];
	heu[0]        = W;
	Literal* p    = lits_->lits_;
	Solver&  s    = *ctx.master();
	bool needData       = lits_->hasWeights();
	bound_[WL::FFB_BTB]	= (sumWeights-bound)+1;// ffb-btb
	bound_[WL::FTB_BFB]	= bound;               // ftb-bfb
	up_                 = undoStart();         // undo stack is initially empty
	undo_[0].data       = 0;
	undo_[up_].data     = 0;
	setBpIndex(1);                             // where to start back propagation
	*p++                = ~W;                  // store constraint literal
	if (hasWeights) *p++= Literal::fromRep(1); // and weight if necessary
	ctx.setFrozen(W.var(), true);              // exempt from variable elimination
	if (needData) ctx.requestData(W.var());
	active_ = s.value(W.var())==value_free     // if unassigned
		? NOT_ACTIVE                             // both constraints are relevant
		: WL::FFB_BTB+s.isFalse(W);              // otherwise only one is relevant
	for (uint32 i = 0, end = (uint32)lits.size(); i != end; ++i) {
		*p++ = heu[i+1] = lits[i].first;         // store constraint literal
		if (hasWeights) {                        // followed by weight 
			*p++= Literal::fromRep(lits[i].second);// if necessary
		}
		addWatch(s, i+1, WL::FTB_BFB);           // watches  lits[idx]
		addWatch(s, i+1, WL::FFB_BTB);           // watches ~lits[idx]
		ctx.setFrozen(lits[i].first.var(), true);// exempt from variable elimination
		if (needData) ctx.requestData(lits[i].first.var());
		undo_[up_+i+1].data = 0;
	}	
	// Initialize heuristic with literals (no weights) in constraint.
	s.strategies().heuristic->newConstraint(s, heu, 1+lits.size(), Constraint_t::static_constraint);
	delete [] heu;
	if (active_ == NOT_ACTIVE) {
		addWatch(s, 0, WL::FTB_BFB);             // watch con in both phases
		addWatch(s, 0, WL::FFB_BTB);             // in order to allow for backpropagation
	}
	else {
		uint32 d = active_;                      // propagate con
		WeightConstraint::propagate(s, ~lits_->lit(0, (ActiveConstraint)active_), d);
	}
}

WeightConstraint::WeightConstraint(Solver& s, const WeightConstraint& other) {
	lits_           = other.lits_->share();
	WL::Iterator it = lits_->begin();
	Literal* heu    = new Literal[lits_->size()];
	heu[0]          = ~*it;
	bound_[0]	      = other.bound_[0];
	bound_[1]	      = other.bound_[1];
	active_         = other.active_;
	if (active_ == NOT_ACTIVE && s.value(heu[0].var()) == value_free) {
		addWatch(s, 0, WL::FTB_BFB);  // watch con in both phases
		addWatch(s, 0, WL::FFB_BTB);  // in order to allow for backpropagation
	}
	for (uint32 i = 1; i < lits_->size(); ++i) {
		heu[i]      = *++it;
		if (s.value(heu[i].var()) == value_free) {
			addWatch(s, i, WL::FTB_BFB);  // watches  lits[i]
			addWatch(s, i, WL::FFB_BTB);  // watches ~lits[i]
		}
	}
	// Initialize heuristic with literals (no weights) in constraint.
	s.strategies().heuristic->newConstraint(s, heu, lits_->size(), Constraint_t::static_constraint);
	delete [] heu;
	std::memcpy(undo_, other.undo_, sizeof(UndoInfo)*(lits_->size()+lits_->hasWeights()));
	up_ = other.up_;
}

WeightConstraint::~WeightConstraint() {
}

Constraint* WeightConstraint::cloneAttach(Solver& other) { 
	void* m = ::operator new(sizeof(WeightConstraint) + (lits_->size()+lits_->hasWeights())*sizeof(UndoInfo));
	return new (m) WeightConstraint(other, *this);

}

void WeightConstraint::addWatch(Solver& s, uint32 idx, ActiveConstraint c) {
	// Add watch only if c is relevant.
	if (uint32(c^1) != active_) {
		// Use LSB to store the constraint that watches the literal.
		s.addWatch(~lits_->lit(idx, c), this, (idx<<1)+c);
	}
}

void WeightConstraint::destroy(Solver*, bool) {
	lits_->release();
	void* mem = static_cast<Constraint*>(this);
	this->~WeightConstraint();
	::operator delete(mem);
}

// Returns the numerical highest decision level watched by this constraint.
uint32 WeightConstraint::highestUndoLevel(Solver& s) const {
	return up_ != undoStart() 
		? s.level(lits_->var(undoTop().idx()))
		: 0;
}

// Updates the bound of sub-constraint c and adds the literal at index idx to the 
// undo stack. If the current decision level is not watched, an undo watch is added
// so that the bound can be adjusted once the solver backtracks.
void WeightConstraint::updateConstraint(Solver& s, uint32 idx, ActiveConstraint c) {
	bound_[c] -= lits_->weight(idx);
	if (highestUndoLevel(s) != s.decisionLevel()) {
		s.addUndoWatch(s.decisionLevel(), this);
	}
	undo_[up_].data = (idx<<2) + (c<<1) + (undo_[up_].data & 1);
	++up_;
	assert(!litSeen(idx));
	toggleLitSeen(idx);
}

// Since clasp uses an eager assignment strategy where literals are assigned as soon
// as they are added to the propagation queue, we distinguish processed from unprocessed literals.
// Processed literals are those for which propagate was already called and the corresponding bound 
// was updated; they are flagged in updateConstraint(). 
// Unprocessed literals are either free or were not yet propagated. During propagation
// we treat all unprocessed literals as free. This way, conflicts are detected early.
// Consider: x :- 3 [a=3, b=2, c=1,d=1] and PropQ: b, ~Body, c. 
// Initially b, ~Body, c are unprocessed and the bound is 3.
// Step 1: propagate(b)    : b is marked as processed and bound is reduced to 1.
//   Now, although we already know that the body is false, we do not backpropagate yet
//   because the body is unprocessed. Deferring backpropagation until the body is processed
//   makes reason computation easier.
// Step 2: propagate(~Body): ~body is marked as processed and bound is reduced to 0.
//   Since the body is now part of our reason set, we can start backpropagation.
//   First we assign the unprocessed and free literal ~a. Literal ~b is skipped, because
//   its complementary literal was already successfully processed. Finally, we force 
//   the unprocessed but false literal ~c to true. This will generate a conflict and 
//   propagation is stopped. Without the distinction between processed and unprocessed
//   lits we would have to skip ~c. We would then have to manually trigger the conflict
//   {b, ~Body, c} in step 3, when propagate(c) sets the bound to -1.
Constraint::PropResult WeightConstraint::propagate(Solver& s, Literal, uint32& d) {
	// determine the affected constraint and its body literal
	ActiveConstraint c  = (ActiveConstraint)(d&1);
	Literal body        = lits_->lit(0, c);
	if ( uint32(c^1) == active_ || s.isTrue(body) ) {
		// the other constraint is active or this constraint is already satisfied; 
		// nothing to do
		return PropResult(true, true);		
	}
	// the constraint is not yet satisfied; update it and
	// check if we can now propagate any literals.
	updateConstraint(s, d >> 1, c);
	if (bound_[c] <= 0 || (lits_->hasWeights() && litSeen(0))) {
		uint32 reasonData = !lits_->hasWeights() ? UINT32_MAX : up_;
		if (!litSeen(0)) {
			// forward propagate constraint to true
			active_ = c;
			return PropResult(s.force(body, this, reasonData), true);
		}
		else {
			// backward propagate false constraint
			uint32 n = getBpIndex();
			for (const uint32 end = lits_->size(); n != end && (bound_[c] - lits_->weight(n)) < 0; ++n) {
				if (!litSeen(n)) {
					active_   = c;
					Literal x = lits_->lit(n, c);
					if (!s.force(x, this, reasonData)) {
						return PropResult(false, true);
					}
				}
			}
			assert(n == 1 || n == lits_->size() || lits_->hasWeights());
			setBpIndex(n);
		}
	}
	return PropResult(true, true);
}

// Builds the reason for p from the undo stack of this constraint.
// The reason will only contain literals that were processed by the 
// active sub-constraint.
void WeightConstraint::reason(Solver& s, Literal p, LitVec& r) {
	assert(active_ != NOT_ACTIVE);
	Literal x;
	uint32 stop = !lits_->hasWeights() ? up_ : s.reasonData(p);
	assert(stop <= up_);
	for (uint32 i = undoStart(); i != stop; ++i) {
		UndoInfo u = undo_[i];
		// Consider only lits that are relevant to the active constraint
		if (u.constraint() == (ActiveConstraint)active_) {
			x = lits_->lit(u.idx(), u.constraint());
			r.push_back( ~x );
		}
	}
}

bool WeightConstraint::minimize(Solver& s, Literal p, CCMinRecursive* rec) {
	assert(active_ != NOT_ACTIVE);
	Literal x;
	uint32 stop = !lits_->hasWeights() ? up_ : s.reasonData(p);
	assert(stop <= up_);
	for (uint32 i = undoStart(); i != stop; ++i) {
		UndoInfo u = undo_[i];
		// Consider only lits that are relevant to the active constraint
		if (u.constraint() == (ActiveConstraint)active_) {
			x = lits_->lit(u.idx(), u.constraint());
			if (!s.ccMinimize(~x, rec)) {
				return false;
			}
		}
	}
	return true;
}

// undoes processed assignments 
void WeightConstraint::undoLevel(Solver& s) {
	setBpIndex(1);
	for (UndoInfo u; up_ != undoStart() && s.value(lits_->var((u=undoTop()).idx())) == value_free;) {
		assert(litSeen(u.idx()));
		toggleLitSeen(u.idx());
		bound_[u.constraint()] += lits_->weight(u.idx());
		--up_;
	}
	if (!litSeen(0)) active_ = NOT_ACTIVE;
}

bool WeightConstraint::simplify(Solver& s, bool) {
	typedef WL::LiteralIterator Iter;
	if (bound_[0] <= 0 || bound_[1] <= 0) {
		for (Iter it = lits_->begin(), end = lits_->end(); it != end; ++it) {
			s.removeWatch( (*it), this );
			s.removeWatch(~(*it), this );
		}
		return true;
	}
	else if (s.value(lits_->var(0)) != value_free && (active_ == NOT_ACTIVE || lits_->hasWeights())) {
		if (active_ == NOT_ACTIVE) {
			Literal W = ~*lits_->begin();
			active_   = WL::FFB_BTB+s.isFalse(W);	
		}
		for (Iter it = lits_->begin(), end = lits_->end(); it != end; ++it) {
			s.removeWatch(active_ == WL::FTB_BFB ? ~(*it) : *it, this);
		}
	}
	if (lits_->unique() && lits_->size() > 4 && (up_ - undoStart()) > lits_->size()/2) {
		Literal*     lits = lits_->lits_;
		const uint32 inc  = 1 + lits_->hasWeights();
		uint32       end  = lits_->size_*inc;
		uint32 i, j, idx  = 1;
		// find first assigned literal - must be there otherwise undo stack would be empty
		for (i = inc; s.value(lits[i].var()) == value_free; i += inc) {
			assert(!litSeen(idx));
			++idx;
		}
		// move unassigned literals down
		// update watches because indices have changed
		for (j = i, i += inc; i != end; i += inc) {
			if (s.value(lits[i].var()) == value_free) {
				lits[j] = lits[i];
				if (lits_->hasWeights()) { lits[j+1] = lits[i+1]; }
				undo_[idx].data = 0;
				assert(!litSeen(idx));
				if (Clasp::GenericWatch* w = s.getWatch(lits[i], this)) {
					w->data = (idx<<1) + 1;
				}
				if (Clasp::GenericWatch* w = s.getWatch(~lits[i], this)) {
					w->data = (idx<<1) + 0;
				}
				j += inc;
				++idx;
			}
			else {
				s.removeWatch(lits[i], this);
				s.removeWatch(~lits[i], this);
			}
		}
		// clear undo stack & update to new size
		up_ = undoStart();
		setBpIndex(1);
		lits_->size_ = idx;
	}
	return false;
}

uint32 WeightConstraint::estimateComplexity(const Solver& s) const {
	weight_t w1 = bound_[0];
	weight_t w2 = bound_[1];
	uint32 r    = 2;
	for (uint32 i = 1; i != lits_->size(); ++i) {
		if (s.value(lits_->var(i)) == value_free) {
			++r;
			if ( (w1 -= lits_->weight(i)) <= 0 ) break;
			if ( (w2 -= lits_->weight(i)) <= 0 ) break;
		}
	}
	return r;
}
}
