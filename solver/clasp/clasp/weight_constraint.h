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
#ifndef CLASP_SMODELS_CONSTRAINTS_H_INCLUDED
#define CLASP_SMODELS_CONSTRAINTS_H_INCLUDED

#ifdef _MSC_VER
#pragma warning (disable : 4200) // nonstandard extension used : zero-sized array
#pragma once
#endif

#include <clasp/constraint.h>

namespace Clasp { namespace mt {

struct SharedWeightConstraintLits;

} // namespace mt

class WeightConstraint;	

//! Class storing the literals of a weight constraint
/*!
 * \see WeightConstraint
 */
class WeightConstraintLits {
public:
	/*!
	 * Logically, we distinguish two constraints: 
	 * FFB_BTB for handling forward false body and backward true body and
	 * FTB_BFB for handling forward true body and backward false body
	 * Physically, we store the literals in one array: ~W=1, l1=w1,...,ln=wn
	 */
	enum ActiveConstraint {
		FFB_BTB   = 0, /**< (SumW-bound)+1 [~W=1, l1=w1,...,ln=wn]; */
		FTB_BFB   = 1, /**< bound [W=1, ~l1=w1,...,~ln=wn] */
	};
	
	/*!
	 * Returns the i'th literal of constraint c, i.e.
	 *  li, iff c == FFB_BTB
	 * ~li, iff c == FTB_BFB
	 */
	Literal	 lit(uint32 i, ActiveConstraint c) const {
		return Literal::fromIndex( lits_[(i<<wc_)].index() ^ c );
	}
	//! Returns the i'th variable of this constraint
	Var var(uint32 i) const { return lits_[(i<<wc_)].var(); }

	//! Returns true for weight- and false for cardinality constraints
	bool hasWeights() const { return wc_ != 0; }

	//! Returns the weight of the i'th literal
	weight_t weight(uint32 i)	const {
		return wc_ == 0 ? weight_t(1) : (weight_t)lits_[(i<<1)+1].asUint();
	}

	//! number of literals (counting W)
	uint32 size() const { return size_;  }

	class LiteralIterator {
	public:
		explicit LiteralIterator(const Literal* p = 0, bool weights = false) : pos(p), inc(1+weights) {}
		LiteralIterator& operator++()   { pos += inc; return *this; }
		LiteralIterator operator++(int) { LiteralIterator t(*this); ++*this; return t; } 
		LiteralIterator& operator--()   { pos -= inc; return *this; }
		LiteralIterator operator--(int) { LiteralIterator t(*this); --*this; return t; } 
		Literal         operator*()const{ return *pos; }
		weight_t weight() const         { return inc == 1 ? 1 : static_cast<weight_t>((pos+1)->asUint()); }
		bool operator==(const LiteralIterator& rhs) const { return pos == rhs.pos; }
		bool operator!=(const LiteralIterator& rhs) const { return pos != rhs.pos; }
	private:
		LiteralIterator& operator=(const LiteralIterator&);
		const Literal* pos;
		const uint32   inc;
	};
	typedef LiteralIterator Iterator;
	Iterator begin() const { return LiteralIterator(lits_, wc_); }
	Iterator end()   const { return LiteralIterator(lits_+(size_<<wc_), wc_); }

	//! Returns true if this object can be shared between solvers
	bool sharable()   const { return shared_ != 0; }

	//! Returns true if this object is currently not shared
	bool unique() const;

	//! Shares this object
	/*!
	 * \pre sharable() == true
	 */
	WeightConstraintLits* share();
	
	//! Releases a reference to this object
	void release();
private:
	WeightConstraintLits() {}
	WeightConstraintLits(uint32 s, bool shared, bool w) : size_(s), shared_(shared), wc_(w) { }
	WeightConstraintLits(const WeightConstraintLits&);
	WeightConstraintLits& operator=(const WeightConstraintLits&);
	friend class  WeightConstraint;
	friend struct mt::SharedWeightConstraintLits;
	uint32   size_   : 30; // number of lits in constraint (counting the literal associated with the constraint)
	uint32   shared_ :  1; // if 1  this object is actually part of a SharedWeightConstraintLits object
	uint32   wc_     :  1; // 1 if this is a weight constraint, otherwise 0 (in that case no weights are stored)
	Literal  lits_[0];     // Literals of constraint: ~B [Bw], l1 [w1], ..., ln-1 [Wn-1]
};

//! Class implementing smodels-like cardinality- and weight constraints.
/*!
 * \ingroup constraint
 * This class represents a constraint of type W == w1*x1 ... wn*xn >= B,
 * where W and each xi are literals and B and each wi are strictly positive integers.
 * The class is used to represent smodels-like weight constraint, i.e.
 * the body of a basic weight rule. In this case W is the literal associated with the body.
 * A cardinality constraint is handled like a weight constraint where all weights are equal to 1.
 *
 * The class implements the following four inference rules:
 * Let L be the set of literals of the constraint,
 * let sumTrue be the sum of the weights of all literals l in L that are currently true,
 * let sumReach be the sum of the weights of all literals l in L that are currently not false,
 * let U = {l in L | value(l.var()) == value_free}
 * - FTB: If sumTrue >= bound: set W to true.
 * - BFB: If W is false: set false all literals l in U for which sumTrue + weight(l) >= bound.
 * - FFB: If sumReach < bound: set W to false.
 * - BTB: If W is true: set true all literals l in U for which sumReach - weight(l) < bound.
 */
class WeightConstraint : public Constraint {
public:
	typedef WeightConstraintLits WL;

	//! Creates a new weight constraint from the given weight literals
	/*!
	 * If the right hand side of the weight constraint is initially true/false (FTB/FFB),
	 * W is assigned appropriately but no constraint is created. Otherwise
	 * the new weight constraint is added to the context.
	 * \param ctx context in which the new constraint is to be used.
	 * \param W the literal that is associated with the constraint
	 * \param lits the literals of the weight constraint
	 * \param bound the lower bound of the weight constraint.
	 * \return false if the constraint is initially conflicting w.r.t the current assignment.
	 * \note Cardinality constraint are represented as weight constraints with all weights equal
	 * to 1.
	 */
	static bool newWeightConstraint(SharedContext& ctx, Literal W, WeightLitVec& lits, weight_t bound, Constraint** out = 0);
	
	// constraint interface
	Constraint* cloneAttach(Solver&);
	bool simplify(Solver& s, bool = false);
	void destroy(Solver*, bool);
	PropResult propagate(Solver& s, Literal p, uint32& data);
	void reason(Solver&, Literal p, LitVec& lits);
	bool minimize(Solver& s, Literal p, CCMinRecursive* r);
	void undoLevel(Solver& s);
	uint32 estimateComplexity(const Solver& s) const;
private:
	WeightConstraint(SharedContext& ctx, Literal W, const WeightLitVec& lits, uint32 bound, uint32 sumW);
	WeightConstraint(Solver& s, const WeightConstraint& other);
	~WeightConstraint();
	static weight_t canonicalize(Solver& s, WeightLitVec& lits, weight_t& bound);
	typedef WL::ActiveConstraint ActiveConstraint;
	
	static const uint32 NOT_ACTIVE = 3u;
	
	// Represents a literal on the undo stack
	// idx()        returns the index of the literal
	// constraint() returns the constraint that added the literal to the undo stack
	// Note: Only 31-bits are used for undo info.
	// The remaining bit is used as a flag for marking processed literals.
	struct UndoInfo {
		explicit UndoInfo(uint32 d = 0) : data(d) {}
		uint32           idx()        const { return data >> 2; }
		ActiveConstraint constraint() const { return static_cast<ActiveConstraint>((data&2) != 0); }
		uint32 data; 
	};
	// is literal idx contained as reason lit in the undo stack?
	bool litSeen(uint32 idx) const { return (undo_[idx].data & 1) != 0; }
	// mark/unmark literal idx
	void toggleLitSeen(uint32 idx) { undo_[idx].data ^= 1; }
	
	void addWatch(Solver& s, uint32 idx, ActiveConstraint c);
	// Updates bound_[c] and adds an undo watch to the solver if necessary.
	// Then adds the literal at position idx to the reason set (and the undo stack).
	void updateConstraint(Solver& s, uint32 idx, ActiveConstraint c);
	// Returns the starting index of the undo stack.
	uint32   undoStart()       const { return lits_->hasWeights(); }
	UndoInfo undoTop()         const { assert(up_ != undoStart()); return undo_[up_-1]; }
	// Returns the decision level of the last assigned literal
	// or 0 if no literal was assigned yet.
	inline uint32	highestUndoLevel(Solver&) const;
	// Returns the index of next literal to look at during backward propagation.
	uint32   getBpIndex() const  { return !lits_->hasWeights() ? 1 : undo_[0].data>>1; }
	void     setBpIndex(uint32 n){ if (lits_->hasWeights()) undo_[0].data = (n<<1)+(undo_[0].data&1); }
	
	WL*      lits_;        // literals of constraint
	uint32   up_     : 30; // undo position; [undoStart(), up_] is the undo stack
	uint32   active_ :  2; // which of the two sub-constraints is currently unit?
	weight_t bound_[2];    // FFB_BTB: (sumW-bound)+1 / FTB_BFB: bound
	UndoInfo undo_[0];     // undo stack + seen flag for each literal
};
}

#endif
