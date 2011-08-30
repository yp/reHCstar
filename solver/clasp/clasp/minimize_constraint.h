// 
// Copyright (c) 2010, Benjamin Kaufmann
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
#ifndef CLASP_MINIMIZE_CONSTRAINT_H_INCLUDED
#define CLASP_MINIMIZE_CONSTRAINT_H_INCLUDED

#ifdef _MSC_VER
#pragma warning (disable : 4200) // nonstandard extension used : zero-sized array
#pragma once
#endif

#include <clasp/constraint.h>
#include <clasp/util/atomic.h>
#include <cassert>

namespace Clasp {
class MinimizeConstraint;
//! Supported minimization modes
/*! 
 * Defines the possible minimization modes used during solving.
 * If optimize is used, an optimal solution is one that is
 * strictly smaller than all previous solutions. Otherwise
 * an optimal solution is one that is no greater than an initially 
 * given optimum.
 */
struct MinimizeMode_t {
	enum Mode {
		optimize,  /**< optimize using < */
		enumerate  /**< enumerate w.r.t given optimum */
	};
};
typedef MinimizeMode_t::Mode MinimizeMode;

//! A type holding data (possibly) shared between a set of minimize constraints
/*!
 * \ingroup shared
 */
struct SharedMinimizeData {
	//! A type to represent a weight at a certain level.
	/*!
	 * Objects of this type are used to create sparse vectors of weights. E.g.
	 * a weight vector (w1@L1, w2@L3, w3@L5) is represented as [<L1,1,w1><L3,1,w2><L5,0,w3>], 
	 * where each <level, next, weight>-tuple is an object of type LevelWeight.
	 */
	struct LevelWeight {
		LevelWeight(uint32 l, weight_t w) : level(l), next(0), weight(w) {}
		uint32   level : 31; /**< the level of this weight */
		uint32   next  :  1; /**< does this weight belong to a sparse vector of weights? */
		weight_t weight;     /**< the weight at this level */
	};
	//! A type for holding sparse vectors of level weights of a multi-level constraint
	typedef PodVector<LevelWeight>::type WeightVec;
	//! A dense vector of weights
	typedef PodVector<wsum_t>::type      SumVec;
	//! A type for holding an optimum and the information for how to proceed
	struct Optimum {
		typedef std::pair<uint32, wsum_t> Step;
		typedef std::atomic<uint32>       Atomic;
		Atomic seq;     /**< sequence number associated with this optimum */
		SumVec opt;     /**< current optimum / one value per level      */
		Step   step;    /**< (opt - step) gives the next bound to check */
	};
	explicit SharedMinimizeData(const SumVec& initSum, const SumVec& initOpt, MinimizeMode m = MinimizeMode_t::optimize);
	~SharedMinimizeData();
	
	void destroy() const;
	
	//! sets the enumeration mode used during solving
	/*!
	 * \note if m is MinimizeMode::enumerate, the caller must explicitly
	 * set an optimum. Otherwise, *all* solutions are considered optimal.
	 */
	void setMode(MinimizeMode m, int hierarchical);
	
	//! number of minimize statements contained in this constraint
	uint32         numRules()      const { return static_cast<uint32>(sum_.size()); }
	uint32         maxLevel()      const { return numRules()-1; }
	//! returns the active minimization mode
	MinimizeMode   mode()          const { return mode_; }
	//! hierarchical optimization enabled?
	bool           hierarchical()  const { return hierarch_ > 0; }
	//! returns the initial sum
	const SumVec&  sum()           const { return sum_; }
	//! returns the current optimum
	const Optimum* optimum()       const { return opt_; }
	//! returns the level of the literal with the given idx
	uint32         getLevel(uint32 litIdx) const {
		return numRules() == 1 ? 0 : weights[lits[litIdx].second].level;
	}
	//! set optOut to (opt - (applyStep * opt->step))
	bool           convert(const Optimum* opt, bool applyStep, wsum_t* optOut) const;
	
	/*!
	 * \name interface for optimization
	 * The following functions shall not be called concurrently.
	 */
	//@{
	
	//! makes opt the new optimum
	/*!
	 * \pre opt is a pointer to an array of size numRules()
	 */
	const Optimum* setOptimum(const wsum_t* opt);
	
	//! sets the next level to optimize
	/*!
	 * \return 
	 *   - true if there is a next level to optimize
	 *   - false otherwise
	 */
	bool optimizeNext() { return optimizeNext(opt_); }
	//@}
private:
	typedef std::atomic<Optimum*> OptPtr;
	SumVec       sum_;     // initial sum
	Optimum      opts_[2]; // buffers for update via "double buffering"
	OptPtr       opt_;     // active optimum
	Optimum*     next_;    // used during update of optimum
	wsum_t       low_;     // current lower bound
	uint32       hierarch_;// hierarchical optimization?
	MinimizeMode mode_;    // how to compare assignments?
public:
	WeightVec       weights;  // sparse vectors of weights - only used for multi-level constraints
	WeightLiteral   lits[0];  // (shared) literals - terminated with posLit(0)
private: 
	SharedMinimizeData(const SharedMinimizeData&);
	SharedMinimizeData& operator=(const SharedMinimizeData&);
	void setStep(Optimum* opt);
	bool optimizeNext(Optimum* opt);
};

//! Helper class for creating minimize constraints
class MinimizeBuilder {
public:
	typedef SharedMinimizeData SharedData;
	MinimizeBuilder();
	~MinimizeBuilder();
	
	bool                hasRules() const { return !initial_.empty(); }
	uint32              numRules() const { return (uint32)initial_.size(); }

	//! adds a minimize statement
	/*!
	 * \param lits the literals of the minimize statement
	 * \param adjustSum the initial sum of the minimize statement
	 */
	MinimizeBuilder&    addRule(const WeightLitVec& lits, wsum_t adjustSum = 0);

	//! sets the initial optimum of level lev to opt
	MinimizeBuilder&    setOptimum(uint32 lev, wsum_t opt);
	
	//! creates a new data object from previously added minimize statements
	/*!
	 * The function creates a new minimize data object from 
	 * the previously added minimize statements. The returned
	 * object can be used to attach one or more MinimizeConstraints.
	 * \param ctx A ctx object used to simplify minimize statements
	 * \param ma  An assumption that should be used for tagging reasons
	 */
	SharedData*         build(SharedContext& ctx, Literal ma = posLit(0));

	//! attaches a minimize constraint to a previously built data object 
	static MinimizeConstraint* attach(Solver& s, SharedData* data, int heuristic = 0);

	/*!
	 * Effect:
	 *  SharedData* d = build(ctx);
	 *  d->setMode(m);
	 *  return attach(*ctx.master(), d, heuristic);
	 */
	MinimizeConstraint* buildAndAttach(SharedContext& ctx, MinimizeMode m = MinimizeMode_t::optimize, int heuristic = 0);
private:
	struct Weight {
		Weight(uint32 lev, weight_t w) : level(lev), weight(w), next(0) {}
		uint32   level;
		weight_t weight;
		Weight*  next;
		static void free(Weight*& w);
	};
	typedef std::pair<Literal, Weight*> LitRep;
	typedef PodVector<LitRep>::type     LitRepVec;
	typedef SharedData::SumVec          SumVec;
	struct CmpByLit {
		bool operator()(const LitRep& lhs, const LitRep& rhs) const;
	};
	struct CmpByWeight {
		bool operator()(const LitRep& lhs, const LitRep& rhs) const;
		int  compare   (const LitRep& lhs, const LitRep& rhs) const;
	};
	void     prepare(SharedContext& ctx);
	void     addAssumption(Literal ma);
	void     adjustSum(LitRep l);
	void     mergeReduceWeight(LitRep& x, LitRep& by);
	weight_t addFlattened(SharedData::WeightVec& x, const Weight& w);
	bool     eqWeight(const SharedData::LevelWeight* lhs, const Weight& rhs);
	LitRepVec lits_;
	SumVec    initial_;
	SumVec    opts_;
	bool      ready_;
};

//! Implements a (meta-)constraint for supporting Smodels-like minimize statements
/*!
 * \ingroup constraint
 * A solver contains at most one minimize constraint, but a minimize constraint
 * may represent several minimize statements. In that case, each minimize statement
 * has a unique level (starting at 0) and minimize statements with a lower level
 * have higher priority. Priorities are handled like in smodels, i.e. statements
 * with lower priority become relevant only if all statements with higher priority
 * currently have an optimal assignment.
 *
 * MinimizeConstraint supports two modes of operation: if mode is set to
 * optimize, solutions are considered optimal only if they are strictly smaller 
 * than previous solutions. Otherwise, if mode is set to enumerate a
 * solution is valid if it is not greater than the initially set optimum.
 * Example: 
 *  m0: {a, b}
 *  m1: {c, d}
 *  All models: {a, c,...}, {a, d,...} {b, c,...}, {b, d,...} {a, b,...} 
 *  Mode = optimize: {a, c, ...} (m0 = 1, m1 = 1}
 *  Mode = enumerate and initial opt=1,1: {a, c, ...}, {a, d,...}, {b, c,...}, {b, d,...} 
 * 
 */
class MinimizeConstraint : public Constraint {
public:
	friend class MinimizeBuilder;
	typedef MinimizeBuilder::SharedData SharedData;
	typedef PodVector<wsum_t>::type SumVec;

	bool attach(Solver& s);

	//! number of minimize statements contained in this constraint
	uint32 numRules() const { return shared()->numRules(); }
	
	//! resets the local optimum to the one stored in the shared data object
	bool restoreOptimum();
	
	//! integrates the next tentative optimum into this constraint
	/*!
	 * Starting from the current optimum stored in the shared data object,
	 * the function tries to integrate the next candidate optimum into
	 * this constraint.
	 *
	 * \pre s is in a non-conflicting fix-point w.r.t propagation, i.e.
	 *      !s.hasConflict() && s.queueSize() == 0
	 *
	 * \return The function returns true if integration succeeded. Otherwise
	 * false is returned and s.hasConflict() is true
	 *
	 * \note If integrateNext() failed, the optimum of this constraint
	 *       is *not* changed. The caller has to resolve the conflict first
	 *       and then integrateNext() shall be called again. 
	 */
	bool integrateNext(Solver& s);

	//! enqueues all literals that are implied by a newly set optimum
	/*!
	 * \pre Neither the assignment nor this constraint are conflicting
	 * \note
	 *   The function removes decision levels until either all
	 *   implied literals where successfully added or the root
	 *   level is hit. In the latter case, false is returned.
	 * \note
	 *   The function *does not* call Solver::propagate(). Hence, the caller
	 *   is responsible for calling Solver::propagate() to further
	 *   propagate any implied literals.
	 */
	bool propagateNewOptimum(Solver& s);
	
	
	//! notifies the constraint about a model
	/*!
	 * setModel must be called whenever the solver finds a model.
	 * The found model is recorded as new optimum in the shared data
	 * object and the decision level on which search should continue is returned.
	 * \return The decision level on which the search should continue or
	 *  UINT32_MAX if search space is exhausted w.r.t this minimize-constraint.
	 *
	 * \note
	 *   Once decision levels are removed, propagateNewOptimum() must be called
	 *   in order to continue enumeration.
	 */
	uint32 setModel(Solver& s);
	
	const SharedData* shared() const { return shared_; }
	SharedData*       shared_unsafe(){ return shared_; }
	
	//! tries to assign as many literals from the last model to false
	bool modelHeuristic(Solver& s);

	// base interface
	PropResult  propagate(Solver& s, Literal p, uint32& data);
	void        undoLevel(Solver& s);
	void        reason(Solver& s, Literal p, LitVec& lits);
	bool        minimize(Solver& s, Literal p, CCMinRecursive* r);
	Constraint* cloneAttach(Solver& other);
	void        destroy(Solver*, bool);

	// FOR TESTING ONLY!
	wsum_t sum(uint32 i) const { return sum_[i]; }
protected:
	MinimizeConstraint(SharedData* d, unsigned heuristic);
	~MinimizeConstraint();
	union UndoInfo {
		UndoInfo() : rep(0) {}
		struct {
			uint32 idx    : 30; // index of literal on stack
			uint32 newDL  :  1; // first literal of new DL?
			uint32 idxSeen:  1; // literal with idx already propagated?
		}      data;
		uint32 rep;
		uint32 index()      const { return data.idx; }
		bool   newDL()      const { return data.newDL != 0u; }
	};
	typedef const WeightLiteral*         Iter;
	uint32    lastUndoLevel(const Solver& s) const;
	void      pushUndo(Solver& s, uint32 litIdx);
	enum PropMode { propagate_new_sum, propagate_new_opt };
	bool      litSeen(uint32 i) const { return undo_[i].data.idxSeen != 0; }
	uint32    updateOpt(Solver* s, bool applyStep, bool model);
	bool      propagateImpl(Solver& s, PropMode m);
	uint32    computeImplicationSet(Solver& s, const WeightLiteral& it, uint32& undoPos);
	uint32&   active(wsum_t* s) { return active_[s == temp_]; }
	bool      greater(wsum_t* lhs, wsum_t* rhs, uint32 len) const {
		while (*lhs == *rhs && --len) { ++lhs, ++rhs; }
		return *lhs > *rhs;
	}
	// Arithmetic operations
	enum   ArithType { SINGLE_LEVEL_ARITH = 0, MULTI_LEVEL_ARITH = 1};
	void   assign(wsum_t* lhs, wsum_t* rhs);
	void   add(weight_t wOrIdx);
	void   subtract(wsum_t* lhs, weight_t wOrIdx);
	bool   implied(wsum_t* lhs, weight_t wOrIdx);
	uint32 convert(wsum_t*)         { return numRules(); }
	const  wsum_t* sumToOpt() const { return sum_; }
	
	wsum_t*     sum_;       // current sum
	wsum_t*     opt_;       // current (local) optimum
	wsum_t*     temp_;      // temporary sum; used in propagateImpl() to compute implication level
	wsum_t*     tempOpt_;   // temporary opt; used in integrateOptimum()
	SharedData* shared_;    // pointer to (read-only) shared data
	Iter        pos_;       // position of literal to look at next
	UndoInfo*   undo_;      // one "seen" flag for each literal +
	uint32      undoTop_;   // undo stack holding assigned literals
	uint32      posTop_;    // stack of saved "look at" positions
	uint32      active_[2]; // first level to look at (one for sum_ and one for temp_)
	const ArithType type_;  // type of arithmetic operations
	bool        signHeu_;   // set preferred sign of literals?
	bool        modelHeu_;  // use model heuristic?
	bool        ownsShared_;// true if this constraints "owns" the shared data
};
} // end namespace Clasp

#endif
