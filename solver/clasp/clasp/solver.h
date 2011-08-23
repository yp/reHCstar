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
#ifndef CLASP_SOLVER_H_INCLUDED
#define CLASP_SOLVER_H_INCLUDED
#ifdef _MSC_VER
#pragma once
#pragma warning (disable : 4996) // 'std::copy': Function call with parameters that may be unsafe
#endif

#include <clasp/solver_types.h>
#include <clasp/shared_context.h>
#include <string>
#include <utility>

namespace Clasp { 
///////////////////////////////////////////////////////////////////////////////
// classes used by the Solver
///////////////////////////////////////////////////////////////////////////////
class   DecisionHeuristic;
class   PostPropagator;
class   SharedContext;
class   SatPreprocessor;

/**
 * \defgroup solver Solver and related classes
 */
//@{

//! Parameter-Object for configuring a solver.
struct SolverStrategies {
	typedef std::auto_ptr<DecisionHeuristic>  Heuristic;
	//! Clasp's two general search strategies
	enum SearchStrategy {
		use_learning, /*!< Analyze conflicts and learn First-1-UIP-clause */
		no_learning   /*!< Don't analyze conflicts - chronological backtracking */
	};
	//! Antecedents to consider during conflict clause minimization
	enum CflMinAntes {
		no_antes              = 0,  /*!< Don't minimize first-uip-clauses */
		binary_antes          = 2,  /*!< Consider only binary antecedents         */
		binary_ternary_antes  = 6,  /*!< Consider binary- and ternary antecedents */
		all_antes             = 7   /*!< Consider all antecedents                 */
	};
	enum ReduceAlgo {
		reduce_linear   = 0,
		reduce_in_place = 1,
		reduce_stable   = 2
	};
	enum ReduceScore {
		score_act,
		score_lbd,
		score_both
	};
	Heuristic           heuristic;
	SearchStrategy      search;
	RNG                 rng;
	int                 saveProgress;
	int                 reverseArcs;
	int                 otfs;
	CflMinAntes         cflMinAntes;
	int                 reduceAlgo;          /*!< shall be one of ReduceAlgo */
	int                 reduceScore;         /*!< shall be one of ReduceScore */
	uint32              reduceGlue;          /*!< clauses with lbd <= reduceGlue are exempt from deletion   */
	bool                strengthenRecursive; /*!< If true, use more expensive recursive nogood minimization */
	bool                randomWatches;
	bool                updateLbd;
	uint32              compress()    const { return compress_; }
	void setCompressionStrategy(uint32 len) {
		compress_ = len == 0 ? static_cast<uint32>(-1) : len;
	}
	int                 compareScore(const LearntConstraint::Activity& lhs, const LearntConstraint::Activity& rhs) const {
		switch (reduceScore) {
			case score_act: 
				if (lhs.scoreAct() != rhs.scoreAct()) { 
					return ((int)lhs.scoreAct()) - ((int)rhs.scoreAct()); 
				}
				// activity scores are equal - fall through
			case score_lbd:  
				if (lhs.scoreLbd() != rhs.scoreLbd()) {
					return ((int)lhs.scoreLbd()) - ((int)rhs.scoreLbd()); 
				}
				// lbd scores are equal - fall through
			case score_both: 
			default: break;
		}
		return ((int)lhs.scoreBoth()) - ((int)rhs.scoreBoth());
	}
	uint32              score(const LearntConstraint::Activity& act) const {
		switch (reduceScore) {
			case score_act:  return act.scoreAct();
			case score_lbd:  return act.scoreLbd();
			case score_both: return act.scoreBoth();
		}
		return 0;
	}
private:
	friend class Solver;
	SolverStrategies(const SolverStrategies&);
	SolverStrategies& operator=(const SolverStrategies&);
	//! creates a default-initialized object.
	SolverStrategies();
	uint32 compress_;
};

//! clasp's Solver class
/*!
 * A Solver-object maintains the state and provides the functions 
 * necessary to implement our CDNL-algorithm. 
 *
 * The interface supports incremental solving (search under assumptions) as well as
 * solution enumeration. To make this possible the solver maintains two special  
 * decision levels: the root-level and the backtrack-level.
 *
 * The root-level is the lowest decision level to which a search
 * can return. Conflicts on the root-level are non-resolvable and end a search - the
 * root-level therefore acts as an artificial top-level during search.
 * Incremental solving is then implemented by first adding a list of unit assumptions
 * and second initializing the root-level to the current decision level.
 * Once search terminates assumptions can be undone by calling clearAssumptions
 * and a new a search can be started using different assumptions.
 *
 * For model enumeration the solver maintains a backtrack-level which restricts
 * backjumping in order to prevent repeating already enumerated solutions.
 * The solver will never backjump above that level and conflicts on the backtrack-level 
 * are resolved by backtracking, i.e. flipping the corresponding decision literal.
 *
 * \see "Conflict-Driven Answer Set Enumeration" for a detailed description of this approach. 
 *
 */
class Solver {
	friend class SharedContext;
	void startInit(SharedContext* ctx, uint32 constraintGuess);
	void setEnumerationConstraint(Constraint* c) {
		if (enum_) enum_->destroy(this, true);
		enum_ = c;
	}
	bool endInit();
public:
	typedef PodVector<Constraint*>::type      ConstraintDB;
	/*!
	 * \name construction/destruction/setup
	 */
	//@{
	//! creates an empty solver object with all strategies set to their default value.
	Solver();
	
	//! destroys the solver object and all contained constraints.
	~Solver();
	
	//! sets the solver id - default is 0
	void   setId(uint32 id) { id_ = id; }
	uint32 id() const       { return id_; }
	
	//! resets a solver object to the state it had after default-construction.
	void reset();

	//! shuffle constraints upon next simplification
	void shuffleOnNextSimplify() { shuffle_ = true; }
	
	//! returns the strategies used by this solver-object. 
	SolverStrategies& strategies() { return strategy_; }

	/*!
	 * \overload SolverStrategies& Solver::strategies()
	 */
	const SolverStrategies& strategies() const { return strategy_; }

	//! returns a pointer to the shared context object of this solver
	const SharedContext* sharedContext() const { return shared_; }

	SatPreprocessor* satPrepro() const;

	//@}	

	/*!
	 * \name static state inspection
	 * Functions for inspecting the static state of the solver
	 * \note validVar(v) is a precondition for all functions that take a variable as 
	 * parameter.
	 */
	//@{
	
	//! returns the number of problem variables.
	/*!
	 * \note The special sentine-var 0 is not counted, i.e. numVars() returns
	 * the number of problem-variables.
	 * To iterate over all problem variables use a loop like:
	 * \code
	 * for (Var i = 1; i <= s.numVars(); ++i) {...}
	 * \endcode
	 */
	uint32 numVars() const { return assign_.numVars() - 1; }

	//! returns true if var represents a valid variable in this solver.
	/*!
	 * \note The range of valid variables is [1..numVars()]. The variable 0
	 * is a special sentinel variable that is always true. 
	 */
	bool validVar(Var var) const { return var <= numVars(); }

	//! number of problem constraints in this solver
	uint32 numConstraints() const;

	//@}
	struct RestartLimit {
		RestartLimit(uint64 c = UINT64_MAX, bool loc = false) : maxConf(c), countLocal(loc) {}
		uint64 maxConf;
		bool   countLocal;
		bool   reached() const { return maxConf == 0; }
	};
	struct LearntLimit {
		LearntLimit(uint64 c = UINT64_MAX, uint32 n = UINT32_MAX) : maxConf(c), maxLearnt(n) {}
		uint64 maxConf;
		uint32 maxLearnt;
		bool   reached() const { return maxConf == 0; }
	};
	/*!
	 * \name DPLL-functions
	 */
	//@{
	//! searches for a model as long as none of the given limits is reached.
	/*!
	 * The search function implements clasp's CDNL-algorithm.
	 * It searches for a model as long as neither of the limits given by r and d
	 * is reached. The limits r and d are updated during search.
	 *
	 * \pre endInit() returned true and !hasConflict()
	 * \param r restart limit  (stop after r.maxConf global or local conflicts)
	 * \param d deletion limit (stop after d.maxConf or if numLearnts() > d.maxLearnt)
	 * \param randProb pick next decision variable randomly with a probability of randProp
	 * \return
	 *  - value_true: if a model was found.
	 *  - value_false: if the problem is unsatisfiable (under assumptions, if any)
	 *  - value_free: if search was stopped because limit was reached. 
	 *  .
	 *
	 * \note search treats the root level as top-level, i.e. it
	 * will never backtrack below that level.
	 */
	ValueRep search(RestartLimit& r, LearntLimit& d, double randProb = 0.0);
	ValueRep search(uint64 maxC, uint32 maxL, bool local = false, double rp  = 0.0) {
		RestartLimit r(maxC, local);
		LearntLimit  d(UINT64_MAX, maxL);
		return search(r, d, rp);
	}

	//! moves the root-level i levels down (i.e. away from the top-level)
	/*!
	 * The root-level is similar to the top-level in that it cannot be
	 * undone during search, i.e. the solver will not resolve conflicts that are on or
	 * above the root-level. 
	 */
	void pushRootLevel(uint32 i = 1) { 
		rootLevel_  = std::min(decisionLevel(), rootLevel_+i); 
		btLevel_    = std::max(btLevel_, rootLevel_);
	}
	//! returns the current root level
	uint32 rootLevel() const { return rootLevel_; }

	//! moves the root-level i levels up (i.e. towards the top-level)
	/*!
	 * The function undos all levels between the new root level and the current decision level,
	 * resets the current backtrack-level, and asserts any implied literals.
	 * \param i number of levels to pop
	 * \post decisionLevel() == rootLevel()
	 * \note The function first calls clearStopConflict() to remove any stop conflicts.
	 * \note If the solver has a root-level conflict prior to the call and resolveConflict
	 *       is true, the function tries to resolve that conflict. 
	 *       If this fails, false is returned.
	 * \note The function *does not* propagate any asserted literals. It is
	 *       the caller's responsibility to call propagate() after the function returned.
	 */
	bool popRootLevel(uint32 i = 1, bool resolveConflict = true);

	//! removes a previously set stop conflict and restores the root level
	void clearStopConflict();

	//! removes any implications made between the top-level and the root-level.
	/*!
	 * The function also resets the current backtrack-level and re-asserts learnt
	 * facts.
	 * \note
	 *   Equivalent to popRootLevel(rootLevel(), false) followed
	 *   by simplify().
	 */
	bool clearAssumptions();

	//! sets the backtracking level to dl
	void setBacktrackLevel(uint32 dl) {
		btLevel_    = std::max(std::min(dl, decisionLevel()), rootLevel_);
	}
	uint32 backtrackLevel() const { return btLevel_; }

	//! returns whether the solver can split-off work from its current guiding path
	bool   splittable() const { return decisionLevel() > rootLevel() && !frozenLevel(rootLevel()+1); }

	//! copies the solver's currrent guiding path to gp
	/*!
	 * \note The solver's guiding path consists of:
	 *   - the decisions from levels [1, rootLevel()]
	 *   - any literals that are implied on a level <= rootLevel() because of newly learnt
	 *     information. This particularly includes literals that were flipped during model enumeration.
	 * 
	 * \param[in,out] gp          where to store the guiding path
	 * \param[in,out] startPos    position in the trail from which copying should start
	 * \param[in,out] numImplied  number of newly implied literals in gp
	 */
	void   updateGuidingPath(LitVec& gp, LitVec::size_type& startPos, uint32& numImplied);

	//! If called on top-level removes SAT-clauses + Constraints for which Constraint::simplify returned true
	/*!
	 * \note if this method is called on a decision-level > 0 it is a noop and will
	 * simply return true.
	 * \return fallse, if a top-level conflict is detected. Otherwise, true.
	 */
	bool simplify();

	//! removes all conditional knowledge, i.e. all previously tagged learnt clauses
	/*!
	 * \see SharedContext::requestTagLiteral()
	 */
	void removeConditional();

	//! resolves all tagged clauses with the tag literal and thereby strengthens the learnt db
	/*!
	 * \see SharedContext::requestTagLiteral()
	 */
	void strengthenConditional();

	//! sets the literal p to true and schedules p for propagation.
	/*!
	 * Setting a literal p to true means assigning the appropriate value to
	 * p's variable. That is: value_false if p is a negative literal and value_true
	 * if p is a positive literal.
	 * \param p the literal that should become true
	 * \param a the reason for the literal to become true or 0 if no reason exists.
	 * 
	 * \return
	 *  - false if p is already false
	 *  - otherwise true.
	 *
	 * \pre hasConflict() == false
	 * \pre a.isNull() == false || decisionLevel() <= rootLevel() || SearchStrategy == no_learning
	 * \post
	 *  p.var() == trueValue(p) || p.var() == falseValue(p) && hasConflict() == true
	 *
	 * \note if setting p to true leads to a conflict, the nogood that caused the
	 * conflict can be requested using the conflict() function.
	 */
	bool force(const Literal& p, const Antecedent& a) {
		assert((!hasConflict() || isTrue(p)) && !shared_->eliminated(p.var()));
		if (assign_.assign(p, decisionLevel(), a)) return true;
		setConflict(p, a, UINT32_MAX);
		return false;
	}
	/*!
	 * \overload bool Solver::force(const Literal&, const Antecedent&)
	 * \pre data == UINT32_MAX || SharedContext::requestData(p.var()) was called during setup
	 */
	bool force(const Literal& p, const Antecedent& a, uint32 data) {
		return data != UINT32_MAX 
			? assign_.assign(p, decisionLevel(), a.constraint(), data) || (setConflict(p, a, data), false)
			: force(p, a);
	}

	//! assigns p at max(dl, backtrackLevel()) because of r
	/*!
	 * \pre dl <= decisionLevel()
	 * \note if dl < backtrackLevel(), p is actually assigned at backtrackLevel() but
	 * the solver stores enough information to reassign p on backtracking
	 */
	bool force(Literal p, uint32 dl, const Antecedent& r, uint32 d = UINT32_MAX) {
		if (dl == decisionLevel() || undoUntil(dl, false) == dl) { return force(p, r, d); }
		if (isTrue(p)) {
			// p is currrently implied by some other constraint but
			// the new constraint is stronger
			setReason(p, r, d);
		}
		// Logically the implication is on level dl but backjumping is
		// bounded by btLevel_ thus p is asserted on the current level.
		// We store enough information so that p can be re-asserted once we backtrack.
		impliedLits_.push_back(ImpliedLiteral( p, dl, r, d ));
		return force(p, r, d);
	}

	//! assumes the literal p if possible
	/*!
	 * If p is currently unassigned, sets p to true and starts a new decision level.
	 * \pre validVar(p.var()) == true
	 * \param p the literal to assume
	 * \return !isFalse(p)
	 */
	bool assume(const Literal& p);

	//! selects and assumes the next branching literal by calling the installed decision heuristic.
	/*!
	 * \pre queueSize() == 0
	 * \note the next decision literal will be selected randomly with a
	 * probability set using the initRandomHeuristic method.
	 * \return 
	 *  - true if the assignment is not total and a literal was assumed (or forced).
	 *  - false otherwise
	 *  .
	 * \see DecisionHeuristic
	 */
	bool decideNextBranch();

	//! Sets a conflict that forces the solver to terminate its search
	/*!
	 * \pre  !hasConflict()
	 * \post hasConflict()
	 *
	 * \note 
	 *   To prevent the solver from resolving the stop conflict, the
	 *   function sets the root level to the current decision level. 
	 *   Call clearStopConflict() to remove the conflict and to restore
	 *   the previous root-level.
	 */
	void setStopConflict();

	/*!
	 * propagates all enqueued literals. If a conflict arises during propagation
	 * propagate returns false and the current conflict (as a set of literals)
	 * is stored in the solver's conflict variable.
	 * \pre !hasConflict()
	 * \see Solver::force
	 * \see Solver::assume
	 * \note shall not be called recursively
	 */
	bool propagate();

	/*!
	 * Does unit propagation and calls x->propagateFixpoint(*this)
	 * for all post propagators up to but not including p.
	 * \note The function is meant to be called only in the context of p
	 */
	bool propagateUntil(PostPropagator* p) { return unitPropagate() && (p == post_.head || post_.propagate(*this, p)); }

	//! executes a one-step lookahead on p.
	/*!
	 * Assumes p and propagates this assumption. If propagations leads to
	 * a conflict false is returned. Otherwise the assumption is undone and 
	 * the function returns true.
	 * \param p the literal to test
	 * \param c the constraint that wants to test p (may be 0)
	 * \pre p is free
	 * \note If c is not null and testing p does not lead to a conflict, 
	         c->undoLevel() is called *before* p is undone. Hence, the
			 range [s.levelStart(s.decisionLevel()), s.assignment().size())
			 contains p followed by all literals that were forced because of p.
	 * \note During propagation of p, only post propagators with priority
	 * in the range [priority_highest, priority_lookahead) are called.
	 */
	bool test(Literal p, Constraint* c);

	//! estimates the number of assignments following from setting p to true.
	/*!
	 * \note for the estimate only binary clauses are considered.
	 */ 
	uint32 estimateBCP(const Literal& p, int maxRecursionDepth = 5) const;
	
	//! removes upto remMax percent of the learnt nogoods.
	/*!
	 * \param remMax percantage of learnt nogoods that should be removed ([0.0f, 1.0f])
	 * \note nogoods that are the reason for a literal to be in the assignment
	 * are said to be locked and won't be removed.
	 */
	void reduceLearnts(float remMax);

	//! resolves the active conflict using the selected strategy
	/*!
	 * If the SearchStrategy is set to learning, resolveConflict implements
	 * First-UIP learning and backjumping. Otherwise it simply applies
	 * chronological backtracking.
	 * \pre hasConflict
	 * \return
	 *  - true if the conflict was successfully resolved
	 *  - false otherwise
	 * \note
	 *  if decisionLevel() == rootLevel() false is returned.
	 */
	bool resolveConflict();

	//! backtracks the last decision and sets the backtrack-level to the resulting decision level.
	/*!
	 * \return
	 *  - true if backtracking was possible
	 *  - false if decisionLevel() == rootLevel()
	 */
	bool backtrack();

	//! undoes all assignments up to (but not including) decision level dl.
	/*!
	 * \pre dl > 0 (assignments made on decision level 0 cannot be undone)
	 * \pre dl <= decisionLevel()
	 * \post decisionLevel == max(dl, max(rootLevel(), btLevel))
	 */
	void undoUntil(uint32 dl);

	/*!
	 * similar to undoUntil but also pops the backtrack-level
	 * to dl if possible
	 * \return the current decision level
	 */
	uint32 undoUntil(uint32 dl, bool popBtLevel);

	//! checks whether there is a model that is symmetric to the current model
	/*!
	 * The function checks for symmetric models, i.e. models that differ only in the 
	 * assignment of variables outside of the solver's assignment. 
	 * Typical example: vars eliminated by the SAT-preprocessor
	 * \param expand if false, any symmetric models are ignored. Otherwise, symmetric models
	 *        are expanded and stored in the solver.
	 * \pre the current assignment is a model
	 */ 
	bool nextSymModel(bool expand);

	//@}  

	/*!
	 * \name learning functions
	 * Functions controling learning
	 */
	//@{
		
	//! wraps SharedContext::add() - only provided for convenience
	void add(Constraint* c);
	
	Constraint* getEnumerationConstraint() const { return enum_; }

	//! adds the unary constraint p to the solver.
	/*!
	 * \note unary constraints are immediately asserted.
	 * \return false if asserting p leads to a conflict.
	 */
	bool addUnary(Literal p, ConstraintType);

	//! wraps SharedContext::addBinary() resp. SharedContext::learnBinary() depending on the type t
	bool addBinary(Literal p, Literal q, ConstraintType t);
	
	//! wraps SharedContext::addTernary() resp. SharedContext::learnTernary() depending on the type t
	bool addTernary(Literal p, Literal q, Literal r, ConstraintType t);
	
	//! adds the learnt constraint c to the solver.
	/*!
	 * \pre endInit() was called.
	 */
	void addLearnt(LearntConstraint* c, uint32 size, bool tagged = false) {
		learnts_.push_back(c);
		stats.addLearnt(size, c->type()); 
		tagged_ += tagged;
	}
	
	//! adds p as post propagator to this solver
	/*!
	 * \pre p was not added previously and is not part of any other solver
	 * \note post propagators are stored in priority order
	 * \see PostPropagator::priority()
	 */
	void addPost(PostPropagator* p)    { post_.add(p); }

	//! removes p from the solver's list of post propagators
	/*!
	 * \note removePost(p) shall only be called during propagation
	 *       of p or if no propagation is currently active.
	 */
	void removePost(PostPropagator* p) { post_.remove(p); }

	//! returns the idx'th learnt constraint
	/*!
	 * \pre idx < numLearntConstraints()
	 */
	LearntConstraint& getLearnt(uint32 idx) const {
		assert(idx < numLearntConstraints());
		return *static_cast<LearntConstraint*>(learnts_[ idx ]);
	}

	bool ccMinimize(Literal p, CCMinRecursive* rec) {
		return seen(p.var()) || 
			(rec && hasLevel(level(p.var())) && rec->checkRecursive(p));
	}

	uint32 computeLbd(Literal p, const Literal* first, const Literal* last);

	//@}

	/*!
	 * \name dynamic state inspection
	 * Functions for inspecting the state of the search.
	 * \note validVar(v) is a precondition for all functions that take a variable as 
	 * parameter.
	 */
	//@{
	
	//! returns the number of assigned variables.
	/*!
	 * \note The special sentinel-var 0 is not counted.
	 */
	uint32 numAssignedVars()  const { return assign_.assigned(); }
	
	//! returns the number of free variables.
	/*!
	 * The number of free variables is the number of vars that are neither
	 * assigned nor eliminated.
	 */
	uint32 numFreeVars()      const { return assign_.free()-1; }

	//! returns the number of constraints that are currently in the solver's learnt database.
	uint32 numLearntConstraints() const { return (uint32)learnts_.size(); }
	
	//! returns the value of v w.r.t the current assignment
	ValueRep value(Var v) const {
		assert(validVar(v));
		return assign_.value(v);
	}
	//! returns the previous value of v
	ValueRep savedValue(Var v) const {
		assert(validVar(v));
		return assign_.saved(v);
	}
	//! returns true if p is true w.r.t the current assignment
	bool isTrue(Literal p) const {
		assert(validVar(p.var()));
		return assign_.value(p.var()) == trueValue(p);
	}
	//! returns true if p is false w.r.t the current assignment
	bool isFalse(Literal p) const {
		assert(validVar(p.var()));
		return assign_.value(p.var()) == falseValue(p);
	}
	//! returns the decision level on which v was assigned.
	/*!
	 * \note the returned value is only meaningful if value(v) != value_free
	 */
	uint32 level(Var v) const {
		assert(validVar(v));
		return assign_.level(v);
	}

	//! returns the reason for p being true
	/*!
	 * \pre p is true w.r.t the current assignment
	 * \note if solver is used in no-learning mode, the reason is always null
	 */
	const Antecedent& reason(Literal p) const {
		assert(isTrue(p));
		return assign_.reason(p.var());
	}
	//! returns the reason for p being true as a set of literals
	void reason(Literal p, LitVec& out) {
		assert(isTrue(p)); out.clear();
		return assign_.reason(p.var()).reason(*this, p, out);
	}
	//! returns the additional reason data associated with p
	uint32 reasonData(Literal p) const {
		return assign_.data(p.var());
	}
	//! updates the reason for p being tue
	/*!
	 * \pre p is true and x is a valid reason for p
	 */
	void setReason(Literal p, const Antecedent& x, uint32 data = UINT32_MAX) {
		assert(isTrue(p));
		assign_.setReason(p.var(), x);
		if (data != UINT32_MAX) { assign_.setData(p.var(), data); }
	}
	//! returns the literal of v being true in the current assignment
	/*!
	 * \pre v is assigned a value in the current assignment
	 */
	Literal trueLit(Var v) const {
		assert(value(v) != value_free);
		return Literal(v, value(v) == value_false);
	}

	//! returns true if v is currently marked as seen.
	/*!
	 * Note: variables assigned on level 0 are always marked.
	 */
	bool seen(Var v) const {
		assert(validVar(v));
		return assign_.seen(v, 3u);
	}
	//! returns true if the literal p is currently marked as seen
	bool seen(Literal p) const {
		assert(validVar(p.var()));
		return assign_.seen(p.var(), uint8(1+p.sign()));
	}
	void markSeen(Var v)    { assert(validVar(v)); assign_.setSeen(v, 3u); }
	void markSeen(Literal p){ assert(validVar(p.var())); assign_.setSeen(p.var(), uint8(1+p.sign())); }
	void clearSeen(Var v)   { assert(validVar(v)); assign_.clearSeen(v);  }

	//! returns the current decision level.
	uint32 decisionLevel() const { return (uint32)levels_.size(); }

	//! returns the starting position of decision level dl in the trail.
	/*!
	 * \pre dl != 0 && dl <= decisionLevel()
	 */
	uint32 levelStart(uint32 dl) const { 
		assert(dl != 0 && dl <= decisionLevel() );
		return levels_[dl-1].trailPos;
	}

	//! returns the decision literal of the decision level dl.
	/*!
	 * \pre dl != 0 && dl <= decisionLevel()
	 */
	Literal decision(uint32 dl) const {
		assert(dl != 0 && dl <= decisionLevel() );
		return assign_.trail[ levels_[dl-1].trailPos ];
	}

	void markLevel(uint32 dl) {
		assert(dl != 0 && dl <= decisionLevel() );
		levels_[dl-1].marked = 1;
	}
	void freezeLevel(uint32 dl) {
		assert(dl != 0 && dl <= decisionLevel() );
		levels_[dl-1].freeze = 1;
	}
	void unmarkLevel(uint32 dl) {
		assert(dl != 0 && dl <= decisionLevel() );
		levels_[dl-1].marked = 0;
	}
	void unfreezeLevel(uint32 dl) {
		assert(dl != 0 && dl <= decisionLevel() );
		levels_[dl-1].freeze = 0;
	}
	bool hasLevel(uint32 dl) const {
		assert(dl != 0 && dl <= decisionLevel() );
		return levels_[dl-1].marked != 0;
	}
	bool frozenLevel(uint32 dl) const {
		assert(dl != 0 && dl <= decisionLevel() );
		return levels_[dl-1].freeze != 0;
	}
	
	//! returns the current (partial) assignment as a set of true literals
	/*!
	 * \note although the special var 0 always has a value it is not considered to be
	 * part of the assignment.
	 */
	const LitVec&     trail()      const { return assign_.trail; }
	const Assignment& assignment() const { return assign_; }
	
	//! returns true, if the current assignment is conflicting
	bool hasConflict() const { return !conflict_.empty(); }
	//! returns the current conflict as a set of literals
	const LitVec& conflict() const { return conflict_; }
	//! returns the most recently derived conflict clause
	const LitVec& conflictClause() const { return cc_; }
	
	//! number of (unprocessed) literals in the propagation queue
	uint32 queueSize()    const { return (uint32) assign_.qSize(); }
	
	uint32 lastSimplify() const { return (uint32) lastSimplify_; }
	void   simplifyDB(ConstraintDB& db);
	bool   minimizeLitRedundant(Literal p);

	SolveStats  stats;  /**< stores some statistics about the solving process */
	//@}

	/*!
	 * \name watch management
	 * Functions for setting/removing watches.
	 * \pre validVar(v)
	 */
	//@{
	//! returns the number of constraints watching the literal p
	uint32 numWatches(Literal p) const;
	//! returns true if the constraint c watches the literal p
	bool    hasWatch(Literal p, Constraint* c) const;
	bool    hasWatch(Literal p, ClauseHead* c) const;
	//! returns c's watch-structure associated with p
	/*!
	 * \note returns 0, if hasWatch(p, c) == false
	 */
	GenericWatch* getWatch(Literal p, Constraint* c) const;

	//! Adds c to the watch-list of p.
	/*!
	 * When p becomes true, c->propagate(p, data, *this) is called.
	 * \post hasWatch(p, c) == true
	 */
	void addWatch(Literal p, Constraint* c, uint32 data = 0) {
		assert(validWatch(p));
		watches_[p.index()].push_right(GenericWatch(c, data));
	}

	//! Adds w to the clause watch-list of p.
	void addWatch(Literal p, const ClauseWatch& w) {
		assert(validWatch(p));
		watches_[p.index()].push_left(w);
	}
	
	//! removes c from p's watch-list.
	/*!
	 * \post hasWatch(p, c) == false
	 */
	void removeWatch(const Literal& p, Constraint* c);
	void removeWatch(const Literal& p, ClauseHead* c);

	//! adds c to the watch-list of decision-level dl
	/*!
	 * Constraints in the watch-list of a decision level are
	 * notified when that decision level is about to be backtracked.
	 * \pre dl != 0 && dl <= decisionLevel()
	 */
	void addUndoWatch(uint32 dl, Constraint* c) {
		assert(dl != 0 && dl <= decisionLevel() );
		if (levels_[dl-1].undo != 0) {
			levels_[dl-1].undo->push_back(c);
		}
		else {
			levels_[dl-1].undo = allocUndo(c);
		}
	}
	
	//! removes c from the watch-list of the decision level dl
	void removeUndoWatch(uint32 dl, Constraint* c);

	void initSavedValue(Var v, ValueRep val) {
		assert(validVar(v));
		if (assign_.saved(v) == value_free) {
			assign_.setSavedValue(v, val);
		}
	}
	
	void* allocSmall()       { return smallAlloc_->allocate(); }
	void  freeSmall(void* m) { smallAlloc_->free(m);    }
	//@}
private:
	struct DLevel {
		explicit DLevel(uint32 pos = 0, ConstraintDB* u = 0)
			: trailPos(pos)
			, marked(0)
			, freeze(0)
			, undo(u) {}
		uint32        trailPos : 30;
		uint32        marked   :  1;
		uint32        freeze   :  1;
		ConstraintDB* undo;
	};
	typedef PodVector<DLevel>::type         DecisionLevels;  
	typedef PodVector<ImpliedLiteral>::type ImpliedLits;
	typedef PodVector<Antecedent>::type     ReasonVec;
	typedef PodVector<WatchList>::type      Watches;
	struct PPList {
		PPList();
		~PPList();
		void add(PostPropagator* p);
		void remove(PostPropagator* p);
		bool propagate(Solver& s, PostPropagator* p);
		void simplify(Solver& s, bool shuf);
		void reset();
		bool isModel(Solver& s);
		bool nextSymModel(Solver& s, bool expand);
		bool init(Solver& s);
		PostPropagator* head;
		PostPropagator* look;
		PostPropagator* saved;
	};
	struct CmpScore {
		CmpScore(const ConstraintDB& learnts, const SolverStrategies& st) : db(learnts), strat(st) {}
		bool operator()(LearntConstraint::Activity lhs, LearntConstraint::Activity rhs) const { 
			return strat.compareScore(lhs, rhs) < 0; 
		}
		bool operator()(uint32 lhsId, uint32 rhsId) const { return (*this)(db[lhsId], db[rhsId]); }
		bool operator()(const Constraint* lhs, const Constraint* rhs) const {
			return strat.compareScore(
				static_cast<const LearntConstraint*>(lhs)->activity(), 
				static_cast<const LearntConstraint*>(rhs)->activity()) < 0;
		}
		const ConstraintDB&     db;
		const SolverStrategies& strat;      
	private: CmpScore& operator=(const CmpScore&);
	};
	inline  bool hasStopConflict() const;
	bool    validWatch(Literal p) const { return p.index() < (uint32)watches_.size(); }
	uint32  mark(uint32 s, uint32 e);
	void    initRandomHeuristic(double randFreq);
	void    freeMem();
	bool    simplifySAT();
	void    simplifyShort(Literal p);
	bool    unitPropagate();
	void    undoLevel(bool sp);
	uint32  analyzeConflict();
	void    otfs(Antecedent& lhs, const Antecedent& rhs, Literal p, bool final);
	ClauseHead* otfsRemove(ClauseHead* c, const LitVec* newC);
	uint32  finalizeConflictClause(ClauseHead* rhs);
	uint32  ccMinimize(LitVec& cc, LitVec& removed, uint32 anteMask, CCMinRecursive* ccMin);
	bool    ccRemovable(Literal p, uint32 anteMask, CCMinRecursive* ccMin);
	bool    ccReverseArc(Literal p, uint32 maxLevel, uint32 maxNew, Antecedent& out);
	void    undoFree(ConstraintDB* x);
	bool    forceImplied();
	void    setConflict(Literal p, const Antecedent& a, uint32 data);
	uint32  reduceLinear(uint32 minR, uint32 maxR);
	uint32  reduceSorted(uint32 minR, uint32 maxR);
	uint32  reduceSortedInPlace(uint32 minR, uint32 maxR);
	ConstraintDB* allocUndo(Constraint* c);
	SolverStrategies  strategy_;    // Strategies used by this solver-object
	Assignment        assign_;      // three-valued assignment
	ConstraintDB      constraints_; // problem constraints
	ConstraintDB      learnts_;     // learnt constraints
	ImpliedLits       impliedLits_; // Lits that were asserted on current dl but are logically implied earlier
	LitVec            conflict_;    // stores conflict-literals for later analysis
	LitVec            cc_;          // temporary: stores conflict clause within analyzeConflict
	LitVec            unconstr_;    // eliminated vars that are unconstraint w.r.t the current model
	DecisionLevels    levels_;      // Stores information (e.g. position in trail) on each decision level
	VarVec            lbdStamp_;    // temporary vector for computing LBD
	Watches           watches_;     // for each literal p: list of constraints watching p
	PPList            post_;        // (possibly empty) list of post propagators
	ClauseInfo        ccInfo_;      // temporary: stores information about conflict clause cc_
	CCMinRecursive*   ccMin_;       // additional data for supporting recursive strengthen
	SmallClauseAlloc* smallAlloc_;  // allocator object for small clauses
	SharedContext*    shared_;      // initialized by master thread - otherwise read-only!
	DecisionHeuristic*randHeuristic_;
	VarVec*           levConflicts_;// For a DL d, levConflicts_[d-1] stores num conflicts when d was started
	ConstraintDB*     undoHead_;    // Free list of undo DBs
	Constraint*       enum_;        // enumeration constraint - set by enumerator
	uint32            id_;          // solver id - only needed with multi-threading
	uint32            units_;       // number of top-level assignments: always marked as seen
	LitVec::size_type lastSimplify_;// number of top-level assignments on last call to simplify
	uint32            rootLevel_;   // dl on which search started.
	uint32            btLevel_;     // When enumerating models: DL of the last unflipped decision from the current model. Can't backjump below this level.
	uint32            tagged_;      // number of conditional clauses
	uint32            lbdTime_;     // temporary counter for computing lbd
	bool              shuffle_;     // shuffle program on next simplify?
};

inline bool isRevLit(const Solver& s, Literal p, uint32 maxL) {
	return s.isFalse(p) && (s.seen(p) || s.level(p.var()) < maxL);
}


//@}

/**
 * \defgroup heuristic Decision Heuristics
 */
//@{

//! Base class for decision heuristics to be used in a Solver.
/*! 
 * During search the decision heuristic is used whenever the DPLL-procedure must pick 
 * a new variable to branch on. Each concrete decision heuristic can implement a
 * different algorithm for selecting the next decision variable.
 */
class DecisionHeuristic {
public:
	DecisionHeuristic() {}
	virtual ~DecisionHeuristic(); 
	
	/*!
	 * Called once after all problem variables are known to the solver.
	 * The default-implementation is a noop.
	 * \param s The solver in which this heuristic is used.
	 */
	virtual void startInit(const Solver& /* s */) {}  

	/*!
	 * Called once after all problem constraints are known to the solver
	 * and the problem was simplified. 
	 * The default-implementation is a noop.
	 * \param s The solver in which this heuristic is used.
	 */
	virtual void endInit(Solver& /* s */) { }  

	/*!
	 * If the heuristic is used in an incremental setting, enable/disable
	 * reinitialization of existing variables.
	 * The default-implementation is a noop. Hence, heuristics will typically
	 * simply (re-)initialize all variables.
	 */
	virtual void reinit(bool /* b */) {}
	
	/*!
	 * Called for each var that changes its state from eliminated back to normal.
	 * The default-implementation is a noop.
	 * \param s Solver in which v is resurrected
	 * \param v The variable that is resurrected
	 */
	virtual void resurrect(const Solver& /* s */, Var /* v */) {}
	
	/*!
	 * Called on decision level 0. Variables that are assigned on this level
	 * may be removed from any decision heuristic.
	 * \note Whenever the solver returns to dl 0, simplify is only called again
	 * if the solver learnt new facts since the last call to simplify.
	 *
	 * \param s The solver that reached decision level 0.
	 * \param st The position in the trail of the first new learnt fact.
	 */
	virtual void simplify(const Solver& /* s */, LitVec::size_type /* st */) { }
	
	/*!
	 * Called whenever the solver backracks.
	 * Literals in the range [s.trail()[st], s.trail().size()) are subject to backtracking.
	 * The default-implementation is a noop.
	 * \param s The solver that is about to backtrack
	 * \param st Position in the trail of the first literal that will be backtracked.
	 */
	virtual void undoUntil(const Solver& /* s */, LitVec::size_type /* st */) {}
	
	/*!
	 * Called whenever a new constraint is added to the solver s.
	 * The default-implementation is a noop.
	 * \param s The solver to which the constraint is added.
	 * \param first First literal of the new constraint.
	 * \param size Size of the new constraint.
	 * \param t Type of the new constraint.
	 * \note first points to an array of size size.
	 */
	virtual void newConstraint(const Solver&, const Literal* /* first */, LitVec::size_type /* size */, ConstraintType /* t */) {}
	
	/*!
	 * Called for each new reason-set that is traversed during conflict analysis.
	 * The default-implementation is a noop.
	 * \param s the solver in which the conflict is analyzed.
	 * \param lits The current reason-set under inspection.
	 * \param resolveLit The literal that is currently resolved.
	 * \note When a conflict is detected the solver passes the conflicting literals
	 * in lits during the first call to updateReason. On that first call resolveLit
	 * is the sentinel-literal.
	 */
	virtual void updateReason(const Solver& /* s */, const LitVec& /* lits */, Literal /* resolveLit */) {}
	
	/*! 
	 * Called whenever the solver must pick a new variable to branch on. 
	 * \param s The solver that needs a new decision variable.
	 * \return
	 *  - true  : if the decision heuristic assumed a literal 
	 *  - false : if no decision could be made because assignment is total or there is a conflict
	 *  .
	 * \post
	 * if true is returned, the heuristic has asserted a literal.
	 */
	bool select(Solver& s) { return s.numFreeVars() != 0 && s.assume(doSelect(s)); }

	//! implements the actual selection process
	/*!
	 * \pre s.numFreeVars() > 0, i.e. there is at least one variable to branch on.
	 * \return 
	 *  - a literal that is currently free or
	 *  - a sentinel literal. In that case, the heuristic shall have asserted a literal!
	 */ 
	virtual Literal doSelect(Solver& /* s */) = 0;

	/*! 
	 * Shall select one of the literals in the range [first, last).
	 * \param s     The solver that needs a new decision variable.
	 * \param first Pointer to first literal in range
	 * \param last  Pointer to the end of the range
	 * \pre [first, last) is not empty and all literals in the range are currently unassigned.
	 * \note The default implementation returns *first
	 */
	virtual Literal selectRange(Solver& /* s */, const Literal* first, const Literal* /* last */) {
		return *first;
	}
protected:
	Literal savedLiteral(const Solver& s, Var var) const {
		Literal r(0,false); ValueRep v;
		if ((v = s.savedValue(var)) != value_free) {
			r = Literal(var, v == value_false);
		}
		return r;
	}
private:
	DecisionHeuristic(const DecisionHeuristic&);
	DecisionHeuristic& operator=(const DecisionHeuristic&);
};

//! selects the first free literal w.r.t to the initial variable order.
class SelectFirst : public DecisionHeuristic {
private:
	Literal doSelect(Solver& s);
};

//@}
}
#endif
