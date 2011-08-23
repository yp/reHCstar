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
#ifndef CLASP_CLAUSE_H_INCLUDED
#define CLASP_CLAUSE_H_INCLUDED

#ifdef _MSC_VER
#pragma warning (disable : 4200) // nonstandard extension used : zero-sized array
#pragma once
#endif

#include <clasp/solver_types.h>

namespace Clasp { 

//! An array of literals that can be shared
class SharedLiterals {
public:
	//! creates a shareable (ref-counted) object containing the literals in lits
	/*!
	 * \note The reference count is set to numRefs
	 */
	static SharedLiterals* newShareable(const LitVec& lits, ConstraintType t, uint32 numRefs = 1) {
		return newShareable(!lits.empty() ? &lits[0]:0, static_cast<uint32>(lits.size()), t, numRefs);
	}
	static SharedLiterals* newShareable(const Literal* lits, uint32 size, ConstraintType t, uint32 numRefs = 1);
	
	//! returns a pointer to the beginning of the literal array
	const Literal* begin() const { return lits_; }
	//! returns a pointer to the end of the literal array
	const Literal* end()   const { return lits_+size(); }
	//! returns the number of literals in the array
	uint32         size()  const { return size_type_ >> 2; }
	//! returns the type of constraint from which the literals originated
	ConstraintType type()  const { return ConstraintType( size_type_ & uint32(3) ); }
	//! simplifies the literals w.r.t to the assignment in s
	/*!
	 * returns the number of non-false literals in this object or 0 if
	 * the array contains a true literal.
	 * \note If the object is currently not shared, simplify() removes
	 * all false literals from the array.
	 */
	uint32 simplify(Solver& s);

	void            release();
  SharedLiterals* share();
	bool            unique()   const { return refCount_ <= 1; }
	uint32          refCount() const { return refCount_; }
private:
	void destroy();
	SharedLiterals(const Literal* lits, uint32 size, ConstraintType t, uint32 numRefs);
	SharedLiterals(const SharedLiterals&);
	SharedLiterals& operator=(const SharedLiterals&);
	std::atomic<int32> refCount_;
	uint32             size_type_;
	Literal            lits_[0];
};

//! A helper-class for creating/adding clauses.
/*!
 * \ingroup constraint
 * This class simplifies clause creation. It hides the special handling of 
 * short, and shared clauses. It also makes sure that learnt clauses watch 
 * the literals from the highest decision levels.
 */
class ClauseCreator {
public:
	//! Creates a new ClauseCreator object
	/*!
	 * \param s the Solver in which to store created clauses.
	 */
	explicit ClauseCreator(Solver* s = 0);
	//! sets the solver in which created clauses are stored.
	void setSolver(Solver& s) { solver_ = &s; }
	//! reserve space for a clause of size s
	void reserve(LitVec::size_type s) { literals_.reserve(s); }
	//! sets the initial activity for the to be created clause
	void setActivity(uint32 a)        { extra_.setActivity(a);  }
	//! sets the initial literal block distance for the to be created clause
	void setLbd(uint32 l1, uint32 l2) { extra_.setLbd(l1, l2); }
	//! discards the current clause
	void clear() { literals_.clear(); }
	
	//! starts the creation of a new clause
	/*!
	 * \pre s.decisionLevel() == 0 || t != Constraint_t::static_constraint
	 */
	ClauseCreator& start(ConstraintType t = Constraint_t::static_constraint);
	
	//! start the creation of an asserting or unit clause.
	/*!
   * All literals except the one passed to startAsserting() must be
	 * false and at least one of the false literals must be assigned on
	 * the current decision level.
	 * \param t The type of the learnt clause: either conflict or loop.
	 * \param assertingLit the literal to be asserted.
	 * \note position 0 is always reserved for the asserting literal.
	 */
	ClauseCreator& startAsserting(ConstraintType t, const Literal& assertingLit = Literal());
	
	//! adds the literal p to the clause.
	/*!
	 * \note p is only added if p is free or the current DL is > 0
	 * \pre clause neither contains p nor ~p.
	 * \pre if clause was started with startAsserting() p is false!
	 */
	ClauseCreator& add(const Literal& p);

	//! returns the literal at position i.
	/*!
	 * \pre i < size()
	 */
	Literal& operator[](LitVec::size_type i)       { return literals_[i]; }
	Literal  operator[](LitVec::size_type i) const { return literals_[i]; }
	//! returns the literals of the clause
	const LitVec&     lits()    const { return literals_; }
	LitVec&           lits()          { return literals_; }
	//! returns the current size of the clause.
	LitVec::size_type size()    const { return literals_.size(); }
	bool              empty()   const { return literals_.empty();}
	//! returns the clause's type.
	ConstraintType    type()    const { return extra_.type(); }
	//! returns the current position of the second watch
	uint32            sw()      const { return 1; }
		
	//! Possible status values of a clause
	enum Status {
		status_subsumed, status_sat, status_empty, status_conflicting, status_unit, status_open
	};

	//! returns the current status of the clause
	/*!
	 * 
	 * \return Given a clause with literals [l1,...,ln] returns
	 *  - status_subsumed   : if exists li, s.isTrue(li) && s.level(li) == 0
	 *  - status_sat        : if exists li, s.isTrue(li)
	 *  - status_empty      : if n == 0 or for all li, (s.isFalse(li) && s.level(li) == 0)
	 *  - status_conflicting: if n != 0 and for all li, s.isFalse(li) and for some lj, s.level(lj) > 0
	 *  - status_unit       : if s.isFree(l1) and for all i > 1, s.isFalse(li)
	 *  - status_open       : else
	 */
	Status status() const;

	//! returns the implication level of the unit-resulting literal
	/*! 
	 * \pre  status() == status_unit
	 */
	uint32 implicationLevel() const { return swLev_; }
	//! returns the (numerical) highest DL on which the clause is conflicting
	/*!
	 * \pre status() == status_conflicting || status_empty
	 */
	uint32 conflictLevel()    const { return fwLev_; }

	//! removes duplicate lits. Marks the clause as sat if it contains p and ~p
	void simplify();
	
	//! A type for storing the result of a clause insertion operation
	/*!
	 * The insertion operations return a tuple r s.th.
	 *   - r.ok is true if the clause was successfully created/integrated.
	 *   - r.ok is false if the clause is conflicting w.r.t the current assignment.
	 *   - r.local is the newly created local representation of the clause
	 */
	struct Result {
		explicit Result(ClauseHead* loc = 0, bool o = true)
			: local(loc)
			, ok(o) {}
		ClauseHead*     local;
		bool            ok;
		operator bool() const { return ok; }
	};

	//! adds the clause to the solver if possible
	/*!
	 * Adds the clause only if it is not conflicting.
	 * \note If the clause to be added is empty, end() fails and s.hasConflict() is set to true.
	 */
	Result end();

	//! reorders lits such that the first two literals are valid watches
	/*!
	 * \note 
	 *   If allFree is true, all literals in lits shall be free w.r.t
	 *   to s. In that case, the watches are picked according to
	 *   the strategy stored in s.
	 */
	static void initWatches(Solver& s, LitVec& lits, bool allFree);

	/*!
	 * \name factory functions
	 * Functions for creating and integrating clauses
	 */
		
	//! Flags controlling clause creation and integration
	enum CreateFlag {
		clause_no_add        = 1, /**< do not add clause to solver db */
		clause_explicit      = 2, /**< force creation of explicit clause even if size <= 3 */
		clause_not_sat       = 4, /**< integrate only if clause is not satisfied w.r.t current assignment */
		clause_not_root_sat  = 8, /**< integrate only if clause is not satisfied w.r.t root level */
		clause_no_release    = 16,/**< do not call release on shared literals */
		clause_known_order   = 32,/**< assume clause is already ordered w.r.t watches */
	};

	//! creates a clause from the literals given in lits
	/*!
	 * \param s     The solver to which the clause should be added
	 * \param lits  The literals of the clause
	 * \param flags flag set controlling creation (see CreateFlag)
	 * \param info  Initial information (e.g. type) for the new clause
	 *
	 * \pre !s.hasConflict() and s.decisionLevel() == 0 or extra.type() != Constraint_t::static_constraint
	 * \pre lits does not contain duplicate/complementary literals
	 * \pre lits is empty or lits[0] is either free or true
	 * \pre if lits[1] is false all lits except for lits[0] are false and no literal was assigned on a higher dl
	 *
	 * \note 
	 *   The function always assumes clause_known_order, i.e. it is the caller's job
	 *   to provide valid watches.
	 *
	 * \note 
	 *   If the given clause is unit (or asserting), the unit-resulting literal is
	 *   asserted on the (numerical) lowest level possible but the new information
	 *   is not immediately propagated, i.e. on return queueSize() may be greater than 0.
	 *   If the given clause is conflicting, the clause is stored as a conflict and 
	 *   resolveConflict() must be called in order to proceed.
	 *
	 * \note 
	 *   The local representation of the clause is always attached to the solver 
	 *   but only added to the solver if clause_no_add is not contained in modeFlags.
	 *   Otherwise, the returned clauses is owned by the caller
	 *   and it is the caller's responsibility to manage it. Furthermore, 
	 *   learnt statistics are *not* updated automatically in that case.
	 *
	 */
	static Result create(Solver& s, LitVec& lits, uint32 flags = 0, const ClauseInfo& info = ClauseInfo());
	
	//! integrates the given clause into the current search of s
	/*!
	 * \pre the assignment in s is not conflicting
	 * \param s           The solver in which the clause should be integrated
	 * \param clause      The clause to be integrated
	 * \param modeFlags   A set of flags controlling integration (see CreateFlag)
	 * \param e           Initial data for the local representation 
	 * 
	 * \note 
	 *   The function behaves similar to ClauseCreator::create() but it does not
	 *   assume clause_known_order, i.e. by default it determines the correct watches
	 *   automatically. Furthermore, the function does not add local representations for 
	 *   implicit clauses (i.e. size <= 3) unless modeFlags contains clause_explicit. 
	 *   In that case, an explicit representation is created. 
	 *   Implicit representations can only be created via ClauseCreator::createClause().
	 *
	 * \note
	 *   The function acts as a sink for the given clause (i.e. it decreases its reference count)
	 *   unless modeFlags contains clause_no_release.
	 *   
	 * \note integrate() is intended to be called in a post propagator. 
	 *   To integrate a set of clauses F, one would use a loop like this:
	 *   \code
	 *   bool MyPostProp::propagate(Solver& s) {
	 *     bool r = true;
	 *     while (!F.empty() && r) {
	 *       SharedLiterals* C = f.pop();
	 *       r = integrate(s, C, ...).ok;
	 *     }
	 *     return r;
	 *   \endcode
	 */
	static Result integrate(Solver& s, SharedLiterals* clause, uint32 modeFlags, const ClauseInfo& e);
	/*!
	 * \overload Result integrate(Solver& s, SharedLiterals* clause, uint32 modeFlags)
	 */
	static Result integrate(Solver& s, SharedLiterals* clause, uint32 modeFlags);
	//@}
private:
	void init(ConstraintType t);
	static inline Status status(const Solver& s, Literal* lits, uint32 flags);
	struct Watches;
	static ClauseHead* newProblemClause(Solver& s, LitVec& lits, const ClauseInfo& e, uint32 flags);
	static ClauseHead* newLearntClause(Solver& s,  LitVec& lits, const ClauseInfo& e, uint32 flags);
	Solver*         solver_;    // solver in which new clauses are stored
	LitVec          literals_;  // literals of the new clause
	ClauseInfo      extra_;     // extra info 
	uint32          fwLev_;     // decision level of the first watched literal
	uint32          swLev_;     // decision level of the second watched literal
	bool            asserting_; // true if clause was started with startAsserting
};

//! class for representing a clause in a solver
class Clause : public ClauseHead {
public:
	typedef Constraint::PropResult PropResult;
	enum { MAX_SHORT_LEN = 5 };
	
	/*!
	 * Creates a new problem clause from the literals in lits.
	 * 
	 * \param s Solver in which the new clause is to be used.
	 * \param lits The literals of the clause.
	 *
	 * \pre lits contains at least two literals and:
	 *  - all literals are currently free.
	 *  - lits does not contain any duplicate literals.
	 *  - if lits contains literal p it must not contain ~p.
	 *  . 
	 * \note The clause must be destroyed using Clause::destroy
	 */ 
	static ClauseHead*  newClause(Solver& s, const LitVec& lits) {
		return newClause(s, &lits[0], static_cast<uint32>(lits.size()), ClauseInfo());
	}
	static ClauseHead*  newClause(Solver& s, const Literal* lits, uint32 size, const ClauseInfo& extra);

	/*!
	 * Creates a new learnt clause from the literals in lits.
	 * 
	 * \param s     Solver in which the new clause is to be used.
	 * \param lits  The literals of the clause.
	 * \param extra Initial extra data of the learnt constraint.
	 *
	 * \pre lits contains at least two literals and:
	 *  - lits does not contain any duplicate/complementary literals.
	 *  - lits[0] and lits[1] are valid watches, i.e. 
	 *    - they are either both not false or
	 *    - lits[1] is false AND
	 *      - for each l != lits[0] in lits, isFalse(l) && level(l) <= level(lits[1])
	 *      - !isFalse(lits[0]) OR level(lits[0]) >= level(lits[1])
	 *  . 
	 * \note The clause must be destroyed using Clause::destroy
	 */ 
	static ClauseHead*  newLearntClause(Solver& s, const LitVec& lits, const ClauseInfo& extra) {
		return newClause(s, &lits[0], static_cast<uint32>(lits.size()), extra);
	}

	/*!
	 * Creates a new contracted clause from the literals contained in lits. A contracted clause
	 * consists of an active head and a (false) tail. Propagation is restricted to the head. 
	 * The tail is only needed to compute reasons from assignments.
	 * 
	 * \param s     Solver in which the new clause is to be used.
	 * \param lits  The literals of the clause.
	 * \param extra Initial extra data of the clause.
	 * \param tail  The starting index of the tail (first literal that should be temporarily removed from the clause).
	 * \param exten Extend head part of clause as tail literals become free?
	 */ 
	static ClauseHead*  newContractedClause(Solver& s, const LitVec& lits, const ClauseInfo& extra, LitVec::size_type tail, bool extend);

	//! creates a new specialized short clause object
	/*!
	 * \pre size <= MAX_SHORT_LEN
	 * \pre mem is either 0 or points to a memory block that was allocated via s.allocSmall();
	 */
	static ClauseHead*  newShortClause(Solver& s, void* mem, const Literal* lits, uint32 size, const ClauseInfo& extra);

	//! Creates a new local surrogate for shared_lits to be used in the given solver
	/*!
	 * \param s The solver in which this clause will be used
	 * \param shared_lits The shared literals of this clause
	 * \param e Initial data for the new (local) clause
	 * \param lits watches and cache literal for the new (local) clause
	 * \param addRef increment ref count of shared_lits
	 */
	static ClauseHead* newShared(Solver& s, SharedLiterals* shared_lits, const ClauseInfo& e, const Literal* lits, bool addRef = true);

	// Constraint-Interface
	
	Constraint* cloneAttach(Solver& other);

	/*!
	 * for a clause [x y p] the reason for p is ~x and ~y. 
	 * \pre *this previously asserted p
	 * \note if the clause is a learnt clause, calling reason increases
	 * the clause's activity.
	 */
	void reason(Solver& s, Literal p, LitVec& lits);

	bool minimize(Solver& m, Literal p, CCMinRecursive* r);

	bool isReverseReason(const Solver& s, Literal p, uint32 maxL, uint32 maxN);

	//! returns true if clause is SAT.
	/*!
	 * removes from the clause all literals that are false.
	 */
	bool simplify(Solver& s, bool = false);

	//! destroys the clause and frees its memory.
	void destroy(Solver* s = 0, bool detach = false);
	
	// LearntConstraint interface

	//! returns true if the clause is currently satisfied.
	/*!
	 * If the clause is currently not satisfied all free literals of
	 * the clause are copied to freeLits.
	 */
	bool isSatisfied(const Solver& s, LitVec& freeLits);
	
	// clause interface
	BoolPair strengthen(Solver& s, Literal p, bool allowToShort);
	void     detach(Solver&);
	uint32   size()                 const;
	void     toLits(LitVec& out)    const;
	bool     contracted()           const { return tail_.data.contracted(); }
	bool     isSmall()              const { return tail_.data.isSmall(); }
private:
	Clause(Solver& s, const ClauseInfo& extra, const Literal* lits, uint32 size, uint32 tail, bool extend);
	Clause(Solver& s, const Clause& other);
	static void* alloc(Solver& s, uint32 lits);
	void         undoLevel(Solver& s);
	bool         updateWatch(Solver& s, uint32 pos);
	typedef std::pair<Literal*, Literal*> LitRange;
	Literal*       longEnd()   { return head_+tail_.data.size(); }
	Literal*       smallEnd()  { return smallLits() + !isSentinel(*smallLits()) + !isSentinel(smallLits()[1]); }
	Literal*       smallLits() { return reinterpret_cast<Literal*>(&tail_); }
	const Literal* smallLits() const { return reinterpret_cast<const Literal*>(&tail_); }
	Literal*     removeFromTail(Solver& s, Literal* it, Literal* end);
	LitRange     tail() {
		return !isSmall()
			? LitRange(head_+ClauseHead::HEAD_LITS, longEnd())
			: LitRange(smallLits(), smallEnd());
	}
};

//! Constraint for Loop-Formulas
/*!
 * \ingroup constraint
 * Special purpose constraint for loop formulas of the form: R -> ~a1, ~a2, ..., ~an
 * where R is a conjunction (B1,...,Bm) of bodies that are false and a1...an are the atoms of
 * an unfounded set.
 * I.e. such a loop formula is equivalent to the following n clauses:
 * ~a1 v B1 v ... v Bm
 * ...
 * ~an v B1 v ... v Bm
 * Representing loop formulas as n clauses is wasteful because each clause
 * contains the same set of bodies. 
 * 
 * The idea behind LoopFormula is to treat the conjunction of atoms as a special
 * "macro-literal" L with the following properties:
 * - isTrue(L), iff for all ai isTrue(~ai) 
 * - isFalse(L), iff for some ai isFalse(~ai) 
 * - L is watchable, iff not isFalse(L)
 * - Watching L means watching all ai.
 * - setting L to true means setting all ai to false.
 * Using this convention the TWL-algo can be implemented as in a clause.
 * 
 * \par Implementation:
 * - The literal-array is divided into two parts, an "active clause" part and an atom part
 * - The "active clause" contains one atom and all bodies: [B1 ... Bj ~ai]
 * - The atom part contains all atoms: [~a1 ... ~an]
 * - Two of the literals of the "active clause" are watched (again: watching an atom means watching all atoms)
 * - If a watched atom becomes true, it is copied into the "active clause" and the TWL-algo starts.
 */
class LoopFormula : public LearntConstraint {
public:
	/*!
	 * Creates a new loop-formula for numAtoms atoms sharing the literals contained in bodyLits.
	 * 
	 * \param s Solver in which the new loop-formula is to be used.
	 * \param bodyLits Pointer to an array of numBodies body-literals
	 * \param numBodies Number of body-literals in bodyLits
	 * \param numAtoms Number of atoms in the loop-formula.
	 *
	 * \pre all body-literals are currently false.
	 */ 
	static LoopFormula* newLoopFormula(Solver& s, Literal* bodyLits, uint32 numBodies, uint32 bodyToWatch, uint32 numAtoms) {
		void* mem = ::operator new( sizeof(LoopFormula) + (numBodies+numAtoms+3) * sizeof(Literal) );
		return new (mem) LoopFormula(s, numBodies+numAtoms, bodyLits, numBodies, bodyToWatch);
	}

	//! Adds an atom to the loop-formula
	/*!
	 * \pre the loop-formula currently contains fewer than numAtoms atoms
	 */
	void addAtom(Literal atom, Solver& s);
	
	//! notifies the installed heuristic about the new constraint.
	void updateHeuristic(Solver& s);
	
	//! returns the size of the loop-formula.
	uint32 size() const;
	
	// Constraint interface
	Constraint* cloneAttach(Solver&) { return 0; }
	PropResult  propagate(Solver& s, Literal p, uint32& data);
	void reason(Solver&, Literal p, LitVec& lits);
	bool minimize(Solver& s, Literal p, CCMinRecursive* ccMin);
	bool simplify(Solver& s, bool = false);
	void destroy(Solver* = 0, bool = false);
	
	// LearntConstraint interface
	bool locked(const Solver& s) const;
	
	//! returns true if the loop-formula is satisfied
	/*!
	 * A loop-formula is satisfied if one body is true or all contained atoms are false.
	 * If the loop-formula is currently not satisfied isSatisfied copies all free literals
	 * to freeLits.
	 */
	bool isSatisfied(const Solver& s, LitVec& freeLits);
	
	//! returns the loop-formula's activity
	/*!
	 * The activity of a loop-formula is increased, whenever reason() is called.
	 */
	Activity activity() const { return makeActivity(activity_, ClauseHead::MAX_LBD); }
	
	//! halves the loop-formula's activity
	void decreaseActivity() { activity_ >>= 1; }
	
	//! returns Constraint_t::learnt_loop
	ConstraintType type() const { return Constraint_t::learnt_loop; }
private:
	LoopFormula(Solver& s, uint32 size, Literal* bodyLits, uint32 numBodies, uint32 bodyToWatch);
	bool watchable(const Solver& s, uint32 idx);
	bool isTrue(const Solver& s, uint32 idx);
	uint32  activity_; // activity of this loop formula
	uint32  end_;      // position of second sentinel
	uint32  size_;     // size of lits_
	uint32  other_;    // stores the position of a literal that was recently true
	Literal lits_[0];  // S B1...Bm ai S a1...an
};

namespace mt {

//! Stores the local part of a shared clause
/*!
 * The local part of a shared clause consists of a
 * clause head and and a pointer to the shared literals.
 * Since the local part is owned by a particular solver 
 * it can be safely modified. Destroying a SharedLitsClause 
 * means destroying the local part and decreasing the
 * shared literals' reference count.
 */
class SharedLitsClause : public ClauseHead {
public:
	//! Creates a new SharedLitsClause to be used in the given solver
	/*!
	 * \param s The solver in which this clause will be used
	 * \param shared_lits The shared literals of this clause
	 * \param e Initial data for the new (local) clause
	 * \param lits watches and cache literal for the new (local) clause
	 * \param addRef increment ref count of shared_lits
	 */
	static ClauseHead* newClause(Solver& s, SharedLiterals* shared_lits, const ClauseInfo& e, const Literal* lits, bool addRef = true);
	
	Constraint*    cloneAttach(Solver& other);
	void           reason(Solver& s, Literal p, LitVec& out);
	bool           minimize(Solver& s, Literal p, CCMinRecursive* rec);
	bool           isReverseReason(const Solver& s, Literal p, uint32 maxL, uint32 maxN);
	bool           simplify(Solver& s, bool);
	void           destroy(Solver* s, bool detach);
	bool           isSatisfied(const Solver& s, LitVec& freeLits);
	uint32         size() const;
	void           toLits(LitVec& out) const;
private:
	SharedLitsClause(Solver& s, SharedLiterals* x, const Literal* lits, const ClauseInfo&,  bool addRef);
	bool     updateWatch(Solver& s, uint32 pos);
	BoolPair strengthen(Solver& s, Literal p, bool allowToShort);
};
}

}
#endif
