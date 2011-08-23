// 
// Copyright (c) 2010-2011, Benjamin Kaufmann
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
#ifndef CLASP_SHARED_CONTEXT_H_INCLUDED
#define CLASP_SHARED_CONTEXT_H_INCLUDED
#ifdef _MSC_VER
#pragma warning (disable : 4200) // nonstandard extension used : zero-sized array
#pragma once
#endif

#include <clasp/literal.h>
#include <clasp/constraint.h>
#include <clasp/util/left_right_sequence.h>
#include <clasp/util/misc_types.h>

/*!
 * \file 
 * Contains some types shared between different solvers
 */
namespace Clasp {
class Solver;
class ClauseInfo;	
class Assignment;
class SharedContext;
class Enumerator;
class SharedLiterals;

/*!
 * \addtogroup solver
 */
//@{
//! Base class for preprocessors working on clauses only
class SatPreprocessor {
public:
	SatPreprocessor() : ctx_(0) {}
	virtual ~SatPreprocessor();
	void setContext(SharedContext& ctx) { ctx_ = &ctx; }
	virtual bool addClause(const LitVec& cl) = 0;
	virtual bool preprocess(bool enumerateModels) = 0;
	virtual void extendModel(Assignment& m, LitVec& open) = 0;
	virtual bool limit(uint32 numCons) const = 0;
	struct Stats {
		Stats() : clRemoved(0), clAdded(0), litsRemoved(0) {}
		uint32 clRemoved;
		uint32 clAdded;
		uint32 litsRemoved;
	} stats;
protected:
	void reportProgress(char t, uint32 min, uint32 max);
	SharedContext*  ctx_;
private:
	SatPreprocessor(const SatPreprocessor&);
	SatPreprocessor& operator=(const SatPreprocessor&);
};
//@}

/**
 * \defgroup shared classes to be shared between solvers
 */
//@{

///////////////////////////////////////////////////////////////////////////////
// Problem statistics
///////////////////////////////////////////////////////////////////////////////
//! A struct for aggregating basic problem statistics
/*!
 * Maintained in SharedContext
 */
struct ProblemStats {
	ProblemStats() { reset(); }
	uint32  vars;
	uint32  vars_eliminated;
	uint32  vars_frozen;
	uint32  constraints;
	uint32  constraints_binary;
	uint32  constraints_ternary;
	void    reset() { std::memset(this, 0, sizeof(*this)); }
	void diff(const ProblemStats& o) {
		vars               = std::max(vars, o.vars)-std::min(vars, o.vars);
		vars_eliminated    = std::max(vars_eliminated, o.vars_eliminated)-std::min(vars_eliminated, o.vars_eliminated);
		vars_frozen        = std::max(vars_frozen, o.vars_frozen)-std::min(vars_frozen, o.vars_frozen);
		constraints        = std::max(constraints, o.constraints) - std::min(constraints, o.constraints);
		constraints_binary = std::max(constraints_binary, o.constraints_binary) - std::min(constraints_binary, o.constraints_binary);
		constraints_ternary= std::max(constraints_ternary, o.constraints_ternary) - std::min(constraints_ternary, o.constraints_ternary);
	} 
};

//! Stores static information about variables
class VarInfo {
public:
	enum FLAGS {
		RESERVED_1 = 0x1u, // reserved for future
		RESERVED_2 = 0x2u, // use
		NANT   = 0x4u, // if this var is an atom, is it in NAnt(P)
		PROJECT= 0x8u, // do we project on this var?
		BODY   = 0x10u,// is this var representing a body?
		EQ     = 0x20u,// is the var representing both a body and an atom?
		ELIM   = 0x40u,// is the variable eliminated?
		FROZEN = 0x80u // is the variable frozen?
	};
	VarInfo() {}
	void  reserve(uint32 maxSize) { info_.reserve(maxSize); }
	void  add(bool body) {
		uint8 m = (!body?0:flag(BODY));
		info_.push_back( m );
	}
	bool      empty()                const { return info_.empty(); }
	uint32    numVars()              const { return (uint32)info_.size(); }
	bool      isSet(Var v, FLAGS f)  const { return (info_[v] & flag(f)) != 0; }
	void      toggle(Var v, FLAGS f)       { info_[v] ^= flag(f); }
	void      clear() { info_.clear(); }
private:
	// Bit:   7     6   5   4    3    2   1     0
	//      frozen elim eq body proj nant reserved
	typedef PodVector<uint8>::type InfoVec;
	static uint8 flag(FLAGS x) { return uint8(x); }
	
	VarInfo(const VarInfo&);
	VarInfo& operator=(const VarInfo&);
	InfoVec info_;
};

//! A class for efficiently storing and propagating binary and ternary clauses 
class ShortImplicationsGraph {
public:
	ShortImplicationsGraph();
	~ShortImplicationsGraph();
	//! make room for nodes number of nodes
	void resize(uint32 nodes);
	//! adds the binary constraint (p, q) to the implication graph
	/*!
	 * \return true iff a new implication was added
	 */
	bool addBinary(Literal p, Literal q, bool learnt, bool shared);
	//! adds the ternary constraint (p, q, r) to the implication graph
	/*!
	 * \return true iff a new implication was added
	 */
	bool addTernary(Literal p, Literal q, Literal r, bool learnt, bool shared);
	
	//! removes p and its implications
	/*!
	 * \pre s.isTrue(p)
	 */
	void removeTrue(Solver& s, Literal p);
	
	//! propagates consequences of p following from binary and ternary clauses
	/*!
	 * \pre s.isTrue(p)
	 */
	bool   propagate(Solver& s, Literal p) const;
	//! propagates immediate consequences of p following from binary clauses only
	bool   propagateBin(Assignment& out, Literal p, uint32 dl) const;
	//! checks whether there is a reverse arc implying p and if so returns it in out
	bool   reverseArc(const Solver& s, Literal p, uint32 maxLev, Antecedent& out) const;
	
	uint32 numBinary() const { return bin_[0]; }
	uint32 numTernary()const { return tern_[0]; }
	uint32 numLearnt() const { return bin_[1] + tern_[1]; }
	uint32 numEdges(Literal p) const;
private:
	ShortImplicationsGraph(const ShortImplicationsGraph&);
	ShortImplicationsGraph& operator=(ShortImplicationsGraph&);
#ifndef DISABLE_MULTI_THREADING
	struct Block;
	typedef std::atomic<Block*> SharedBlockPtr;
	typedef bk_lib::left_right_sequence<Literal, std::pair<Literal,Literal>, 64-sizeof(SharedBlockPtr)> ImpListBase;
	struct ImplicationList : public ImpListBase {
		ImplicationList() : ImpListBase() { learnt = 0; }
		ImplicationList(const ImplicationList& other) : ImpListBase(other), learnt(other.learnt) {}
		~ImplicationList();
		bool hasLearnt(Literal q, Literal r = negLit(0)) const;
		void addLearnt(Literal q, Literal r = negLit(0));
		bool empty() const { return ImpListBase::empty() && learnt == 0; }
		void move(ImplicationList& other);
		void clear(bool b);
		SharedBlockPtr learnt; 
	};
#else
	typedef bk_lib::left_right_sequence<Literal, std::pair<Literal,Literal>, 64> ImplicationList;
#endif
	ImplicationList& getList(Literal p) { return graph_[p.index()]; }
	void remove_bin(ImplicationList& w, Literal p);
	void remove_tern(ImplicationList& w, Literal p);
	typedef PodVector<ImplicationList>::type ImpLists;
	ImpLists   graph_;     // one implication list for each literal
	uint32     bin_[2];    // number of binary constraints (0: problem / 1: learnt)
	uint32     tern_[2];   // number of ternary constraints(0: problem / 1: learnt)
};

//! Base class for distributing learnt knowledge between solvers
class Distributor {
public:
	Distributor(uint32 maxShare, uint32 typesToShare, uint32 maxLbd);
	virtual ~Distributor();
	SharedLiterals* publish(const Solver& source, const Literal* lits, uint32 size, const ClauseInfo& extra);
	virtual uint32  receive(const Solver& in, SharedLiterals** out, uint32 maxOut) = 0;
protected:
	virtual void    doPublish(const Solver& source, SharedLiterals* lits) = 0;
private:
	Distributor(const Distributor&);
	Distributor& operator=(const Distributor&);
	uint32 maxShare_;
	uint32 lbdMax_    : 16;
	uint32 typeMask_  : 16;
};

	
//! Aggregates information to be shared between solver objects
/*!
 * Among other things, SharedContext objects store 
 * static information on variables, the (possibly unused) 
 * symbol table, as well as the binary and ternary 
 * implication graph of the input problem.
 * 
 * Furthermore, a SharedContext object stores a distinguished
 * master solver that is used to store and simplify problem constraints.
 *
 * Once initialization is completed, other solvers s can 
 * attach to this object by calling ctx->attach(s).
 */
class SharedContext {
public:
	typedef std::auto_ptr<SatPreprocessor> SatPrepro;
	typedef ProblemStats                   Stats;
	typedef LitVec::size_type              size_type;
	typedef ShortImplicationsGraph         BTIG;
	
	enum InitMode { init_share_symbols };
	
	//! creates a new object for sharing variables and the binary and ternary implication graph.
	SharedContext();
	//! creates a new object that shares its symbol table with rhs
	SharedContext(const SharedContext& rhs,  InitMode m);
	
	~SharedContext();
	
	//! enables sharing of initial problem constraints
	/*!
	 * If this function is not called, problem constraints are
	 * cloned whenever a solver attaches to this object. Otherwise,
	 * they are shared.
	 */
	void enableConstraintSharing() {  shareConstr_ = true; }

	//! enables sharing of learnt constraints
	/*!
	 * If this function is not called, learnt constraints are
	 * not shared between different solvers. Otherwise,
	 * sharing is possible and controlled by a distribution 
	 * strategy.
	 */
	void enableLearntSharing(Distributor* d) { 
		distributor_.reset(d); 
	}

	//! resets this object to the state after default construction
	void reset();
	
	//! returns true if var represents a valid variable in this object
	/*!
	 * \note The range of valid variables is [1..numVars()]. The variable 0
	 * is a special sentinel variable. 
	 */
	bool validVar(Var var) const { return var <= numVars(); }

	//! returns the number of problem variables.
	/*!
	 * \note The special sentinel-var 0 is not counted, i.e. numVars() returns
	 * the number of problem-variables.
	 * To iterate over all problem variables use a loop like:
	 * \code
	 * for (Var i = 1; i <= numVars(); ++i) {...}
	 * \endcode
	 */
	uint32 numVars() const { return varInfo_.numVars() - 1; }

	//! returns the number of eliminated vars
	uint32 numEliminatedVars()const { return problem_.vars_eliminated; }

	//! reserve space for at least varGuess variables 
	/*!
	 * \param varGuess number of vars to reserve space for
	 * \note if the number of variables is known upfront, passing the correct value
	 * for varGuess avoids repeated regrowing of the state data structures
	 */
	void reserveVars(uint32 varGuess);

	//! adds a new variable of type t.
	/*!
	 * \param t  type of the new variable (either Var_t::atom_var or Var_t::body_var)
	 * \param eq true if var represents both an atom and a body. In that case
	 *           t is the variable's primary type and determines the preferred literal.
	 * \return The index of the new variable
	 * \note The problem variables are numbered from 1 onwards!
	 */
	Var addVar(VarType t, bool eq = false);

	//! requests a special tag literal for tagging conditional knowledge
	/*!
	 * Once a tag literal p is set, newly learnt clauses containing ~p are
	 * tagged as "conditional". Conditional clauses can be removed from a solver
	 * by calling Solver::removeConditional(). Furthermore, calling 
	 * Solver::strengthenConditional() removes ~p from conditional clauses and
	 * transforms them to unconditional knowledge.
	 *
	 * \note Typically, the tag literal is an initial assumption and hence true during 
	 *       the whole search. 
	 */
	void    requestTagLiteral();
	Literal tagLiteral() const   { return tag_; }
	void    removeTagLiteral();

	//! request additional reason data slot for variable v
	void    requestData(Var v);

	//! returns the type of variable v.
	/*!
	 * If v was added with parameter eq=true, the return value
	 * is Var_t::atom_body_var.
	 */
	VarType type(Var v) const {
		assert(validVar(v));
		return varInfo_.isSet(v, VarInfo::EQ)
			? Var_t::atom_body_var
			: VarType(Var_t::atom_var + varInfo_.isSet(v, VarInfo::BODY));
	}

	//! returns the preferred decision literal of variable v w.r.t its type
	/*!
	 * \return 
	 *  - posLit(v) if type(v) == body_var
	 *  - negLit(v) if type(v) == atom_var
	 * \note if type(v) is atom_body_var, the preferred literal is determined
	 *       by v's primary type, i.e. the one that was initially passed to addVar().
	 */
	Literal preferredLiteralByType(Var v) const {
		assert(validVar(v));
		return Literal(v, !varInfo_.isSet(v, VarInfo::BODY));
	}

	//! returns true if v is currently eliminated, i.e. no longer part of the problem
	bool eliminated(Var v)  const     { assert(validVar(v)); return varInfo_.isSet(v, VarInfo::ELIM); }
	//! returns true if v is excluded from variable elimination
	bool frozen(Var v)      const     { assert(validVar(v)); return varInfo_.isSet(v, VarInfo::FROZEN); }
	//! returns true if v is a projection variable
	bool project(Var v)     const     { assert(validVar(v)); return varInfo_.isSet(v, VarInfo::PROJECT);}
	//! returns true if v is contained in a negative loop or head of a choice rule
	bool nant(Var v)        const     { assert(validVar(v)); return varInfo_.isSet(v, VarInfo::NANT);}
	
	//! freezes/defreezes a variable (a frozen var is exempt from SatELite preprocessing)
	void setFrozen(Var v, bool b);
	void setProject(Var v, bool b)    { assert(validVar(v)); if (b != varInfo_.isSet(v, VarInfo::PROJECT)) varInfo_.toggle(v, VarInfo::PROJECT); }
	void setNant(Var v, bool b)       { assert(validVar(v)); if (b != varInfo_.isSet(v, VarInfo::NANT))    varInfo_.toggle(v, VarInfo::NANT);    }
	/*!
	 * \name problem specification
	 * Functions for adding a problem to the master solver.
	 * Problem specification is a four-stage process:
	 * - first, add variables to the SharedContext object 
	 * - second, call startAddConstraints()
	 * - third, add problem constraints
	 * - finally, endInit() shall be called to finished the initialization process
	 * .
	 * \note After endInit() was called, other solvers can be attached to this object
	 * \note In incremental setting, the process must be repeated for each incremental step.
	 * 
	 * \note Problem specification is *not* thread-safe, i.e. during initialization no other thread shall
	 * access the context.
	 *
	 * \note unique() is a precondition for all functions in this group!
	 */
	//@{

	bool unique() const   { return shareCount_ == 1; } 

	//! starts an initialization phase
	/*!
	 * Must be called once before constraints can be added.
	 */
	Solver& startAddConstraints(uint32 constraintGuess = 100);

	//! returns the number of problem constraints
	uint32 numConstraints()   const;

	//! eliminates the variable v
	/*!
	 * \pre v must not occur in any constraint 
	 *  and frozen(v) == false and value(v) == value_free
	 */
	void eliminate(Var v);

	//! adds the constraint c to the master solver
	/*!
	 * \pre endInit() was not called.
	 */
	void add(Constraint* c);

	//! adds the unary constraint p to the master solver.
	/*!
	 * \note unary constraints are immediately asserted.
	 * \return false if asserting p leads to a conflict.
	 */
	bool addUnary(Literal p);

	//! adds the binary constraint (p, q) to the shared implication graph
	void addBinary(Literal p, Literal q) { btig_.addBinary(p, q, false, !unique()); }
	
	//! adds the ternary constraint (p, q, r) to the shared implication graph
	void addTernary(Literal p, Literal q, Literal r) { btig_.addTernary(p, q, r, false, !unique()); }

	//! adds p as post propagator to the master solver
	/*!
	 * \pre p was not added previously and is not part of any other solver
	 * \note post propagators are stored in priority order
	 * \see PostPropagator::priority()
	 */
	void addPost(PostPropagator* p);

	//! attaches the given enumerator to this object
	/*!
	 * \note ownership is transferred
	 * \note In incremental setting, the enumerator must be reattached in
	 *       each incremental step by calling addEnumerator(enumerator());
	 */
	void        addEnumerator(Enumerator* en);
	Enumerator* enumerator() const { return enumerator_.get(); }

	//! finishes initialization of the master solver
	/*!
	 * endInit must be called once before search is started. After endInit()
	 * was called, a number (shareCount-1) of other solvers can be attached to the 
	 * shared context and learnt constraints may be added to solver.
	 * \return 
	 *  - false if the constraints are initially conflicting. True otherwise.
	 * \note
	 * The master solver can't recover from top-level conflicts, i.e. if endInit()
	 * returned false, the solver is in an unusable state.
	 */
	bool endInit(uint32 shareCount = 1);

	//! attaches s to this object
	/*!
	 * \pre other != master()
	 * \note It is safe to attach multiple solvers concurrently
	 * but the master solver shall not change during the whole
	 * operation.
	 */
	bool attach(Solver& other);

	//! detaches s from this object
	/*!
	 * The function removes any enumeration related constraints from s.
	 * Shall be called once after search has stopped.
	 * \note The function is concurrency-safe w.r.t to different solver objects, 
	 *       i.e. in a parallel search different solvers may call detach()
	 *       concurrently.
	 */
	void detach(Solver& s);

	//! estimates the problem complexity
	/*!
	 * \return sum of c->estimateComplexity(*master()) for each problem 
	 *         constraint c.
	 */
	uint32 problemComplexity() const;

	//! size of top-level after last call to endInit()
	size_type topLevelSize() const { return lastTopLevel_; }
	//@}

	/*!
	 * \name learning
	 * Functions for distributing knowledge
	 * 
	 * \note The functions in this group can be safely called 
	 * from multiple threads.
	 */
	//@{

	//! learns the binary clause (p, q)
	/*!
	 * \note 
	 *   The binary clause (p, q) is stored in a shared data-structure but
	 *   threads are not informed about the new clause.
	 *   It is the caller's responsibility to distribute the new information.
	 */
	bool learnBinary(Literal p, Literal q) { return btig_.addBinary(p, q, true, !unique()); }
	
	//! learns the ternary clause (p, q, r)
	/*!
	 * \note 
	 *   The ternary clause (p, q, r) is stored in a shared data-structure but
	 *   threads are not informed about the new clause.
	 *   It is the caller's responsibility to distribute the new information.
	 */
	bool learnTernary(Literal p, Literal q, Literal r) { return btig_.addTernary(p, q, r, true, !unique()); }
	
	//! distributes the clause in lits to other threads
	/*!
	 * The function first calls the distribution strategy 
	 * to decides whether the clause is a good candidate for distribution.
	 * If so, it distributes the clause and returns a handle to the
	 * now shared literals of the clause. Otherwise, it returns 0.
	 *
	 * \param owner The solver that created the clause.
	 * \param lits  The literals of the clause.
	 * \param size  The number of literals in the clause
	 * \param extra Additional information about the clause
	 * \note 
	 *   If the return value is not null, it is the caller's 
	 *   responsibility to release the returned handle (i.e. by calling release()).
	 */
	SharedLiterals* distribute(const Solver& owner, const Literal* lits, uint32 size, const ClauseInfo& extra) const;
	void            distribute(const Solver& owner, Literal p, Literal q, Literal r, ConstraintType type) const;

	//! tries to receive at most maxOut clauses
	/*!
	 * The function queries the distribution object for new clauses to be delivered to
	 * the solver target. Clauses are stored in out.
	 * \return the number of clauses received
	 */
	uint32 receive(const Solver& target, SharedLiterals** out, uint32 maxOut) const;

	//! returns the number of learnt binary and ternary clauses
	uint32       numLearntShort() const { return btig_.numLearnt(); }

	//@}
	//! returns the master solver associated with 
	Solver*      master()    const   { return master_; }	
	bool         shareConstraints() const { return shareConstr_; }
	const Stats& stats()     const   { return problem_; }
	uint32       numBinary() const   { return btig_.numBinary();  }
	uint32       numTernary()const   { return btig_.numTernary(); }
	SymbolTable& symTab()    const   { return symTabPtr_->symTab; }
	void         simplifyBtig(Literal p) { btig_.removeTrue(*master(), p); }
	const BTIG&  shortImplications() const { return btig_; }
	uint32       numShortImplications(Literal p) const { return btig_.numEdges(p); }
	SatPrepro    satPrepro;  // preprocessor
private:
	SharedContext(const SharedContext&);
	SharedContext& operator=(const SharedContext&);
	typedef std::auto_ptr<Distributor> DistPtr;
	typedef std::auto_ptr<Enumerator>  EnumPtr;
	struct SharedSymTab {
		SharedSymTab() : refs(1) {}
		SymbolTable symTab;
		uint32      refs;
	}*           symTabPtr_;   // pointer to shared symbol table
	ProblemStats problem_;     // problem statistics1
	VarInfo      varInfo_;     // varInfo[v] stores info about variable v
	BTIG         btig_;        // binary-/ternary implication graph
	Solver*      master_;      // master solver, responsible for init
	EnumPtr      enumerator_;  // enumerator object
	DistPtr      distributor_; // object for distributing learnt knowledge
	Literal      tag_;         // literal for tagging learnt constraints
	size_type    lastInit_;    // size of master's db after last init
	size_type    lastTopLevel_;// size of master's top-level after last init
	uint32       shareCount_;  // number of objects sharing this object
	bool         shareConstr_; // sharing of problem constraints enabled?
};

//@}
}
#endif
