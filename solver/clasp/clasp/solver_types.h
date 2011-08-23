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
#ifndef CLASP_SOLVER_TYPES_H_INCLUDED
#define CLASP_SOLVER_TYPES_H_INCLUDED
#ifdef _MSC_VER
#pragma once
#endif

#include <clasp/literal.h>
#include <clasp/constraint.h>
#include <clasp/util/left_right_sequence.h>
#include <clasp/util/misc_types.h>
#include <clasp/util/type_manip.h>

/*!
 * \file 
 * Contains some types used by a Solver
 */
namespace Clasp {
class SharedLiterals;

/*!
 * \addtogroup solver
 */
//@{

///////////////////////////////////////////////////////////////////////////////
// Statistics
///////////////////////////////////////////////////////////////////////////////
//! A struct for aggregating statistics relevant for parallel solving
/*!
 * Always associated with one solver (thread)
 */
struct ParallelStats {
	ParallelStats() { reset(); }
	double  cpuTime; /**< (Estimated) cpu time of the current solver */
	uint64  splits;  /**< Number of split requests handled */
	uint64  gps;     /**< Number of guiding paths received */
	uint64  gpLits;  /**< Sum of literals in received guiding paths */
	uint64  shared;  /**< Number of nogoods shared */
	bool    terminated;
	void    reset() { 
		std::memset(this, 0, sizeof(*this)); 
		cpuTime    = 0.0;
		terminated = false;
	}
	void    accu(const ParallelStats& o) {
		cpuTime += o.cpuTime;
		splits  += o.splits;
		gps     += o.gps;
		gpLits  += o.gpLits;
		shared  += o.shared;
		terminated |= o.terminated;
	}

	void newGP(LitVec::size_type length) {
		++gps;
		gpLits += length;
	}
};

//! A struct for jump statistics
struct JumpStats {
	JumpStats() { reset(); }
	void reset(){ std::memset(this, 0, sizeof(*this)); }
	void accu(const JumpStats& o) {
		jumps   += o.jumps;
		bJumps  += o.bJumps;
		jumpSum += o.jumpSum;
		boundSum+= o.boundSum;
		if (o.maxJump   > maxJump)   maxJump   = o.maxJump;
		if (o.maxJumpEx > maxJumpEx) maxJumpEx = o.maxJumpEx;
		if (o.maxBound  > maxBound)  maxBound  = o.maxBound;
	}
	void    update(uint32 dl, uint32 uipLevel, uint32 bLevel) {
		++jumps;
		jumpSum += dl - uipLevel; 
		maxJump = std::max(maxJump, dl - uipLevel);
		if (uipLevel < bLevel) {
			++bJumps;
			boundSum += bLevel - uipLevel;
			maxJumpEx = std::max(maxJumpEx, dl - bLevel);
			maxBound  = std::max(maxBound, bLevel - uipLevel);
		}
		else { maxJumpEx = maxJump; }
	}
	double avgJumpLen()   const { return jumpSum/std::max(1.0,(double)jumps); }
	double avgJumpLenEx() const { return (jumpSum-boundSum)/std::max(1.0,(double)jumps); }
	uint64  jumps;    /**< Number of backjumps (i.e. number of analyzed conflicts) */
	uint64  bJumps;   /**< Number of backjumps that were bounded */
	uint64  jumpSum;  /**< Number of levels that could be skipped w.r.t first-uip */
	uint64  boundSum; /**< Number of levels that could not be skipped because of backtrack-level*/
	uint32  maxJump;  /**< Longest possible backjump */
	uint32  maxJumpEx;/**< Longest executed backjump (< maxJump if longest jump was bounded) */
	uint32  maxBound; /**< Max difference between uip- and backtrack-level */
};

//! A struct for holding core solving statistics used by a solver
struct CoreStats {
	CoreStats() { reset(); }
	void reset(){ std::memset(this, 0, sizeof(*this)); }
	void accu(const CoreStats& o) {
		choices  += o.choices;
		conflicts+= o.conflicts;
		restarts += o.restarts;
		models   += o.models;
		
		binary   += o.binary;
		ternary  += o.ternary;
		deleted  += o.deleted;
		modLits  += o.modLits;
		for (int i = 0; i != Constraint_t::max_value; ++i) {
			learnts[i] += o.learnts[i];
			lits[i]    += o.lits[i];
		}
	}
	void addLearnt(uint32 size, ConstraintType t) {
		assert(t != Constraint_t::static_constraint && t <= Constraint_t::max_value);
		++learnts[t-1];
		lits[t-1] += size;
		binary += (size == 2);
		ternary+= (size == 3);
	}
	void removeLearnt(uint32 size, ConstraintType t) {
		assert(t != Constraint_t::static_constraint && t <= Constraint_t::max_value);
		--learnts[t-1];
		lits[t-1] -= size;
		binary -= (size == 2);
		ternary-= (size == 3);
	}
	void addModel(uint32 size) {
		++models;
		modLits += size;
	}
	void removeModel(uint32 size) {
		--models;
		modLits -= size;
	}
	uint64 choices;   /**< Number of choices performed */
	uint64 conflicts; /**< Number of conflicts found */
	uint64 restarts;  /**< Number of restarts */ 
	uint64 models;    /**< Number of models found */

	uint64 learnts[Constraint_t::max_value]; /**< Number of learnt nogoods of type t-1 */
	uint64 lits[Constraint_t::max_value];    /**< Sum of literals in nogoods of type t-1 */
	uint64 binary;    /**< Number of learnt binary nogoods */
	uint64 ternary;   /**< Number of learnt ternary nogoods*/
	uint64 deleted;   /**< Sum of learnt nogoods removed */
	uint64 modLits;   /**< Sum of decision literals in models */
};

//! A struct for aggregating statistics of one solve operation
struct SolveStats : public CoreStats {
	SolveStats() : jumps(0), parallel(0) { }
	SolveStats(const SolveStats& o) : CoreStats(o), jumps(0), parallel(0) { 
		if (o.jumps)    jumps    = new JumpStats(*o.jumps); 
		if (o.parallel) parallel = new ParallelStats(*o.parallel);
	}
	~SolveStats() { delete jumps; delete parallel; }
	void enableJumpStats()    { if (!jumps)   jumps    = new JumpStats();   }
	void enableParallelStats(){ if (!parallel)parallel = new ParallelStats(); }
	void reset() {  
		CoreStats::reset();
		if (jumps)   jumps->reset();
		if (parallel)parallel->reset();
	}
	void accu(const SolveStats& o) {
		CoreStats::accu(o);
		if (jumps && o.jumps) jumps->accu(*o.jumps);
		if (parallel && o.parallel) parallel->accu(*o.parallel);
	}
	void updateJumps(uint32 dl, uint32 uipLevel, uint32 bLevel) {
		if (!jumps) return;
		jumps->update(dl, uipLevel, bLevel);
	}
	JumpStats*     jumps;    /**< optional jump statistics */
	ParallelStats* parallel; /**< optional parallel statistics */
private:
	SolveStats& operator=(const SolveStats&);
};
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Clauses
///////////////////////////////////////////////////////////////////////////////
//! Type for passing additional clause information 
class ClauseInfo {
public:
	enum { MAX_LBD = (1<<7)-1, MAX_ACTIVITY = (1<<15)-1 }; 
	ClauseInfo() 
		: act_(0), lbd_(MAX_LBD), rootLbd_(MAX_LBD), type_(0), tag_(0) {}
	ConstraintType type()     const { return static_cast<ConstraintType>(type_); }
	uint32         activity() const { return static_cast<uint32>(act_); }
	uint32         lbd()      const { return static_cast<uint32>(lbd_); }
	uint32         rootLbd()  const { return static_cast<uint32>(rootLbd_); }
	bool           tagged()   const { return tag_ != 0; }
	ClauseInfo& setType(ConstraintType t) { type_ = static_cast<uint32>(t); return *this; }
	ClauseInfo& setActivity(uint32 act)   {
		if (act > MAX_ACTIVITY) act = MAX_ACTIVITY;
		act_ = act;
		return *this;
	}
	ClauseInfo& setLbd(uint32 a_lbd, uint32 a_rootLbd) {
		assert(a_lbd >= a_rootLbd);
		if (a_lbd > static_cast<uint32>(MAX_LBD)) {
			a_lbd    = static_cast<uint32>(MAX_LBD);
			a_rootLbd= std::min(a_lbd, a_rootLbd);
		}
		lbd_     = a_lbd;
		rootLbd_ = a_rootLbd;
		return *this;
	}
	
	ClauseInfo& tag() { tag_ = 1; return *this; }
private:
	uint32 act_     : 15; // Initial clause activity
	uint32 lbd_     :  7; // Literal block distance in the range [0, MAX_LBD]
	uint32 rootLbd_ :  7; // LBD w.r.t to current root level ([0, MAX_LBD])
	uint32 type_    :  2; // Type of clause
	uint32 tag_     :  1; // conditional clause?
};

//! (Abstract) base class for clause types
/*!
 * ClauseHead is used to enforce a common memory-layout for all clauses.
 * It contains the two watched literals and a cache literal to improve
 * propagation performance. A virtual call to Constraint::propagate()
 * is only needed if the other watch is not true and the cache literal
 * is false.
 */
class ClauseHead : public LearntConstraint {
public:
	enum { HEAD_LITS = 3, MAX_LBD = (1<<5)-1, TAGGED_CLAUSE = 1023};
	explicit ClauseHead(const ClauseInfo& init);
	// base interface
	//! propagates the head and calls propagateTail() if necessary
	PropResult propagate(Solver& s, Literal, uint32& data);
	//! type of clause
	ConstraintType type() const    { return ConstraintType(info_.data.type); }
	//! true if this clause currently is the antecedent of an assignment
	bool     locked(const Solver& s) const;
	//! returns the activity of this clause
	Activity activity() const       { return makeActivity(info_.data.act, info_.data.lbd); }
	//! halves the activity of this clause
	void     decreaseActivity()     { info_.data.act >>= 1; }
	//! downcast from LearntConstraint
	ClauseHead* clause()           { return this; }
	
	// clause interface
	typedef std::pair<bool, bool> BoolPair;
	//! increase activity
	void bumpActivity()     { info_.data.act += (info_.data.act != ClauseInfo::MAX_ACTIVITY); }
	//! adds watches for first two literals in head to solver
	void attach(Solver& s);
	//! returns true if head is satisfied w.r.t current assignment in s
	bool satisfied(const Solver& s);
	//! conditional clause?
	bool tagged() const     { return info_.data.key == uint32(TAGGED_CLAUSE); }
	uint32 lbd()  const     { return info_.data.lbd; }
	void setLbd(uint32 lbd) { info_.data.lbd = std::min(lbd, (uint32)MAX_LBD); }
	//! removes watches from s
	virtual void     detach(Solver& s);
	//! returns the size of this clause
	virtual uint32   size()              const = 0;
	//! returns the literals of this clause in out
	virtual void     toLits(LitVec& out) const = 0;
	//! returns true if this clause is a valid "reverse antecedent" for p
	virtual bool     isReverseReason(const Solver& s, Literal p, uint32 maxL, uint32 maxN) = 0;
	//! removes p from clause if possible
	/*!
	 * \return
	 *   The first component of the returned pair specifies whether or not
	 *   p was removed from the clause.
	 *   The second component of the returned pair specifies whether
	 *   the clause should be kept (false) or removed (true). 
	 */
	virtual BoolPair strengthen(Solver& s, Literal p, bool allowToShort = true) = 0;
protected:
	friend struct ClauseWatch;
	bool         toImplication(Solver& s);
	void         clearTagged(){ info_.data.key = 0; }
	//! shall replace the watched literal at position pos with a non-false literal
	/*!
	 * \pre pos in [0,1] 
	 * \pre s.isFalse(head_[pos]) && s.isFalse(head_[2])
	 * \pre head_[pos^1] is the other watched literal
	 */
	virtual bool updateWatch(Solver& s, uint32 pos) = 0;
	union TailInfo {
		struct TInfo {
			uint32 sizeExt_; // or literal if small
			uint32 idx_;     // or literal if small
			void     init(bool small)  { sizeExt_ = idx_ = (small ? negLit(0).asUint() : 1); }
			bool     isSmall()   const { return (sizeExt_ & 1u) == 0;  }
			bool     contracted()const { return (sizeExt_ & 3u) == 3u; }
			uint32   size()      const { return sizeExt_ >> 2; }
			uint32   idx()       const { return idx_; }
			void     setSize(uint32 size) { sizeExt_ = (size << 2) | (sizeExt_ & 3u); }
			void     setIdx(uint32  idx)  { idx_ = idx; }
			void     markContracted()     { sizeExt_ |= 2u;  }
			void     clearContracted()    { sizeExt_ &= ~2u; }
		}                   data;
		SharedLiterals*   shared;
	}       tail_;
	union Info { 
		Info() : rep(0) {}
		explicit Info(const ClauseInfo& i);
		struct {
			uint32 act : 15; // activity of clause
			uint32 key : 10; // lru key of clause
			uint32 lbd :  5; // lbd of clause
			uint32 type:  2; // type of clause
		}      data;
		uint32 rep;
	}       info_;
	Literal head_[HEAD_LITS]; // two watched literals and one cache literal
};
//! Allocator for small (at most 32-byte) clauses
class SmallClauseAlloc {
public:
	SmallClauseAlloc();
	~SmallClauseAlloc();
	void* allocate() {
		if(freeList_ == 0) {
			allocBlock();
		}
		Chunk* r   = freeList_;
		freeList_  = r->next;
		return r;
	}
	void   free(void* mem) {
		Chunk* b = reinterpret_cast<Chunk*>(mem);
		b->next  = freeList_;
		freeList_= b;
	}
private:
	SmallClauseAlloc(const SmallClauseAlloc&);
	SmallClauseAlloc& operator=(const SmallClauseAlloc&);
	struct Chunk {
		Chunk*        next; // enforce ptr alignment
		unsigned char mem[32 - sizeof(Chunk*)];
	};
	struct Block {
		enum { num_chunks = 1023 };
		Block* next;
		unsigned char pad[32-sizeof(Block*)];
		Chunk  chunk[num_chunks];
	};
	void allocBlock();
	Block*  blocks_;
	Chunk*  freeList_;
};
///////////////////////////////////////////////////////////////////////////////
// Watches
///////////////////////////////////////////////////////////////////////////////
//! Represents a clause watch in a Solver.
struct ClauseWatch {
	//! clause watch: clause head
	explicit ClauseWatch(ClauseHead* a_head) : head(a_head) { }
	ClauseHead* head;
	struct EqHead {
		explicit EqHead(ClauseHead* h) : head(h) {}
		bool operator()(const ClauseWatch& w) const { return head == w.head; }
		ClauseHead* head;
	};
};

//! Represents a generic watch in a Solver.
struct GenericWatch {
	//! a constraint and some associated data
	explicit GenericWatch(Constraint* a_con, uint32 a_data = 0) : con(a_con), data(a_data) {}
	//! calls propagate on the stored constraint and passes the stored data to that constraint.
	Constraint::PropResult propagate(Solver& s, Literal p) { return con->propagate(s, p, data); }
	
	Constraint* con;    /**< The constraint watching a certain literal */
	uint32      data;   /**< Additional data associated with this watch - passed to constraint on update */

	struct EqConstraint {
		explicit EqConstraint(Constraint* c) : con(c) {}
		bool operator()(const GenericWatch& w) const { return con == w.con; }
		Constraint* con;
	};	
};

//! Watch list type
typedef bk_lib::left_right_sequence<ClauseWatch, GenericWatch, 0> WatchList;
inline void releaseVec(WatchList& w) { w.clear(true); }

///////////////////////////////////////////////////////////////////////////////
// Assignment
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//! Type for storing reasons for variable assignments together with additional data
/*!
 * \note On 32-bit systems additional data is stored in the high-word of antecedents 
 */
struct ReasonStore32 : PodVector<Antecedent>::type {
	uint32  dataSize() const     { return (uint32)size(); }
	void    dataResize(uint32)   {}
	uint32  data(uint32 v) const { return static_cast<uint32>((*this)[v].asUint()>>32);}
	void    setData(uint32 v, uint32 data) { encode((*this)[v], data); }
	static  void encode(Antecedent& a, uint32 data) {
		a.asUint() = (uint64(data)<<32) | static_cast<uint32>(a.asUint());
	}
	struct value_type {
		value_type(const Antecedent& a, uint32 d) : ante_(a) {
			if (d != UINT32_MAX) { encode(ante_, d); }
		}
		const Antecedent& ante() const { return ante_;      }
		      uint32      data() const { return UINT32_MAX; }
		Antecedent ante_;
	};
};

//! Type for storing reasons for variable assignments together with additional data
/*
 * \note On 64-bit systems additional data is stored in a separate container.
 */
struct ReasonStore64 : PodVector<Antecedent>::type {
	uint32  dataSize() const               { return (uint32)data_.size(); }
	void    dataResize(uint32 nv)          { if (nv > dataSize()) data_.resize(nv, UINT32_MAX); }
	uint32  data(uint32 v) const           { return data_[v]; }
	void    setData(uint32 v, uint32 data) { dataResize(v+1); data_[v] = data; }
	VarVec  data_;
	struct  value_type : std::pair<Antecedent, uint32> {
		value_type(const Antecedent& a, uint32 d) : std::pair<Antecedent, uint32>(a, d) {}
		const Antecedent& ante() const { return first;  }
		      uint32      data() const { return second; }
	};
};

//! Stores assignment related information.
/*!
 * For each variable v, the class stores 
 *  - v's current value (value_free if unassigned)
 *  - the decision level on which v was assign (only valid if value(v) != value_free)
 *  - the reason why v is in the assignment (only valid if value(v) != value_free)
 *  - (optionally) some additional data associated with the reason
 *  .
 * Furthermore, the class stores the sequences of assignments as a set of
 * true literals in its trail-member.
 */
class Assignment  {
public:
	typedef PodVector<uint32>::type     AssignVec;
	typedef PodVector<uint8>::type      SavedVec;
	typedef bk_lib::detail::if_then_else<
		sizeof(Constraint*)==sizeof(uint64)
		, ReasonStore64
		, ReasonStore32>::type            ReasonVec;
	typedef ReasonVec::value_type       ReasonWithData;
	Assignment() : front(0), eliminated_(0) { }
	LitVec            trail;   // assignment sequence
	LitVec::size_type front;   // "propagation queue"
	bool              qEmpty() const { return front == trail.size(); }
	uint32            qSize()  const { return (uint32)trail.size() - front; }
	Literal           qPop()         { return trail[front++]; }
	void              qReset()       { front  = trail.size(); }

	//! number of variables in the three-valued assignment
	uint32            numVars()    const { return (uint32)assign_.size(); }
	//! number of assigned variables
	uint32            assigned()   const { return (uint32)trail.size();   }
	//! number of free variables
	uint32            free()       const { return numVars() - (assigned()+eliminated_);   }
	//! returns the largest possible decision level
	uint32            maxLevel()   const { return (1u<<28)-1; }
	//! returns v's value in the three-valued assignment
	ValueRep          value(Var v) const { return ValueRep(assign_[v] & 3u); }
	//! returns v's previously saved value in the three-valued assignment
	ValueRep          saved(Var v) const { return v < saved_.size() ? saved_[v] : value_free; }
	//! returns the decision level on which v was assigned if value(v) != value_free
	uint32            level(Var v) const { return assign_[v] >> 4u; }
	//! returns the reason for v being assigned if value(v) != value_free
	const Antecedent& reason(Var v)const { return reason_[v]; }
	//! returns the number of allocated data slots
	uint32            numData()    const { return reason_.dataSize(); }
	//! returns the reason data associated with v
	uint32            data(Var v)  const { assert(v < reason_.dataSize()); return reason_.data(v); }

	//! resize to nv variables
	void resize(uint32 nv) {
		assign_.resize(nv);
		reason_.resize(nv);
	}
	//! adds var to assignment - initially the new var is unassigned
	Var addVar() {
		assign_.push_back(0);
		reason_.push_back(0);
		return numVars()-1;
	}
	//! allocates data slots for nv variables to be used for storing additional reason data
	void requestData(uint32 nv) {
		reason_.dataResize(nv);
	}
	//! eliminate var from assignment
	void eliminate(Var v) {
		assert(value(v) == value_free && "Can not eliminate assigned var!\n");
		setValue(v, value_true);
		++eliminated_;
	}

	//! assigns p.var() on level lev to the value that makes p true and store x as reason for the assignment
	/*!
	 * \return true if the assignment is consistent. False, otherwise.
	 * \post if true is returned, p is in trail. Otherwise, ~p is.
	 */
	bool assign(Literal p, uint32 lev, const Antecedent& x) {
		const Var      v   = p.var();
		const ValueRep val = value(v);
		if (val == value_free) {
			assign_[v] = (lev<<4) + trueValue(p);
			reason_[v] = x;
			trail.push_back(p);
			return true;
		}
		return val == trueValue(p);
	}
	bool assign(Literal p, uint32 lev, Constraint* c, uint32 data) {
		const Var      v   = p.var();
		const ValueRep val = value(v);
		if (val == value_free) {
			assign_[v] = (lev<<4) + trueValue(p);
			reason_[v] = c;
			reason_.setData(v, data);
			trail.push_back(p);
			return true;
		}
		return val == trueValue(p);
	}
	//! undos all assignments in the range trail[first, last).
	/*!
	 * \param first first assignment to be undone
	 * \param save  if true, previous assignment of a var is saved before it is undone
	 */
	void undoTrail(LitVec::size_type first, bool save) {
		if (!save) { popUntil<&Assignment::clearValue>(trail[first]); }
		else       { saved_.resize(numVars(), 0); popUntil<&Assignment::saveAndClear>(trail[first]); }
		front  = trail.size();
	}
	//! undos the last assignment
	void undoLast() { clearValue(trail.back().var()); trail.pop_back(); }
	//! returns the last assignment as a true literal
	Literal last() const { return trail.back(); }
	Literal&last()       { return trail.back(); }
	//! sets val as "previous value" of v
	void setSavedValue(Var v, ValueRep val) {
		if (saved_.size() <= v) saved_.resize(v+1, 0);
		saved_[v] = val;
	}
	/*!
	 * \name implementation functions
	 * Low-level implementation functions. Use with care and only if you
	 * know what you are doing!
	 */
	//@{
	bool seen(Var v, uint8 m) const { return (assign_[v] & (m<<2)) != 0; }
	void setSeen(Var v, uint8 m)    { assign_[v] |= (m<<2); }
	void clearSeen(Var v)           { assign_[v] &= ~uint32(12); }
	void clearValue(Var v)          { assign_[v] = 0; }
	void setValue(Var v, ValueRep val) {
		assert(value(v) == val || value(v) == value_free);
		assign_[v] = val;
	}
	void setReason(Var v, const Antecedent& a) { reason_[v] = a;  }
	void setData(Var v, uint32 data) { reason_.setData(v, data); }
	void copyAssignment(Assignment& o) const { o.assign_ = assign_; }
	//@}
private:
	Assignment(const Assignment&);
	Assignment& operator=(const Assignment&);
	void    saveAndClear(Var v) { saved_[v] = value(v); clearValue(v); }
	template <void (Assignment::*op)(Var v)>
	void popUntil(Literal stop) {
		Literal p;
		do {
			p = trail.back(); trail.pop_back();
			(this->*op)(p.var());
		} while (p != stop);
	}
	AssignVec assign_; // for each var: three-valued assignment
	ReasonVec reason_; // for each var: reason for being assigned (+ optional data)
	SavedVec  saved_;  // for each var: previous assignment
	uint32    eliminated_;
};

//! Stores information about a literal that is implied on an earlier level than the current decision level.
struct ImpliedLiteral {
	typedef Assignment::ReasonWithData AnteInfo;
	ImpliedLiteral(Literal a_lit, uint32 a_level, const Antecedent& a_ante, uint32 a_data = UINT32_MAX) 
		: lit(a_lit)
		, level(a_level)
		, ante(a_ante, a_data) {
	}
	Literal     lit;    /**< The implied literal */
	uint32      level;  /**< The earliest decision level on which lit is implied */
	AnteInfo    ante;   /**< The reason why lit is implied on decision-level level */
};

struct CCMinRecursive {
	enum State { state_open = 0, state_poison = 1, state_removable = 2 };
	void  init(uint32 numV) { extra.resize(numV,0); }
	State state(Literal p) const { return State(extra[p.var()]); }
	bool  checkRecursive(Literal p) {
		if (state(p) == state_open) { p.clearWatch(); dfsStack.push_back(p); }
		return state(p) != state_poison;
	}
	void  markVisited(Literal p, State st) {
		if (state(p) == state_open) {
			visited.push_back(p.var());
		}
		extra[p.var()] = static_cast<uint8>(st);
	}
	void clear() {
		for (; !visited.empty(); visited.pop_back()) {
			extra[visited.back()] = 0;
		}
	}
	typedef PodVector<uint8>::type DfsState;
	LitVec   dfsStack;
	VarVec   visited;
	DfsState extra;
};
//@}
}
#endif
