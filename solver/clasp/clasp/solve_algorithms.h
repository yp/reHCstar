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
#ifndef CLASP_SOLVE_ALGORITHMS_H_INCLUDED
#define CLASP_SOLVE_ALGORITHMS_H_INCLUDED

#ifdef _MSC_VER
#pragma warning (disable : 4200) // nonstandard extension used : zero-sized array
#pragma once
#endif

#include <clasp/literal.h>
#include <clasp/constraint.h>

/*!
 * \file 
 * Defines top-level functions for solving problems.
 */
namespace Clasp { 

class  Solver;
class  SharedContext;
struct SharedMinimizeData;
class  MinimizeConstraint;
class  Enumerator;


// Implements clasp's configurable schedule-strategies.
// Note: Currently all schedule-strategies can be easily implemented using one class.
// In the future, as new strategies emerge, the class should be replaced with a strategy hierarchy
struct ScheduleStrategy {
public:
	ScheduleStrategy(uint32 base = 0, double grow = 0.0, uint64 outer = 0) {
		init(base, grow, outer);
	}
	//! configure the schedule
	/*!
	 * \param base  initial interval or run-length
	 * \param inc   grow factor
	 * \param outer max restart interval, repeat sequence if reached
	 * \note
	 *  if base is equal to 0, schedule is disabled.
	 *  if inc is equal to 0, luby-schulde is used and base is interpreted as run-length
	 */
	void   init(uint32 base, double grow, uint64 outer);
	void   reset()  { idx_ = 0; }
	uint64 next();
	uint64 current() const;
private:
	uint64 geomR() const;
	uint64 lubyR() const;
	double  grow_;
	uint64  outer_;
	uint32  base_;
	uint32  idx_;
};

///////////////////////////////////////////////////////////////////////////////
// Parameter-object for the solve function
// 
///////////////////////////////////////////////////////////////////////////////
//! Aggregates restart-parameters to configure restarts during search.
/*!
 * clasp currently supports four different restart-strategies:
 *  - fixed-interval restarts: restart every n conflicts
 *  - geometric-restarts: restart every n1 * n2^k conflict (k >= 0)
 *  - inner-outer-geometric: similar to geometric but sequence is repeated once bound outer is reached. Then, outer = outer*n2
 *  - luby's restarts: see: Luby et al. "Optimal speedup of las vegas algorithms."
 *  .
 */
struct RestartParams {
	RestartParams() : sched(100, 1.5), local(false), bounded(false), resetOnModel(false) {}
	ScheduleStrategy sched; /**< restart schedule to use */
	bool local;             /**< local restarts, i.e. restart if number of conflicts in *one* branch exceed threshold */
	bool bounded;           /**< allow (bounded) restarts after first solution was found */
	bool resetOnModel;      /**< repeat restart strategy after each solution */
};

//! Aggregates parameters for the nogood deletion heuristic used during search
class ReduceParams {
private:
	double frac_;
	double inc_;
	double max_;
	float  rf_;
	uint32 base_;
	uint32 iMin_, iMax_;
public:
	ReduceParams() : frac_(3.0), inc_(1.1), max_(3.0), rf_(.75f)
		, base_(5000), iMin_(10), iMax_(UINT32_MAX)
		, sched(0,0,0), reduceOnRestart(false), estimate(false), disable(false) {}
	
	//! sets the initial problem size used to compute initial db size
	void setProblemSize(uint32 base)       { base_ = base; }
	//! sets initial minimum and maximum size of db
	void setInit(uint32 iMin, uint32 iMax) { iMin_ = iMin; iMax_ = iMax; }
	//! sets fraction of nogoods to delete on reduction
	void setReduceFraction(double rf) {
		rf_ = (float)std::max(0.01, std::min(rf, 1.0));
	}
	uint32 initMin()const { return iMin_; }
	uint32 initMax()const { return iMax_; }
	uint32 init()   const;
	double inc()    const { return inc_; }
	double frac()   const { return frac_; }
	uint32 bound()  const { return (uint32)std::min(std::max(base_, iMin_)*max_, double(std::numeric_limits<uint32>::max())); }
	float  reduceFrac() const { return rf_; }
	
	ScheduleStrategy sched; /**< secondary deletion schedule */
	bool reduceOnRestart;   /**< delete some nogoods on each restart */
	bool estimate;          /**< use estimate of problem complexity to init problem size */
	bool disable;           /**< do not delete any nogoods */
	//! configure reduce strategy
	/*!
	 * \param frac     init db size to problemSize/frac
	 * \param inc      grow factor applied after each restart
	 * \param maxF     stop growth once db size is > problemSize*maxF
	 */
	void setStrategy(double frac, double inc, double maxF) {
		frac_ = std::max(0.0001, frac);
		inc_  = std::max(1.0   , inc);
		max_  = std::max(0.0001, maxF);
	}
};

struct SearchLimits {
	SearchLimits(uint64 conf = UINT64_MAX, uint64 r = UINT64_MAX) 
		: conflicts(conf)
		, restarts(r) {
	}
	uint64 conflicts;
	uint64 restarts;
};

//! Parameter-Object for configuring search-parameters
/*!
 * \ingroup solver
 */
struct SolveParams {
	//! creates a default-initialized object.
	/*!
	 * The following parameters are used:
	 * restart      : quadratic: 100*1.5^k / no restarts after first solution
	 * shuffle      : disabled
	 * deletion     : initial size: vars()/3, grow factor: 1.1, max factor: 3.0, do not reduce on restart
	 * randomization: disabled
	 * randomProp   : 0.0 (disabled)
	 * enumerator   : no
	 */
	SolveParams();
	
	RestartParams restart;
	ReduceParams  reduce;
	mutable SearchLimits  limits;

	//! sets the shuffle-parameters to use during search.
	/*!
	 * \param first   Shuffle program after first restarts
	 * \param next    Re-Shuffle program every next restarts
	 * \note
	 *  if first is equal to 0, shuffling is disabled.
	 */
	void setShuffleParams(uint32 first, uint32 next) {
		shuffleFirst_ = first;
		shuffleNext_  = next;
	}

	//! sets the randomization-parameters to use during search.
	/*!
	 * \param runs number of initial randomized-runs
	 * \param cfl number of conflicts comprising one randomized-run
	 */
	void setRandomizeParams(uint32 runs, uint32 cfls) {
		if (runs > 0 && cfls > 0) {
			randRuns_ = runs;
			randConflicts_ = cfls;
		}
	}

	//! sets the probability with which choices are made randomly instead of with the installed heuristic.
	void setRandomProbability(double p) {
		if (p >= 0.0 && p <= 1.0) {
			randFreq_ = p;
		}
	}
	// accessors
	uint32  randRuns()      const { return randRuns_; }
	uint32  randConflicts() const { return randConflicts_; }
	double  randomProbability() const { return randFreq_; }
	
	uint32      shuffleBase() const { return shuffleFirst_; }
	uint32      shuffleNext() const { return shuffleNext_; }
private:
	double        randFreq_;
	uint32        randRuns_;
	uint32        randConflicts_;
	uint32        shuffleFirst_;
	uint32        shuffleNext_;
};


///////////////////////////////////////////////////////////////////////////////
// Basic solve functions
///////////////////////////////////////////////////////////////////////////////

//! Basic sequential search
/*!
 * \ingroup solver
 * \relates Solver
 * \param ctx The context containing the problem.
 * \param p   The solve parameters to use.
 *
 * \return
 *  - true: if the search stopped before the search-space was exceeded.
 *  - false: if the search-space was completely examined.
 * 
 */
bool solve(SharedContext& ctx, const SolveParams& p);

//! Basic sequential search under assumptions
/*!
 * \ingroup solver
 * \relates Solver
 * The use of assumptions allows for incremental solving. Literals contained
 * in assumptions are assumed to be true during search but are undone before solve returns.
 *
 * \param ctx The context containing the problem.
 * \param p   The solve parameters to use.
 * \param assumptions The list of initial unit-assumptions
 *
 * \return
 *  - true: if the search stopped before the search-space was exceeded.
 *  - false: if the search-space was completely examined.
 * 
 */
bool solve(SharedContext& ctx, const SolveParams& p, const LitVec& assumptions);

///////////////////////////////////////////////////////////////////////////////
// General solve
///////////////////////////////////////////////////////////////////////////////

//! Interface for solve algorithms
/*!
 * \ingroup solver
 * \relates Solver
 * SolveAlgorithm objects wrap an enumerator and
 * implement concrete solve algorithms
 */
class SolveAlgorithm {
public:
	explicit SolveAlgorithm();
	virtual ~SolveAlgorithm();
	
	//! force termination of current solve process
	/*!
	 * shall return true if termination is supported, otherwise false
	 */
	virtual bool   terminate() = 0;

	//! Runs the solve algorithm
	/*!
	 * \param ctx    A fully initialized context object containing the problem.
	 * \param p      The solve parameters for the master solver.
	 * \param assume A list of initial unit-assumptions
	 *
	 * \return
	 *  - true: if the search stopped before the search-space was exceeded.
	 *  - false: if the search-space was completely examined.
	 *
	 * \note 
	 * The use of assumptions allows for incremental solving. Literals contained
	 * in assumptions are assumed to be true during search but are undone before solve returns.
	 */
	bool solve(SharedContext& ctx, const SolveParams& p, LitVec assume);
protected:
	//! The default implementation simply forwards the call to the enumerator
	virtual bool backtrackFromModel(Solver& s);
	//! The default implementation simply forwards the call to the enumerator
	virtual void reportProgress(int type, Solver& s, uint64 maxCfl, uint32 maxLearnt);
	virtual bool doSolve(Solver& s, const SolveParams& p, const LitVec& assume) = 0;
	bool     initPath(Solver& s,  const LitVec& gp);
	ValueRep solvePath(Solver& s, const SolveParams& p);
private:
	SolveAlgorithm(const SolveAlgorithm&);
	SolveAlgorithm& operator=(const SolveAlgorithm&);
};

//! A basic algorithm for single-threaded sequential solving without search-space splitting
class SimpleSolve : public SolveAlgorithm {
public:
	SimpleSolve() : SolveAlgorithm() {}
	bool   terminate();
private:
	bool   doSolve(Solver& s, const SolveParams& p, const LitVec& assume);
};

}
#endif
