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
#include <clasp/solve_algorithms.h>
#include <clasp/solver.h>
#include <clasp/minimize_constraint.h>
#include <clasp/enumerator.h>
#include <cmath>
using std::log;
namespace Clasp { 
/////////////////////////////////////////////////////////////////////////////////////////
// SolveParams
/////////////////////////////////////////////////////////////////////////////////////////
SolveParams::SolveParams() 
	: randFreq_(0.0)
	, randRuns_(0), randConflicts_(0)
	, shuffleFirst_(0), shuffleNext_(0) {
}

uint32 ReduceParams::init() const {
	if (disable) return UINT32_MAX;
	uint32 ret = static_cast<uint32>(base_/frac_);
	if      (ret < iMin_) ret = iMin_;
	else if (ret > iMax_) ret = iMax_;
	return ret;
}
/////////////////////////////////////////////////////////////////////////////////////////
// Schedule
/////////////////////////////////////////////////////////////////////////////////////////
void ScheduleStrategy::init(uint32 base, double grow, uint64 outer) {
	grow_ = grow == 0.0 || grow >= 1.0 ? grow : 1.0;
	outer_= (outer ? outer : UINT64_MAX);
	base_ = base;
	idx_  = 0;
}
uint64 ScheduleStrategy::current() const {
	uint64 x;
	if      (base_ == 0)      x = UINT64_MAX;
	else if (grow_ == 0)      x = lubyR();
	else                      x = geomR();
	return x;
}
uint64 ScheduleStrategy::next() {
	++idx_;
	uint64 x = current();
	if (base_ != 0 && grow_ != 0 && x > outer_) {
		idx_    = 0;
		outer_  = static_cast<uint64>(outer_*grow_);
		x       = geomR();
	}
	return x;
}
uint64 ScheduleStrategy::geomR() const { 
	return static_cast<uint64>(base_ * pow(grow_, (double)idx_)); 
}
uint64 ScheduleStrategy::lubyR() const {
	uint32 k = idx_+1;
	while (k) {
		uint32 nk = static_cast<uint32>(log((double)k) / log(2.0)) + 1;
		if (k == ((uint32(1) << nk) - 1)) {
			return base_ * (uint32(1) << (nk-1));
		}
		k -= uint32(1) << (nk-1);
		++k;
	}
	return base_;
}
/////////////////////////////////////////////////////////////////////////////////////////
// solve
/////////////////////////////////////////////////////////////////////////////////////////
bool solve(SharedContext& ctx, const SolveParams& p) {
	return SimpleSolve().solve(ctx, p, LitVec());
}

bool solve(SharedContext& ctx, const SolveParams& p, const LitVec& assumptions) {
	return SimpleSolve().solve(ctx, p, assumptions);
}

/////////////////////////////////////////////////////////////////////////////////////////
// SolveAlgorithm
/////////////////////////////////////////////////////////////////////////////////////////
SolveAlgorithm::SolveAlgorithm()  {
}
SolveAlgorithm::~SolveAlgorithm() {}

bool SolveAlgorithm::backtrackFromModel(Solver& s) { 
	return s.sharedContext()->enumerator()->backtrackFromModel(s) == Enumerator::enumerate_continue;
}

void SolveAlgorithm::reportProgress(int t, Solver& s, uint64 maxCfl, uint32 maxL) {
	return s.sharedContext()->enumerator()->reportProgress(Enumerator::ProgressType(t), s, maxCfl, maxL);
}
bool SolveAlgorithm::solve(SharedContext& ctx, const SolveParams& p, LitVec assume) {
	assert(ctx.master() && "SharedContext not initialized!\n");
	if (!isSentinel(ctx.tagLiteral())) {
		assume.push_back(ctx.tagLiteral());
	}
	bool r = doSolve(*ctx.master(), p, assume);
	ctx.detach(*ctx.master());
	return r;
}
bool SolveAlgorithm::initPath(Solver& s, const LitVec& path) {
	assert(!s.hasConflict() && s.decisionLevel() == 0);
	for (LitVec::size_type i = 0, end = path.size(); i != end; ++i) {
		Literal p = path[i];
		if (s.value(p.var()) == value_free) {
			s.assume(p); --s.stats.choices;
			// increase root level - assumption can't be undone during search
			s.pushRootLevel();
			if (!s.propagate())  return false;
		}
		else if (s.isFalse(p)) return false;
	}
	return true;
}

ValueRep SolveAlgorithm::solvePath(Solver& s, const SolveParams& p) {
	if (s.hasConflict()) return false;
	double maxLearnts         = p.reduce.init();
	const double boundLearnts = p.reduce.bound();
	if (maxLearnts < s.numLearntConstraints()) {
		maxLearnts = static_cast<double>(s.numLearntConstraints()) + p.reduce.initMin();
		maxLearnts = std::min(maxLearnts, (double)UINT32_MAX);
	}
	ScheduleStrategy rs = p.restart.sched;
	ScheduleStrategy ds = p.reduce.sched;
	ValueRep result = value_free;
	uint32 randRuns = p.randRuns();
	double randFreq = randRuns == 0 ? p.randomProbability() : 1.0;
	uint32 shuffle  = p.shuffleBase();
	SearchLimits lim= p.limits;
	uint64 off      = s.stats.conflicts;
	if (lim.restarts==0) lim.restarts = 1;
	Solver::RestartLimit rlimit(randRuns == 0 ? rs.current() : p.randConflicts(), p.restart.local);
	Solver::LearntLimit  dlimit(ds.current(), (uint32)maxLearnts);
	Enumerator::ProgressType t = Enumerator::progress_restart;
	while (result == value_free && lim.conflicts > 0 && lim.restarts > 0) {
		if (rlimit.maxConf > lim.conflicts) { rlimit.maxConf = lim.conflicts; }
		reportProgress(t, s, rlimit.maxConf, dlimit.maxLearnt);
		result         = s.search(rlimit, dlimit, randFreq);
		lim.conflicts -= std::min((s.stats.conflicts-off), lim.conflicts);
		off            = s.stats.conflicts;
		if (result == value_true) {
			if (!backtrackFromModel(s)) {
				break; // No more models requested
			}
			else {
				result   = value_free; // continue enumeration
				randRuns = 0;          // but cancel remaining probings
				randFreq = p.randomProbability();
				t        = Enumerator::progress_model;
				if (p.restart.resetOnModel) {
					rs.reset();
				}
				// After the first solution was found, we allow further restarts only if this
				// is compatible with the enumerator used. 
				rlimit.maxConf  = !p.restart.bounded && s.backtrackLevel() > s.rootLevel()
					? static_cast<uint64>(-1)
					: rs.current();
				dlimit.maxConf  = ds.current();
			}
		}
		else if (result == value_free){  // limit reached
			if (!rlimit.reached()) {
				s.reduceLearnts(p.reduce.reduceFrac());
				if (s.numLearntConstraints() >= dlimit.maxLearnt) {
					maxLearnts       += s.numLearntConstraints();
					dlimit.maxLearnt  = (uint32)maxLearnts;
				}
				dlimit.maxConf   = dlimit.reached() ? ds.next() : ds.current();
				t                = Enumerator::progress_reduce;
				continue;
			}
			// rlimit reached - do restart
			s.undoUntil(0);
			t = Enumerator::progress_restart;
			if (randRuns == 0) {
				rlimit.maxConf = rs.next();
				if (p.reduce.reduceOnRestart) { s.reduceLearnts(.33f); }
				if (maxLearnts != (double)UINT32_MAX && maxLearnts < boundLearnts && (s.numLearntConstraints()+rlimit.maxConf) > maxLearnts) {
					maxLearnts = std::min(maxLearnts*p.reduce.inc(), (double)UINT32_MAX);
					dlimit.maxLearnt = (uint32)maxLearnts;
				}
				if (++s.stats.restarts == shuffle) {
					shuffle += p.shuffleNext();
					s.shuffleOnNextSimplify();
				}
			}
			else if (--randRuns != 0) {
				rlimit.maxConf = p.randConflicts();
			}
			else {
				rlimit.maxConf = rs.current();
				randFreq       = p.randomProbability();
			}
			--lim.restarts;
		}
	}
	p.limits = lim;
	return result;
}
/////////////////////////////////////////////////////////////////////////////////////////
// SimpleSolve
/////////////////////////////////////////////////////////////////////////////////////////
bool SimpleSolve::terminate() { return false; }
bool SimpleSolve::doSolve(Solver& s, const SolveParams& p, const LitVec& assume) {
	s.stats.reset();
	MinimizeConstraint* min = 0;
	Enumerator*  enumerator = s.sharedContext()->enumerator();
	bool complete = true;
	// Remove any existing assumptions and simplify problem.
	// If this fails, the problem is unsat, even under no assumptions.
	bool hasWork  = s.clearAssumptions();
	while (hasWork) {
		// Add assumptions.
		// If this fails, the problem is unsat under the current assumptions
		// but not necessarily unsat.
		if (initPath(s, assume) && (!min || min->integrateNext(s))) {
			complete = (solvePath(s, p) != value_free && s.decisionLevel() == s.rootLevel());
		}
		// finished current work
		hasWork    = false;
		if (enumerator->optimizeHierarchical() && enumerator->enumerated > 0 && complete) {
			// we are done with the current level			
			min = enumerator->constraint(s)->minimize();
			// restore previous optimum - this relaxes the constraint 
			// to the last known good value
			min->restoreOptimum();
			// check if there are more levels to minize
			hasWork  = min->shared_unsafe()->optimizeNext();
		}
		// Finally, remove the assumptions again and restore
		// the solver to a usable state if possible.
		if (!s.clearAssumptions()) {
			complete = true;
			hasWork  = false;
		}
	}
	enumerator->reportResult(complete);
	return !complete;
}

}
