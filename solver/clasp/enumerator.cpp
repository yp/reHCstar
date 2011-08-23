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
#include <clasp/enumerator.h>
#include <clasp/minimize_constraint.h>
#include <clasp/solver.h>
namespace Clasp { 

/////////////////////////////////////////////////////////////////////////////////////////
// Enumerator::SolveConstraint
/////////////////////////////////////////////////////////////////////////////////////////
Enumerator::SolveConstraint::SolveConstraint() : mini_(0) {}
Constraint::PropResult Enumerator::SolveConstraint::propagate(Solver&, Literal, uint32&) { return PropResult(true, false); }
void Enumerator::SolveConstraint::reason(Solver&, Literal, LitVec&) {}
void Enumerator::SolveConstraint::destroy(Solver* s, bool detach) {
	if (mini_) mini_->destroy(s, detach);
	Constraint::destroy(s, detach);
}
MinimizeConstraint* Enumerator::SolveConstraint::cloneMini(Solver& s) const {
	MinimizeConstraint* ret = 0;
	if (mini_) {
		ret = MinimizeBuilder::attach(s, mini_->shared_unsafe(), (s.id() & 1) == 0);
	}
	return ret;
}
/////////////////////////////////////////////////////////////////////////////////////////
// Enumerator
/////////////////////////////////////////////////////////////////////////////////////////
Enumerator::ProgressReport::ProgressReport()  {}
Enumerator::ProgressReport::~ProgressReport() {}
Enumerator::Report::Report()  {}
Enumerator::Enumerator(Report* r) : enumerated(0), numModels_(1), report_(r), progress_(0), mini_(0), restartOnModel_(false)  {}
Enumerator::~Enumerator()                               { if (mini_) mini_->destroy(); }
void Enumerator::setReport(Report* r)                   {  report_   = r; }
void Enumerator::enableProgressReport(ProgressReport* r){  progress_ = r; }
void Enumerator::setRestartOnModel(bool r)              { restartOnModel_ = r; }
void Enumerator::updateModel(Solver&)                   {}
bool Enumerator::ignoreSymmetric() const                { return optimize(); }
bool Enumerator::optimize()        const                { return mini_ && mini_->mode() == MinimizeMode_t::optimize; }
bool Enumerator::optimizeHierarchical() const           { return mini_ && mini_->hierarchical(); }
void Enumerator::setMinimize(SharedMinimizeData* min, int h) {  
	mini_         = min;   
	minHeuristic_ = h;
}
void Enumerator::startInit(SharedContext& ctx) { 
	enumerated = 0;
	doInit(ctx, 0, true);
}
void Enumerator::enumerate(uint64 m) {
	numModels_ = m; 
}
Enumerator::SolveConstraint* Enumerator::endInit(SharedContext& ctx, uint32 t) { 
	SolveConstraint* c = doInit(ctx, t, false);
	if (mini_) {
		if (!c) { 
			struct MinimizeHolder : public SolveConstraint {
				Constraint* cloneAttach(Solver&s) { 
					MinimizeHolder* r = new MinimizeHolder();
					r->setMinimize(cloneMini(s));
					return r;
				}
			};
			c = new MinimizeHolder(); 
		}
		c->setMinimize(MinimizeBuilder::attach(*ctx.master(), mini_, minHeuristic_));
	}
	return c;
}
Enumerator::SolveConstraint* Enumerator::constraint(const Solver& s) const {
	return static_cast<Enumerator::SolveConstraint*>(s.getEnumerationConstraint());
}
bool Enumerator::continueFromModel(Solver& s) {
	if (restartOnModel_) {
		s.undoUntil(0);
	}
	if (mini_) {
		if (optimizeHierarchical()) {
			s.strengthenConditional();
		}
		MinimizeConstraint* m = constraint(s)->minimize();
		if (mini_->mode() == MinimizeMode_t::optimize && !m->propagateNewOptimum(s)) {
			return false;
		}
		if (restartOnModel_ && s.queueSize() == 0) {
			return m->modelHeuristic(s);
		}
	}
	return true;
}
Enumerator::Result Enumerator::backtrackFromModel(Solver& s) {
	assert(s.numFreeVars() == 0 && !s.hasConflict());
	bool expandSym = !ignoreSymmetric();
	do {
		if (!onModel(s)) { // enough models enumerated?
			return enumerate_stop_enough;
		}
		// Process symmetric models, i.e. models that differ only in the 
		// assignment of atoms outside of the solver's assignment. 
		// Typical example: vars eliminated by the SAT-preprocessor
	} while (s.nextSymModel(expandSym));
	if (activeLevel_ <= s.rootLevel() || !(backtrack(s) && continueFromModel(s))) {
		s.undoUntil(0, true);
		return enumerate_stop_complete;
	}
	return enumerate_continue;
}
bool Enumerator::trivial(bool disjoint) {
	return !optimize() && isTrivial(disjoint);
}
bool Enumerator::update(Solver& s, bool disjoint) {
	bool ret = true;
	if (optimize()) {
		ret      = constraint(s)->minimize()->integrateNext(s);
		disjoint = true; // enforced by minimize constraint
	}
	return ret && updateConstraint(s, disjoint);
}

bool Enumerator::onModel(Solver& s) {
	updateModel(s); // Hook for derived classes
	SolveConstraint* c = constraint(s);
	activeLevel_ = mini_ ? c->minimize()->setModel(s)+1 : s.decisionLevel();
	++enumerated;
	if (report_) {
		report_->reportModel(s, *this);
	}
	return numModels_ == 0 || --numModels_ != 0;
}

}
