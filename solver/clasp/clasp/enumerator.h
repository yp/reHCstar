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
#ifndef CLASP_ENUMERATOR_H_INCLUDED
#define CLASP_ENUMERATOR_H_INCLUDED

#ifdef _MSC_VER
#pragma once
#endif
#include <clasp/literal.h>
#include <clasp/constraint.h>

namespace Clasp { 
class  Solver;
class  MinimizeConstraint;
struct SharedMinimizeData;
class  SharedContext;
class  SatPreprocessor;

/**
 * \defgroup enumerator Enumerators and related classes
 */
//@{

//! Interface for enumerating models
/*!
 * Enumerators are global w.r.t to one search operation, that is
 * even if the search operation itself is a parallel search, there
 * shall be only one enumerator and it is the user's responsibility
 * to protect the enumerator with appropriate locking.
 *
 * Concrete enumerators may create and install an enumeration constraint 
 * in each solver to support (parallel) enumeration
 */
class Enumerator {
public:
	enum ProgressType    { progress_restart, progress_reduce, progress_model };
	//! Interface used by the enumerator to report restart events
	class ProgressReport {
	public:
		ProgressReport();
		virtual ~ProgressReport();
		//! preprocessing progress
		virtual void reportPreProgress(char /* type */, const SatPreprocessor& /* p */, uint32 /* min */, uint32 /* max */) {}
		//! solving progress
		virtual void reportProgress(ProgressType /* t */, const Solver& /* s */, uint64 /* maxCfl */, uint32 /* maxL */) {}
	private:
		ProgressReport(const ProgressReport&);
		ProgressReport& operator=(const ProgressReport&);
	};
	//! Interface used by the enumerator to report model events
	class Report : public ProgressReport {
	public:
		Report();
		//! the solver has found a new model
		virtual void reportModel(const Solver& /* s */, const Enumerator& /* self */) {}
		//! enumeration has terminated
		virtual void reportSolution(const Enumerator& /* self */, bool /* complete */) {}
	};
	//! A solver-local (i.e. thread-local) constraint to support enumeration
	class SolveConstraint : public Constraint {
	public:
		SolveConstraint();
		//! detaches and destroys minimize constraint if any
		void destroy(Solver* s, bool detach);
		void setMinimize(MinimizeConstraint* m) { mini_ = m;    } 
		MinimizeConstraint* minimize()  const   { return mini_; }
		MinimizeConstraint* cloneMini(Solver& s) const;
	protected:
		//! returns (true, false)
		PropResult propagate(Solver&, Literal, uint32&);
		//! noop
		void       reason(Solver&, Literal, LitVec&);
	private:
		MinimizeConstraint* mini_;
	};
	explicit Enumerator(Report* r = 0);
	virtual ~Enumerator();
	
	//! How to continue after a model was found
	enum Result {
		enumerate_continue,      /**< continue enumeration */
		enumerate_stop_complete, /**< stop because search space is exceeded  */
		enumerate_stop_enough    /**< stop because enough models have been enumerated */
	};

	//! sets the report callback to be used during enumeration
	/*!
	 * \note Ownership is *not* transferred and r must be valid
	 * during complete search
	 */
	void setReport(Report* r);
	//! enables progress reporting via the given report callback
	void enableProgressReport(ProgressReport* r);

	//! if true, does a search-restart every time a model is found
	void setRestartOnModel(bool r);
	
	//! starts initialization of this enumerator
	/*!
	 * Must be called once before search is started and before solvers
	 * can be attached. Shall freeze relevant variables.
	 * \note In the incremental setting, startInit() must be called once for each incremental step
	 */
	void startInit(SharedContext& ctx);

	//! sets the (shared) minimize object to be used during enumeration
	/*!
	 * \note Ownership is transferred.
	 */
	void setMinimize(SharedMinimizeData* min, int useHeuristic);
	const SharedMinimizeData* minimize() const { return mini_; }

	//! sets the maximum number of models to enumerate
	/*!
	 * \param numModels number of models to enumerate (0 means all models)
	 */
	void enumerate(uint64 numModels);

	bool enumerate() const { return numModels_ != 1; }

	//! completes the initialization of this enumerator
	/*!
	 * Sets the number of solvers sharing this enumerator.
	 */
	SolveConstraint* endInit(SharedContext& ctx, uint32 shareCount);

	uint64 enumerated; /**< Number of models enumerated so far */

	//! returns the enumeration constraint attached to s
	SolveConstraint* constraint(const Solver& s) const;

	//! called whenever a solver has found a model.
	/*!
	 * The return value determines how search should proceed.
	 * If enumerate_continue is returned, the enumerator has
	 * removed at least one decision level and search should
	 * continue from the new level.
	 * Otherwise, the search shall be stopped.
	 * \pre The solver contains a model and DL = s.decisionLevel()
	 * \post If enumerate_continue is returned:
	 *  - s.decisionLevel() < DL and s.decisionLevel() >= s.rootLevel()
	 *
	 * \note If enumerate_continue is returned, the caller must call s.propagate()
	 * \note The function is not concurrency-safe, i.e. in a parallel search
	 *       at most one solver shall call the function at any one time.
	 */
	Result backtrackFromModel(Solver& s);

	//! returns whether or not this enumerator supports full restarts once a model was found
	virtual bool supportsRestarts() const { return true; }
	
	//! returns whether enumerator is trivial w.r.t search scheme
	/*!
	 * If trivial() returns false, this signals that other solvers
	 * must be updated whenever any solver finds a new model.
	 */
	bool trivial(bool disjointPath);

	//! updates solver s with new model-related information
	/*!
	 * The function is used to transfer information in
	 * parallel search. Whenever a solver s1 found a new model
	 * (i.e. called backtrackFromModel()), all other solvers
	 * shall eventually call update() to incorporate any new information.
	 */
	bool update(Solver& s, bool disjointPath);

	//! called once after search has stopped
	/*!
	 * The function notifies the installed Report object
	 */
	void reportResult(bool complete) {
		if (report_) report_->reportSolution(*this, complete);
	}

	void reportPreProgress(char type, const SatPreprocessor& p, uint32 min, uint32 max) {
		if (progress_) progress_->reportPreProgress(type, p, min, max); 
	}

	void reportProgress(ProgressType t, const Solver& s, uint64 maxCfl, uint32 maxL) {
		if (progress_) progress_->reportProgress(t, s, maxCfl, maxL); 
	}
	
	//! returns true if optimization is active
	bool   optimize() const;
	bool   optimizeHierarchical() const;
protected:
	uint32 getHighestActiveLevel() const { return activeLevel_; }
	bool   continueFromModel(Solver& s);
	virtual SolveConstraint* doInit(SharedContext& s, uint32 numThreads, bool start) = 0;
	virtual bool backtrack(Solver& s) = 0;
	virtual bool isTrivial(bool disjoint) const = 0;
	virtual bool updateConstraint(Solver& s, bool disjoint) = 0;
	virtual void updateModel(Solver& s);
	virtual bool ignoreSymmetric() const;
private:
	Enumerator(const Enumerator&);
	Enumerator& operator=(const Enumerator&);
	bool onModel(Solver& s);
	uint64              numModels_;
	Report*             report_;
	ProgressReport*     progress_;
	SharedMinimizeData* mini_;
	uint32              activeLevel_;
	int                 minHeuristic_;
	bool                restartOnModel_;
};

class NullEnumerator : public Enumerator {
public:
	NullEnumerator() {}
private:
	SolveConstraint* doInit(SharedContext&, uint32, bool) { return 0; }
	bool backtrack(Solver&)              { return false; }
	bool isTrivial(bool) const           { return true; }
	bool updateConstraint(Solver&, bool) { return true; }

};

//@}

}

#endif
