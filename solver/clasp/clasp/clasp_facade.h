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
#ifndef CLASP_CLASP_FACADE_H_INCLUDED
#define CLASP_CLASP_FACADE_H_INCLUDED

#ifdef _MSC_VER
#pragma warning (disable : 4200) // nonstandard extension used : zero-sized array
#pragma once
#endif

#if !defined(CLASP_VERSION)
#define CLASP_VERSION "2.0.2"
#endif
#if !defined(CLASP_LEGAL)
#define CLASP_LEGAL \
"Copyright (C) Benjamin Kaufmann\n"\
"License GPLv2+: GNU GPL version 2 or later <http://gnu.org/licenses/gpl.html>\n"\
"clasp is free software: you are free to change and redistribute it.\n"\
"There is NO WARRANTY, to the extent permitted by law."
#endif

#if defined(WITH_THREADS)
#undef WITH_THREADS
#endif
#ifndef DISABLE_MULTI_THREADING
#define WITH_THREADS 1
#else
#define WITH_THREADS 0
#endif

#include <clasp/literal.h>
#include <clasp/solver.h>
#include <clasp/enumerator.h>
#include <clasp/solve_algorithms.h>
#include <clasp/heuristics.h>
#include <clasp/lookahead.h>
#include <clasp/program_builder.h>
#include <clasp/unfounded_check.h>
#include <clasp/reader.h>
#include <clasp/util/misc_types.h>
#include <string>

/*!
 * \file 
 * This file provides a facade around the clasp library. 
 * I.e. a simplified interface for (incrementally) solving a problem using
 * some configuration (set of parameters).
 */
namespace Clasp {

/////////////////////////////////////////////////////////////////////////////////////////
// Parameter configuration
/////////////////////////////////////////////////////////////////////////////////////////
//! Options that control the decision heuristic
struct HeuristicOptions {
	HeuristicOptions();
	typedef Lookahead::Type LookaheadType;
	DecisionHeuristic* createHeuristic() const;
	std::string   name;        /**< name of decision heuristic */
	LookaheadType lookahead;   /**< type of lookahead */
	int           lookaheadNum;/**< number of times lookahead is used */
	int           loops;       /**< consider loops in heuristic (0: no, 1: yes, -1: let heuristic decide) */
	union Extra {
	int      berkMax;          /**< only for Berkmin */
	int      vmtfMtf;          /**< only for Vmtf    */
	} extra;
	bool     berkMoms;         /**< only for Berkmin */
	bool     berkHuang;	       /**< only for Berkmin */
	bool     berkOnce;         /**< only for Berkmin */
	bool     nant;             /**< only for unit    */
};

//! Local options holding options for one solver instance
class LocalOptions {
public:
	typedef DefaultUnfoundedCheck::ReasonStrategy LoopMode;
	explicit LocalOptions(Solver* s = 0);
	~LocalOptions();
	void setSolver(Solver* s) { solver_ = s; }
	bool validate(std::string& err);
	HeuristicOptions& heuristic();
	void              initFrom(const LocalOptions& other);
	void              applyHeuristic();
	Solver&           solver() { return *solver_; }
	SolveParams       solve;   /**< strategies used during solving */ 
	LoopMode          loopRep; /**< how to represent loops? */
private:
	LocalOptions(const LocalOptions&);
	LocalOptions& operator=(const LocalOptions&);
	Solver*           solver_;
	HeuristicOptions* heuristic_;
};

//! Options for configuring parallel (multi-threaded) solving
struct ThreadOptions {  
	ThreadOptions() : forceGP(false), genTemplate(false) {}
	struct Distribution { /**< Nogood distribution options */
		Distribution();
		enum Filter { filter_no = 0, filter_gp = 1, filter_sat = 2, filter_heuristic = 3 };
		uint32 grace;       /**< lower bound on number of shared nogoods to keep */
		uint8  lbd;         /**< upper bound on lbd for sharing nogoods */
		uint8  types;       /**< types of nogoods to share */
		uint8  filter;      /**< filter for integrating shared nogoods */
		bool   copyProblem; /**< copy problem - default=false: share problem */
	} distribute;         /**< nogood distribution parameters */
	struct GRestarts {    /**< Options for configuring global restarts */
		GRestarts():maxR(0) {}
		uint32           maxR;
		ScheduleStrategy sched;
	} restarts;           /**< global restart strategy */
	std::string portfolio;/**< portfolio file */
	bool   forceGP;        /**< force guiding path scheme */
	bool   genTemplate;    /**< generate a template portfolio file and exit */
};


//! Global options controlling global solving algorithms
struct GlobalOptions {
public:
	typedef ProgramBuilder::EqOptions        EqOptions;
	typedef Enumerator::ProgressReport       Progress;
	typedef std::auto_ptr<Enumerator> EnumPtr;
	GlobalOptions();
	Enumerator*     initEnumerator(Enumerator::Report* r = 0);
	bool consequences() const { return enumerate.brave || enumerate.cautious; }
	const char* cbType()const { return consequences() ? (enumerate.brave ? "Brave" : "Cautious") : "none"; }
	SharedContext ctx;    /**< Context-object used by all solvers */
	ThreadOptions thread; /**< Options for parallel solving */
	EqOptions     eq;     /**< options for equivalence preprocessing */
	struct Optimize {     /**< Optimization options */
		Optimize() : hierarch(0), heu(0), no(false), all(false), sat(false) {}
		WeightVec vals;     /**< initial values for optimize statements */
		int    hierarch;    /**< use hierarchical optimization scheme */
		int    heu;         /**< consider optimize statements in heuristics */
		bool   no;          /**< ignore optimize statements */
		bool   all;         /**< compute all models <= vals */
		bool   sat;         /**< force max sat */
	}    opt;
	struct EnumOptions {  /**< Enumeration options */
		EnumOptions() : progress(0), numModels(-1), projectOpts(7), project(false)
		              , record(false), restartOnModel(false), brave(false)
		              , cautious(false), onlyPre(false)  {}
		Progress* progress; /**< enable progress reporting? */
		int  numModels;     /**< number of models to compute */
		int  projectOpts;   /**< options for projection */
		bool project;       /**< enable projection */
		bool record;        /**< enable solution recording */
		bool restartOnModel;/**< restart after each model */
		bool brave;         /**< enable brave reasoning */
		bool cautious;      /**< enable cautious reasoning */
		bool onlyPre;       /**< stop after preprocessing step? */
	} enumerate;
};

class ClaspFacade;

//! Interface for controling incremental solving
class IncrementalControl {
public:
	IncrementalControl();
	virtual ~IncrementalControl(); 
	//! Called before an incremental step is started
	virtual void initStep(ClaspFacade& f)  = 0;
	//! Called after an incremental step finished
	/*!
	 * \return
	 *  - true to signal that solving should continue with next step
	 *  - false to terminate the incremental solving loop
	 */
	virtual bool nextStep(ClaspFacade& f)  = 0;
private:
	IncrementalControl(const IncrementalControl&);
	IncrementalControl& operator=(const IncrementalControl&);
};

//! Parameter-object that groups & validates options
class ClaspConfig : public GlobalOptions {
public:
	ClaspConfig();
	~ClaspConfig();
	uint32        numThreads() const;
	LocalOptions* threadConfig(uint32 i) const;
	LocalOptions* master() const { return threadConfig(0); }
	bool          validate(std::string& err);
	void          applyHeuristic();
	void          reset();
	void          setThreads(uint32 i);
private:
	ClaspConfig(const ClaspConfig&);
	ClaspConfig& operator=(const ClaspConfig&);
	typedef PodVector<LocalOptions*>::type SolverOptions;
	SolverOptions          solvers;
};
/////////////////////////////////////////////////////////////////////////////////////////
// ClaspFacade
/////////////////////////////////////////////////////////////////////////////////////////
//! Provides a simplified interface for (incrementally) solving a given problem
class ClaspFacade : public Enumerator::Report {
public:
	//! Defines the possible solving states
	enum State { 
		state_not_started,
		state_read,        /*!< problem is read from input */
		state_preprocess,  /*!< problem is prepared */
		state_solve,       /*!< search is active */
		num_states
	};
	//! Defines important event types
	enum Event { 
		event_state_enter, /*!< a new state was entered */
		event_state_exit,  /*!< about to exit from the active state */
		event_p_prepared,  /*!< problem was transformed to nogoods */
		event_model        /*!< a model was found */
	};
	//! Defines possible solving results
	enum Result { result_unsat, result_sat, result_unknown }; 
	//! Callback interface for notifying about important steps in solve process
	class Callback {
	public:
		virtual ~Callback() {}
		//! State transition. Called on entering/exiting a state
		/*!
		 * \param e is either event_state_enter or event_state_exit
		 * \note f.state() returns the active state
		 */
		virtual void state(Event e, ClaspFacade& f) = 0;
		//! Some operation triggered an important event
		/*!
		 * \param s the solver that triggered the event
		 * \param e an event that is neither event_state_enter nor event_state_exit 
		 */
		virtual void event(const Solver& s, Event e, ClaspFacade& f) = 0;
		//! Some configuration option is unsafe/unreasonable w.r.t the current problem
		virtual void warning(const char* msg)       = 0;
	};
	ClaspFacade();
	/*!
	 * Returns the number of hardware threads available on this system
	 * or 1 if this information is not available.
	 * 
	 * \note The number of hardware threads typically depends on the number of CPUs or cores 
	 * or hyperthreading units. 
	 */
	static uint32 hardware_concurrency();

	/*!
	 * Solves the problem given in problem using the given configuration.
	 * \pre config is valid, i.e. config.valid() returned true
	 * \note Once solve() returned, the result of the computation can be
	 * queried via the function result().
	 * \note if config.onlyPre is true, solve() returns after
	 * the preprocessing step (i.e. once the solver is prepared) and does not start a search.
	 */
	void solve(Input& problem, ClaspConfig& config, Callback* c);

	/*!
	 * Incrementally solves the problem given in problem using the given configuration.
	 * \pre config is valid, i.e. config.valid() returned true
	 * \note Call result() to get the result of the computation
	 * \note config.onlyPre is ignored in incremental setting!
	 *
	 * solveIncremental() runs a simple loop that is controlled by the
	 * given IncrementalControl object inc.
	 * \code
	 * do {
	 *   inc.initStep(*this);
	 *   read();
	 *   preprocess();
	 *   solve();
	 * } while (inc.nextStep(*this));
	 * \endcode
	 * 
	 */
	void solveIncremental(Input& problem, ClaspConfig& config, IncrementalControl& inc, Callback* c);

	//! returns the result of a computation
	Result result() const { return result_; }
	//! returns true if search-space was completed. Otherwise false.
	bool   more()   const { return more_; }
	//! returns the active state
	State  state()  const { return state_; }
	//! returns the current incremental step (starts at 0)
	int    step()   const { return step_; }
	//! returns the current input problem
	Input* input() const { return input_; }
	//! tries to terminate an active search
	bool   terminate() const {  return ctrl_ && ctrl_->terminate(); }
	
	const ClaspConfig* config() const { return config_; }

	//! returns the ProgramBuilder-object that was used to transform a logic program into nogoods
	/*!
	 * \note A ProgramBuilder-object is only created if input()->format() == Input::SMODELS
	 * \note The ProgramBuilder-object is destroyed after the event
	 *       event_p_prepared was fired. Call releaseApi to disable auto-deletion of api.
	 *       In that case you must later manually delete it!
	 */
	ProgramBuilder* api() const  { return api_.get();     }
	ProgramBuilder* releaseApi() { return api_.release(); }

	void warning(const char* w) const { if (cb_) cb_->warning(w); }
private:
	ClaspFacade(const ClaspFacade&);
	ClaspFacade& operator=(const ClaspFacade&);
	struct AutoState {
		AutoState(ClaspFacade* f, State s) : self_(f) { self_->setState(s, event_state_enter); }
		~AutoState() { self_->setState(self_->state(), event_state_exit); delete self_->ctrl_; self_->ctrl_ = 0; }
		ClaspFacade* self_;
	};
	typedef SingleOwnerPtr<ProgramBuilder> Api;
	typedef SingleOwnerPtr<SharedDependencyGraph> GraphPtr;
	// -------------------------------------------------------------------------------------------  
	// Status information
	void setState(State s, Event e)          { state_ = s; if (cb_) cb_->state(e, *this); }
	void fireEvent(const Solver& s, Event e) { if (cb_) cb_->event(s, e, *this); }
	// -------------------------------------------------------------------------------------------
	// Enumerator::Report interface
	void reportModel(const Solver& s, const Enumerator&) {
		result_ = result_sat;
		fireEvent(s, event_model);
	}
	void reportSolution(const Enumerator& e, bool complete) {
		more_ = !complete;
		if (!more_ && e.enumerated == 0) {
			result_ = result_unsat;
		}
	}
	// -------------------------------------------------------------------------------------------
	// Internal setup functions
	void   validateWeak();
	void   validateWeak(ClaspConfig& cfg);
	void   init(Input&, ClaspConfig&, IncrementalControl*, Callback* c);
	bool   read();
	bool   preprocess();
	void   setProblemSize() const;
	bool   configureMinimize(SharedMinimizeData* min) const;
	bool   initEnumeration(SharedMinimizeData* min);
	void   setGraph();
	void   initSolveObject(ClaspConfig& config);
	// -------------------------------------------------------------------------------------------
	ClaspConfig*           config_;
	IncrementalControl*    inc_;
	Callback*              cb_;
	Input*                 input_;
	SolveAlgorithm*        ctrl_;
	GraphPtr               graph_;
	Api                    api_;
	Result                 result_;
	State                  state_;
	int                    step_;
	bool                   more_;
};

}
#endif
