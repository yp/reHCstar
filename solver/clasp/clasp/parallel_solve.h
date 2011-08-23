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
#ifndef CLASP_PARALLEL_SOLVE_H_INCLUDED
#define CLASP_PARALLEL_SOLVE_H_INCLUDED

#ifdef _MSC_VER
#pragma once
#endif

#ifndef DISABLE_MULTI_THREADING

#include <clasp/solve_algorithms.h>
#include <clasp/constraint.h>
#include <clasp/shared_context.h>
#include <clasp/util/thread.h>
#include <clasp/util/multi_queue.h>
#include <clasp/solver_types.h>

/*!
 * \file 
 * Defines classes controlling multi-threaded parallel solving
 * 
 */
namespace Clasp { namespace mt {

class ParallelHandler;

//! A parallel algorithm for multi-threaded solving with and without search-space splitting
/*!
 * The class adapts clasp's basic solve algorithm
 * to a parallel solve algorithm that solves
 * a problem using a given number of threads.
 * It supports guiding path based solving, portfolio based solving, as well
 * as a combination of these two approaches.
 */
class ParallelSolve : public SolveAlgorithm {
public:
	ParallelSolve();
	~ParallelSolve();
	// "setup" interface
	
	//! returns number of threads that can run concurrently on the current hardware
	static uint32  hardwareThreads();
	//! used to signal an error in a client thread
	static const   uint32 invalidId = UINT32_MAX;
	
	//! initializes this object for solving with at most maxThreads
	/*!
	 * \pre numThreads >= 1
	 * \pre ctx.shareCount() >= maxThreads
	 */
	SharedContext* init(int maxThreads, SharedContext& ctx);
	//! enables/disables search-space splitting
	void   setForceGP(bool gp);
	//! configures nogood distribution
	void   setIntegrate(uint32 grace, uint8 filter);
	//! configures "global" restarts (only used if search-space splitting is active)
	void   setRestarts(uint32 maxR, const ScheduleStrategy& rs);
	//! must be called once for each solver except the master
	/*!
	 * \pre s.id() != ParallelSolve::invalidId
	 */
	void   addThread(Solver& s, SolveParams& p);
	
	// solve interface 
	
	//! returns the number of active threads
	uint32 numThreads()            const;
	bool   integrateUseHeuristic() const { return intHeuristic_; }
	uint32 integrateGrace()        const { return intGrace_; }
	uint32 integrateFlags()        const { return intFlags_; }
	//! terminates current solving process
	bool   terminate();
	void   requestRestart();
	// interface for controlling the solve process
	struct Message_t {
		enum Type {
			  /* control */ quit_message = 0, terminate_message = 1, unsat_message = 2, restart_message = 3
			, /* other */   not_a_control_message = 4, split_message = 8, update_message = 16
		};
		struct Local {
			Local() : any(0), update(0), doUpdate(false) {}
			uint32 any;
			uint32 update;
			bool   doUpdate;
		};
	};
	typedef Message_t::Type MessageType;	
	bool   handleMessages(Solver& s, Message_t::Local& mq);
	bool   checkUpdates(Message_t::Local& mq) const;
	void   pushWork(LitVec& gp);
	void   postMessage(MessageType);
private:
	ParallelSolve(const ParallelSolve&);
	ParallelSolve& operator=(const ParallelSolve&);
	typedef SingleOwnerPtr<const LitVec> PathPtr;
	enum ErrorCode { error_none = 0, error_oom = 1, error_runtime = 2, error_other = 4 };
	// -------------------------------------------------------------------------------------------
	// Thread setup 
	struct EntryPoint;
	void   reserveThreads();
	void   destroyThread(uint32 id);
	void   joinThreads();
	// -------------------------------------------------------------------------------------------
	// Algorithm steps
	void   startSolve(Solver& s, const SolveParams& p, const LitVec& assume);
	bool   endSolve(Solver& s); 
	bool   doSolve(Solver& s, const SolveParams& p, const LitVec& assume);
	void   initQueue();
	bool   requestWork(Solver& s, PathPtr& out);
	bool   backtrackFromModel(Solver& s);
	void   terminateUnsat(Solver& s, bool complete);
	void   terminateComplete(bool complete);
	bool   waitOnRestart(Solver& s);
	void   exception(Solver& s, PathPtr& path, ErrorCode e, const char* what);
	// -------------------------------------------------------------------------------------------
	struct SharedData;
	SharedData*         shared_;      // Shared control data
	ParallelHandler**   thread_;      // Thread-locl control data
	struct ThreadControlFlags {
	  ThreadControlFlags();
		enum Flag { 
			  quit_flag = 1, terminate_flag = 2, unsat_flag = 4, restart_flag = 8
			, ctrl_message = quit_flag|terminate_flag|unsat_flag|restart_flag
			, inhibit_restart_flag = 16, complete_flag = 32, allow_gp_flag = 64 
		};
		void reset();
		bool quit()         const { return (control & uint32(quit_flag)) != 0; }
		bool terminate()    const { return (control & uint32(quit_flag|terminate_flag)) != 0; }
		bool ok()           const { return (control & uint32(terminate_flag)) != 0; }
		bool restart()      const { return (control & uint32(unsat_flag|restart_flag))  != 0; }
		bool unsatRestart() const { return (control & uint32(unsat_flag)) != 0; }
		bool hasControlMsg()const { return (control & uint32(ctrl_message)) != 0; }
		bool complete()     const { return (control & uint32(complete_flag)) != 0; }
		bool allowRestart() const { return (control & uint32(inhibit_restart_flag)) == 0; }
		bool allowSplit()   const { return (control & uint32(allow_gp_flag)) != 0; }
		bool setFlag(Flag flag);
		std::atomic<int>    workReq;   // > 0: someone needs work
		std::atomic<uint32> restartReq;// == numThreads(): restart
		std::atomic<uint32> control;   // set of active control flags
		std::atomic<uint32> messages;  // global message counter
		std::atomic<uint32> updates;   // global number of "update"  messages
	}                    msg_;
	ScheduleStrategy     globalR_;   // global restart strategy
	uint32               maxThreads_;   // number of threads alloacated 
	uint32               maxRestarts_;  // disable global restarts once reached 
	uint32               intGrace_;     // grace period for clauses to integrate
	uint32               intFlags_;     // bitset controlling clause integration
	bool                 intHeuristic_; // use heuristic in clause integration
	bool                 forceGP_;	    // force guiding path mode even if portfolio is used
};


//! A per-solver (i.e. thread) class that implements message handling and knowledge integration
/*!
 * The class adds itself as a post propagator to the given solver. Each time
 * propagateFixpoint() is called (i.e. on each new decision level), it checks
 * for new lemmas to integrate and synchronizes the search with any new models.
 * Furthermore, it adds a second (high-priority) post propagator for message handling.
 */
class ParallelHandler : public PostPropagator {
public:
	//! creates a new parallel handler to be used in the given solve group
	/*!
	 * \param ctrl The object controlling the parallel solve operation
	 */
	explicit ParallelHandler(ParallelSolve& ctrl);
	~ParallelHandler();

	//! attaches this parallel handler to the given solver
	/*!
	 * \param s The solver in which this object is used
	 * \param p The solving parameters under which the solver operates
	 */
	void attach(Solver& s, const SolveParams& p);

	//! detaches this object from its solver
	void detach();

	//! executes F in a new thread
	template <class F>
	void run(F f, Solver& s, SolveParams& p) {
		assert(!joinable() && solver_ == 0);
		std::thread(f, &s, &p).swap(thread_);
		assert(joinable());
	}
	
	void setError(int e) { error_ = e; }
	int  error() const   { return error_; }
	
	//! true if *this has an associated thread of execution, false otherwise. 
	bool joinable() const { return thread_.joinable(); }
	//! waits for the thread of execution associated with *this to finish.
	/*!
	 * \note the function is a noop of !joinable()
	 */
	int join() { if (joinable()) { thread_.join(); } return error_; }
	
	// overridden methods
	
	//! returns a priority suited for a post propagators that is non-deterministic
	uint32 priority() const { return priority_general + 100; }

	//! integrates new information
	bool propagateFixpoint(Solver& s);
	bool propagate(Solver& s) { return ParallelHandler::propagateFixpoint(s); }
	
	//! checks whether new information has invalidated current model
	bool isModel(Solver& s);

	// own interface
	// TODO: make functions virtual once necessary 
	
	//! returns true if handler's guiding path is disjoint from all others
	bool disjointPath() const { return gp_.split; }
	//! returns true if the handler can currently split off a new guiding path
	bool splittable()   const;
	
	//! called before solver starts to solve given guiding path
	/*!
	 * \param gp      the new guiding path
	 * \param restart request restart after restart number of conflicts
	 * \param isSplit true if gp resulted from a split
	 */
	void prepareForGP(const LitVec& gp, uint64 restart, bool isSplit);

	/*!
	 * \name Message handlers
	 * \note 
	 *   Message handlers are intended as callbacks for ParallelSolve::handleMessages().
	 *   They shall not change the assignment of the solver object.
	 */
	//@{
	
	//! algorithm is about to terminate
	/*!
	 * removes this object from the solver's list of post propagators
	 */
	void handleTerminateMessage();

	//! request for split
	/*!
	 * Splits off a new guiding path and adds it to the control object.
	 * \pre the guiding path of this object is "splittable"
	 */
	void handleSplitMessage();

	//! request for (global) restart
	/*!
	 * \return true if restart is valid, else false
	 */
	bool handleRestartMessage();

	//! mark all update messages <= up as handled
	void setUpdate(uint32 up) { messageHandler_.msg.update = up; }

	SolveStats aggStats;  // aggregated statistics over all gps
	//@}  
private:
	static void threadMain(ParallelHandler* h);
	bool simplify(Solver& s, bool re);
	bool integrateClauses(Solver& s, SharedLiterals** arr, uint32 num);
	void add(ClauseHead* h);
	ParallelSolve* ctrl() const { return messageHandler_.ctrl; }
	typedef LitVec::size_type size_type;
	typedef PodVector<Constraint*>::type ClauseDB;
	std::thread        thread_;     // active thread or empty for master
	ClauseDB           integrated_; // my integrated clauses
	Solver*            solver_;     // my solver
	const SolveParams* params_;     // my solving params
	size_type          intTail_;    // where to put next clause
	int                error_;      // error code or 0 if ok
	struct GP {
		LitVec      path;     // current guiding path
		uint64      restart;  // don't give up before restart number of conflicts
		size_type   pos;      // pos in trail
		uint32      impl;     // number of additional implied literals
		bool        split;    // does gp result from a split?
		void reset(uint64 r = UINT64_MAX, bool sp = false) {
			path.clear();
			restart = r;
			pos     = 0;
			impl    = 0;
			split   = sp;
		}
	} gp_;
	struct MessageHandler : PostPropagator {
		typedef ParallelSolve::Message_t::Local MessageQ;
		explicit MessageHandler(ParallelSolve* c) : ctrl(c) {}
		uint32 priority() const { return PostPropagator::priority_highest; }
		bool   propagateFixpoint(Solver& s) { return ctrl->handleMessages(s, msg); }
		bool   propagate(Solver& s)         { return MessageHandler::propagateFixpoint(s); }
		bool   pendingUpdate() const        { return msg.doUpdate; }
		bool   doUpdate(Solver& s, const GP& gp);
		ParallelSolve* ctrl; // get messages from here
		MessageQ       msg;  // local messages
	}    messageHandler_;
};

class GlobalQueue : public Distributor {
public:
	GlobalQueue(uint32 maxShare, uint32 typesToShare, uint32 maxLbd);
	~GlobalQueue();
	uint32  receive(const Solver& in, SharedLiterals** out, uint32 maxOut);
protected:
	void    doPublish(const Solver& source, SharedLiterals* lits);
private:
	void release();
	struct ClauseNode {
		ClauseNode()
			: targetMask(0), lits(0) {}
		uint64          targetMask;
		SharedLiterals* lits;
	};
	class Queue : public MultiQueue<ClauseNode> {
	public:
		typedef MultiQueue<ClauseNode> base_type;
		using base_type::publish;
		Queue(uint32 m) : base_type(m) {}
	};
	struct ThreadInfo {
		Queue::ThreadId id;
		char pad[64 - sizeof(Queue::ThreadId)];
	};
	Queue::ThreadId& getThreadId(uint32 sId) const {
		return threadId_[sId].id;
	}
	Queue*               queue_;
	ThreadInfo*          threadId_;
};

} }
#endif

#endif
