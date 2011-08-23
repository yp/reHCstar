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
#ifndef DISABLE_MULTI_THREADING
#include <clasp/parallel_solve.h>
#include <clasp/solver.h>
#include <clasp/clause.h>
#include <clasp/clasp_facade.h>
#include <clasp/util/timer.h>
#include <clasp/minimize_constraint.h>
#include <tbb/concurrent_queue.h>

#if !defined(NDEBUG) && !defined(ENABLE_LOG)
#define ENABLE_LOG 1
#endif

#if defined(ENABLE_LOG) && ENABLE_LOG == 1
#include <stdio.h>
#define LOG(id, fmt, ...) ( fprintf(stderr, "THREAD[%u]: " fmt "\n", (id), ##__VA_ARGS__), fflush(stderr) )
#else
#define LOG(id, fmt, ...)
#endif

namespace Clasp { namespace mt {
/////////////////////////////////////////////////////////////////////////////////////////
// BarrierSemaphore
/////////////////////////////////////////////////////////////////////////////////////////
// A combination of a barrier and a semaphore
class BarrierSemaphore {
public:
	explicit BarrierSemaphore(int counter = 0, int maxParties = 1) : counter_(counter), active_(maxParties) {}
	// Initializes this object
	// PRE: no thread is blocked on the semaphore
	//      (i.e. internal counter is >= 0)
	// NOTE: not thread-safe
	void unsafe_init(int counter = 0, int maxParties = 1) {
		counter_ = counter;
		active_  = maxParties;
	}
	// Returns the current semaphore counter.
	int  counter()   { std::unique_lock<tbb::mutex> lock(semMutex_); return counter_; }
	// Returns the number of parties required to trip this barrier.
	int  parties()   { std::unique_lock<tbb::mutex> lock(semMutex_); return active_;  } 
	// Returns true if all parties are waiting at the barrier
	bool active()    { std::unique_lock<tbb::mutex> lock(semMutex_); return unsafe_active(); }
	
	// barrier interface
	
	// Increases the barrier count, i.e. the number of 
	// parties required to trip this barrier.
	void addParty() {
		std::unique_lock<tbb::mutex> lock(semMutex_);
		++active_;
	}
	// Decreases the barrier count and resets the barrier
	// if reset is true. 
	// PRE: the thread does not itself wait on the barrier
	void removeParty(bool reset) {
		std::unique_lock<tbb::mutex> lock(semMutex_);
		assert(active_ > 0);
		--active_;
		if      (reset)           { unsafe_reset(0); }
		else if (unsafe_active()) { counter_ = -active_; lock.unlock(); semCond_.notify_one(); }
	}
	// Waits until all parties have arrived, i.e. called wait.
	// Exactly one of the parties will receive a return value of true,
	// the others will receive a value of false.
	// Applications shall use this value to designate one thread as the
	// leader that will eventually reset the barrier thereby unblocking the other threads.
	bool wait() {
		std::unique_lock<tbb::mutex> lock(semMutex_);
		if (--counter_ >= 0) { counter_ = -1; }
		return unsafe_wait(lock);
	}
	// Resets the barrier and unblocks any waiting threads.
	void reset(uint32 semCount = 0) {
		std::unique_lock<tbb::mutex> lock(semMutex_);
		unsafe_reset(semCount);
	}
	// semaphore interface
	
	// Decrement the semaphore's counter.
	// If the counter is zero or less prior to the call
	// the calling thread is suspended.
	// Returns false to signal that all but the calling thread
	// are currently blocked.
	bool down() {
		std::unique_lock<tbb::mutex> lock(semMutex_);
		if (--counter_ >= 0) { return true; }
		return !unsafe_wait(lock);
	}
	// Increments the semaphore's counter and resumes 
	// one thread which has called down() if the counter 
	// was less than zero prior to the call.
	void up() {
		bool notify;
		{
			std::unique_lock<tbb::mutex> lock(semMutex_);
			notify    = (++counter_ < 1);
		}
		if (notify) semCond_.notify_one();
	}
private:
	BarrierSemaphore(const BarrierSemaphore&);
	BarrierSemaphore& operator=(const BarrierSemaphore&);
	typedef std::condition_variable cv;
	bool    unsafe_active() const { return -counter_ >= active_; }
	void    unsafe_reset(uint32 semCount) {
		int prev = counter_;
		counter_ = semCount;
		if (prev < 0) { semCond_.notify_all(); }
	}
	// Returns true for the leader, else false
	bool    unsafe_wait(std::unique_lock<tbb::mutex>& lock) {
		assert(counter_ < 0);
		// don't put the last thread to sleep!
		if (!unsafe_active()) {
			semCond_.wait(lock);
		}
		return unsafe_active();
	}	
	cv         semCond_;  // waiting threads
	tbb::mutex semMutex_; // mutex for updating counter
	int        counter_;  // semaphore's counter
	int        active_;   // number of active threads
};
/////////////////////////////////////////////////////////////////////////////////////////
// ParallelSolve::Impl
/////////////////////////////////////////////////////////////////////////////////////////
struct ParallelSolve::SharedData {
	typedef tbb::concurrent_queue<const LitVec*> queue;
	SharedData() : path(0) { reset(0, 1); }
	void reset(SharedContext* a_ctx, uint32 numT) {
		clearQueue();
		init.reset();
		workSem.unsafe_init(0, numT);
		ctx    = a_ctx;
		path   = 0;
		nextId = 1;
	}
	void clearQueue() {
		for (const LitVec* a = 0; workQ.try_pop(a); ) {
			if (a != path) { delete a; }
		}
	}
	Enumerator* enumerator() const { return ctx->enumerator(); }
	SharedContext*    ctx;     // shared context object
	const LitVec*     path;    // initial guiding path - typically empty
	Timer<RealTime>   init;    // thread init time
	tbb::mutex        modelM;  // model-mutex 
	BarrierSemaphore  workSem; // work-semaphore
	queue             workQ;   // work-queue
	uint32            nextId;  // next solver id to use
};

// Callable function object to be used as "main" function
// for threads. Simply forwards to ParallelSolve::doSolve()
struct ParallelSolve::EntryPoint {
	explicit EntryPoint(ParallelSolve* self) : self_(self) {}
	void operator()(Solver* s, const SolveParams* p) { self_->doSolve(*s, *p, LitVec()); }
	ParallelSolve* self_;
};
/////////////////////////////////////////////////////////////////////////////////////////
// ParallelSolve
/////////////////////////////////////////////////////////////////////////////////////////
ParallelSolve::ThreadControlFlags::ThreadControlFlags() {
	reset();
}
void ParallelSolve::ThreadControlFlags::reset() {
	workReq    = 0;
	restartReq = 0;
	control    = 0;
	messages   = 0;
	updates    = 0;
}

// sets flag flag if not yet set
// return: true if flag was set, false if flag is already set
bool ParallelSolve::ThreadControlFlags::setFlag(Flag flag) {
	for (uint32 x = control, y = flag; (x & y) == 0;) {
		if (control.compare_and_swap(x|y, x) == x) {
			return true;
		}
		x = control;
	}
	return false;
}

ParallelSolve::ParallelSolve()
	: SolveAlgorithm()
	, shared_(new SharedData)
	, thread_(0)
	, globalR_(0,0.0,0)
	, maxThreads_(0)
	, maxRestarts_(0)
	, intGrace_(1024)
	, intFlags_(ClauseCreator::clause_not_root_sat | ClauseCreator::clause_no_add)
	, intHeuristic_(false)
	, forceGP_(true) {
}

ParallelSolve::~ParallelSolve() {
	if (maxThreads_ > 1) {
		// algorithm was not started but there may be active threads -
		// force orderly shutdown
		ParallelSolve::terminate();
		shared_->workSem.removeParty(true);
		joinThreads();
		assert(numThreads() == 0);
	}
	destroyThread(0);
	delete shared_;
}

uint32 ParallelSolve::hardwareThreads() {
	return std::thread::hardware_concurrency();
}

void ParallelSolve::setForceGP(bool gp) {
	forceGP_ = gp;
}

void ParallelSolve::setIntegrate(uint32 grace, uint8 filter) {
	typedef ThreadOptions::Distribution Dist;
	intGrace_     = grace;
	intFlags_     = ClauseCreator::clause_not_root_sat | ClauseCreator::clause_no_add;
	intHeuristic_ = filter == Dist::filter_heuristic;
	if (filter == Dist::filter_sat) {
		intFlags_ |= (ClauseCreator::clause_not_root_sat | ClauseCreator::clause_not_sat);
	}
	else if (filter == Dist::filter_gp) {
		intFlags_ |= ClauseCreator::clause_not_root_sat;
	}
}

void ParallelSolve::setRestarts(uint32 maxR, const ScheduleStrategy& rs) {
	maxRestarts_ = maxR;
	globalR_     = rs;
}

// prepare parallel solving
SharedContext* ParallelSolve::init(int numT, SharedContext& ctx) {
	assert(numT >= 1 && "Illegal number of threads");
	shared_->reset(&ctx, static_cast<uint32>(numT));
	msg_.reset();
	reserveThreads();
	shared_->init.start();
	return &ctx;
}

uint32 ParallelSolve::numThreads() const {
	return shared_->workSem.parties();
}

void ParallelSolve::reserveThreads() {
	destroyThread(0);
	uint32 numT = numThreads();
	thread_     = new ParallelHandler*[numT];
	maxThreads_ = 0;
	for (uint32 i = 0; i != numT; ++i) {
#pragma message TODO("replace with CACHE_LINE_ALIGNED alloc")
		uint32 b  = ((sizeof(ParallelHandler)+63) / 64) * 64;
		thread_[i]= new (::operator new( b )) ParallelHandler(*this);
		++maxThreads_;
	}
}

void ParallelSolve::destroyThread(uint32 id) {
	assert(id <= maxThreads_);
	if (id < maxThreads_) {
		assert(!thread_[id]->joinable() && "Shutdown not completed!");
		thread_[id]->~ParallelHandler();
		::operator delete(thread_[id]);
		if (id == 0) {
			assert(maxThreads_ == 1 && "Shutdown not completed!");
			delete [] thread_;
			thread_     = 0;
			maxThreads_ = 0;
		}
	}
}

// joins with and destroys all active threads
void ParallelSolve::joinThreads() {
	shared_->init.reset();
	shared_->init.start();
	LOG(0, "STARTING  shutdown");
	int error = thread_[0]->error();
	for (uint32 i = 1, end = maxThreads_; i != end; ++i) {
		if (thread_[i]->join() > error) {
			error = thread_[i]->error();
		}
		destroyThread(i);
	}
	thread_[0]->setError(msg_.ok() ? thread_[0]->error() : error);
	maxThreads_ = 1;
	shared_->init.stop();
	LOG(0, "COMPLETED shutdown after %.3fs", shared_->init.total());
}

// attach new client thread
void ParallelSolve::addThread(Solver& s, SolveParams& p) {
	uint32 id = shared_->nextId++;
	assert(id < maxThreads_ && (s.id() == 0 || s.id() == id));
	assert(thread_ && thread_[id]  && "ParallelSolve::init() not called!");
	assert(!thread_[id]->joinable()&& "IDs are not unique or thread already initialized!");
	s.setId(id);
	thread_[id]->run( EntryPoint(this), s, p ); // start in a new thread
}

// called by all threads once before solving begins
// establishes solver<->handler connection and attaches client
// threads to shared context
void ParallelSolve::startSolve(Solver& s, const SolveParams& p, const LitVec& assume) {
	uint32 id = s.id();
	LOG(id, "ENTERING  algorithm");
	if (!msg_.terminate()) {
		if (shared_->ctx->master()->stats.jumps) {
			s.stats.enableJumpStats();
		}
		thread_[id]->attach(s, p); // connect solver<->handler
		if (id == 0) {
			// explicity init parallel handler because Solver::endInit() was
			// already called.
			thread_[0]->init(s);
			shared_->path = &assume;
		}
		else if (!shared_->ctx->attach(s) && !msg_.quit()) {
			// problem is UNSAT - terminate algorithm
			LOG(s.id(), "TRIVIALLY UNSAT!");
			terminateComplete(true);
			return;
		}
		// wait for other threads
		if (shared_->workSem.wait()) {
			assert(msg_.workReq == 0 && shared_->workQ.empty());
			initQueue(); // add initial path to work queue
			shared_->init.stop();
			LOG(id, "COMPLETED initialization after %.3fs", shared_->init.total());
			shared_->workSem.reset();
		}
	}
}

// main solve loop executed by all threads
bool ParallelSolve::doSolve(Solver& s, const SolveParams& p, const LitVec& assume) {
	PathPtr a(0);
	try {
		startSolve(s, p, assume);
		MinimizeConstraint* min = 0;
		for (; requestWork(s, a);) {
			thread_[s.id()]->prepareForGP(*a, globalR_.current(), a.is_owner());
			s.stats.reset();
			if (initPath(s, *a) && (!min || min->integrateNext(s))) {
				if (solvePath(s, p) == value_free) {
					terminateComplete(false);
					s.stats.parallel->terminated = true;
				}
				s.clearStopConflict();
			}
			if (min || shared_->enumerator()->optimizeHierarchical()) {
				// path is unsat and hierarchical optimization is active;
				// relax minimize constraint to avoid problems during simplification
				min  = shared_->enumerator()->constraint(s)->minimize();
				min->restoreOptimum();
			}
		}
	}
	catch (const std::bad_alloc&)      { exception(s,a,error_oom, "std::bad_alloc" ); }
	catch (const std::runtime_error& e){ exception(s,a,error_runtime, e.what()); }
	catch (...)                        { exception(s,a,error_other, "unknown");  }
	return endSolve(s);
}

// called once by each solver just after solving process has stopped
// detaches the solver<->handler connection and joins with
// all client threads
bool ParallelSolve::endSolve(Solver& s) {
	assert(msg_.terminate() || thread_[s.id()]->error() != error_none);
	LOG(s.id(), "LEAVING   algorithm");
	uint32 id    = s.id();
	int    error = thread_[id]->error();
	thread_[id]->detach();
	if (!msg_.quit() && s.sharedContext() && error == error_none) {
		s.clearAssumptions();
		if (s.hasConflict()) {
			msg_.setFlag(ThreadControlFlags::complete_flag);
		}
	}
	// this thread is leaving
	shared_->workSem.removeParty(msg_.terminate());
	s.stats.reset();
	if (error != error_none) {
		// solver terminated because of an exception
		s.reset();
		s.setId(ParallelSolve::invalidId);
		std::swap(s.stats.jumps, thread_[id]->aggStats.jumps);
		std::swap(s.stats.parallel, thread_[id]->aggStats.parallel);
	}
	// copy aggregated statistics for this solver to solver
	s.stats.accu(thread_[id]->aggStats);
	id == 0
		? joinThreads()            // join with the other solvers
		: shared_->ctx->detach(s); // detach client threads from shared context object
	if (s.stats.parallel) {
		s.stats.parallel->cpuTime = ThreadTime::getTime();
	}
	if (id == 0) {
		error = thread_[0]->error();
		if (error != error_none) {
			switch (error) {
				case error_oom    : throw std::bad_alloc();
				case error_runtime: throw std::runtime_error("RUNTIME ERROR!");
				default:            throw std::runtime_error("UNKNOWN ERROR!");
			}
		}
		shared_->enumerator()->reportResult(msg_.complete());
	}
	return !msg_.complete();
}

void ParallelSolve::exception(Solver& s, PathPtr& path, ErrorCode e, const char* what) {
	try {
		LOG(s.id(), "EXCEPTION %s", what); (void) what;
		thread_[s.id()]->setError(e);
		if (s.id() == 0 || shared_->workSem.active()) { 
			terminate();
			return;
		}
		else if (msg_.allowSplit() && path.get()) {
			shared_->workQ.push(path.release());
			shared_->workSem.up();
		}
	}
	catch(...) { terminate(); }
}

// forced termination from outside
bool ParallelSolve::terminate() {
	// do not use postMessage() to avoid possible
	// deadlock because of unblockAll()!
	msg_.setFlag(ThreadControlFlags::quit_flag);
	++msg_.messages;
	return true;
}

// terminated from inside of algorithm
void ParallelSolve::terminateComplete(bool complete) {
	postMessage(Message_t::terminate_message);
	if (complete) {
		msg_.setFlag(ThreadControlFlags::complete_flag);
	}
}

// terminated because of unsat, check if there is more to do
void ParallelSolve::terminateUnsat(Solver& s, bool complete) {
	if (!msg_.terminate()) {
		if (complete && shared_->enumerator()->optimizeHierarchical() && s.popRootLevel(s.rootLevel(), false)) {
			// Problem is unsat and hierarchical optimization is active.
			// The active level is at its optimum, but lower-priority levels might
			// still be non-optimal.
			// Notify other threads to prepare for solving next priority level.
			LOG(s.id(), "COMPLETED current optimization step"); 
			postMessage(Message_t::unsat_message);
			assert(s.decisionLevel() == 0);
			return;
		}
		// nothing more to do - post terminate message
		terminateComplete(complete);
	}
}

// post message to all threads
void ParallelSolve::postMessage(MessageType m) {
	if (m < Message_t::not_a_control_message) {
		// control message - notify all if new
		ThreadControlFlags::Flag cf = static_cast<ThreadControlFlags::Flag>((uint32(1)<<m));
		if (msg_.setFlag(cf)) {
			++msg_.messages;
			shared_->workSem.reset();
		}
	}
	else if (m == Message_t::split_message) {
		++msg_.workReq;
		++msg_.messages;
	}
	else if (m == Message_t::update_message) {
		++msg_.updates;
		++msg_.messages;
	}
	else {
		assert("ERROR: Unknown message type!\n");
	}
}


// tries to get new work for the given solver
bool ParallelSolve::requestWork(Solver& s, PathPtr& out) { 
	const LitVec* a;
	while (!msg_.terminate()) {
		// only clear path and stop conflict - we don't propagate() here
		// because we would then have to handle any eventual conflicts
		if (!s.popRootLevel(s.rootLevel(), false)) {
			// s has a real top-level conflict - 
			// there is nothing we can do to continue
			terminateComplete(true); 
		}	
		else if (shared_->workQ.try_pop(a)) {
			assert(s.decisionLevel() == 0);
			// got new work from work-queue
			out = a;
			// do not take over ownership of initial gp!
			if (a == shared_->path) { out.release(); } 
			// propagate any new facts now that we have new work
			if (s.simplify())      { return true; }
			// s now has a conflict - either an artifical stop conflict
			// or a real conflict - we'll handle it in the next iteration
			// via the call to popRootLevel()
		}
		else if (msg_.restart()) {
			// a restart request is active - we are fine with
			// this but did not yet had a chance to react on it
			waitOnRestart(s);
		}
		else if (msg_.allowSplit()) {
			// gp mode is active - request a split	
			// and wait until someone has work for us
			postMessage(Message_t::split_message);
			if (!shared_->workSem.down() && !msg_.restart()) {
				// we are the last man standing, there is no
				// work left - quitting time?
				terminateUnsat(s, true);
			}
		}
		else {
			// portfolio mode is active - no splitting, no work left
			// quitting time? 
			terminateUnsat(s, true);
		}
	}
	return false;
}

// handles an active restart request
// returns true to signal that s should restart otherwise false
bool ParallelSolve::waitOnRestart(Solver& s) {
	assert(msg_.restart());
	// react once
	if (!thread_[s.id()]->handleRestartMessage()) {
		msg_.setFlag(ThreadControlFlags::inhibit_restart_flag);
	}
	LOG(s.id(), "WAITING   for restart");
	if (shared_->workSem.wait()) {
		// last man standing - complete restart request
		msg_.workReq   = 0;
		msg_.restartReq= 0;
		bool hasUnsat  = msg_.unsatRestart();
		msg_.control  -= (msg_.control & (uint32(ThreadControlFlags::unsat_flag|ThreadControlFlags::restart_flag)));
		if (!hasUnsat) {
			globalR_.next();
			if ( (maxRestarts_ -= maxRestarts_ > 0) == 0 ) {
				globalR_.init(0, 0.0, 0);
			}
		}
		if (hasUnsat || msg_.allowRestart()) {
			initQueue();
		}
		if (hasUnsat && shared_->enumerator()->optimizeHierarchical()) {
			Enumerator*         en = shared_->enumerator();
			MinimizeConstraint* m  = en->constraint(s)->minimize();
			if (en->enumerated == 0 || !m->shared_unsafe()->optimizeNext()) {
				terminateComplete(true);
			}
		}
		LOG(s.id(), "COMPLETED restart");
		// wake up all blocked threads
		shared_->workSem.reset();
	}
	return msg_.terminate() || (!msg_.restart() && msg_.allowRestart());
}

// If guiding path scheme is active only one
// thread can start with gp (typically empty) and this
// thread must then split up search-space dynamically.
// Otherwise, all threads can start with initial gp
// TODO:
//  heuristic for initial splits?
void ParallelSolve::initQueue() {
	shared_->clearQueue();
	int end = numThreads();
	if (forceGP_) {
		msg_.setFlag(ThreadControlFlags::allow_gp_flag);
		end = 1;
	}
	else {
		msg_.setFlag(ThreadControlFlags::inhibit_restart_flag);
	}
	for (int i = 0; i != end; ++i) {
		shared_->workQ.push(shared_->path);
	}
}

// adds work to the work-queue
void ParallelSolve::pushWork(LitVec& work) { 
	LitVec* v = new LitVec;
	v->swap(work);
	shared_->workQ.push(v);
	shared_->workSem.up();
}

// called whenever some solver has found a model
bool ParallelSolve::backtrackFromModel(Solver& s) { 
	Enumerator::Result r;
	{
		// grab lock - models must be processed sequentially
		// in order to simplify printing and to avoid duplicates
		// in all non-trivial enumeration modes
		std::unique_lock<tbb::mutex> lock(shared_->modelM);
		// first check if the model is still valid once all
		// information is integrated into the solver
		uint32 dl = s.decisionLevel();
		if (!thread_[s.id()]->isModel(s)) {
			// model no longer a (unique) model - continue search
			s.stats.removeModel(dl);
			return true;
		}
		r = shared_->enumerator()->backtrackFromModel(s);
		if (r == Enumerator::enumerate_stop_enough || (r == Enumerator::enumerate_stop_complete && s.decisionLevel() == 0)) {
			// must be called while holding the lock - otherwise
			// we have a race condition with solvers that
			// are currently blocking on the mutex and we could enumerate 
			// more models than requested by the user
			terminateComplete(s.decisionLevel() == 0);
		}
		else if (!shared_->enumerator()->trivial(thread_[s.id()]->disjointPath())) {
			// enumerator is not trivial w.r.t current search scheme
			// force update of other solvers
			postMessage(Message_t::update_message);
			// ignore message in the active solver
			thread_[s.id()]->setUpdate(msg_.updates);
		}
		if (shared_->enumerator()->enumerated == 1 && !shared_->enumerator()->supportsRestarts()) {
			msg_.setFlag(ThreadControlFlags::inhibit_restart_flag);
			msg_.setFlag(ThreadControlFlags::allow_gp_flag);
		}
	}
	return r == Enumerator::enumerate_continue && !msg_.terminate();
}

// checks for new update messages
// POST: m.doUpdate if there is an active update message for m
// return: m.doUpdate
bool ParallelSolve::checkUpdates(Message_t::Local& m) const {
	uint32 updates     = msg_.updates;
	uint32 seen        = m.update;
	m.update           = updates;
	return (m.doUpdate = (m.doUpdate || seen != updates));
}

// updates s with new messages and uses s to create a new guiding path
// if necessary and possible
bool ParallelSolve::handleMessages(Solver& s, Message_t::Local& m) {
	// check if there are new messages for s
	uint32 hId         = s.id();
	uint32 gt          = msg_.messages;
	ParallelHandler* h = thread_[hId];
	if (m.any != gt) {
		if (msg_.hasControlMsg()) {
			if (msg_.terminate()) {
				LOG(s.id(), "TERMINATE message received!");
				h->handleTerminateMessage();
				s.stats.parallel->terminated = true;
				s.setStopConflict();
				return false;
			}
			if (msg_.restart() && waitOnRestart(s)) {
				s.setStopConflict();
				return false;
			}
			return true;
		}
		checkUpdates(m);
		uint32 mt  = msg_.workReq; 
		if (mt <= 0) {
			m.any    = gt;
		}
		else if (h->splittable()) {
			// First declare split request as handled
			// and only then do the actual split.
			// This way, we minimize the chance for 
			// "over"-splitting, i.e. one split request handled
			// by more than one thread.
			if (--msg_.workReq == 0) {
				m.any  = gt;
			}
			h->handleSplitMessage();
		}
	}
	return true;
}

void ParallelSolve::requestRestart() {
	if (msg_.allowRestart() && ++msg_.restartReq == numThreads()) {
		postMessage(Message_t::restart_message);
	}
}

////////////////////////////////////////////////////////////////////////////////////
// ParallelHandler
/////////////////////////////////////////////////////////////////////////////////////////
ParallelHandler::ParallelHandler(ParallelSolve& ctrl) 
	: solver_(0)
	, params_(0)
	, intTail_(0)
	, error_(0)
	, messageHandler_(&ctrl) {
}

ParallelHandler::~ParallelHandler() {
	detach();
}

// adds this as post propagator to s
void ParallelHandler::attach(Solver& s, const SolveParams& p) {
	detach();
	aggStats.reset();
	s.stats.enableParallelStats();
	aggStats.enableParallelStats();
	if (s.stats.jumps) { aggStats.enableJumpStats(); }
	messageHandler_.msg = MessageHandler::MessageQ();
	gp_.reset();
	intTail_= 0;
	error_  = 0;
	s.addPost(&messageHandler_);
	s.addPost(this);
	solver_ = &s;
	params_ = &p;
}

// detaches this from its solver
void ParallelHandler::detach() {
	if (solver_) {
		handleTerminateMessage();
		aggStats.accu(solver_->stats);
		while (!integrated_.empty()) {
			ClauseHead* c = (ClauseHead*)integrated_.back();
			integrated_.pop_back();
			if (error() == 0) {
				if (c->locked(*solver_)) {
					solver_->addLearnt(c, c->size());
				}
				else {
					c->destroy(solver_, true);
				}
			}
			else { c->destroy(0, false); }
		}
		integrated_.clear();
		intTail_= 0;
		solver_ = 0;
		params_ = 0;
	}
}

// true if solver can be used to split-off a new guiding path
bool ParallelHandler::splittable() const {
	return solver_->splittable() && !messageHandler_.pendingUpdate();
}

void ParallelHandler::prepareForGP(const LitVec& out, uint64 restart, bool fromSplit) {
	gp_.reset(restart, fromSplit);
	aggStats.parallel->newGP(out.size());
	aggStats.accu(solver_->stats);
}

// detach from solver, i.e. ignore any further messages 
void ParallelHandler::handleTerminateMessage() {
	if (solver_ && this->next != this) {
		// mark removed propagators by creating "self-loop"
		solver_->removePost(&messageHandler_);
		messageHandler_.next = &messageHandler_;
		solver_->removePost(this);
		this->next = this;
	}
}

// split-off new guiding path and add it to solve object
void ParallelHandler::handleSplitMessage() {
	assert(solver_ && "ParallelHandler::handleSplitMessage(): not attached!");
	Solver& s = *solver_;
	s.updateGuidingPath(gp_.path, gp_.pos, gp_.impl);
	LitVec newGP(gp_.path);
	s.pushRootLevel();
	newGP.push_back(~s.decision(s.rootLevel()));
	++s.stats.parallel->splits;
	gp_.split = true;
	ctrl()->pushWork(newGP);
}

bool ParallelHandler::handleRestartMessage() {
	// TODO
	// we may want to implement some heuristic, like
	// computing a local var order. 
	return true;
}

bool ParallelHandler::simplify(Solver& s, bool sh) {
	ClauseDB::size_type i, j, end = integrated_.size();
	for (i = j = 0; i != end; ++i) {
		Constraint* c = integrated_[i];
		if (c->simplify(s, sh)) { 
			c->destroy(&s, false); 
			intTail_ -= (i < intTail_); 
		}
		else                    { 
			integrated_[j++] = c;  
		}
	}
	shrinkVecTo(integrated_, j);
	if (intTail_ > integrated_.size()) intTail_ = integrated_.size();
	return false;
}

// integrate new information from models
bool ParallelHandler::MessageHandler::doUpdate(Solver& s, const GP& gp) {
	// Skip update of model while assumption literal is not yet assigned.
	// This is necessary during hierarchical optimization because otherwise
	// integrating a too strong optimum might irrevocably force the assumption literal
	// which would defeat its purpose.
	if (s.isTrue(s.sharedContext()->tagLiteral())) {
		if (!s.sharedContext()->enumerator()->update(s, gp.split)) {
			return false;
		}
		msg.doUpdate = false;
	}
	return true;	
}

bool ParallelHandler::propagateFixpoint(Solver& s) {
	if (messageHandler_.pendingUpdate() && !messageHandler_.doUpdate(s, gp_)) {
		return false;
	}
	SharedLiterals* temp[5];
	uint32 rec;
	do {
		rec = s.sharedContext()->receive(s, temp, 5);
		if (!integrateClauses(s, temp, rec)) {
			return false;
		}
		if (s.queueSize() > 0 && !s.propagateUntil(this)) {
			return false;
		}
	} while (rec == 5);
	if (s.stats.conflicts >= gp_.restart) {
		ctrl()->requestRestart();
		gp_.restart *= 2;
	}
	return true;
}

// checks whether s still has a model once all 
// information from previously found models were integrated 
bool ParallelHandler::isModel(Solver& s) {
	assert(s.numFreeVars() == 0);
	// either no unprocessed updates or still a model after
	// updates was integrated
	return !ctrl()->checkUpdates(messageHandler_.msg)
		|| (messageHandler_.doUpdate(s, gp_) && s.numFreeVars() == 0);
}

bool ParallelHandler::integrateClauses(Solver& s, SharedLiterals** arr, uint32 num) {
	assert(!s.hasConflict() && &s == solver_);
	ClauseCreator::Result ret(0, true);
	ClauseInfo e; e.setType(Constraint_t::learnt_other);
	uint32 intFlags = ctrl()->integrateFlags();
	for (uint32 i = 0; i != num && ret.ok; ++i) {
		ret = ClauseCreator::integrate(s, arr[i], intFlags, e);
		if (ret.local) {
			add(ret.local);
		}
	}
	return ret.ok;
}

void ParallelHandler::add(ClauseHead* h) {
	if (intTail_ < integrated_.size()) {
		ClauseHead* o = (ClauseHead*)integrated_[intTail_];
		integrated_[intTail_] = h;
		assert(o);
		if (!ctrl()->integrateUseHeuristic() || o->locked(*solver_) || o->activity().score() > 0) {
			solver_->addLearnt(o, o->size());
		}
		else {
			o->destroy(solver_, true);
		}
	}
	else {
		integrated_.push_back(h);
	}
	if (++intTail_ >= ctrl()->integrateGrace()) {
		intTail_ = 0;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// GlobalQueue
/////////////////////////////////////////////////////////////////////////////////////////
GlobalQueue::GlobalQueue(uint32 maxT, uint32 types, uint32 lbd) : Distributor(maxT, types, lbd), queue_(0) {
	assert(maxT < 64);
	queue_     = new Queue(maxT);
	threadId_  = new ThreadInfo[maxT];
	for (uint32 i = 0; i != maxT; ++i) {
		threadId_[i].id = queue_->addThread();
	}
}
GlobalQueue::~GlobalQueue() {
	release();
}
void GlobalQueue::release() {
	if (queue_) {
		for (uint32 i = 0; i != queue_->maxThreads(); ++i) {
			Queue::ThreadId& id = getThreadId(i);
			uint64 mask = uint64(1) << i;
			for (ClauseNode n; queue_->tryConsume(id, n); ) { 
				if ((n.targetMask & mask) != 0) {
					n.lits->release();
				}
			}
		}
		delete queue_;
		queue_ = 0;
		delete [] threadId_;
	}
}
void GlobalQueue::doPublish(const Solver& s, SharedLiterals* lits) {
	assert(lits->refCount() == queue_->maxThreads());
	ClauseNode n;
	n.targetMask = ((uint64(1) << queue_->maxThreads()) - 1) ^ (uint64(1) << s.id());
	n.lits       = lits;
	queue_->publish(n, getThreadId(s.id()));
}
uint32 GlobalQueue::receive(const Solver& in, SharedLiterals** out, uint32 maxn) {
	uint32 r = 0;
	Queue::ThreadId& id = getThreadId(in.id());
	uint64 mask = uint64(1) << in.id();
	for (ClauseNode n; r != maxn && queue_->tryConsume(id, n); ) {
		if ( (n.targetMask & mask) != 0 ) {
			out[r++] = n.lits;
		}
	}
	return r;
}
} } // namespace Clasp::mt
#endif
