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
#include <clasp/clasp_facade.h>
#include <clasp/model_enumerators.h>
#include <clasp/cb_enumerator.h>
#include <clasp/weight_constraint.h>
#include <clasp/minimize_constraint.h>
#include <clasp/parallel_solve.h>
namespace Clasp {
/////////////////////////////////////////////////////////////////////////////////////////
// HeuristicOptions
/////////////////////////////////////////////////////////////////////////////////////////
HeuristicOptions::HeuristicOptions() 
	: name("berkmin")
	, lookahead(Lookahead::no_lookahead)
	, lookaheadNum(-1)
	, loops(-1)
	, berkMoms(true)
	, berkHuang(false)
	, berkOnce(false) 
	, nant(false) {
		extra.berkMax = -1;
}

DecisionHeuristic* HeuristicOptions::createHeuristic() const {
	DecisionHeuristic* heu = 0;
	if (name == "berkmin") {
		bool   l = loops == -1 || loops == 1;
		uint32 m = extra.berkMax < 0 ? 0 : extra.berkMax;
		heu = new ClaspBerkmin(m, l, berkMoms, berkHuang, berkOnce);
	}
	else if (name == "vmtf") {
		uint32 m = extra.vmtfMtf < 0 ? 8 : extra.vmtfMtf;
		heu = new ClaspVmtf( m, loops == 1);
	}
	else if (name == "vsids") {
		heu = new ClaspVsids(loops == 1);
	}
	else if (name == "none") {
		heu = new SelectFirst();
	}
	if (lookahead != Lookahead::no_lookahead || lookaheadNum != -1) {
		return new UnitHeuristic(lookahead, nant, heu, lookaheadNum);
	}
	return heu;
}

/////////////////////////////////////////////////////////////////////////////////////////
// LocalOptions
/////////////////////////////////////////////////////////////////////////////////////////
LocalOptions::LocalOptions(Solver* s) : loopRep(DefaultUnfoundedCheck::common_reason), solver_(s), heuristic_(new HeuristicOptions()) {
	assert(s);
}
LocalOptions::~LocalOptions() {
	delete heuristic_;
}
HeuristicOptions& LocalOptions::heuristic() {
	if (heuristic_) return *heuristic_;
	heuristic_ = new HeuristicOptions();
	return *heuristic_;
}

void LocalOptions::applyHeuristic() {
	if (heuristic_) {
		DecisionHeuristic* h = heuristic_->createHeuristic();
		delete heuristic_;
		heuristic_ = 0;		
		solver_->strategies().heuristic.reset(h);
	}
}

bool LocalOptions::validate(std::string& err) {
	if (solver_->strategies().search == SolverStrategies::no_learning) {
		if (heuristic_ && heuristic_->name != "unit" && heuristic_->name != "none") {
			err  = "Heuristic '";
			err += heuristic_->name;
			err += "' requires lookback strategy!";
			return false;
		}
		SolverStrategies* s = &solver_->strategies();
		s->cflMinAntes = SolverStrategies::no_antes;
		s->setCompressionStrategy(0);
		s->saveProgress = 0;
		solve.restart.sched.init(0, 0.0, 0);
		solve.reduce.setStrategy(-1.0, 0.0, 0.0);
		solve.setRandomizeParams(0,0);
		solve.setShuffleParams(0,0);
		solve.restart.local = solve.restart.bounded = solve.reduce.reduceOnRestart = false;
	}
	return true;
}

void LocalOptions::initFrom(const LocalOptions& other) {
	if (other.heuristic_) {
		heuristic()      = *other.heuristic_;
	}
	SolverStrategies* st = &solver_->strategies();
	st->rng          = other.solver_->strategies().rng;
	st->search       = other.solver_->strategies().search;
	st->saveProgress = other.solver_->strategies().saveProgress;
	st->cflMinAntes  = other.solver_->strategies().cflMinAntes;
	st->strengthenRecursive = other.solver_->strategies().strengthenRecursive;
	st->randomWatches= other.solver_->strategies().randomWatches;
	st->setCompressionStrategy(other.solver_->strategies().compress());
	
	solve   = other.solve;
	loopRep = other.loopRep;

}

/////////////////////////////////////////////////////////////////////////////////////////
// GlobalOptions & ClaspConfig
/////////////////////////////////////////////////////////////////////////////////////////
GlobalOptions::GlobalOptions() { }

ThreadOptions::Distribution::Distribution() 
	: grace(1024), lbd(0), types(0), filter(filter_gp), copyProblem(false)  {} 

Enumerator* GlobalOptions::initEnumerator(Enumerator::Report* r) {
	ModelEnumerator* e = 0;
	Enumerator* ret    = 0;
	if (consequences()) {
		ret = new CBConsequences(enumerate.brave ? CBConsequences::brave_consequences : CBConsequences::cautious_consequences);
	}
	else if (enumerate.record) {
		ret = (e = new RecordEnumerator());
	}
	else {
		ret = (e = new BacktrackEnumerator(enumerate.projectOpts));
	}
	if (e) { e->setEnableProjection(enumerate.project); }
	ret->setRestartOnModel(enumerate.restartOnModel);
	ret->setReport(r);
	if (enumerate.progress) {
		ret->enableProgressReport(enumerate.progress);
	}
	return ret;
}

ClaspConfig::ClaspConfig() {
	solvers.push_back(new LocalOptions(ctx.master()));
}

ClaspConfig::~ClaspConfig() {
	setThreads(1);
	delete solvers.back();
	solvers.clear();
}

bool ClaspConfig::validate(std::string& err) {
	for (uint32 i = 0; i != numThreads(); ++i) {
		if (!threadConfig(i)->validate(err)) {
			return false;
		}
	}
	if (enumerate.brave && enumerate.cautious) {
		err = "Options 'brave' and 'cautious' are mutually exclusive!";
		return false;
	}
	if (enumerate.restartOnModel) { enumerate.record  = true; }
	return true;
}

uint32 ClaspConfig::numThreads() const { return (uint32)solvers.size(); }

LocalOptions* ClaspConfig::threadConfig(uint32 i) const {
	assert(i < numThreads());
	return solvers[i];
}

void ClaspConfig::applyHeuristic() {
	for (uint32 i = 0; i != numThreads(); ++i) {
		threadConfig(i)->applyHeuristic();
	}
}

void ClaspConfig::reset() {
	ctx.reset();
	thread   = ThreadOptions();
	eq       = EqOptions();
	opt      = Optimize();
	enumerate= EnumOptions();
	master()->setSolver(ctx.master());
	for (uint32 i = 1; i < numThreads(); ++i) {
		threadConfig(i)->solver().reset();
	}
}

void ClaspConfig::setThreads(uint32 i) {
	if (i == 0) { i = 1; }
	while (i < numThreads()) {
		LocalOptions* x = solvers.back();
		Solver*       s = &x->solver();
		solvers.pop_back();
		delete x;
		delete s;
	}
	while (i > numThreads()) {
		Solver* x = new Solver();
		solvers.push_back(new LocalOptions(x));
	}
}
IncrementalControl::IncrementalControl()  {}
IncrementalControl::~IncrementalControl() {}
/////////////////////////////////////////////////////////////////////////////////////////
// ClaspFacade
/////////////////////////////////////////////////////////////////////////////////////////
ClaspFacade::ClaspFacade() 
	: config_(0)
	, inc_(0)
	, cb_(0)
	, input_(0)
	, ctrl_(0)
	, graph_(0)
	, api_(0)
	, result_(result_unknown)
	, state_(state_not_started)
	, step_(0)
	, more_(true) {
}

void ClaspFacade::init(Input& problem, ClaspConfig& config, IncrementalControl* inc, Callback* c) {
	config_ = &config;
	inc_    = inc;
	cb_     = c;
	input_  = &problem;
	ctrl_   = 0;
	graph_  = 0;
	api_    = 0;
	result_ = result_unknown;
	state_  = state_not_started;
	step_   = 0;
	more_   = true;
	validateWeak();
	config.applyHeuristic();
	if (config.numThreads() > 1 && !config.thread.distribute.copyProblem) {
		config.ctx.enableConstraintSharing();
	}
}

void ClaspFacade::validateWeak(ClaspConfig& cfg) {
	bool warnUnit = true;
	bool warnInit = true;
	for (uint32 i = 0; i != cfg.numThreads(); ++i) {
		HeuristicOptions& h = cfg.threadConfig(i)->heuristic();
		if (h.name == "unit") {
			if (h.lookahead == Lookahead::no_lookahead) {
				if (warnUnit) {
					warning("Unit-heuristic implies lookahead. Forcing atom-lookahead!");
					warnUnit = false;
				}
				h.lookahead = Lookahead::atom_lookahead;
			}
			else if (h.lookaheadNum != -1) {
				if (warnInit) {
					warning("Unit-heuristic implies lookahead. Ignoring 'initial-lookahead'!");
					warnInit = false;
				}
				h.lookaheadNum = -1;
			}
		}
		else if (h.lookaheadNum != -1 && h.lookahead == Lookahead::no_lookahead) {
			h.lookahead = Lookahead::atom_lookahead;
		}
	}
}

void ClaspFacade::validateWeak() {
	if (inc_) {
		if (config_->enumerate.onlyPre) {
			warning("Option 'onlyPre' is ignored in incremental setting!"); 
			config_->enumerate.onlyPre = false;
		}
	}
	if (config_->eq.noSCC && config_->eq.iters != 0) {
		warning("Supported models requires --eq=0. Disabling eq-preprocessor!");
		config_->eq.noEq();
	}
	if (config_->opt.all || config_->opt.no) {
		config_->opt.hierarch = 0;
	}
	validateWeak(*config_);
}

// Non-incremental solving...
void ClaspFacade::solve(Input& problem, ClaspConfig& config, Callback* c) {
	init(problem, config, 0, c);
	LocalOptions* master = config.master();
	if (!read() || !preprocess()) {
		result_ = result_unsat;
		more_   = false;
		reportSolution(*config.ctx.enumerator(), true);
	}
	else if (!config.enumerate.onlyPre) {
		setProblemSize();
		AutoState state(this, state_solve);
		initSolveObject(config); // deleted by AutoState
		more_   = ctrl_->solve(config.ctx, master->solve, LitVec());
	}
}

// Incremental solving...
void ClaspFacade::solveIncremental(Input& problem, ClaspConfig& config, IncrementalControl& inc, Callback* c) {
	init(problem, config, &inc, c);
	LitVec assume;
	LocalOptions* master = config_->master();
	do {
		inc.initStep(*this);
		result_   = result_unknown;
		more_     = true;
		if (!read() || !preprocess()) {
			result_ = result_unsat;
			more_   = false;
			reportSolution(*config.ctx.enumerator(), true);
			break;
		}
		else {
			setProblemSize();
			AutoState state(this, state_solve);
			assume.clear();
			problem.getAssumptions(assume);
			initSolveObject(config);
			more_    = ctrl_->solve(config.ctx, master->solve, assume);
			if (result_ == result_unknown && !more_) {
				// initial assumptions are unsat
				result_ = result_unsat;
			}
		}
	} while (inc.nextStep(*this) && ++step_);
}

// Creates a ProgramBuilder-object if necessary and reads
// the input by calling input_->read().
// Returns false, if the problem is trivially UNSAT.
bool ClaspFacade::read() {
	AutoState state(this, state_read);
	Input::ApiPtr ptr(&config_->ctx);
	if (input_->format() == Input::SMODELS) {
		if (step_ == 0) {
			api_   = new ProgramBuilder();
			api_->startProgram(config_->ctx, config_->eq);
		}
		if (inc_) {
			api_->updateProgram();
		}
		ptr.api= api_.get();
	}
	if (config_->opt.hierarch > 0 && !config_->opt.no) {
		config_->ctx.requestTagLiteral();
	}
	if (!input_->read(ptr, config_->enumerate.numModels)) {
		return false;
	}
	return true;
}

// Prepare the solving state:
//  - if necessary, transforms the input to nogoods by calling ProgramBuilder::endProgram()
//  - fires event_p_prepared after input was transformed to nogoods
//  - adds any minimize statements to the solver and initializes the enumerator
//  - calls Solver::endInit().
// Returns false, if the problem is trivially UNSAT.
bool ClaspFacade::preprocess() {
	AutoState state(this, state_preprocess);
	SharedContext& ctx = config_->ctx;
	SharedMinimizeData* m = 0;
	Input::ApiPtr ptr(&ctx);
	if (api_.get()) {
		if (!api_->endProgram()) {
			fireEvent(*ctx.master(), event_p_prepared);
			return false;
		}
		setGraph();
		ptr.api = api_.get();
	}
	if (!config_->opt.no && step_ == 0) {
		MinimizeBuilder builder;
		input_->addMinimize(builder, ptr);
		if (builder.hasRules()) {
			if (!config_->opt.vals.empty()) {
				const WeightVec& opt = config_->opt.vals;
				for (uint32 i = 0; i != opt.size(); ++i) {
					builder.setOptimum(i, opt[i]);
				}
			}
			m = builder.build(ctx, config_->ctx.tagLiteral());
		}
		if (!builder.hasRules() || (builder.numRules() == 1 && config_->opt.hierarch < 2)) {
			config_->ctx.removeTagLiteral();
		}
	}
	fireEvent(*ctx.master(), event_p_prepared);
	if (!inc_ && api_.is_owner()) {
		api_ = 0;
	}
	if (!initEnumeration(m) || !config_->ctx.endInit(config_->numThreads())) {
		return false;
	}
	return true;
}

// Configures the given minimize constraint and adds it to the enumerator.
// Optimize values that are given in config are added to min.
bool ClaspFacade::configureMinimize(SharedMinimizeData* min) const {
	min->setMode(config_->opt.all ? MinimizeMode_t::enumerate : MinimizeMode_t::optimize, config_->opt.hierarch);
	if (config_->consequences()) {
		warning("Minimize statements: Consequences may depend on enumeration order!");
	}
	if (config_->enumerate.project) {
		for (const WeightLiteral* it = min->lits; !isSentinel(it->first); ++it) {
			if ( !config_->ctx.project(it->first.var()) ) {
				warning("Projection: Optimization values may depend on enumeration order!");
				break;
			}
		}
	}
	config_->ctx.enumerator()->setMinimize(min, config_->opt.heu);
	return true;
}

// Finalizes initialization of model enumeration.
// Configures and adds an eventual minimize constraint,
// sts the number of models to compute and adds warnings
// if this number conflicts with the preferred number of the enumerator.
bool ClaspFacade::initEnumeration(SharedMinimizeData* min)  {
	Enumerator* e = step_ == 0 ? config_->initEnumerator(this) : config_->ctx.enumerator();
	config_->ctx.addEnumerator(e);
	GlobalOptions::EnumOptions& opts = config_->enumerate;
	if (step_ == 0) {
		if (min && !configureMinimize(min)) {
			return false;
		}
		uint32 defM = !e->minimize() && !config_->consequences();
		if (opts.numModels == -1) { 
			opts.numModels = defM; 
		}
		else if (opts.numModels > 0 && defM == 0) {
			if (config_->consequences()) {
				warning("Option '--number' not 0: last model may not cover consequences!");
			}
			if (e->minimize()) {
				warning("Option '--number' not 0: Optimality of last model not guaranteed!");
			}
		}
		if (config_->numThreads() > 1) {
			bool forceST = opts.project && !opts.record;
			if (forceST) {
				warning("Multi-Threading disabled - Selected reasoning mode not supported!");
				config_->setThreads(1);
			}
			else if (config_->numThreads() > ClaspFacade::hardware_concurrency()) {
				warning("Oversubscription: Number of threads exceeds hardware concurrency!");
			}
		}
	}
	config_->ctx.enumerator()->enumerate(opts.numModels);
	return true;
}

void ClaspFacade::setGraph() {
	if (!graph_.get() && api_->dependencyGraph() && api_->dependencyGraph()->nodes() > 0) {
		graph_ = api_->dependencyGraph(true);
		DefaultUnfoundedCheck* ufs = new DefaultUnfoundedCheck(config_->master()->loopRep);
		ufs->attachTo(*config_->ctx.master(), graph_.get());
	}
}

// Computes a value that represents the problem size.
// The value is then used by the reduce-heuristic
// to determine the initial learnt db size.
void ClaspFacade::setProblemSize() const {
	const SharedContext& ctx  = config_->ctx;
	uint32 estimate = 0;
	uint32 ps;
	for (uint32 i = 0; i != config_->numThreads(); ++i) {
		LocalOptions* o = config_->threadConfig(i);
		if (o->solve.reduce.estimate) {
			if (estimate == 0) {
				estimate = ctx.problemComplexity();
			}
			ps = estimate;
		}
		else if (input_->format() != Input::DIMACS) {
			double r = ctx.numVars() / std::max(1.0, double(ctx.numConstraints()));
			if (r < 0.1 || r > 10.0) {
				ps = std::max(ctx.numVars(), ctx.numConstraints());
			}
			else {
				ps = std::min(ctx.numVars(), ctx.numConstraints());
			}
		}
		else {
			ps = ctx.numConstraints();
		}
		o->solve.reduce.setProblemSize(ps);
	}
}

#ifndef DISABLE_MULTI_THREADING
uint32 ClaspFacade::hardware_concurrency() { return mt::ParallelSolve::hardwareThreads(); }
void ClaspFacade::initSolveObject(ClaspConfig& config) {
	ctrl_    = 0;
	int numT = config.numThreads();
	if (step_ != 0) {
		assert(config.master()->solver().id() != mt::ParallelSolve::invalidId);
		int invalid = 0;
		for (uint32 i = 1; i != config.numThreads(); ++i) {
			invalid += (config.threadConfig(i)->solver().id() == mt::ParallelSolve::invalidId);
		}
		numT -= invalid;
	}
	if (numT == 1) {
		ctrl_ = new SimpleSolve();
		return;
	}
	// parallel solving
	mt::ParallelSolve* ctrl = new mt::ParallelSolve();
	ctrl_ = ctrl;
	std::auto_ptr<mt::GlobalQueue> dist;
	if (config.thread.distribute.types != 0) {
		dist.reset(new mt::GlobalQueue(numT, config.thread.distribute.types, config.thread.distribute.lbd));
		config.ctx.enableLearntSharing(dist.release());
	}
	ctrl->setForceGP(config.thread.forceGP);
	ctrl->setIntegrate(config.thread.distribute.grace, config.thread.distribute.filter);
	ctrl->setRestarts(config.thread.restarts.maxR, config.thread.restarts.sched);
	ctrl->init(numT, config.ctx);
	for (int i = 1; i != numT; ++i) {
		LocalOptions* ti = config.threadConfig(i);
		ti->solver().stats.reset();
		if (ti->solver().id() != mt::ParallelSolve::invalidId) {
			if (ti->solver().sharedContext() == 0 && graph_.get()) {
				DefaultUnfoundedCheck* u = new DefaultUnfoundedCheck(ti->loopRep);
				u->attachTo(ti->solver(), graph_.get());
			}
			ctrl->addThread(ti->solver(), ti->solve);
		}
	}
}
#else
uint32 ClaspFacade::hardware_concurrency() { return 1; }
void ClaspFacade::initSolveObject(ClaspConfig&) {
	ctrl_ = new SimpleSolve();
}
#endif

}
