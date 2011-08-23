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
#include <clasp/shared_context.h>
#include <clasp/solver.h>
#include <clasp/clause.h>
#include <clasp/enumerator.h>
#include <clasp/util/thread.h>
namespace Clasp {
/////////////////////////////////////////////////////////////////////////////////////////
// ShortImplicationsGraph::ImplicationList
/////////////////////////////////////////////////////////////////////////////////////////
namespace {
	inline bool isBinary(const Literal* imp) { return imp->watched(); }
}
#ifndef DISABLE_MULTI_THREADING
struct ShortImplicationsGraph::Block {
	typedef std::atomic<uint32> atomic_size;
	enum { block_cap = (64 - (sizeof(atomic_size)+sizeof(SharedBlockPtr)))/sizeof(Literal) };
	explicit Block() {
		for (int i = 0; i != block_cap; ++i) { data[i] = posLit(0); }
		size_lock = 0;
		next      = 0;
	}
	const Literal* begin() const { return data; }
	const Literal* end()   const { return data+size(); }
	Literal*       end()         { return data+size(); }
	uint32         size()  const { return size_lock >> 1; }
	void add(const Literal* x, uint32 a_size, uint32 cs) {
		std::copy(x, x+a_size, data+cs);
		if (a_size == 1) {
			data[cs].watch();
		}
	}
	bool tryLock(uint32& size) {
		uint32 s = size_lock;
		if ((s & 1) == 0 && size_lock.compare_and_swap(s | 1, s) == s) {
			size = s >> 1;
			return true;
		}
		return false;
	}
	void unlock(uint32 newSize) {
		size_lock = (newSize << 1);
	}
	SharedBlockPtr next;
	atomic_size    size_lock;
	Literal        data[block_cap];
};

#define FOR_EACH_LEARNT(x, Y) \
	if (0) ; else \
	for (Block* b = (x).learnt; b ; b = b->next) \
		for (const Literal* Y = b->begin(), *endof = b->end(); Y != endof; Y += 2 - Y->watched())

ShortImplicationsGraph::ImplicationList::~ImplicationList() {
	clear(true);
}

void ShortImplicationsGraph::ImplicationList::clear(bool b) {
	ImpListBase::clear(b);
	for (Block* x = learnt; x; ) {
		Block* t = x;
		x = x->next;
		delete t;
	}
	learnt = 0;
}

void ShortImplicationsGraph::ImplicationList::addLearnt(Literal p, Literal q) {
	Literal nc[2] = {p, q};
	uint32 size   = 1 + !isSentinel(q);
	for (Block* x;;) {
		x = learnt;
		if (x) {
			uint32 cs;
			if (x->tryLock(cs)) {
				if ( (cs + size) <=  Block::block_cap ) {
					x->add(nc, size, cs);
					x->unlock(cs + size);
				}
				else {
					Block* t = new Block();
					t->add(nc, size, 0);
					t->unlock(size);
					t->next   = x; // x is full and remains locked forever
					x = learnt= t; // publish new block - unlocks x and learnt
				}
				return;
			}
			else { 
				std::this_thread::yield();
			}
		}
		else {
			x = new Block();
			if (learnt.compare_and_swap(x, 0) != 0) {
				delete x;
			}
		}
	}
}

bool ShortImplicationsGraph::ImplicationList::hasLearnt(Literal q, Literal r) const {
	const bool binary = isSentinel(r);
	FOR_EACH_LEARNT(*this, imp) {
		if (imp[0] == q || imp[0] == r) {
			// binary clause subsumes new bin/tern clause
			if (isBinary(imp))                           { return true; }
			// existing ternary clause subsumes new tern clause
			if (!binary && (imp[1] == q || imp[1] == r)) { return true; }
		}
	}
	return false;
}

void ShortImplicationsGraph::ImplicationList::move(ImplicationList& other) {
	ImpListBase::move(other);
	delete learnt;
	learnt       = other.learnt;
	other.learnt = 0;
}
#else
#define FOR_EACH_LEARNT(x, Y) if (const Literal* Y = 0) ; else while(0)
#endif
/////////////////////////////////////////////////////////////////////////////////////////
// ShortImplicationsGraph
/////////////////////////////////////////////////////////////////////////////////////////
ShortImplicationsGraph::ShortImplicationsGraph() {
	bin_[0]  = bin_[1]  = 0;
	tern_[0] = tern_[1] = 0;
}
ShortImplicationsGraph::~ShortImplicationsGraph() {
	PodVector<ImplicationList>::destruct(graph_);
}
void ShortImplicationsGraph::resize(uint32 nodes) {
	if (graph_.capacity() >= nodes) {
		graph_.resize(nodes);
	}
	else {
		ImpLists temp; temp.resize(nodes);
		for (ImpLists::size_type i = 0; i != graph_.size(); ++i) {
			temp[i].move(graph_[i]);
		}
		graph_.swap(temp);
	}
}

uint32 ShortImplicationsGraph::numEdges(Literal p) const { return graph_[p.index()].size(); }

bool ShortImplicationsGraph::addBinary(Literal p, Literal q, bool learnt, bool shared) {
	uint32 added = 0;
	if (!shared) {
		if (!learnt) { p.clearWatch(); q.clearWatch(); }
		else         { p.watch(); q.watch(); }
		getList(~p).push_left(q);
		getList(~q).push_left(p);
		added = 1;
	}
#ifndef DISABLE_MULTI_THREADING
	else if (learnt && !getList(~p).hasLearnt(q)) {
		getList(~p).addLearnt(q);
		getList(~q).addLearnt(p);
		added = 1;
	}
#endif
	bin_[learnt] += added;
	return added > 0;
}

bool ShortImplicationsGraph::addTernary(Literal p, Literal q, Literal r, bool learnt, bool shared) {
	uint32 added = 0;
	if (!shared) {
		if (!learnt) { p.clearWatch(); q.clearWatch(); r.clearWatch(); }
		else         { p.watch(); q.watch(); r.watch(); }
		getList(~p).push_right(std::make_pair(q, r));
		getList(~q).push_right(std::make_pair(p, r));
		getList(~r).push_right(std::make_pair(p, q));
		added = 1;
	}
#ifndef DISABLE_MULTI_THREADING
	else if (learnt && !getList(~p).hasLearnt(q, r)) {
		getList(~p).addLearnt(q, r);
		getList(~q).addLearnt(p, r);
		getList(~r).addLearnt(p, q);
		added = 1;
	}
#endif
	tern_[learnt] += added;
	return added > 0;
}

void ShortImplicationsGraph::remove_bin(ImplicationList& w, Literal p) {
	w.erase_left_unordered(std::find(w.left_begin(), w.left_end(), p)); 
	w.try_shrink(); 
}
void ShortImplicationsGraph::remove_tern(ImplicationList& w, Literal p) {
	w.erase_right_unordered(std::find_if(w.right_begin(), w.right_end(), PairContains<Literal>(p))); 
	w.try_shrink();	
}

void ShortImplicationsGraph::removeTrue(Solver& s, Literal p) {
	typedef ImplicationList SWL;
	SWL& negPList = graph_[(~p).index()];
	SWL& pList    = graph_[ (p).index()];
	// remove every binary clause containing p -> clause is satisfied
	for (SWL::left_iterator it = negPList.left_begin(), end = negPList.left_end(); it != end; ++it) {
		--bin_[it->watched()];
		remove_bin(graph_[(~*it).index()], p);
	}
	// remove every ternary clause containing p -> clause is satisfied
	for (SWL::right_iterator it = negPList.right_begin(), end = negPList.right_end(); it != end; ++it) {
		--tern_[it->first.watched()];
		remove_tern(graph_[ (~it->first).index() ], p);
		remove_tern(graph_[ (~it->second).index() ], p);
	}
	// transform ternary clauses containing ~p to binary clause
	for (SWL::right_iterator it = pList.right_begin(), end = pList.right_end(); it != end; ++it) {
		Literal q = it->first;
		Literal r = it->second;
		--tern_[q.watched()];
		remove_tern(graph_[(~q).index()], ~p);
		remove_tern(graph_[(~r).index()], ~p);
		if (s.value(q.var()) == value_free && s.value(r.var()) == value_free) {
			// clause is binary on dl 0
			addBinary(q, r, false, false);
		}
		// else: clause is SAT and removed when the satisfied literal is processed
	}
	graph_[(~p).index()].clear(true);
	graph_[ (p).index()].clear(true);
}

bool ShortImplicationsGraph::propagate(Solver& s, Literal p) const {
	const ImplicationList& x = graph_[p.index()];
	if (x.empty()) return true;
	ImplicationList::const_right_iterator rEnd = x.right_end(); // prefetch
	// Binary BCP
	for (ImplicationList::const_left_iterator it = x.left_begin(), end = x.left_end(); it != end; ++it) {
		if (!s.isTrue(*it) && !s.force(*it, Antecedent(p))) { return false; }
	}
	// Ternary BCP
	Literal q, r;
	for (ImplicationList::const_right_iterator it = x.right_begin(); it != rEnd; ++it) {
		if (!s.isTrue(q = it->first) && !s.isTrue(r = it->second)) {
			if (s.isFalse(q)) { q = r; r = it->first; }
			if (s.isFalse(r) && !s.force(q, Antecedent(p, ~r))) {
				return false;
			}
		}
	}
	// combined short BCP on learnt clauses
	FOR_EACH_LEARNT(x, imp) {
		if (!s.isTrue(q = imp[0]) && (isBinary(imp) || !s.isTrue(r = imp[1]))) {
			if (!isBinary(imp)) {
				if (s.isFalse(q)) { q = r; r = imp[0]; }
				if (s.isFalse(r) && !s.force(q, Antecedent(p, ~r))) {
					return false;
				}
			}
			else if (!s.force(q, Antecedent(p))) {
				return false;
			}
		}
	}
	return true;
}
bool ShortImplicationsGraph::reverseArc(const Solver& s, Literal p, uint32 maxLev, Antecedent& out) const {
	const ImplicationList& x = graph_[p.index()];
	if (x.empty()) return false;
	ImplicationList::const_right_iterator rEnd = x.right_end(); // prefetch
	for (ImplicationList::const_left_iterator it = x.left_begin(), end = x.left_end(); it != end; ++it) {
		if (isRevLit(s, *it, maxLev)) { 
			out = Antecedent(~*it);
			return true;
		}	
	}
	for (ImplicationList::const_right_iterator it = x.right_begin(); it != rEnd; ++it) {
		if (isRevLit(s, it->first, maxLev) && isRevLit(s, it->second, maxLev)) {
			out = Antecedent(~it->first, ~it->second);
			return true;
		}
	}
	FOR_EACH_LEARNT(x, imp) {
		if (isRevLit(s, imp[0], maxLev) && (isBinary(imp) || isRevLit(s, imp[1], maxLev))) {
			out = !isBinary(imp) ? Antecedent(~imp[0], ~imp[1]) : Antecedent(~imp[0]);
			return true;
		}
	}
	return false;
}
bool ShortImplicationsGraph::propagateBin(Assignment& out, Literal p, uint32 level) const {
	const ImplicationList& x = graph_[p.index()];
	Antecedent ante(p);
	for (ImplicationList::const_left_iterator it = x.left_begin(), end = x.left_end(); it != end; ++it) {
		if (!out.assign(*it, level, p)) {
			return false;
		}
	}
	return true;
}
/////////////////////////////////////////////////////////////////////////////////////////
// SatPreprocessor
/////////////////////////////////////////////////////////////////////////////////////////
SatPreprocessor::~SatPreprocessor() {}
void SatPreprocessor::reportProgress(char t, uint32 min, uint32 max) {
	ctx_->enumerator()->reportPreProgress(t, *this, min, max);
}

/////////////////////////////////////////////////////////////////////////////////////////
// SharedContext
/////////////////////////////////////////////////////////////////////////////////////////
SharedContext::SharedContext() 
	: symTabPtr_(new SharedSymTab()), master_(new Solver), distributor_(0), lastInit_(0), lastTopLevel_(0), shareCount_(1), shareConstr_(false) {
	addVar(Var_t::atom_body_var);
	problem_.vars = 0;
	addEnumerator(0);
	Antecedent::checkPlatformAssumptions();
}

SharedContext::SharedContext(const SharedContext& rhs,  InitMode) 
	: master_(new Solver), distributor_(0), lastInit_(0), lastTopLevel_(0), shareCount_(1), shareConstr_(false) {
	symTabPtr_ = rhs.symTabPtr_;
	++symTabPtr_->refs;
	addVar(Var_t::atom_body_var);
	problem_.vars = 0;
	addEnumerator(0);
}

SharedContext::~SharedContext() {
	delete master_;
	if (--symTabPtr_->refs == 0) delete symTabPtr_;
}

void SharedContext::reset() {
	this->~SharedContext();
	new (this) SharedContext();	
}

void SharedContext::reserveVars(uint32 varGuess) {
	varInfo_.reserve(varGuess);
}

Var SharedContext::addVar(VarType t, bool eq) {
	Var v = varInfo_.numVars();
	varInfo_.add(t == Var_t::body_var);
	if (eq) varInfo_.toggle(v, VarInfo::EQ);
	++problem_.vars;
	return v;
}

void SharedContext::requestTagLiteral() {
	if (tag_ == posLit(0)) {
		tag_ = negLit(0);
	}
}

void SharedContext::removeTagLiteral() {
	assert(master());
	if (!isSentinel(tag_)) {
		master()->force(tag_, 0);
	}
	tag_ = posLit(0);
}

void SharedContext::requestData(Var v) {
	master()->assign_.requestData(v + 1);
}

void SharedContext::setFrozen(Var v, bool b) {
	assert(validVar(v)); 
	if (b != varInfo_.isSet(v, VarInfo::FROZEN)) {
		varInfo_.toggle(v, VarInfo::FROZEN);
		b ? ++problem_.vars_frozen : --problem_.vars_frozen;
	}
}
void SharedContext::eliminate(Var v) {
	assert(validVar(v) && unique() && master()->decisionLevel() == 0); 
	if (!eliminated(v)) {
		varInfo_.toggle(v, VarInfo::ELIM);
		++problem_.vars_eliminated;
		// assign var to true - no longer a decision variable!
		master()->assign_.eliminate(v);
	}
}


Solver& SharedContext::startAddConstraints(uint32 constraintGuess) {
	assert(unique());
	if (tag_ == negLit(0)) {
		// add aux var for tag literal
		tag_ = posLit(addVar(Var_t::atom_var));
		setFrozen(tag_.var(), true);
		--problem_.vars;
	}
	btig_.resize((numVars()+1)<<1);
	if (satPrepro.get()) {
		satPrepro->setContext(*this);
	}
	master_->startInit(this, constraintGuess);
	lastInit_ = master_->constraints_.size();
	return *master_;
}

void SharedContext::add(Constraint* c) {
	assert(unique());
	++problem_.constraints;
	master()->constraints_.push_back(c);
}

bool SharedContext::addUnary(Literal p) {
	assert(unique());
	return master()->addUnary(p, Constraint_t::static_constraint);
}

void SharedContext::addPost(PostPropagator* p) {
	assert(unique());
	return master()->addPost(p);
}

void SharedContext::addEnumerator(Enumerator* en) {
	if (en == 0) {
		en = new NullEnumerator();
	}
	enumerator_.reset(en);
	enumerator_->startInit(*this);
}
 
uint32 SharedContext::numConstraints() const {
	return numBinary() + numTernary() + master()->constraints_.size();
}

bool SharedContext::endInit(uint32 shareCount) {
	assert(unique());
	if (master()->hasConflict()) return false;
	if (!master()->post_.init(*master())) {
		return false;
	}
	struct Holder {
		~Holder() { if (con) con->destroy(master, true); }
		Solver*     master;
		Constraint* con;
	} enumC = { master(), enumerator()->endInit(*this, shareCount) };
	if (satPrepro.get()) {
		SatPrepro temp(satPrepro.release());
		bool r = temp->preprocess(enumerator()->enumerate());
		satPrepro.reset(temp.release());
		if (!r) return false;
	}
	master_->setEnumerationConstraint(enumC.con);
	enumC.con = 0;
	if (!master()->endInit()) return false;
	lastTopLevel_ = master()->units_;
	shareCount_   = shareCount;
	problem_.constraints_binary = btig_.numBinary();
	problem_.constraints_ternary= btig_.numTernary();
	return true;
}

bool SharedContext::attach(Solver& other) {
	assert(!unique() && master_ != &other);
	Var oldV = other.numVars();
	other.startInit(this, static_cast<uint32>(master_->constraints_.size()-lastInit_));
	// 1. clone assignment
	other.assign_.requestData(master()->assign_.numData());
	LitVec::size_type prevTs = other.trail().size();
	const LitVec& trail      = master()->trail();
	Antecedent null;
	for (LitVec::size_type i = other.units_; i < trail.size(); ++i) {
		if (!other.force(trail[i], null)) {
			return false;
		}
		other.markSeen(trail[i].var());
	}
	other.units_        = static_cast<uint32>(trail.size());
	other.lastSimplify_ = other.constraints_.empty() ? trail.size() : prevTs;
	if (satPrepro.get()) {
		for (Var v = oldV+1; v <= other.numVars(); ++v) {
			if (eliminated(v) && other.value(v) == value_free) {
				other.assign_.eliminate(v);
			}
		}
	}
	// 2. clone & attach constraints
	const Solver::ConstraintDB& db = master()->constraints_;
	for (LitVec::size_type i = lastInit_; i < db.size(); ++i) {
		if (Constraint* c = db[i]->cloneAttach(other)) {
			other.constraints_.push_back(c);
		}
		if (other.hasConflict()) return false;
	}
	Constraint* c = master_->getEnumerationConstraint();
	other.setEnumerationConstraint( c ? c->cloneAttach(other) : 0 );
	// 3. endInit
	return (other.post_.init(other) && other.endInit())
		|| (detach(other), false);
}

void SharedContext::detach(Solver& s) {
	s.setEnumerationConstraint(0);
	if (master() == &s) {
		shareCount_ = 1;
	}
}

uint32 SharedContext::problemComplexity() const {
	uint32 r = numBinary() + numTernary();
	for (uint32 i = 0; i != master()->constraints_.size(); ++i) {
		r += master()->constraints_[i]->estimateComplexity(*master());
	}
	return r;
}

SharedLiterals* SharedContext::distribute(const Solver& s, const Literal* lits, uint32 size, const ClauseInfo& extra) const {
	SharedLiterals* ret = distributor_.get() ? distributor_->publish(s, lits, size, extra) : 0;
	if (ret && s.stats.parallel) ++s.stats.parallel->shared;
	return ret;
}

void SharedContext::distribute(const Solver& s, Literal p, Literal q, Literal r, ConstraintType t) const {
	if (distributor_.get()) {
		Literal lits[3] = {p, q, r};
		SharedLiterals* shared = distributor_->publish(s, lits, 1 + !isSentinel(q) + !isSentinel(r), ClauseInfo().setType(t));
		if (shared)           { shared->release(); }
		if (s.stats.parallel) { ++s.stats.parallel->shared; }
	}
}

uint32 SharedContext::receive(const Solver& target, SharedLiterals** out, uint32 maxOut) const {
	if (distributor_.get()) {
		return distributor_->receive(target, out, maxOut);
	}
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////
// Distributor
/////////////////////////////////////////////////////////////////////////////////////////
Distributor::Distributor(uint32 maxShare, uint32 typesToShare, uint32 maxLbd) : maxShare_(maxShare), lbdMax_(maxLbd), typeMask_(typesToShare) {}
Distributor::~Distributor() {}

SharedLiterals* Distributor::publish(const Solver& source, const Literal* lits, uint32 size, const ClauseInfo& extra) {
	SharedLiterals* res = 0;
	if (maxShare_ > 1 && !extra.tagged() && (size <= 3 || (extra.rootLbd() <= lbdMax_ && (extra.type() & typeMask_) != 0))) {
		res = SharedLiterals::newShareable(lits, size, extra.type(), maxShare_);
		doPublish(source, res);
	}
	return res;
}

}
