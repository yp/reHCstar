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
#include <clasp/preprocessor.h>
#include <clasp/program_builder.h>
namespace Clasp { namespace {

struct LessBodySize {
	LessBodySize(const BodyList& bl) : bodies_(&bl) {}
	bool operator()(Var b1, Var b2 ) const {
		return (*bodies_)[b1]->size() < (*bodies_)[b2]->size()
			|| ((*bodies_)[b1]->size() == (*bodies_)[b2]->size() && (*bodies_)[b1]->type() < (*bodies_)[b2]->type());
	}
private:
	const BodyList* bodies_;
};
}

Literal Preprocessor::newBodyLit(PrgBodyNode* n) {
	if      (!backprop_||n->value()!=value_false) { 
		                      return posLit(prg_->vars_.add(Var_t::body_var)); }
	else if (false_ == 0) { return posLit(false_ = prg_->vars_.add(Var_t::body_var)); }
	else                  { return posLit(false_); }
}

// If the program is defined incrementally, marks atoms from previous steps as supported.
// Logically, the functions acts as if the program contains a choice rule over
// all atoms defined in the previous steps minus those to be unfrozen and hence
// defined in this step.	
// Pre: Atoms that are true/false have their value-member set
void Preprocessor::updatePreviouslyDefinedAtoms(Var startAtom, bool strong) {
	follow_.clear();
	if (prg_->incData_) {
		// mark all atom to be unfrozen in this step
		for (VarVec::const_iterator it = prg_->incData_->unfreeze_.begin(), end = prg_->incData_->unfreeze_.end(); it != end; ++it) {
			assert(*it < startAtom);
			PrgAtomNode* a = prg_->atoms_[*it];
			Literal      x = a->literal(); x.watch();
			a->setLiteral(x);
		}
	}
	LitVec::size_type old = prg_->compute_.size();
	for (Var i = 1; i < startAtom; ++i) {
		PrgAtomNode* a   = prg_->atoms_[i];
		if (a->literal().watched()) { 
			Literal x = a->literal(); x.clearWatch(); 
			a->setLiteral(x); 
			if (!strong) { a->preds.clear(); }
			else         { setRootAtom(a->literal(), i); };
		}
		else if (!a->eq() && a->value() != value_false) {
			if (strong) {
				setRootAtom(a->literal(), i);
				nodes_[i].aSeen = 1;
				propagateAtomVar((Var)i, a);
			}
			else {
				for (VarVec::size_type b = 0; b != a->posDep.size(); ++b) {
					PrgBodyNode* body = prg_->bodies_[a->posDep[b]];
					bool isSupp = body->isSupported();
					body->onPosPredSupported(Var(i));
					if (!isSupp && body->isSupported()) {
						prg_->initialSupp_.push_back(a->posDep[b]);
					}
				}
			}
		}
		else {
			assert(!a->eq() || (a->posDep.empty() && a->negDep.empty() && "EQ-Atom in new Program - rule not simplified!"));
			if (a->posDep.size() + a->negDep.size() + a->preds.size() > 0) {
				prg_->setCompute(i, false);
			}
		}
	}
	applyCompute(prg_->compute_, old, strong);
}

// Back-propagates atoms from the compute statement.
bool Preprocessor::applyCompute(LitVec& compute, uint32 start, bool strong) {
	uint32 newStart = (uint32)compute.size();
	for (LitVec::size_type i = start; i < compute.size(); ++i) {
		Var        atomId = prg_->getEqAtom(compute[i].var());
		PrgAtomNode* atom = prg_->atoms_[atomId];
		if (strong) {
			// We can safely remove atom from bodies where it occurs negatively.
			// If atom is weak_true, subgoal can't be true.
			// If atom is false, subgoal is true.
			for (VarVec::size_type j = 0; j != atom->negDep.size(); ++j) {
				setSimplifyBody(atom->negDep[j]);
			}	
		}
		if (compute[i].sign()){// compute false
			if (!atom->setValue(value_false)) return false;
			// atom is false, thus all its normal bodies must be false
			for (HeadVec::iterator it = atom->preds.begin(), end = atom->preds.end(); it != end; ++it) {
				PrgBodyNode* b = prg_->bodies_[it->node()];
				b->removeHead(atomId);
				if (it->normal()) {
					if (!b->propagateFalse(it->node(), prg_->atoms_)) return false;
					// backpropagate false body
					if (backprop_ && !b->backpropagate(*prg_, compute)) { return false; }
				}
				else if (!b->hasHeads() && b->value() != value_false) {
					b->setIgnore(true);
				}
			}
			atom->preds.clear();
		}
		else if (!atom->setValue(value_weak_true)) return false;
		else if (backprop_ && atom->preds.size() == 1) {
			PrgBodyNode* b = prg_->bodies_[atom->preds[0].node()];
			if (!b->setValue(value_weak_true) || !b->backpropagate(*prg_, compute)) return false;
		}
	}
	compute.erase(compute.begin()+newStart, compute.end());
	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////
// simple preprocessing
//
// Simplifies the program by computing max consequences.
// If bodyEq == false, adds a variable for each supported atom and body
// Otherwise: Adds a variable for each supported atom and each body not equivalent to 
// some atom.
/////////////////////////////////////////////////////////////////////////////////////////
bool Preprocessor::preprocessSimple(bool bodyEq) {
	if (!applyCompute(prg_->compute_, 0, false)) return false;
	if (prg_->incData_) { updatePreviouslyDefinedAtoms(prg_->incData_->startAtom_, false); }
	Literal marked = negLit(0); marked.watch();
	if (bodyEq) {
		std::stable_sort(prg_->initialSupp_.begin(), prg_->initialSupp_.end(), LessBodySize(prg_->bodies_));
		for (VarVec::size_type i = 0; i != prg_->initialSupp_.size(); ++i) {
			prg_->bodies_[prg_->initialSupp_[i]]->setLiteral(marked);
		}
	}
	for (VarVec::size_type i = 0; i < prg_->initialSupp_.size(); ++i) {
		PrgBodyNode* b = prg_->bodies_[prg_->initialSupp_[i]];
		b->buildHeadSet();
		if (!bodyEq) b->setLiteral(newBodyLit(b));
		for (HeadVec::const_iterator h = b->heads_begin(); h != b->heads_end(); ++h) {
			bool m = bodyEq && b->size() == 0 && h->normal();
			if ( !prg_->atoms_[h->node()]->hasVar() ) {
				prg_->atoms_[h->node()]->preds.clear();
				if (m) {
					prg_->atoms_[h->node()]->setLiteral(posLit(0));
					prg_->atoms_[h->node()]->setValue(value_true);
				}
				else {
					prg_->atoms_[h->node()]->setLiteral(posLit(prg_->vars_.add(Var_t::atom_var)));
				}
			}
			if (prg_->atoms_[h->node()]->preds.empty()) {
				const VarVec& bl = prg_->atoms_[h->node()]->posDep;
				for (VarVec::const_iterator it = bl.begin(); it != bl.end(); ++it) {
					PrgBodyNode* nb = prg_->bodies_[*it];
					if (m) { nb->setLiteral(marked); }
					if (!nb->isSupported() && nb->onPosPredSupported(h->node())) { prg_->initialSupp_.push_back( *it ); }
				}
			}
			prg_->atoms_[h->node()]->preds.push_back(HeadEdge(prg_->initialSupp_[i], h->type()));
		}
	}
	if (bodyEq) {
		std::pair<uint32, uint32> ignore;
		for (VarVec::size_type i = 0; i < prg_->initialSupp_.size(); ++i) {
			uint32 bodyId   = prg_->initialSupp_[i];
			PrgBodyNode* b  = prg_->bodies_[bodyId];
			if (b->literal().watched() && !b->simplifyBody(*prg_, bodyId, ignore, *this, true)) {
				return false;
			}
			Literal l; l.watch();
			if      (b->ignore())     { l = negLit(0); }  
			else if (b->size() == 0)  { l = posLit(0); }
			else if (b->size() == 1)  { 
				l = b->posSize() == 1 ? prg_->atoms_[b->pos(0)]->literal() : ~prg_->atoms_[b->neg(0)]->literal(); 
				prg_->vars_.setAtomBody(l.var());
				prg_->incEqs(Var_t::atom_body_var);
			}
			else {
				for (HeadVec::const_iterator it = b->heads_begin(); it != b->heads_end(); ++it) {
					PrgAtomNode* a = prg_->atoms_[it->node()];
					if (a->preds.size() == 1 && it->normal()) {
						l = a->literal();
						prg_->incEqs(Var_t::atom_body_var);
						break;
					}
				}
				if (l.watched()) l = newBodyLit(b);
			}
			b->setLiteral(l);
		}
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////
// equivalence preprocessing
//
// Computes max consequences and minimizes the number of necessary variables 
// by computing equivalence-classes.
/////////////////////////////////////////////////////////////////////////////////////////
bool Preprocessor::preprocessEq(uint32 maxIters) {
	nodes_.resize( prg_->bodies_.size() >= prg_->atoms_.size() ? prg_->bodies_.size() : prg_->atoms_.size() );
	Var startAtom   = prg_->incData_?prg_->incData_->startAtom_ : 0;
	Var stopAtom    = startAtom;
	Var startVar    = prg_->incData_?prg_->incData_->startVar_  : 1;
	pass_           = 0;
	maxPass_        = maxIters;
	if (!applyCompute(prg_->compute_, 0, true)) return false;
	do {
		++pass_;
		prg_->vars_.shrink(startVar);
		if (prg_->vars_.size() <= false_) false_ = 0;
		litToNode_.clear();
		for (Var i = startAtom; i < stopAtom; ++i) { 
			prg_->atoms_[i]->clearLiteral(false); 
		} 
		if (!classifyProgram(startAtom, stopAtom)) return false;
	} while (stopAtom != (uint32)prg_->atoms_.size());
	// replace equivalent atoms in minimize rules
	for (ProgramBuilder::MinimizeRule* r = prg_->minimize_; r; r = r->next_) {
		for (WeightLitVec::iterator it = r->lits_.begin(); it != r->lits_.end(); ++it) {
			it->first = Literal(prg_->getEqAtom(it->first.var()), it->first.sign());
		}
	}
	// replace equivalent atoms in compute statement
	for (LitVec::iterator it = prg_->compute_.begin(), end = prg_->compute_.end(); it != end; ++it) {
		if (prg_->getEqAtom(it->var()) != it->var()) {
			*it = Literal(prg_->getEqAtom(it->var()), it->sign());
		}
	}
	return true;
}

// Computes necessary equivalence-classes starting from the supported bodies
// of a program.
bool Preprocessor::classifyProgram(uint32 startAtom, uint32& stopAtom) {
	assert(ok_);
	Var bodyId, bodyEqId; PrgBodyNode* body;
	Var atomId, atomEqId; PrgAtomNode* atom;
	VarVec::size_type index = 0;
	follow_.clear();
	updatePreviouslyDefinedAtoms(startAtom, true);
	std::stable_sort(prg_->initialSupp_.begin(), prg_->initialSupp_.end(), LessBodySize(prg_->bodies_));
	for (VarVec::size_type i = 0;;) {
		while ( (bodyId = nextBodyId(index)) != varMax ) {
			body        = prg_->bodies_[bodyId];
			bodyEqId    = addBodyVar(bodyId, body);
			if (!ok_) return false;
			for (HeadVec::const_iterator it = body->heads_begin(), end = body->heads_end(); it != end; ++it) {
				atomId    = it->node();
				atom      = prg_->atoms_[atomId];
				atomEqId  = addAtomVar(*it, atom, body);
				if (!ok_) return false;
				if (atomId != atomEqId) {
					atom->clearLiteral(true); // equivalent atoms don't need vars
					// remove atom from head of body - equivalent atoms don't need definitions
					setSimplifyHeads(bodyEqId);
				}
			}
			setHasBody(body->literal());
			if (bodyId != bodyEqId) {
				// body is replaced with bodyEqId - move heads to bodyEqId
				body->replace(*prg_->bodies_[bodyEqId], bodyEqId,*this);
				prg_->incEqs(Var_t::body_var);
				assert(body->ignore());
			}
		}
		follow_.clear();
		index = 0;
		// select next unclassified supported body
		for (; i < prg_->initialSupp_.size(); ++i) {
			bodyId  = prg_->initialSupp_[i];
			body    = prg_->bodies_[bodyId];
			if (!body->visited() && !body->ignore()) {
				follow_.push_back(bodyId);
				break;
			}
			else if (body->ignore() && body->hasVar()) {
				body->clearLiteral(false);
			}
		}
		if (follow_.empty()) break;
	}
	return simplifyClassifiedProgram(startAtom, stopAtom);
}

bool Preprocessor::simplifyClassifiedProgram(uint32 startAtom, uint32& stopAtom) {
	stopAtom = (uint32)prg_->atoms_.size();
	prg_->initialSupp_.clear();
	LitVec::size_type compute = prg_->compute_.size();
	for (uint32 i = 0; i != (uint32)prg_->bodies_.size(); ++i) {
		PrgBodyNode* b = prg_->bodies_[i];
		if (!b->visited() || !b->hasVar()) {
			// !b->visited(): body is unsupported
			// !b->hasVar() : body is eq to other body or was derived to false
			// In either case, body is no longer relevant and can be ignored.
			b->clearLiteral(true); b->setIgnore(true);
		}
		else { 
			assert(!b->ignore());
			bool hadHeads = b->hasHeads();
			std::pair<uint32, uint32> hash;
			if (nodes_[i].sBody == 1 && !b->simplifyBody(*prg_, i, hash, *this, true)) {
				return false;
			}
			if (nodes_[i].sHead == 1 && !b->simplifyHeads(*prg_, *this, true)) {
				return false;
			}
			b->setVisited(false);
			if (hadHeads && b->value() == value_false && pass_ != maxPass_) {
				// New false body. If it was derived to false, we can ignore the body.
				// If it was forced to false (i.e. it is a selfblocker), reclassify 
				// as soon as it becomes supported.
				if      (b->ignore())         { newFalseBody(b, hash.first); }
				else if (b->resetSupported()) { prg_->initialSupp_.push_back(i); }
				// Reclassify, because an atom may have lost its source
				stopAtom = 0;
			}
			else if (!b->hasHeads() && b->value() != value_false) {
				// Body is no longer needed. All heads are either superfluous or equivalent
				// to other atoms. 
				if (pass_ != maxPass_) {
					stopAtom = (getRootAtom(b->literal()) != varMax) * stopAtom;
					b->clearLiteral(true);
					b->setIgnore(true);
				}
			}
			else if (b->value() == value_true && b->var() != 0 && pass_ != maxPass_) {
				// New fact body. Merge with existing fact body if any.
				// Reclassify, if we derived at least one new fact.
				if (newFactBody(b, i, hash.first)) {
					stopAtom = 0;
				}
			}
			else if (backprop_ && pass_ != maxPass_ && !b->ignore() && b->value() != value_free && nodes_[i].sBody == 1) {
				if (!b->backpropagate(*prg_, prg_->compute_)) { return false; }
				if (prg_->compute_.size()>compute)            { stopAtom = 0; }
				if (b->resetSupported())                      { prg_->initialSupp_.push_back(i); }
			}
			else if (b->resetSupported()) {
				prg_->initialSupp_.push_back(i);
			}
			nodes_[i].sBody = 0; nodes_[i].sHead = 0;
			nodes_[i].known = 0; 
			if (pass_ != maxPass_ && hash.first != hash.second && !b->ignore()) {
				// The body has changed - remove old entry and mark
				// so that we check for an equivalent body in the next round
				removeFromIndex(b, hash.first);
				setSimplifyBody(i);
			}
		}
	}
	if (VarVec* unfreeze = prg_->incData_ ? &prg_->incData_->unfreeze_ : 0) {
		for (VarVec::const_iterator it = unfreeze->begin(), end = unfreeze->end(); it != end; ++it) {
			ValueRep x = simplifyAtom(*it, false);
			if (x != value_true) {
				if (x == value_false)  { return false; }
				if (pass_ != maxPass_) { stopAtom = 0; }
				stopAtom = 0;
			}
		}
	}
	bool more = (pass_ != maxPass_ && stopAtom == prg_->atoms_.size());
	bool clear= stopAtom != prg_->atoms_.size();
	for (uint32 i = startAtom; i != (uint32)prg_->atoms_.size(); ++i) {
		ValueRep x = simplifyAtom(Var(i), more);
		if (x != value_true) {
			if (x == value_false) { return false; }
			if (more) {
				more     = false;
				clear    = true;
				stopAtom = i;
			}
		}
		if (clear) { prg_->atoms_[i]->clearLiteral(false); }
	}
	return applyCompute(prg_->compute_,(uint32)compute,true);
}

// Derived a new fact body. The body is eq to True, therefore does not need a separate variable.
bool Preprocessor::newFactBody(PrgBodyNode* body, uint32 id, uint32 oldHash) {
	removeFromIndex(body, oldHash);
	uint32 otherId = findEqBody(body, 0);
	if (otherId != varMax) {
		PrgBodyNode* other = prg_->bodies_[otherId];
		// Found an equivalent body, merge heads... 
		uint32 normal = body->replace(*other, otherId, *this);
		if (!other->simplifyHeads(*prg_, *this, true)) {
			return false;
		}
		nodes_[otherId].sHead = 0;
		return normal > 0;
	}
	// No equivalent body found. 
	prg_->bodyIndex_.insert(ProgramBuilder::BodyIndex::value_type(0, id));
	body->resetSupported();
	prg_->initialSupp_.push_back(id);
	// Only reclassify if this body derives new facts 
	for (HeadVec::const_iterator it = body->heads_begin(), end = body->heads_end(); it != end; ++it) {
		if (it->normal()) return true;
	}
	return false;
}

// body became false after it was used to derived its heads.
void Preprocessor::newFalseBody(PrgBodyNode* body, uint32 oldHash) {
	body->clearLiteral(true); body->setIgnore(true);
	for (HeadVec::const_iterator it = body->heads_begin(), end = body->heads_end(); it != end; ++it) {
		setSimplifyBodies(it->node());
	}
	removeFromIndex(body, oldHash);
}

// Simplify the classified atom with the given id.
// Update list of bodies defining this atom and check
// if atom has a distinct var although it is eq to some body.
// Return:
//  value_false    : conflict
//  value_true     : ok
//  value_weak_true: ok but atom should be reclassified
ValueRep Preprocessor::simplifyAtom(Var id, bool reclass) {
	if (PrgAtomNode* a = prg_->atoms_[id]->hasVar() ? prg_->atoms_[id] : 0) {
		ValueRep v   = a->value();
		ValueRep ret = value_true;
		PrgAtomNode::SimpRes r(true, UINT32_MAX);
		if (nodes_[id].asBody == 1 && (r = a->simplifyBodies(id, *prg_, true)).first == false) {
			return value_false;
		}
		nodes_[id].asBody = 0; nodes_[id].aSeen = 0;
		if (a->value() != v) {
			if (a->value() == value_false) {
				prg_->compute_.push_back(negLit(id));
				ret = value_weak_true;
			}
			else if (a->value() == value_weak_true) {
				prg_->compute_.push_back(posLit(id));
			}
		}
		if (reclass && reclassify(a, id, r.second)) {
			ret = value_weak_true;
		}
		return ok_ ? ret : value_false;
	}
	return value_true;
}

// check if atom has a distinct var although it is eq to some body
bool Preprocessor::reclassify(PrgAtomNode* a, uint32 atomId, uint32 diffLits) {
	if ((a->preds.empty() && a->hasVar()) || 
		(a->preds.size() == 1
		&& a->preds[0].normal()
		&& prg_->bodies_[a->preds[0].node()]->literal() != a->literal())) {
		return true;
	}
	else if (a->preds.size() > 1 && diffLits == 1 && prg_->bodies_[a->preds[0].node()]->literal() != a->literal() && getRootAtom(a->literal()) != varMax) {
		// a is equivalent to eq
		Literal x       = prg_->bodies_[a->preds[0].node()]->literal();
		PrgAtomNode* eq = prg_->atoms_[getRootAtom(x)];
		bool stableTruth= eq->value() == value_true || eq->value() == value_false;
		a->setLiteral(x);
		if (!prg_->mergeEqAtoms(atomId, getRootAtom(x))) {
			return ok_ = false;
		}
		// update bodies containing a	
		LitVec temp; temp.reserve(a->posDep.size() + a->negDep.size());
		for (VarVec::size_type i = 0; i != a->posDep.size(); ++i) { temp.push_back(posLit(a->posDep[i])); }
		for (VarVec::size_type i = 0; i != a->negDep.size(); ++i) { temp.push_back(negLit(a->negDep[i])); }
		a->posDep.clear(); a->negDep.clear();
		for (VarVec::size_type i = 0; i != temp.size(); ++i) {
			Var bodyId      = temp[i].var();
			PrgBodyNode* bn = prg_->bodies_[ bodyId ];
			stableTruth     = stableTruth || (temp[i].sign() && eq->value() == value_weak_true);
			if (!bn->ignore()) {
				if (!stableTruth) {
					(temp[i].sign() ? eq->negDep : eq->posDep).push_back(bodyId);
				}
				bool wasSup = temp[i].sign() || bn->isSupported();
				std::pair<uint32, uint32> hash;
				if (!bn->simplifyBody(*prg_, bodyId, hash, *this, true)) {
					return ok_ = false;
				}
				if (!bn->simplifyHeads(*prg_, *this, true)) {
					return ok_ = false;
				}
				removeFromIndex(bn, hash.first);
				uint32 otherId = findEqBody(bn, hash.second);
				if (otherId != varMax) {
					mergeBodies(bodyId, otherId);
				}
				else {
					prg_->bodyIndex_.insert(ProgramBuilder::BodyIndex::value_type(hash.second, bodyId));
					if (!wasSup && bn->resetSupported()) {
						prg_->initialSupp_.push_back(bodyId);
					}
				}
			}
		}
		// remove a from heads of defining bodies
		for (VarVec::size_type i = 0; i != a->preds.size(); ++i) {
			PrgBodyNode* bn = prg_->bodies_[a->preds[i].node()];
			bn->removeHead(atomId);
			if (!bn->hasHeads()) {
				bn->setIgnore(true);
			}
		}
		return true;
	}
	return false;
}

// associates a variable with the body if necessary
uint32 Preprocessor::addBodyVar(Var bodyId, PrgBodyNode* body) {
	// make sure we don't add an irrelevant body
	assert(body->value() == value_false || body->hasHeads());
	assert(body->isSupported() && !body->eq());
	body->clearLiteral(false);        // clear var in case we are iterating
	body->setVisited(true);           // mark as seen, so we don't classify the body again
	bool changed  = nodes_[bodyId].sBody == 1;
	bool known    = nodes_[bodyId].known == body->size();
	std::pair<uint32, uint32> hashes;
	if (changed && !body->simplifyBody(*prg_, bodyId, hashes, *this, known)) {
		ok_ = false;
		return bodyId;
	}
	nodes_[bodyId].sHead |= (pass_ == 1 && (body->heads_end()-body->heads_begin())>1);
	if (nodes_[bodyId].sHead == 1){ // remove any duplicates from head
		ok_ = body->simplifyHeads(*prg_, *this, false);
		nodes_[bodyId].sHead = 0;
		assert(ok_);
		if (!body->hasHeads() && body->value() == value_free) {
			body->setIgnore(true);
			return bodyId;
		}
	}
	if (body->ignore())    { return bodyId; }
	if (changed) {
		// body has changed - check for equivalent body
		removeFromIndex(body, hashes.first);
		uint32 otherId = findEqBody(body, hashes.second);
		if (otherId != varMax) {
			// found an equivalent body - try to merge them
			ok_ = mergeBodies(bodyId, otherId);
			return body->eqNode();
		}
		// The body is still unique, remember it under its new hash
		prg_->bodyIndex_.insert(ProgramBuilder::BodyIndex::value_type(hashes.second, bodyId));
		nodes_[bodyId].sBody = backprop_ && body->value() != value_free;
	}
	if (!known) {
		body->setLiteral(body->size() > 0 
			? newBodyLit(body) 
			: posLit(0));
		setSimplifyBody(bodyId);        // simplify strongly later
		return bodyId;
	}
	if      (body->size() == 0) { 
		body->setLiteral( posLit(0) ); 
	}
	else if (body->size() == 1) { // body is equivalent to an atom or its negation
		// We know that body is equivalent to an atom. Now check if the atom is itself
		// equivalent to a body. If so, the body is equivalent to the atom's body.
		PrgAtomNode* a = 0; // eq-Atom
		PrgBodyNode* r = 0; // possible eq-body
		if (body->posSize() == 1) {
			a = prg_->atoms_[body->pos(0)];
			body->setLiteral( a->literal() );
		}
		else {
			body->setLiteral( ~prg_->atoms_[body->neg(0)]->literal() );
			Var dualAtom = getRootAtom(body->literal());
			a = dualAtom != varMax ? prg_->atoms_[dualAtom] : 0;
		}
		if (a && a->preds.size() == 1) r = prg_->bodies_[a->preds[0].node()];
		if (r && allowMerge(bodyId, a->preds[0].node())) {
			ok_ = mergeBodies(bodyId, a->preds[0].node());
			return body->eqNode();
		}
		prg_->incEqs(Var_t::atom_body_var);
		prg_->vars_.setAtomBody(body->var());
	}
	else { // body is a conjunction - start new eq class
		body->setLiteral( newBodyLit(body) );
	}
	return !body->eq() ? bodyId : body->eqNode();
}

// associates a variable with the atom if necessary
// b is the supported body that brings a into the closure
uint32 Preprocessor::addAtomVar(HeadEdge it, PrgAtomNode* a, PrgBodyNode* b) {
	Var atomId = it.node();
	if (b->value() == value_true) {
		// Since b is true, it is always a valid support for a, a can never become unfounded. 
		// So ignore it during SCC check and unfounded set computation.
		a->setIgnore(true);
		if (it.normal() && !a->setValue(value_true)) {
			ok_ = false;
			return atomId;
		}
	}
	if (nodes_[atomId].aSeen == 0) {
		// first time we visit this atom
		nodes_[atomId].aSeen = 1;
		// if set of defining bodies of this atom has changed, update it...
		if (nodes_[atomId].asBody == 1 && !a->simplifyBodies(atomId, *prg_, false).first) {
			ok_ = false;
			return atomId;
		}
		if (pass_ == 1) { nodes_[atomId].asBody = 1; }
		if (!a->hasVar()) {
			if (it.normal() && (b->value() == value_true || a->preds.size() == 1)) {
				// Atom is equivalent to b
				prg_->incEqs(Var_t::atom_body_var);
				a->setLiteral( b->literal() );
				prg_->vars_.setAtomBody(a->var());
			}
			else {
				a->setLiteral( posLit(prg_->vars_.add(Var_t::atom_var)) ); 
			}
			Var eqAtomId = getRootAtom(a->literal());
			if ( eqAtomId == varMax ) {
				setRootAtom(a->literal(), atomId);
			}
			else if (!prg_->mergeEqAtoms(atomId, eqAtomId)) {
				ok_ = false;
				return atomId;
			}
		}
		assert(getRootAtom(a->literal()) != varMax);
		return propagateAtomVar(atomId, a);
	}
	else if (hasBody(b->literal()) && a->preds.size() > 1) {
		// The eq-class of b already contains a body. Check 
		// if it is one of the other bodies defining a.
		nodes_[atomId].asBody = 1;
	}
	return prg_->getEqAtom(atomId);
}

uint32 Preprocessor::propagateAtomVar(Var atomId, PrgAtomNode* a) {
	// If atom a has a truth-value or is eq to a', we'll remove
	// it from all bodies. If there is an atom x, s.th. a.lit == ~x.lit, we mark all
	// bodies containing both a and x for simplification in order to detect
	// duplicates/contradictory body-literals.
	// In case that a == a', we also mark all bodies containing a
	// for head simplification in order to detect rules like: a' :- a,B. and a' :- B,not a.
	uint32 eqId       = prg_->getEqAtom(atomId);
	PrgAtomNode* eqAt = prg_->atoms_[eqId];
	bool stableTruth  = a->value() == value_true || a->value() == value_false;
	bool removeAtom   = stableTruth || eqAt != a;
	bool removeNeg    = removeAtom  || a->value() == value_weak_true;
	if (!removeAtom && getRootAtom(~a->literal()) != varMax) {
		PrgAtomNode* comp = prg_->atoms_[getRootAtom(~a->literal())];
		assert(a->literal() == ~comp->literal());
		// propagate any truth-value to complementary eq-class
		ValueRep cv   = value_free;
		uint32   mark = 0;
		if (a->value() != value_free && (cv = (value_false | (a->value()^value_true))) != comp->value()) {
			mark        = 1;
			ok_        &= comp->setValue(cv);
		}
		// temporarily mark all bodies containing ~a
		VarVec* deps[2] = {&comp->posDep, &comp->negDep};
		for (uint32 d = 0; d != 2; ++d) {
			for (VarVec::iterator it = deps[d]->begin(); it != deps[d]->end(); ++it) { 
				nodes_[*it].mBody  = 1;
				nodes_[*it].sBody |= mark;
			}
		}
	}
	for (VarVec::size_type i = 0; i != a->posDep.size();) {
		Var bodyId      = a->posDep[i];
		PrgBodyNode* bn = prg_->bodies_[ bodyId ];
		if (bn->ignore()) { a->posDep[i] = a->posDep.back(); a->posDep.pop_back(); continue;}
		bool wasSup = bn->isSupported();
		bool isSup  = wasSup || (a->value() != value_false && bn->onPosPredSupported(atomId));
		if (++nodes_[bodyId].known == bn->size() && !bn->visited() && isSup) {
			follow_.push_back( bodyId );
		}
		else if (isSup && !wasSup) {
			prg_->initialSupp_.push_back(bodyId);
		}
		if (removeAtom || nodes_[bodyId].mBody == 1) { 
			setSimplifyBody(bodyId);
			if (eqAt != a) {
				setSimplifyHeads(bodyId);
				if (!stableTruth) { eqAt->posDep.push_back( bodyId ); }
			}
		}
		++i;
	}
	for (VarVec::size_type i = 0; i != a->negDep.size();) {
		Var bodyId      = a->negDep[i];
		PrgBodyNode* bn = prg_->bodies_[ bodyId ];
		if (bn->ignore()) { a->negDep[i] = a->negDep.back(); a->negDep.pop_back(); continue;}
		if (++nodes_[bodyId].known == bn->size() && !bn->visited()) {
			follow_.push_back( bodyId );
		}
		if (removeNeg || nodes_[bodyId].mBody == 1) { 
			setSimplifyBody(bodyId);
			if (eqAt != a) {
				setSimplifyHeads(bodyId);
				if (!stableTruth && a->value() != value_weak_true) { eqAt->negDep.push_back( bodyId ); }
			}
		}
		++i;
	}
	if (removeNeg)  { a->negDep.clear(); }
	if (removeAtom) { a->posDep.clear(); }
	else  if  (getRootAtom(~a->literal()) != varMax) {
		// unmark  bodies containing ~a
		PrgAtomNode* comp = prg_->atoms_[getRootAtom(~a->literal())];
		VarVec* deps[2] = {&comp->posDep, &comp->negDep};
		for (uint32 d = 0; d != 2; ++d) {
			for (VarVec::iterator it = deps[d]->begin(); it != deps[d]->end(); ++it) { 
				nodes_[*it].mBody = 0;
			}
		}
	}
	return eqId;
}

// body b has changed - remove old entry from body node index
void Preprocessor::removeFromIndex(PrgBodyNode* b, uint32 hash) {
	ProgramBuilder::BodyRange ra = prg_->bodyIndex_.equal_range(hash);
	for (; ra.first != ra.second && prg_->bodies_[ra.first->second] != b; ++ra.first) {;}
	if (ra.first != ra.second) prg_->bodyIndex_.erase(ra.first);
}

// search for a body with given hash that is equal to b
uint32 Preprocessor::findEqBody(PrgBodyNode* b, uint32 hash) {
	WeightLitVec lits;
	uint32 ret = varMax;
	bool marked= false;
	ProgramBuilder::BodyRange ra = prg_->bodyIndex_.equal_range(hash);
	for (bool sorted = false; ra.first != ra.second;) {
		PrgBodyNode& other = *prg_->bodies_[ra.first->second];
		if (other.ignore()) {
			prg_->bodyIndex_.erase(ra.first++);
			continue;
		}
		if (b->type() == other.type() && b->size() == other.size() && b->posSize() == other.posSize() && b->bound() == other.bound()) {
			if (!marked) {	
				// mark all literals in body
				for (uint32 i = 0, end = b->size(); i != end; ++i) { 
					prg_->ruleState_.addToBody(b->goal(i)); 
					if (b->type() == PrgBodyNode::SUM_BODY) { lits.push_back(WeightLiteral(b->goal(i), b->weight(i))); }
				}
				marked = true;
			}
			if (prg_->allLitsMarked(other) && (b->type() != PrgBodyNode::SUM_BODY || prg_->eqWeights(other, lits, sorted))) {
				ret = ra.first->second;
				break;
			}			
		}
		++ra.first;
	}
	if (marked) {
		for (uint32 i = 0, end = b->size(); i != end; ++i) { 
			prg_->ruleState_.popFromRule(b->goal(i).var()); 
		}
	}
	return ret;
}

// check whether we can replace body with root
// for this, 
//  - the literals of the bodies must be equal
//  - the bodies must be compatible w.r.t choice semantic
//  - the merge must be safe w.r.t positive loops
bool Preprocessor::allowMerge(uint32 bodyId, uint32 rootId) {
	PrgBodyNode* body = prg_->bodies_[bodyId];
	PrgBodyNode* root = prg_->bodies_[rootId];
	if (body == root) {
		return false;
	}
	if (body->literal() == root->literal()) {
		if (root->posSize() <= body->posSize()) {
			return true;
		}
		// can't merge bodies - merge values
		if (!body->mergeValue(root)) {
			return ok_ = false;
		}
		assert(body->value() == root->value());
		if (root->value() == value_false) {
			root->propagateFalse(rootId, prg_->atoms_);
			body->propagateFalse(bodyId, prg_->atoms_);
		}
	}
	return false;
}

// body is equivalent to rootId and about to be replaced with rootId
// if one of them is a selfblocker, after the merge both are
bool Preprocessor::mergeBodies(Var bodyId, Var rootId) {
	PrgBodyNode* body = prg_->bodies_[bodyId];
	PrgBodyNode* root = prg_->bodies_[rootId];
	assert(!root->ignore() && root != body);
	assert(body->value() != value_false || !body->hasHeads());
	assert(root->value() != value_false || !root->hasHeads());
	body->setEq(rootId);
	if (body->value() != value_free && root->value() != body->value()) {
		if (!root->mergeValue(body)) return false;
		setSimplifyBody(rootId);
	}
	if (root->value() == value_false && (root->hasHeads()+body->hasHeads()) != 0) {
		PrgBodyNode* newFalse = !root->hasHeads() ? body   : root;
		Var          falseId  = newFalse == root  ? rootId : bodyId;
		newFalse->propagateFalse(falseId, prg_->atoms_);
		setSimplifyBody(rootId);
	}
	if (root->visited()) {
		body->setLiteral( root->literal() );  
	}
	else if (body->hasHeads()) {
		// root is not yet classified. Move body's heads to root so that they are added 
		// when root is eventually classified. After the move, we can ignore body.
		body->replace(*root, rootId, *this);
	}
	return true;
}
uint32 Preprocessor::replaceComp(uint32 id) const {
	PrgAtomNode* a = prg_->atoms_[id];
	if (a->eq()) {
		uint32 eqId      = a->eqNode();
		PrgAtomNode* eq  = prg_->atoms_[eqId];
		if (eq->hasVar() && getRootAtom(~eq->literal()) != varMax) {
			for (HeadVec::iterator it = a->preds.begin(), end = a->preds.end(); it != end; ++it) {
				PrgBodyNode* bn  = prg_->bodies_[it->node()];
				if (bn->eq()) {
					// predecessor is actually equivalent to some other body 
					// but the link was not yet replaced - do so now
					*it = HeadEdge(bn->eqNode(), it->type());
					bn = prg_->bodies_[bn->eqNode()];
				}
				if (bn->literal() == eq->literal() && bn->size() == 1 && bn->negSize() == 1) {
					return getRootAtom(~eq->literal());
				}
			}
		}
	}
	return varMax;
}

}
