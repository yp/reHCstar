/**
 *
 *                              reHC-*
 * Haplotyping with Recombinations, Errors, and Missing Genotypes
 *
 * Copyright (C) 2010,2011  Yuri Pirola <yuri.pirola(-at-)gmail.com>
 *
 * Distributed under the terms of the GNU General Public License (GPL)
 *
 *
 * This file is part of reHC-* (reHCstar),
 * previously known as ZRHC-* (ZRHCstar).
 *
 * reHC-* is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reHC-* is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with reHC-*.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 *
 * clasp_sat_solver_interface.cpp
 *
 * Structures to represent a light-weight interface to an internal SAT solver.
 * It only wants to narrow down the interface to the SAT solver.
 *
 * This file is specific for clasp (2.0.1)
 *
 **/
// Please include "sat_solver_interface.hpp" and not "xxx_sat_solver_interface.hpp"
#include "sat_solver_interface.hpp"

#ifdef INTERNAL_SAT_SOLVER
#ifdef USE_CLASP

#include <boost/foreach.hpp>

#include "clasp/solve_algorithms.h"
#include "clasp/enumerator.h"

using namespace Clasp;

void
SAT_solver_iface_t::add_clause(const std::set<lit_t>& clause) {
  L_TRACE("Adding a clause with " << clause.size() << " literals to the solver...");
  if (_solved) {
	 L_ERROR("Impossible to add a clause to an already solved instance! Abort.");
	 MY_FAIL;
	 return;
  }
  _clauses.push_back(clause.size());
  BOOST_FOREACH( lit_t lit, clause ) {
	 unsigned int var= std::abs(lit);
	 _n_var= std::max(var, _n_var);
	 _clauses.push_back(lit);
  }
};

bool
SAT_solver_iface_t::solve() {
  if (_solved) {
	 L_DEBUG("Instance already solved. Skipping...");
  } else {
	 L_DEBUG("Solving the SAT instance using the internal SAT solver '"
				<< BOOST_PP_STRINGIZE(SAT_SOLVER) << "'...");
	 ctx.reserveVars(_n_var+1);
	 while (_n_var >= ctx.numVars()) ctx.addVar(Var_t::atom_var);
	 ctx.symTab().startInit();
	 ctx.symTab().endInit(SymbolTable::map_direct, _n_var+1);
	 L_TRACE("No. of variables: " << ctx.numVars());

	 ctx.enumerator()->enumerate(1);
	 ctx.startAddConstraints();

	 ClauseCreator _cc(ctx.master());
	 while (!_solved && !_clauses.empty()) {
		lit_t next_size= _clauses.front();
		_clauses.pop_front();
		L_TRACE("Adding a new clause...")
		_cc.start();
		while (next_size>0) {
		  lit_t lit= _clauses.front();
		  _clauses.pop_front();
		  L_TRACE("   literal " << lit);
		  unsigned int var= std::abs(lit);
		  _cc.add(Literal(var, lit<0));
		  --next_size;
		}
		const bool res= _cc.end();
		if (!res) {
		  _solved= true;
		  _sat= false;
		}
	 }
	 ctx.endInit();
	 if (_solved) return _sat;

	 L_DEBUG("Solving...");
	 SolveParams params;
	 const bool ret= Clasp::solve(ctx, params);
	 L_DEBUG("...finished!");
	 L_DEBUG("The internal SAT solver gave: " << (ret ? "TRUE" : "FALSE"));

	 if ( ret ) {
		_sat= true;
		_solver= ctx.master();
		for (var_t var= 1; var <= _n_var; ++var) {
		  if (_solver->value(var) == value_free) {
			 L_TRACE("Setting free variable " << var <<
						" as false (and propagating).");
			 _solver->assume(negLit(var));
			 _solver->propagate();
		  }
		}
	 } else {
		_sat= false;
	 }
	 _solved= true;
  }
  return _sat;
};


bool
SAT_solver_iface_t::model(const var_t _var) {
  const var_t var= _var+1;
  if (!_solved || !_sat) {
	 L_ERROR("The instance has not yet been solved or it is not satisfiable!");
	 MY_FAIL;
  }
  if (var > _n_var) {
	 L_ERROR("Asking for a not-existent variable " << var << ". "
				"Maximum variables: " << _n_var);
	 MY_FAIL;
  }
  const ValueRep val= ctx.master()->value(var);
  if (val ==  value_true)
	 return true;
  else if (val == value_false)
	 return false;
  else {
	 L_ERROR("Variable " << var << " is (still) free!");
	 MY_FAIL;
	 return false;
  }
};

#endif
#endif
