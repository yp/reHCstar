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
 * ms_sat_solver_interface.cpp
 *
 * Structures to represent a light-weight interface to an internal SAT solver.
 * It only wants to narrow down the interface to the SAT solver.
 *
 * This file is specific for MiniSat (2.2.0)
 *
 **/
// Please include "sat_solver_interface.hpp" and not "xxx_sat_solver_interface.hpp"
#include "sat_solver_interface.hpp"

#ifdef INTERNAL_SAT_SOLVER
#ifdef USE_MINISAT

#include "Vec.h"
#include "SolverTypes.h"

#include <boost/foreach.hpp>

using namespace Minisat;

void
SAT_solver_iface_t::add_clause(const std::set<lit_t>& clause) {
  L_TRACE("Adding a clause to the solver...");
  if (_solved) {
	 L_ERROR("Impossible to add a clause to an already solved instance! Abort.");
	 MY_FAIL;
	 return;
  }
  vec<Lit> sc;
  BOOST_FOREACH( lit_t lit, clause ) {
	 unsigned int var= std::abs(lit)-1;
	 while (var >= _solver->nVars()) _solver->newVar();
	 sc.push( mkLit( var, lit<0 ) );
  }
  _solver->addClause(sc);
};

bool
SAT_solver_iface_t::solve() {
  if (_solved) {
	 L_DEBUG("Instance already solved. Skipping...");
  } else {
	 lbool ret;
	 L_DEBUG("Solving the SAT instance using the internal SAT solver '"
				<< BOOST_PP_STRINGIZE(SAT_SOLVER) << "'...");
	 L_DEBUG("Simplification...");
	 _solver->eliminate(true);
	 L_DEBUG("Simplification terminated.");
	 if (!_solver->okay()) {
		L_DEBUG("Unsatisfiable by simplification!");
		ret= l_False;
	 } else {
		L_DEBUG("Solving...");
		vec<Lit> dummy;
		ret= _solver->solveLimited(dummy);
		L_DEBUG("...finished!");
	 }
	 L_DEBUG("The internal SAT solver gave: " <<
				( (ret == l_True) ?
				  "TRUE" :
				  ((ret == l_False) ?
					"FALSE" :
					((ret == l_Undef) ?
					 "UNDEF" :
					 "???") ) ) );
	 if ( ret == l_True ) {
		_sat= true;
	 } else {
		_sat= false;
	 }
	 _solved= true;
  }
  return _sat;
};


bool
SAT_solver_iface_t::model(const var_t var) const {
  if (!_solved || !_sat) {
	 L_ERROR("The instance has not yet solved or it is not satisfiable!");
	 MY_FAIL;
  }
  if (var >= _solver->nVars()) {
	 L_ERROR("Asking for a not-existent variable " << var << ". "
				"Maximum variables: " << _solver->nVars());
	 MY_FAIL;
  }
  MY_ASSERT(_solver->model[var] != l_Undef);
  return _solver->model[var] == l_True;
};

#endif
#endif
