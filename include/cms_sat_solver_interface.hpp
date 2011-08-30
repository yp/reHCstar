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
 * cms_sat_solver_interface.hpp
 *
 * Structures to represent a light-weight interface to an internal SAT solver.
 * It only wants to narrow down the interface to the SAT solver.
 *
 * This file is specific for CryptoMiniSat (>=2.7)
 *
 **/

#ifndef __CMS_SAT_SOLVER_INTERFACE_HPP__
#define __CMS_SAT_SOLVER_INTERFACE_HPP__

#include "log.hpp"
#include "utility.hpp"

#include <set>

#include "Solver.h"

using namespace CMSat;


class SAT_solver_iface_t
  :
  public log_able_t<SAT_solver_iface_t>
{
private:
  Solver* _solver;

  bool _solved;
  bool _sat;

public:
  SAT_solver_iface_t()
		:_solver(new Solver()), _solved(false), _sat(false)
  {
	 SolverConf& conf= _solver->conf;
	 conf.verbosity= 3;
	 conf.doFindXors= false;
	 conf.doFindEqLits= false;
	 conf.doRegFindEqLits= false;
	 conf.libraryUsage= true;
	 conf.restrictPickBranch= 0;
#ifdef TRY_LOW_MEMORY_USAGE
	 conf.doSatELite= false;
#endif
  };

  ~SAT_solver_iface_t() {
	 delete _solver;
  };

  bool is_solved() const {
	 return _solved;
  };

  void add_clause(const std::set<lit_t>& clause);

#ifndef AVOID_XOR_CLAUSES
  void add_xor_clause(const std::set<lit_t>& clause);
#endif

  bool solve();

  bool model(const var_t var) const;

  size_t no_of_vars() const {
	 return _solver->nVars();
  };

};

#endif // __CMS_SAT_SOLVER_INTERFACE_HPP__
