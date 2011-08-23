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
 * clasp_sat_solver_interface.hpp
 *
 * Structures to represent a light-weight interface to an internal SAT solver.
 * It only wants to narrow down the interface to the SAT solver.
 *
 * This file is specific for clasp (2.0.1)
 *
 **/

#ifndef __CLASP_SAT_SOLVER_INTERFACE_HPP__
#define __CLASP_SAT_SOLVER_INTERFACE_HPP__

#include "log.hpp"
#include "utility.hpp"

#include <set>
#include <deque>

#include "clasp/solver.h"
#include "clasp/clause.h"
#include "clasp/satelite.h"


class SAT_solver_iface_t
  :
  public log_able_t<SAT_solver_iface_t>
{
private:

  bool _solved;
  bool _sat;
  bool _started;

  unsigned int _n_var;

  std::deque<lit_t> _clauses;
  Clasp::SharedContext ctx;
  Clasp::Solver* _solver;

public:
  SAT_solver_iface_t()
		:_solved(false), _sat(false), _started(false), _n_var(0),
		 _solver(NULL)
  {
	 Clasp::SatElite::SatElite* pre = new Clasp::SatElite::SatElite();
pre->options.maxIters = 20;
pre->options.maxOcc   = 25;
pre->options.maxTime  = 120;
pre->options.maxFrozen= 1.0;
ctx.satPrepro.reset(pre);

  };

  ~SAT_solver_iface_t() {
  };

  bool is_solved() const {
	 return _solved;
  };

  void add_clause(const std::set<lit_t>& clause);

  bool solve();

  bool model(const var_t var);

  size_t no_of_vars() const {
	 return _n_var;
  };

};

#endif // __CLASP_SAT_SOLVER_INTERFACE_HPP__
