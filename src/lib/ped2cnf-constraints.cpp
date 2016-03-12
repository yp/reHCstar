/**
 *
 *                              reHC-*
 * Haplotyping with Recombinations, Errors, and Missing Genotypes
 *
 * Copyright (C) 2011  Yuri Pirola <yuri.pirola(-at-)gmail.com>
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
 * ped2cnf-constraints.cpp
 *
 * Classes that represent constraint handling.
 *
 **/


#include "ped2cnf-constraints.hpp"

#include <cmath>

#include <boost/foreach.hpp>


void
constraint_handler_t::handle_constraints(pedcnf_t& cnf,
													  const individuals_variables_t& variables) const {
  _handle_constraints(cnf, variables);
};

void
all_false_constraints_t::_handle_constraints(pedcnf_t& cnf,
															const individuals_variables_t& variables) const {
  INFO("No " << _description << " are allowed.");
  BOOST_FOREACH(const individual_variables_t& ivar, variables) {
	 BOOST_FOREACH(const var_t& var, ivar) {
		cnf.add_clause(pedcnf_t::clause_t{ -var });
	 }
  }
};

void
at_most_individual_constraints_t::_handle_constraints(pedcnf_t& cnf,
																		const individuals_variables_t& variables) const {
  INFO("For each individual, number of " << _description << " <= "
		 << _rate << " * no. of possible " << _description);
  BOOST_FOREACH(const individual_variables_t& ivar, variables) {
	 const size_t k= std::ceil(_rate*ivar.size());
	 add_card_constraint_less_or_equal_than(cnf, ivar, k);
  }
};

void
at_most_global_constraints_t::_handle_constraints(pedcnf_t& cnf,
																  const individuals_variables_t& variables) const {
  individual_variables_t all_vars;
  BOOST_FOREACH(const individual_variables_t& ivar, variables) {
	 all_vars.insert(all_vars.end(), ivar.begin(), ivar.end());
  }
  const size_t k= std::ceil(_rate*all_vars.size());
  INFO("Globally, number of " << _description << " <= " << k <<
		 " (over " << all_vars.size() << ")");
  add_card_constraint_less_or_equal_than(cnf, all_vars, k);
};

void
at_most_individual_constraints_abs_t::_handle_constraints(pedcnf_t& cnf,
																			 const individuals_variables_t& variables) const {
  INFO("For each individual, number of " << _description << " <= " << _limit);
  BOOST_FOREACH(const individual_variables_t& ivar, variables) {
	 add_card_constraint_less_or_equal_than(cnf, ivar, _limit);
  }
};

void
at_most_global_constraints_abs_t::_handle_constraints(pedcnf_t& cnf,
																		const individuals_variables_t& variables) const {
  individual_variables_t all_vars;
  BOOST_FOREACH(const individual_variables_t& ivar, variables) {
	 all_vars.insert(all_vars.end(), ivar.begin(), ivar.end());
  }
  INFO("Globally, number of " << _description << " <= " << _limit <<
		 " (over " << all_vars.size() << ")");
  add_card_constraint_less_or_equal_than(cnf, all_vars, _limit);
};

void
interval_global_constraints_abs_t::_handle_constraints(pedcnf_t& cnf,
																		 const individuals_variables_t& variables) const {
  individual_variables_t all_vars;
  BOOST_FOREACH(const individual_variables_t& ivar, variables) {
	 all_vars.insert(all_vars.end(), ivar.begin(), ivar.end());
  }
  INFO("Globally, number of " << _description << " between " << _vmin <<
		 " and " << _vmax << " (over " << all_vars.size() << ")");
  add_card_constraint_between(cnf, all_vars, _vmin, _vmax);
};

void
at_most_windowed_constraints_t::_handle_constraints(pedcnf_t& cnf,
																	 const individuals_variables_t& variables) const {
  BOOST_FOREACH(const individual_variables_t& ivar, variables) {
	 INFO("For each window of length " << _window_length << ", "
			<< "number of " << _description << " <= " << _max_true_in_window);
	 add_uniform_card_constraint_less_or_equal_than(cnf, ivar,
																	_window_length,
																	_max_true_in_window);
  }
};

void
composite_constraints_t::_handle_constraints(pedcnf_t& cnf,
															const individuals_variables_t& variables) const {
  BOOST_FOREACH(const constraint_handler_t& handler, _handlers) {
	 handler.handle_constraints(cnf, variables);
  }
}
