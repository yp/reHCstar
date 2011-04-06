/**
 *
 *                              ZRHC-*
 * Zero-Recombinant Haplotype Configuration with missing genotypes
 *
 * Copyright (C) 2010  Yuri Pirola <yuri.pirola(-at-)gmail.com>
 *
 * Distributed under the terms of the GNU General Public License (GPL)
 *
 *
 * This file is part of ZRHC-* (ZRHCstar).
 *
 * ZRHC-* is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ZRHC-* is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ZRHC-*.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 *
 * ped2cnf-errors.cpp
 *
 * Classes that represent error handling.
 *
 **/


#include "ped2cnf-errors.hpp"

#include <cmath>

#include <boost/foreach.hpp>


void
error_handler_t::prepare_errors(const pedcnf_t& cnf,
										  const size_t pedigree_size,
										  const size_t genotype_length,
										  individuals_errors_t& errors) const {
  for (size_t i= 0; i<pedigree_size; ++i) {
	 {
		individual_errors_t ies;
		errors.push_back(ies);
	 }
	 individual_errors_t& ies= errors.back();
	 for (size_t l= 0; l<genotype_length; ++l) {
		if (cnf.has_e(i, l)) {
		  ies.push_back(cnf.get_e(i, l));
		}
	 }
  }
};

void
error_handler_t::handle_errors(pedcnf_t& cnf,
										 const individuals_errors_t& errors) const {
  _handle_errors(cnf, errors);
};

void
no_errors_handler_t::_handle_errors(pedcnf_t& cnf,
												const individuals_errors_t& errors) const {
  BOOST_FOREACH(const individual_errors_t& ierr, errors) {
	 BOOST_FOREACH(const var_t& serr, ierr) {
		cnf.add_clause<1>((lit_t[]){ -serr });
	 }
  }
};

void
whole_individual_genotype_error_handler_t::_handle_errors(pedcnf_t& cnf,
																			 const individuals_errors_t& errors) const {
  BOOST_FOREACH(const individual_errors_t& ierr, errors) {
	 const size_t k= std::ceil(_error_rate*ierr.size());
	 DEBUG("Generating cardinality constraints for " << ierr.size() <<
			 " error variables (<= " << k << ")...");
	 add_card_constraint_less_or_equal_than(cnf, ierr, k);
  }
};

void
whole_pedigree_genotype_error_handler_t::_handle_errors(pedcnf_t& cnf,
																		  const individuals_errors_t& errors) const {
  individual_errors_t all_errs;
  BOOST_FOREACH(const individual_errors_t& ierr, errors) {
	 all_errs.insert(all_errs.end(), ierr.begin(), ierr.end());
  }
  const size_t k= std::ceil(_error_rate*all_errs.size());
  DEBUG("Generating cardinality constraints for " << all_errs.size() <<
		  " error variables (<= " << k << ")...");
  add_card_constraint_less_or_equal_than(cnf, all_errs, k);
};

void
windowed_error_handler_t::_handle_errors(pedcnf_t& cnf,
													  const individuals_errors_t& errors) const {
  BOOST_FOREACH(const individual_errors_t& ierr, errors) {
	 DEBUG("Generating uniform cardinality constraints for " << ierr.size() <<
			 " error variables in windows of length " << _window_length <<
			 " requiring at most " << _max_error_in_window << " errors...");
	 add_uniform_card_constraint_less_or_equal_than(cnf, ierr,
																	_window_length,
																	_max_error_in_window);
  }
};

void
composite_error_handler_t::_handle_errors(pedcnf_t& cnf,
														const individuals_errors_t& errors) const {
  BOOST_FOREACH(const error_handler_t& handler, _handlers) {
	 handler.handle_errors(cnf, errors);
  }
}
