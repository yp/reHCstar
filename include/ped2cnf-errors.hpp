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
 * ped2cnf-errors.hpp
 *
 * Classes that represent error handling.
 *
 **/

#ifndef __PED2CNF_ERRORS_HPP__
#define __PED2CNF_ERRORS_HPP__

#include "log.hpp"
#include "utility.hpp"
#include "pedcnf.hpp"

#include <vector>

#include <boost/utility.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

class error_handler_t
  :boost::noncopyable {
public:

  typedef std::vector< var_t > individual_errors_t;

  typedef std::vector< individual_errors_t > individuals_errors_t;

private:

  virtual void
  _handle_errors(pedcnf_t& cnf,
					  const individuals_errors_t& errors) const = 0;

protected:
  my_logger logger;

public:

  error_handler_t()
		:logger(get_my_logger("error-handler"))
  {};

  virtual
  ~error_handler_t()
  {};

  void
  handle_errors(pedcnf_t& cnf,
					 const individuals_errors_t& errors) const;

  void
  prepare_errors(const pedcnf_t& cnf,
					  const size_t pedigree_size,
					  const size_t genotype_length,
					  individuals_errors_t& errors) const;

};

class no_errors_handler_t
  :public error_handler_t
{
private:

  virtual void
  _handle_errors(pedcnf_t& cnf, const individuals_errors_t& errors) const;

};

class whole_individual_genotype_error_handler_t
  :public error_handler_t
{
private:
  const double _error_rate;

  virtual void
  _handle_errors(pedcnf_t& cnf, const individuals_errors_t& errors) const;

public:
  whole_individual_genotype_error_handler_t(const double error_rate)
		:_error_rate(error_rate)
  {};

};

class whole_pedigree_genotype_error_handler_t
  :public error_handler_t
{
private:
  const double _error_rate;

  virtual void
  _handle_errors(pedcnf_t& cnf, const individuals_errors_t& errors) const;

public:
  whole_pedigree_genotype_error_handler_t(const double error_rate)
		:_error_rate(error_rate)
  {};

};


class windowed_error_handler_t
  :public error_handler_t
{
private:
  const unsigned int _max_error_in_window;
  const unsigned int _window_length;

  virtual void
  _handle_errors(pedcnf_t& cnf, const individuals_errors_t& errors) const;

public:
  windowed_error_handler_t(const int max_error_in_window,
									const int window_length)
		:_max_error_in_window(max_error_in_window),
		 _window_length(window_length)
  {};
};


class composite_error_handler_t
  :public error_handler_t
{
private:
  boost::ptr_vector<error_handler_t> _handlers;

  virtual void
  _handle_errors(pedcnf_t& cnf, const individuals_errors_t& errors) const;

public:
  void add(error_handler_t* handler) {
	 _handlers.push_back(handler);
  }

};


#endif // __PED2CNF_ERRORS_HPP__
