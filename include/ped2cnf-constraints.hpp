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
 * ped2cnf-constraints.cpp
 *
 * Classes that represent constraint handling.
 *
 **/

#ifndef __PED2CNF_CONSTRAINTS_HPP__
#define __PED2CNF_CONSTRAINTS_HPP__

#include "log.hpp"
#include "utility.hpp"
#include "pedcnf.hpp"

#include <vector>

#include <boost/utility.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

class constraint_handler_t
  :boost::noncopyable {
public:

  typedef std::vector< var_t > individual_variables_t;

  typedef std::vector< individual_variables_t > individuals_variables_t;

private:

  virtual void
  _handle_constraints(pedcnf_t& cnf,
							 const individuals_variables_t& variables) const = 0;

protected:
  my_logger logger;

  const std::string _description;

public:

  constraint_handler_t(const std::string& description="true variables")
		:logger(get_my_logger("constraint-handler")), _description(description)
  {};

  virtual
  ~constraint_handler_t()
  {};

  void
  handle_constraints(pedcnf_t& cnf,
							const individuals_variables_t& variables) const;

};

class all_false_constraints_t
  :public constraint_handler_t
{
private:

  virtual void
  _handle_constraints(pedcnf_t& cnf, const individuals_variables_t& variables) const;

public:
  all_false_constraints_t(const std::string& description)
		:constraint_handler_t(description)
  {};
};

class at_most_individual_constraints_t
  :public constraint_handler_t
{
private:
  const double _rate;

  virtual void
  _handle_constraints(pedcnf_t& cnf, const individuals_variables_t& variables) const;

public:
  at_most_individual_constraints_t(const double rate,
											  const std::string& description="true variables")
		:constraint_handler_t(description), _rate(rate)
  {};

};

class at_most_global_constraints_t
  :public constraint_handler_t
{
private:
  const double _rate;

  virtual void
  _handle_constraints(pedcnf_t& cnf, const individuals_variables_t& variables) const;

public:
  at_most_global_constraints_t(const double rate,
										 const std::string& description="true variables")
		:constraint_handler_t(description), _rate(rate)
  {};

};

class at_most_individual_constraints_abs_t
  :public constraint_handler_t
{
private:
  const unsigned int _limit;

  virtual void
  _handle_constraints(pedcnf_t& cnf, const individuals_variables_t& variables) const;

public:
  at_most_individual_constraints_abs_t(const unsigned int limit,
													const std::string& description="true variables")
		:constraint_handler_t(description), _limit(limit)
  {};

};

class at_most_global_constraints_abs_t
  :public constraint_handler_t
{
private:
  const unsigned int _limit;

  virtual void
  _handle_constraints(pedcnf_t& cnf, const individuals_variables_t& variables) const;

public:
  at_most_global_constraints_abs_t(const unsigned int limit,
											  const std::string& description="true variables")
		:constraint_handler_t(description), _limit(limit)
  {};

};

class interval_global_constraints_abs_t
  :public constraint_handler_t
{
private:
  const unsigned int _vmin;
  const unsigned int _vmax;

  virtual void
  _handle_constraints(pedcnf_t& cnf, const individuals_variables_t& variables) const;

public:
  interval_global_constraints_abs_t(const unsigned int vmin, const unsigned int vmax,
												const std::string& description="true variables")
		:constraint_handler_t(description), _vmin(vmin), _vmax(vmax)
  {};

};


class at_most_windowed_constraints_t
  :public constraint_handler_t
{
private:
  const unsigned int _max_true_in_window;
  const unsigned int _window_length;

  virtual void
  _handle_constraints(pedcnf_t& cnf, const individuals_variables_t& variables) const;

public:
  at_most_windowed_constraints_t(const int max_true_in_window,
											const int window_length,
											const std::string& description="true variables")
		:constraint_handler_t(description),
		 _max_true_in_window(max_true_in_window),
		 _window_length(window_length)
  {};
};


class composite_constraints_t
  :public constraint_handler_t
{
private:
  boost::ptr_vector<constraint_handler_t> _handlers;

  virtual void
  _handle_constraints(pedcnf_t& cnf, const individuals_variables_t& variables) const;

public:
  composite_constraints_t()
		:constraint_handler_t("")
  {};

  void add(constraint_handler_t* handler) {
	 _handlers.push_back(handler);
  }

};

class error_handler_t {
private:

  typedef constraint_handler_t::individuals_variables_t individuals_variables_t;
  typedef constraint_handler_t::individual_variables_t individual_variables_t;

  const constraint_handler_t& _handler;

public:
  error_handler_t(const constraint_handler_t& handler)
		:_handler(handler)
  {};

  void handle_errors(pedcnf_t& cnf,
							const size_t pedigree_size,
							const size_t genotype_length) const {
	 individuals_variables_t variables;
	 for (size_t i= 0; i<pedigree_size; ++i) {
		{
		  individual_variables_t ivs;
		  variables.push_back(ivs);
		}
		individual_variables_t& ivs= variables.back();
		for (size_t l= 0; l<genotype_length; ++l) {
		  if (cnf.has_e(i, l)) {
			 ivs.push_back(cnf.get_e(i, l));
		  }
		}
	 }
	 _handler.handle_constraints(cnf, variables);
  };
};

class global_recombination_handler_t {
private:

  typedef constraint_handler_t::individuals_variables_t individuals_variables_t;
  typedef constraint_handler_t::individual_variables_t individual_variables_t;

  const constraint_handler_t& _handler;

public:
  global_recombination_handler_t(const constraint_handler_t& handler)
		:_handler(handler)
  {};

  void handle_recombinations(pedcnf_t& cnf,
									  const size_t pedigree_size,
									  const size_t genotype_length) const {
	 individuals_variables_t variables;
// Recombinations in both haplotypes
	 for (size_t i= 0; i<pedigree_size; ++i) {
		{
		  individual_variables_t ivs;
		  variables.push_back(ivs);
		}
		individual_variables_t& ivs= variables.back();
		for (size_t l= 1; l<genotype_length; ++l) {
		  if (cnf.has_rp(i, l)) {
			 ivs.push_back(cnf.get_rp(i, l));
		  }
		  if (cnf.has_rm(i, l)) {
			 ivs.push_back(cnf.get_rm(i, l));
		  }
		}
	 }
	 _handler.handle_constraints(cnf, variables);
  };
};

class separated_recombination_handler_t {
private:

  typedef constraint_handler_t::individuals_variables_t individuals_variables_t;
  typedef constraint_handler_t::individual_variables_t individual_variables_t;

  const constraint_handler_t& _handler;

public:
  separated_recombination_handler_t(const constraint_handler_t& handler)
		:_handler(handler)
  {};

  void handle_recombinations(pedcnf_t& cnf,
									  const size_t pedigree_size,
									  const size_t genotype_length) const {
	 individuals_variables_t variables;
// Recombinations in paternal haplotypes
	 for (size_t i= 0; i<pedigree_size; ++i) {
		{
		  individual_variables_t ivs;
		  variables.push_back(ivs);
		}
		individual_variables_t& ivs= variables.back();
		for (size_t l= 1; l<genotype_length; ++l) {
		  if (cnf.has_rp(i, l)) {
			 ivs.push_back(cnf.get_rp(i, l));
		  }
		}
	 }
	 _handler.handle_constraints(cnf, variables);
// Recombinations in maternal haplotypes
	 for (size_t i= 0; i<pedigree_size; ++i) {
		individual_variables_t& ivs= variables[i];
		ivs.clear();
		for (size_t l= 1; l<genotype_length; ++l) {
		  if (cnf.has_rm(i, l)) {
			 ivs.push_back(cnf.get_rm(i, l));
		  }
		}
	 }
	 _handler.handle_constraints(cnf, variables);
  };
};

#endif // __PED2CNF_CONSTRAINTS_HPP__
