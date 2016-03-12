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
 * pedcnf.hpp
 *
 * Structures to represent SAT instances derived from pedigrees.
 *
 **/

#ifndef __PEDCNF_HPP__
#define __PEDCNF_HPP__


#include "data.hpp"
#include "log.hpp"
#include "utility.hpp"

// Include the SAT solver interface (if asked to do so)
#include "sat_solver_interface.hpp"

#include <map>
#include <set>
#include <vector>
#include <deque>
#include <ostream>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/logic/tribool.hpp>


class ped_var_kind
  :public enum_like_t<ped_var_kind, 10, 10>
{
private:

  typedef enum_like_t<ped_var_kind, 10, 10> base;

  ped_var_kind(const int val)
		:base(val)
  {};

public:

  ped_var_kind(const ped_var_kind& pvk)
		:base(pvk)
  {};

  static const ped_var_kind SP; // Grand-parental source (from father)
  static const ped_var_kind SM; // Grand-parental source (from mother)
  static const ped_var_kind P; // Paternal allele
  static const ped_var_kind M; // Maternal allele
  static const ped_var_kind RP; // Recombination events (from father)
  static const ped_var_kind RM; // Recombination events (from mother)
  static const ped_var_kind E; // Errors
  static const ped_var_kind PM; // Paternal multi allele
  static const ped_var_kind MM; // Maternal multi allele
  static const ped_var_kind DUMMY; // Dummy

  static const int int_values[];
  static const std::string str_values[];
  static const ped_var_kind enum_values[];

};


class pedcnf_t
  :
  public log_able_t<pedcnf_t>
{

// Types
private:

  typedef boost::tuple<size_t, size_t, size_t> index_var_t;

public:

  typedef boost::tuple<ped_var_kind, size_t, size_t, size_t> pedvar_t;
  typedef std::map<index_var_t, lit_t> varmap_t;
  typedef std::vector<pedvar_t> varvec_t;
  typedef std::vector<bool> valvec_t;
  typedef std::set<lit_t> clause_t;
#ifndef AVOID_XOR_CLAUSES
  typedef std::set<lit_t> xor_clause_t;
#endif

#ifndef ONLY_INTERNAL_SAT_SOLVER
  typedef std::deque< clause_t > clauses_t;
#ifndef AVOID_XOR_CLAUSES
  typedef std::deque< xor_clause_t > xor_clauses_t;
#endif
#endif // ONLY_INTERNAL_SAT_SOLVER

// Data members
private:

  varmap_t _sp; // Grand-parental source (from father)
  varmap_t _sm; // Grand-parental source (from mother)
  varmap_t _p; // Paternal allele
  varmap_t _m; // Maternal allele
  varmap_t _rp; // Recombination events (from father)
  varmap_t _rm; // Recombination events (from mother)
  varmap_t _e; // Errors
  varmap_t _pm; // Paternal allele
  varmap_t _mm; // Maternal allele
  size_t _next_dummy;

  varvec_t _vars;
  valvec_t _vals;

  size_t _no_of_clauses;
#ifndef AVOID_XOR_CLAUSES
  size_t _no_of_xor_clauses;
#endif

#ifndef ONLY_INTERNAL_SAT_SOLVER
  clauses_t _clauses;
#ifndef AVOID_XOR_CLAUSES
  xor_clauses_t _xor_clauses;
#endif
#endif // ONLY_INTERNAL_SAT_SOLVER


#ifdef INTERNAL_SAT_SOLVER

protected:
  SAT_solver_iface_t _solver;

#endif // INTERNAL_SAT_SOLVER

public:


// Methods
private:

  lit_t get_var3(varmap_t& map,
					  const ped_var_kind& var_kind,
					  const size_t i1, const size_t i2, const size_t i3);

  lit_t get_var3(const varmap_t& map,
					  const size_t i1, const size_t i2, const size_t i3) const;

  bool has_var3(const varmap_t& map,
					 const size_t i1, const size_t i2, const size_t i3) const;

  bool get_val3(const varmap_t& map,
					 const size_t i1, const size_t i2, const size_t i3) const;

  lit_t get_var(varmap_t& map,
					 const ped_var_kind& var_kind,
					 const size_t i1, const size_t i2) {
	 return get_var3(map, var_kind, i1, i2, 0);
  };

  lit_t get_var(const varmap_t& map,
					 const size_t i1, const size_t i2) const {
	 return get_var3(map, i1, i2, 0);
  };

  bool has_var(const varmap_t& map,
					const size_t i1, const size_t i2) const {
	 return has_var3(map, i1, i2, 0);
  };

  bool get_val(const varmap_t& map,
					const size_t i1, const size_t i2) const {
	 return get_val3(map, i1, i2, 0);
  };


public:

  pedcnf_t()
		:_next_dummy(0), _no_of_clauses(0)
#ifndef AVOID_XOR_CLAUSES
		, _no_of_xor_clauses(0)
#endif
  {};

  ~pedcnf_t() {
  };

  lit_t get_sp(const size_t i, const size_t l);

  lit_t get_sm(const size_t i, const size_t l);

  lit_t get_p(const size_t i, const size_t l);

  lit_t get_m(const size_t i, const size_t l);

  lit_t get_rp(const size_t i, const size_t l);

  lit_t get_rm(const size_t i, const size_t l);

  lit_t get_pm(const size_t i, const size_t l, const size_t j);

  lit_t get_mm(const size_t i, const size_t l, const size_t j);

  lit_t get_e(const size_t i, const size_t l);

  bool has_sp(const size_t i, const size_t l) const;

  bool has_sm(const size_t i, const size_t l) const;

  bool has_p(const size_t i, const size_t l) const;

  bool has_m(const size_t i, const size_t l) const;

  bool has_rp(const size_t i, const size_t l) const;

  bool has_rm(const size_t i, const size_t l) const;

  bool has_pm(const size_t i, const size_t l, const size_t j) const;

  bool has_mm(const size_t i, const size_t l, const size_t j) const;

  bool has_e(const size_t i, const size_t l) const;

  lit_t generate_dummy();

  lit_t get_sp(const size_t i, const size_t l) const;

  lit_t get_sm(const size_t i, const size_t l) const;

  lit_t get_p(const size_t i, const size_t l) const;

  lit_t get_m(const size_t i, const size_t l) const;

  lit_t get_rp(const size_t i, const size_t l) const;

  lit_t get_rm(const size_t i, const size_t l) const;

  lit_t get_pm(const size_t i, const size_t l, const size_t j) const;

  lit_t get_mm(const size_t i, const size_t l, const size_t j) const;

  lit_t get_e(const size_t i, const size_t l) const;

  bool sp(const size_t i, const size_t l) const {
	 return get_val(_sp, i, l);
  };

  bool sm(const size_t i, const size_t l) const {
	 return get_val(_sp, i, l);
  };

  bool p(const size_t i, const size_t l) const {
	 return get_val(_p, i, l);
  };

  bool m(const size_t i, const size_t l) const {
	 return get_val(_m, i, l);
  };

  bool rp(const size_t i, const size_t l) const {
	 return get_val(_rp, i, l);
  };

  bool rm(const size_t i, const size_t l) const {
	 return get_val(_rm, i, l);
  };

  bool pm(const size_t i, const size_t l, const size_t j) const {
	 return get_val3(_pm, i, l, j);
  };

  bool mm(const size_t i, const size_t l, const size_t j) const {
	 return get_val3(_mm, i, l, j);
  };

  bool e(const size_t i, const size_t l) const {
	 return get_val(_e, i, l);
  };

  const varvec_t& vars() const {
	 return _vars;
  };

  const valvec_t& vals() const {
	 return _vals;
  };

#ifndef ONLY_INTERNAL_SAT_SOLVER

  const clauses_t& clauses() const {
	 return _clauses;
  };

#ifndef AVOID_XOR_CLAUSES
  const xor_clauses_t& xor_clauses() const {
	 return _xor_clauses;
  };
#endif

#endif // ONLY_INTERNAL_SAT_SOLVER

  size_t no_of_clauses() const {
#ifndef AVOID_XOR_CLAUSES
	 return _no_of_clauses + _no_of_xor_clauses;
#else
	 return _no_of_clauses;
#endif
  };

  void add_clause(const clause_t& clause);

#ifndef AVOID_XOR_CLAUSES
  void add_xor_clause(const xor_clause_t& clause);
#endif

  void add_clause(const lit_t* const clause, const size_t LEN) {
	 add_clause(clause_t(clause, clause+LEN));
  };

#ifndef ONLY_INTERNAL_SAT_SOLVER
  bool is_satisfying_assignment() const;

  std::ostream& clauses_to_dimacs_format(std::ostream& out) const;
  std::ostream& clauses_to_dimacs_format(std::ostream& out,
																 const std::string& note) const;
  std::ostream& clauses_to_dimacs_format(std::ostream& out,
																 const std::vector< std::string >& notes) const;

  std::string clauses_to_dimacs_format() const {
	 std::ostringstream out;
	 clauses_to_dimacs_format(out);
	 return out.str();
  };
  std::string clauses_to_dimacs_format(const std::string& note) const {
	 std::ostringstream out;
	 clauses_to_dimacs_format(out, note);
	 return out.str();
  };
  std::string clauses_to_dimacs_format(const std::vector< std::string > notes) const {
	 std::ostringstream out;
	 clauses_to_dimacs_format(out, notes);
	 return out.str();
  };
#endif // ONLY_INTERNAL_SAT_SOLVER

// Read the assignment from a file like the following one:
// SAT/UNSAT
// 1 -2 3 4 0
  boost::tribool assignment_from_minisat_format(std::istream& in);


#ifdef INTERNAL_SAT_SOLVER
  bool solve() {
	 const bool ret= _solver.solve();
	 if (ret) {
		for (unsigned int var = 0; var != _solver.no_of_vars(); var++) {
		  _vals[var]= _solver.model(var);
		}
	 } else {
	 }
	 return ret;
  };
#endif // INTERNAL_SAT_SOLVER

};



std::ostream&
operator<<(std::ostream& out, const pedcnf_t::pedvar_t& var);

std::ostream&
operator<<(std::ostream& out, const pedcnf_t::clause_t& clause);


void
add_card_constraint_less_or_equal_than(pedcnf_t& cnf,
													const std::vector<var_t>& in_vars,
													const size_t k);

void
add_card_constraint_between(pedcnf_t& cnf,
									 const std::vector<var_t>& in_vars,
									 const size_t k1, const size_t k2);

void
add_uniform_card_constraint_less_or_equal_than(pedcnf_t& cnf,
															  const std::vector<var_t>& in_vars,
															  const size_t window_size,
															  const size_t k);


#endif // __PEDCNF_HPP__
