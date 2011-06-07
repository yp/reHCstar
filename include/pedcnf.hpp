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
#include <ostream>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>


class ped_var_kind_reHCstar
  :public enum_like_t<ped_var_kind_reHCstar, 8, 8>
{
private:

  typedef enum_like_t<ped_var_kind_reHCstar, 8, 8> base;

  ped_var_kind_reHCstar(const int val)
		:base(val)
  {};

public:

  ped_var_kind_reHCstar(const ped_var_kind_reHCstar& pvk)
		:base(pvk)
  {};

  static const ped_var_kind_reHCstar SP; // Grand-parental source (from father)
  static const ped_var_kind_reHCstar SM; // Grand-parental source (from mother)
  static const ped_var_kind_reHCstar P; // Paternal allele
  static const ped_var_kind_reHCstar M; // Maternal allele
  static const ped_var_kind_reHCstar RP; // Recombination events (from father)
  static const ped_var_kind_reHCstar RM; // Recombination events (from mother)
  static const ped_var_kind_reHCstar E; // Errors
  static const ped_var_kind_reHCstar DUMMY; // Dummy

  static const int int_values[];
  static const std::string str_values[];
  static const ped_var_kind_reHCstar enum_values[];

};

/*
  class ped_var_kind_r0HCstar
  :public enum_like_t<ped_var_kind_r0HCstar, 8, 8>
{
private:

  typedef enum_like_t<ped_var_kind_r0HCstar, 7, 7> base;

  ped_var_kind_r0HCstar(const int val)
		:base(val)
  {};

public:

  ped_var_kind_r0HCstar(const ped_var_kind_r0HCstar& pvk)
		:base(pvk)
  {};

  static const ped_var_kind_r0HCstar SP; // Grand-parental source (from father)
  static const ped_var_kind_r0HCstar SM; // Grand-parental source (from mother)
  static const ped_var_kind_r0HCstar H; // Paternal haplotype
  static const ped_var_kind_r0HCstar W; // Homozygosity status
  static const ped_var_kind_r0HCstar RP; // Recombination events (from father)
  static const ped_var_kind_r0HCstar RM; // Recombination events (from mother)
  static const ped_var_kind_r0HCstar DUMMY; // Dummy

  static const int int_values[];
  static const std::string str_values[];
  static const ped_var_kind_r0HCstar enum_values[];

};
*/

template <class KIND>
class pedcnf_t
  :
  public log_able_t<pedcnf_t>
{

// Types
private:

  typedef boost::tuple<size_t, size_t> index_var_t;

public:

  typedef boost::tuple<KIND, size_t, size_t> pedvar_t;
  typedef std::map<index_var_t, lit_t> varmap_t;
  typedef std::vector<pedvar_t> varvec_t;
  typedef std::vector<bool> valvec_t;
  typedef std::set<lit_t> clause_t;
#ifndef AVOID_XOR_CLAUSES
  typedef std::set<lit_t> xor_clause_t;
#endif

#ifndef ONLY_INTERNAL_SAT_SOLVER
  typedef std::set< clause_t > clauses_t;
#ifndef AVOID_XOR_CLAUSES
  typedef std::set< xor_clause_t > xor_clauses_t;
#endif
#endif // ONLY_INTERNAL_SAT_SOLVER

// Data members
private:

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
protected:

  lit_t get_var(varmap_t& map,
					 const KIND& var_kind,
					 const size_t i1, const size_t i2);

  lit_t get_var(const varmap_t& map,
					 const size_t i1, const size_t i2) const;

  bool has_var(const varmap_t& map,
					const size_t i1, const size_t i2) const;

  bool get_val(const varmap_t& map,
					const size_t i1, const size_t i2) const;


public:

  pedcnf_t()
		:_next_dummy(0), _no_of_clauses(0)
#ifndef AVOID_XOR_CLAUSES
		, _no_of_xor_clauses(0)
#endif
  {};

  ~pedcnf_t() {
  };

  lit_t generate_dummy();

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

  template <int LEN>
  void add_clause(const lit_t* const clause) {
	 add_clause(clause_t(clause, clause+LEN));
  };

#ifndef AVOID_XOR_CLAUSES
  template <int LEN>
  void add_xor_clause(const lit_t* const clause) {
	 add_xor_clause(xor_clause_t(clause, clause+LEN));
  };
#endif

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
  bool assignment_from_minisat_format(std::istream& in);


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




template <class KIND>
std::ostream&
operator<<(std::ostream& out, const pedcnf_t<KIND>::pedvar_t& var) {
  out << var.get<0>() << "_" << var.get<1>() << "_" << var.get<2>();
  return out;
}


template <class KIND>
lit_t
pedcnf_t::generate_dummy() {
  _vars.push_back(boost::make_tuple(KIND::DUMMY, _next_dummy, 0));
  _vals.push_back(false);
  ++_next_dummy;
  return _vars.size();
};


template <class KIND>
inline lit_t
pedcnf_t<KIND>::get_var(varmap_t& map,
								const KIND& var_kind,
								const size_t i1, const size_t i2) {
  varmap_t::iterator it= map.find(boost::make_tuple(i1, i2));
  if (it == map.end()) {
	 _vars.push_back(boost::make_tuple(var_kind, i1, i2));
	 _vals.push_back(false);
	 std::pair< varmap_t::iterator, bool> ret= map.insert(std::make_pair(boost::make_tuple(i1, i2), _vars.size()));
	 MY_ASSERT_DBG(ret.second);
	 it= ret.first;
  }
  return it->second;
};

inline lit_t
pedcnf_t::get_var(const varmap_t& map,
						const size_t i1, const size_t i2) const {
  varmap_t::const_iterator it= map.find(boost::make_tuple(i1, i2));
  if (it == map.end()) {
	 return -1;
  } else {
	 return it->second;
  }
};

template <class KIND>
inline bool
pedcnf_t<KIND>::has_var(const varmap_t& map,
								const size_t i1, const size_t i2) const {
  return get_var(map, i1, i2) != -1;
};

template <class KIND>
bool
pedcnf_t<KIND>::get_val(const varmap_t& map,
								const size_t i1, const size_t i2) const {
  lit_t var= get_var(map, i1, i2);
  MY_ASSERT_DBG( (0 < var) && ((size_t)var <= _vals.size()) );
  return _vals[var-1];
};


#ifndef ONLY_INTERNAL_SAT_SOLVER

template <class KIND>
bool
pedcnf_t<KIND>::is_satisfying_assignment() const {
  L_DEBUG("Checking if value assignment satisfies the or-clauses...");
  MY_ASSERT_DBG( _vars.size() == _vals.size() );
  bool ris= true;
  BOOST_FOREACH(const clause_t& clause, _clauses) {
	 if (ris) {
		bool intris= false;
		BOOST_FOREACH(const lit_t& var, clause) {
		  MY_ASSERT( var != 0 );
		  MY_ASSERT( (size_t)abs(var) <= _vars.size() );
		  if (var>0) {
			 intris= intris || _vals[var-1];
		  } else {
			 intris= intris || !_vals[-var-1];
		  }
		}
		if (!intris) {
		  L_DEBUG("Clause " << tostr(clause) << " is not satisfied.");
		}
		ris= ris && intris;
	 }
  }
  if (!ris) {
	 L_DEBUG("The assignment does not satisfy all the or-clauses.");
  } else {
	 L_DEBUG("The assignment satisfies all the or-clauses.");
#ifndef AVOID_XOR_CLAUSES
	 L_DEBUG("Checking if value assignment satisfies the xor-clauses...");
	 BOOST_FOREACH(const xor_clause_t& clause, _xor_clauses) {
		if (ris) {
		  bool intris= false;
		  BOOST_FOREACH(const lit_t& var, clause) {
			 MY_ASSERT( var != 0 );
			 MY_ASSERT( (size_t)abs(var) <= vars().size() );
			 bool val_var= vals()[abs(var)-1];
			 if (var<0) val_var= !val_var;
			 intris= (intris && !val_var) || (!intris && val_var);
		  }
		  if (!intris) {
			 L_DEBUG("Clause " << tostr(clause) << " is not satisfied.");
		  }
		  ris= ris && intris;
		}
	 }
	 if (ris) {
		L_DEBUG("The assignment satisfies all the or- and xor-clauses.");
	 } else {
		L_DEBUG("The assignment does not satisfy the xor-clauses.");
	 }
#endif
  }
  return ris;
};


template <class KIND>
std::ostream&
pedcnf_t<KIND>::clauses_to_dimacs_format(std::ostream& out) const {
  return clauses_to_dimacs_format(out, "SAT instance");
}

template <class KIND>
std::ostream&
pedcnf_t<KIND>::clauses_to_dimacs_format(std::ostream& out,
													  const std::string& note) const {
  return clauses_to_dimacs_format(out, std::vector< std::string >(1, note));
}

template <class KIND>
std::ostream&
pedcnf_t<KIND>::clauses_to_dimacs_format(std::ostream& out,
													  const std::vector< std::string >& notes) const {
  BOOST_FOREACH(const std::string& s, notes) {
	 out << "c " << s << std::endl;
  }
  if (_no_of_xor_clauses > 0)
	 out << "c extended syntax: or- and xor-clauses" << std::endl;
  out << "c" << std::endl;
  size_t i= 1;
  for (pedcnf_t<KIND>::varvec_t::const_iterator it= _vars.begin();
		 it != _vars.end();
		 ++it, ++i) {
	 out << "c v " << std::setw(7) << i << " " << *it << std::endl;
  }
  out << "c" << std::endl;
  out << "p cnf " << _vars.size() << " " << (_no_of_clauses + _no_of_xor_clauses) << std::endl;
  for (clauses_t::const_iterator it= _clauses.begin();
		 it != _clauses.end();
		 ++it) {
	 out << *it << std::endl;
  }
  for (xor_clauses_t::const_iterator it= _xor_clauses.begin();
		 it != _xor_clauses.end();
		 ++it) {
	 out << "x" << *it << std::endl;
  }
  return out;
};

#endif // ONLY_INTERNAL_SAT_SOLVER


// Read the assignment from a file like the following one:
// SAT/UNSAT
// 1 -2 3 4 0
template <class KIND>
bool
pedcnf_t<KIND>::assignment_from_minisat_format(std::istream& in) {
  std::string line;
  std::getline(in, line);
  boost::trim(line);
  MY_ASSERT_DBG( (line == "SAT") || (line == "UNSAT") );
  if (line == "UNSAT") {
	 L_DEBUG("UNSATISFIABLE clauses.");
	 return false;
  } else if (line == "SAT") {
	 L_DEBUG("Satisfiable clauses. Reading assignment...");
	 std::getline(in, line);
	 boost::trim(line);
	 std::istringstream is(line);
	 lit_t value;
	 while ((is >> value) && (value != 0)) {
		MY_ASSERT( (size_t)abs(value) <= _vals.size() );
		_vals[(size_t)abs(value)-1]= (value>0);
	 }
	 MY_ASSERT( ! (is >> value) ); // zero must be the last element
#ifndef ONLY_INTERNAL_SAT_SOLVER
	 MY_ASSERT_DBG( is_satisfying_assignment() );
#endif // ONLY_INTERNAL_SAT_SOLVER
	 return true;
  } else {
	 MY_FAIL;
  }
  return false;
};


template <class KIND>
std::ostream&
operator<<(std::ostream& out, const pedcnf_t<KIND>::clause_t& clause) {
  for (pedcnf_t::clause_t::const_iterator it= clause.begin();
		 it != clause.end();
		 ++it) {
	 out << std::setw(10) << *it << " ";
  }
  out << "     0";
  return out;
};

template <class KIND>
void
pedcnf_t<KIND>::add_clause(const clause_t& clause) {
#ifndef ONLY_INTERNAL_SAT_SOLVER
  _clauses.insert(clause);
#endif
#ifdef INTERNAL_SAT_SOLVER
  _solver.add_clause(clause);
#endif
  ++_no_of_clauses;
};

#ifndef AVOID_XOR_CLAUSES
template <class KIND>
void
pedcnf_t<KIND>::add_xor_clause(const xor_clause_t& clause) {
#ifndef ONLY_INTERNAL_SAT_SOLVER
  _xor_clauses.insert(clause);
#endif
#ifdef INTERNAL_SAT_SOLVER
  _solver.add_xor_clause(clause);
#endif
  ++_no_of_xor_clauses;
};
#endif








class pedcnf_reHCstar_t
  :
  public pedcnf_t<ped_var_kind_reHCstar>
{


// Data members
private:

  varmap_t _sp; // Grand-parental source (from father)
  varmap_t _sm; // Grand-parental source (from mother)
  varmap_t _p; // Paternal allele
  varmap_t _m; // Maternal allele
  varmap_t _rp; // Recombination events (from father)
  varmap_t _rm; // Recombination events (from mother)
  varmap_t _e; // Errors

// Methods
public:

  pedcnf_reHCstar_t() {};

  ~pedcnf_reHCstar_t() {
  };

  lit_t get_sp(const size_t i, const size_t l);

  lit_t get_sm(const size_t i, const size_t l);

  lit_t get_p(const size_t i, const size_t l);

  lit_t get_m(const size_t i, const size_t l);

  lit_t get_rp(const size_t i, const size_t l);

  lit_t get_rm(const size_t i, const size_t l);

  lit_t get_e(const size_t i, const size_t l);

  bool has_sp(const size_t i, const size_t l) const;

  bool has_sm(const size_t i, const size_t l) const;

  bool has_p(const size_t i, const size_t l) const;

  bool has_m(const size_t i, const size_t l) const;

  bool has_rp(const size_t i, const size_t l) const;

  bool has_rm(const size_t i, const size_t l) const;

  bool has_e(const size_t i, const size_t l) const;

  lit_t generate_dummy();

  lit_t get_sp(const size_t i, const size_t l) const;

  lit_t get_sm(const size_t i, const size_t l) const;

  lit_t get_p(const size_t i, const size_t l) const;

  lit_t get_m(const size_t i, const size_t l) const;

  lit_t get_rp(const size_t i, const size_t l) const;

  lit_t get_rm(const size_t i, const size_t l) const;

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

  bool e(const size_t i, const size_t l) const {
	 return get_val(_e, i, l);
  };

};

bool
pedcnf_reHCstar_t::has_sp(const size_t i, const size_t l) const {
  return has_var(_sp, i, l);
};

bool
pedcnf_reHCstar_t::has_sm(const size_t i, const size_t l) const {
  return has_var(_sm, i, l);
};

bool
pedcnf_reHCstar_t::has_p(const size_t i, const size_t l) const {
  return has_var(_p, i, l);
}

bool
pedcnf_reHCstar_t::has_m(const size_t i, const size_t l) const {
  return has_var(_m, i, l);
};

bool
pedcnf_reHCstar_t::has_rp(const size_t i, const size_t l) const {
  return has_var(_rp, i, l);
};

bool
pedcnf_reHCstar_t::has_rm(const size_t i, const size_t l) const {
  return has_var(_rm, i, l);
};

bool
pedcnf_reHCstar_t::has_e(const size_t i, const size_t l) const {
  return has_var(_e, i, l);
};


lit_t
pedcnf_reHCstar_t::get_sp(const size_t i, const size_t l) {
  return get_var(_sp, ped_var_kind_reHCstar::SP, i, l);
};

lit_t
pedcnf_reHCstar_t::get_sm(const size_t i, const size_t l) {
  return get_var(_sm, ped_var_kind_reHCstar::SM, i, l);
};

lit_t
pedcnf_reHCstar_t::get_p(const size_t i, const size_t l) {
  return get_var(_p, ped_var_kind_reHCstar::P, i, l);
};

lit_t
pedcnf_reHCstar_t::get_m(const size_t i, const size_t l) {
  return get_var(_m, ped_var_kind_reHCstar::M, i, l);
};

lit_t
pedcnf_reHCstar_t::get_rp(const size_t i, const size_t l) {
  return get_var(_rp, ped_var_kind_reHCstar::RP, i, l);
};

lit_t
pedcnf_reHCstar_t::get_rm(const size_t i, const size_t l) {
  return get_var(_rm, ped_var_kind_reHCstar::RM, i, l);
};

lit_t
pedcnf_reHCstar_t::get_e(const size_t i, const size_t l) {
  return get_var(_e, ped_var_kind_reHCstar::E, i, l);
};

lit_t
pedcnf_reHCstar_t::get_sp(const size_t i, const size_t l) const {
  return get_var(_sp, i, l);
};

lit_t
pedcnf_reHCstar_t::get_sm(const size_t i, const size_t l) const {
  return get_var(_sm, i, l);
};

lit_t
pedcnf_reHCstar_t::get_p(const size_t i, const size_t l) const {
  return get_var(_p, i, l);
};

lit_t
pedcnf_reHCstar_t::get_m(const size_t i, const size_t l) const {
  return get_var(_m, i, l);
};

lit_t
pedcnf_reHCstar_t::get_rp(const size_t i, const size_t l) const {
  return get_var(_rp, i, l);
};

lit_t
pedcnf_reHCstar_t::get_rm(const size_t i, const size_t l) const {
  return get_var(_rm, i, l);
};

lit_t
pedcnf_reHCstar_t::get_e(const size_t i, const size_t l) const {
  return get_var(_e, i, l);
};


template <class BASE>
class instance_t
{
private:
  BASE& b;
public:
  instance_t(BASE& _b)
		:b(_b)
  {};

  lit_t generate_dummy() {
	 return b.generate_dummy();
  };

  void add_clause(const clause_t& clause) {
	 b.add_clause(clause);
  };

#ifndef AVOID_XOR_CLAUSES
  void add_xor_clause(const xor_clause_t& clause) {
	 b.add_xor_clause(clause);
  };
#endif

  template <int LEN>
  void add_clause(const lit_t* const clause) {
	 b.add_clause<LEN>(clause);
  };

#ifndef AVOID_XOR_CLAUSES
  template <int LEN>
  void add_xor_clause(const lit_t* const clause) {
	 b.add_xor_clause<LEN>(clause);
  };
#endif

};


template <class PEDCNF>
void
add_card_constraint_less_or_equal_than(instance_t<PEDCNF>& cnf,
													const std::vector<var_t>& in_vars,
													const size_t k);

template <class PEDCNF>
void
add_uniform_card_constraint_less_or_equal_than(instance_t<PEDCNF>& cnf,
															  const std::vector<var_t>& in_vars,
															  const size_t window_size,
															  const size_t k);


#endif // __PEDCNF_HPP__
