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
 * pedcnf.cpp
 *
 * Structures to represent SAT instances derived from pedigrees.
 *
 **/

#include "pedcnf.hpp"
#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/trim.hpp>

const ped_var_kind ped_var_kind::S(0);
const ped_var_kind ped_var_kind::P(1);
const ped_var_kind ped_var_kind::M(2);
const ped_var_kind ped_var_kind::E(3);
const ped_var_kind ped_var_kind::DUMMY(4);
const int ped_var_kind::int_values[]={0, 1, 2, 3, 4};
const std::string ped_var_kind::str_values[]={"s", "p", "m", "e", "dummy"};
const ped_var_kind ped_var_kind::enum_values[]={S, P, M, E, DUMMY};

std::ostream&
operator<<(std::ostream& out, const pedcnf_t::pedvar_t& var) {
  out << var.get<0>() << "_" << var.get<1>() << "_" << var.get<2>();
  return out;
}


lit_t
pedcnf_t::generate_dummy() {
  _vars.push_back(boost::make_tuple(ped_var_kind::DUMMY, _next_dummy, 0));
  _vals.push_back(false);
  ++_next_dummy;
  return _vars.size();
};

inline lit_t
pedcnf_t::get_var(varmap_t& map,
						const ped_var_kind& var_kind,
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

bool
pedcnf_t::get_val(const varmap_t& map,
						const size_t i1, const size_t i2) const {
  lit_t var= get_var(map, i1, i2);
  MY_ASSERT_DBG( (0 < var) && ((size_t)var <= _vals.size()) );
  return _vals[var-1];
};

lit_t
pedcnf_t::get_s(const size_t p, const size_t i) {
  return get_var(_s, ped_var_kind::S, p, i);
};

lit_t
pedcnf_t::get_p(const size_t i, const size_t l) {
  return get_var(_p, ped_var_kind::P, i, l);
};

lit_t
pedcnf_t::get_m(const size_t i, const size_t l) {
  return get_var(_m, ped_var_kind::M, i, l);
};

lit_t
pedcnf_t::get_e(const size_t i, const size_t l) {
  return get_var(_e, ped_var_kind::E, i, l);
};

lit_t
pedcnf_t::get_s(const size_t p, const size_t i) const {
  return get_var(_s, p, i);
};

lit_t
pedcnf_t::get_p(const size_t i, const size_t l) const {
  return get_var(_p, i, l);
};

lit_t
pedcnf_t::get_m(const size_t i, const size_t l) const {
  return get_var(_m, i, l);
};

lit_t
pedcnf_t::get_e(const size_t i, const size_t l) const {
  return get_var(_e, i, l);
};


#ifndef ONLY_INTERNAL_SAT_SOLVER

bool
pedcnf_t::is_satisfying_assignment() const {
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
  }
  return ris;
};


std::ostream&
pedcnf_t::clauses_to_dimacs_format(std::ostream& out) const {
  return clauses_to_dimacs_format(out, "SAT instance");
}

std::ostream&
pedcnf_t::clauses_to_dimacs_format(std::ostream& out,
											  const std::string& note) const {
  return clauses_to_dimacs_format(out, std::vector< std::string >(1, note));
}

std::ostream&
pedcnf_t::clauses_to_dimacs_format(std::ostream& out,
											  const std::vector< std::string >& notes) const {
  BOOST_FOREACH(const std::string& s, notes) {
	 out << "c " << s << std::endl;
  }
  if (_no_of_xor_clauses > 0)
	 out << "c extended syntax: or- and xor-clauses" << std::endl;
  out << "c" << std::endl;
  size_t i= 1;
  for (pedcnf_t::varvec_t::const_iterator it= _vars.begin();
		 it != _vars.end();
		 ++it, ++i) {
	 out << "c v " << std::setw(5) << i << " " << *it << std::endl;
  }
  out << "c" << std::endl;
  out << "p cnf " << _vars.size() << " " << _clauses.size() << std::endl;
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
bool
pedcnf_t::assignment_from_minisat_format(std::istream& in) {
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


std::ostream&
operator<<(std::ostream& out, const pedcnf_t::clause_t& clause) {
  for (pedcnf_t::clause_t::const_iterator it= clause.begin();
		 it != clause.end();
		 ++it) {
	 out << std::setw(6) << *it << " ";
  }
  out << "     0";
  return out;
};

void
pedcnf_t::add_clause(const clause_t& clause) {
#ifndef ONLY_INTERNAL_SAT_SOLVER
  _clauses.insert(clause);
#endif
#ifdef INTERNAL_SAT_SOLVER
  _solver.add_clause(clause);
#endif
  ++_no_of_clauses;
};

void
pedcnf_t::add_xor_clause(const xor_clause_t& clause) {
#ifndef ONLY_INTERNAL_SAT_SOLVER
  _xor_clauses.insert(clause);
#endif
#ifdef INTERNAL_SAT_SOLVER
  _solver.add_xor_clause(clause);
#endif
  ++_no_of_xor_clauses;
};


