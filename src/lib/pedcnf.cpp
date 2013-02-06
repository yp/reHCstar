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
 * pedcnf.cpp
 *
 * Structures to represent SAT instances derived from pedigrees.
 *
 **/

#include "pedcnf.hpp"
#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/trim.hpp>

const ped_var_kind ped_var_kind::SP(0);
const ped_var_kind ped_var_kind::SM(1);
const ped_var_kind ped_var_kind::P(2);
const ped_var_kind ped_var_kind::M(3);
const ped_var_kind ped_var_kind::RP(4);
const ped_var_kind ped_var_kind::RM(5);
const ped_var_kind ped_var_kind::E(6);
const ped_var_kind ped_var_kind::PM(7);
const ped_var_kind ped_var_kind::MM(8);
const ped_var_kind ped_var_kind::DUMMY(9);
const int ped_var_kind::int_values[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
const std::string ped_var_kind::str_values[]={"sp", "sm", "p", "m", "rp", "rm", "e", "pm", "mm", "dummy"};
const ped_var_kind ped_var_kind::enum_values[]={SP, SM, P, M, RP, RM, E, PM, MM, DUMMY};

std::ostream&
operator<<(std::ostream& out, const pedcnf_t::pedvar_t& var) {
  out << var.get<0>() << "_" << var.get<1>() << "_" << var.get<2>() << "_" << var.get<3>();
  return out;
}


lit_t
pedcnf_t::generate_dummy() {
  _vars.push_back(boost::make_tuple(ped_var_kind::DUMMY, _next_dummy, 0, 0));
  _vals.push_back(false);
  ++_next_dummy;
  return _vars.size();
};

inline lit_t
pedcnf_t::get_var3(varmap_t& map,
						 const ped_var_kind& var_kind,
						 const size_t i1, const size_t i2, const size_t i3) {
  varmap_t::iterator it= map.find(boost::make_tuple(i1, i2, i3));
  if (it == map.end()) {
	 _vars.push_back(boost::make_tuple(var_kind, i1, i2, i3));
	 _vals.push_back(false);
	 std::pair< varmap_t::iterator, bool> ret= map.insert(std::make_pair(boost::make_tuple(i1, i2, i3), _vars.size()));
	 MY_ASSERT_DBG(ret.second);
	 it= ret.first;
  }
  return it->second;
};

inline lit_t
pedcnf_t::get_var3(const varmap_t& map,
						 const size_t i1, const size_t i2, const size_t i3) const {
  varmap_t::const_iterator it= map.find(boost::make_tuple(i1, i2, i3));
  if (it == map.end()) {
	 return -1;
  } else {
	 return it->second;
  }
};

inline bool
pedcnf_t::has_var3(const varmap_t& map,
						 const size_t i1, const size_t i2, const size_t i3) const {
  return get_var3(map, i1, i2, i3) != -1;
};

bool
pedcnf_t::has_sp(const size_t i, const size_t l) const {
  return has_var(_sp, i, l);
};

bool
pedcnf_t::has_sm(const size_t i, const size_t l) const {
  return has_var(_sm, i, l);
};

bool
pedcnf_t::has_p(const size_t i, const size_t l) const {
  return has_var(_p, i, l);
}

bool
pedcnf_t::has_m(const size_t i, const size_t l) const {
  return has_var(_m, i, l);
};

bool
pedcnf_t::has_pm(const size_t i, const size_t l, const size_t j) const {
  return has_var3(_pm, i, l, j);
}

bool
pedcnf_t::has_mm(const size_t i, const size_t l, const size_t j) const {
  return has_var3(_mm, i, l, j);
};

bool
pedcnf_t::has_rp(const size_t i, const size_t l) const {
  return has_var(_rp, i, l);
};

bool
pedcnf_t::has_rm(const size_t i, const size_t l) const {
  return has_var(_rm, i, l);
};

bool
pedcnf_t::has_e(const size_t i, const size_t l) const {
  return has_var(_e, i, l);
};


bool
pedcnf_t::get_val3(const varmap_t& map,
						 const size_t i1, const size_t i2, const size_t i3) const {
  lit_t var= get_var3(map, i1, i2, i3);
  MY_ASSERT_DBG( (0 < var) && ((size_t)var <= _vals.size()) );
  return _vals[var-1];
};

lit_t
pedcnf_t::get_sp(const size_t i, const size_t l) {
  return get_var(_sp, ped_var_kind::SP, i, l);
};

lit_t
pedcnf_t::get_sm(const size_t i, const size_t l) {
  return get_var(_sm, ped_var_kind::SM, i, l);
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
pedcnf_t::get_pm(const size_t i, const size_t l, const size_t j) {
  return get_var3(_pm, ped_var_kind::PM, i, l, j);
};

lit_t
pedcnf_t::get_mm(const size_t i, const size_t l, const size_t j) {
  return get_var3(_mm, ped_var_kind::MM, i, l, j);
};

lit_t
pedcnf_t::get_rp(const size_t i, const size_t l) {
  return get_var(_rp, ped_var_kind::RP, i, l);
};

lit_t
pedcnf_t::get_rm(const size_t i, const size_t l) {
  return get_var(_rm, ped_var_kind::RM, i, l);
};

lit_t
pedcnf_t::get_e(const size_t i, const size_t l) {
  return get_var(_e, ped_var_kind::E, i, l);
};

lit_t
pedcnf_t::get_sp(const size_t i, const size_t l) const {
  return get_var(_sp, i, l);
};

lit_t
pedcnf_t::get_sm(const size_t i, const size_t l) const {
  return get_var(_sm, i, l);
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
pedcnf_t::get_pm(const size_t i, const size_t l, const size_t j) const {
  return get_var3(_pm, i, l, j);
};

lit_t
pedcnf_t::get_mm(const size_t i, const size_t l, const size_t j) const {
  return get_var3(_mm, i, l, j);
};

lit_t
pedcnf_t::get_rp(const size_t i, const size_t l) const {
  return get_var(_rp, i, l);
};

lit_t
pedcnf_t::get_rm(const size_t i, const size_t l) const {
  return get_var(_rm, i, l);
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
#ifndef AVOID_XOR_CLAUSES
  if (_no_of_xor_clauses > 0)
	 out << "c extended syntax: or- and xor-clauses" << std::endl;
#endif // AVOID_XOR_CLAUSES
  // out << "c" << std::endl;
  // size_t i= 1;
  // for (pedcnf_t::varvec_t::const_iterator it= _vars.begin();
  // 		 it != _vars.end();
  // 		 ++it, ++i) {
  // 	 out << "c v " << std::setw(7) << i << " " << *it << std::endl;
  // }
  // out << "c" << std::endl;
  size_t n_clauses= _no_of_clauses;
#ifndef AVOID_XOR_CLAUSES
  n_clauses += _no_of_xor_clauses;
#endif // AVOID_XOR_CLAUSES
  out << "p cnf " << _vars.size() << " " << n_clauses << std::endl;
  for (clauses_t::const_iterator it= _clauses.begin();
		 it != _clauses.end();
		 ++it) {
	 out << *it << '\n';
  }
#ifndef AVOID_XOR_CLAUSES
  for (xor_clauses_t::const_iterator it= _xor_clauses.begin();
		 it != _xor_clauses.end();
		 ++it) {
	 out << "x" << *it << '\n';
  }
#endif // AVOID_XOR_CLAUSES
  out.flush();
  return out;
};

#endif // ONLY_INTERNAL_SAT_SOLVER

#ifdef USE_PLAIN_CNF_FORMAT
// *********************************************************
//      !!!  OLD/PLAIN FORMAT  !!!
// Read the assignment from a file like the following one:
// SAT/UNSAT
// 1 -2 3 4 0
// *********************************************************
boost::tribool
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

#else // USE_PLAIN_CNF_FORMAT

// *********************************************************
//      !!!  DIMACS FORMAT  !!!
// Read the assignment from a file like the following one:
// s SATISFIABLE/UNSATISFIABLE
// v 1 -2 3 4
// v 5 -6 -7 0
// *********************************************************
boost::tribool
pedcnf_t::assignment_from_minisat_format(std::istream& in) {
  int status= 0; // 0 = undecided, 1 = SAT, -1 = UNSAT
  std::string line;
  while (std::getline(in, line)) {
	 boost::trim(line);
	 if ((line.length()==0) || boost::starts_with(line, "c")) {
// The line is a comment, discard.
	 } else if (boost::starts_with(line, "s ")) {
		if (line == "s SATISFIABLE") {
		  L_DEBUG("The instance is SATISFIABLE");
		  status= 1;
		} else if (line == "s UNSATISFIABLE") {
		  L_DEBUG("The instance is UNSATISFIABLE");
		  status= -1;
		} else {
		  L_FATAL("Read a wrong status line: >" << line << "<");
		  MY_FAIL;
		}
	 } else if (boost::starts_with(line, "v ")) {
		std::istringstream is(line);
// Remove the initial "v"
		char c;
		is >> c;
		lit_t value;
		while ((is >> value) && (value != 0)) {
		  MY_ASSERT( (size_t)abs(value) <= _vals.size() );
		  _vals[(size_t)abs(value)-1]= (value>0);
		}
	 } else {
		L_FATAL("Read an unparsable line: >" << line << "<");
		MY_FAIL;
	 }
  }
  if (status == 1) { // SAT
#ifndef ONLY_INTERNAL_SAT_SOLVER
	 MY_ASSERT_DBG( is_satisfying_assignment() );
#endif // ONLY_INTERNAL_SAT_SOLVER
	 return true;
  } else if (status == -1) { // UNSAT
	 return false;
  } else { // undecided
	 return boost::indeterminate;
  }
};

#endif // USE_PLAIN_CNF_FORMAT

std::ostream&
operator<<(std::ostream& out, const pedcnf_t::clause_t& clause) {
  for (pedcnf_t::clause_t::const_iterator it= clause.begin();
		 it != clause.end();
		 ++it) {
//	 out << std::setw(10) << *it << " ";
	 out << *it << " ";
  }
  out << "0";
  return out;
};

void
pedcnf_t::add_clause(const clause_t& clause) {
#ifndef ONLY_INTERNAL_SAT_SOLVER
  _clauses.push_back(clause);
#endif
#ifdef INTERNAL_SAT_SOLVER
  _solver.add_clause(clause);
#endif
  ++_no_of_clauses;
};

#ifndef AVOID_XOR_CLAUSES
void
pedcnf_t::add_xor_clause(const xor_clause_t& clause) {
#ifndef ONLY_INTERNAL_SAT_SOLVER
  _xor_clauses.push_back(clause);
#endif
#ifdef INTERNAL_SAT_SOLVER
  _solver.add_xor_clause(clause);
#endif
  ++_no_of_xor_clauses;
};
#endif

