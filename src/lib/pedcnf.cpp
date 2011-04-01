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


inline static void
add_half_adder(pedcnf_t& cnf,
					const var_t x, const var_t y,
					var_t& s, var_t& cout) {
  s= cnf.generate_dummy();
  cout= cnf.generate_dummy();
  cnf.add_clause<3>((lit_t[]){ x, -y,  s});
  cnf.add_clause<3>((lit_t[]){-x,  y,  s});
  cnf.add_clause<3>((lit_t[]){-x, -y,  cout});
};

inline static void
add_full_adder(pedcnf_t& cnf,
					const var_t x, const var_t y, const var_t cin,
					var_t& s, var_t& cout) {
  s= cnf.generate_dummy();
  cout= cnf.generate_dummy();
  cnf.add_clause<4>((lit_t[]){ x,  y, -cin,  s});
  cnf.add_clause<4>((lit_t[]){ x, -y,  cin,  s});
  cnf.add_clause<4>((lit_t[]){-x,  y,  cin,  s});
  cnf.add_clause<4>((lit_t[]){-x, -y, -cin,  s});
  cnf.add_clause<3>((lit_t[]){-x, -y,   cout});
  cnf.add_clause<3>((lit_t[]){-x, -cin, cout});
  cnf.add_clause<3>((lit_t[]){-y, -cin, cout});
};

inline static size_t
pow2_of_floor_log2(const size_t n) {
  size_t p= 1;
  while ((p << 1) <= n) {
	 p= p << 1;
  }
  return p;
};

inline static size_t
pow2_of_ceiling_log2(const size_t n) {
  size_t p= 1;
  while (p < n) {
	 p= p << 1;
  }
  return p;
};


static void
generate_counter(pedcnf_t& cnf,
					  const std::vector<var_t>& x,
					  const size_t first, const size_t n,
					  std::vector<var_t>& d,
					  my_logger& logger) {
  const size_t len= n-first;
  d.clear();
  TRACE(" --> recursively generating counter for variables [" << first <<
			 ", " << n-1 << "] ...");
  if (len == 1) {
	 d.push_back(x[first]);
  } else if (len == 2) {
	 size_t d1, d2;
	 add_half_adder(cnf, x[first], x[first+1], d1, d2);
	 d.push_back(d1);
	 d.push_back(d2);
  } else {
	 const size_t tprime= first + pow2_of_floor_log2(len) - 1;
	 std::vector<var_t> d1, d2;
	 generate_counter(cnf, x, first, tprime, d1, logger);
	 if (tprime < n-1) {
		generate_counter(cnf, x, tprime, n-1, d2, logger);
		TRACE("     counter [" << first << ", " << n-1 << "] = "
				"[" << first << ", " << tprime-1 << "] + "
				"[" << tprime << ", " << n-2 << "] + " << n-1);
	 } else {
		TRACE("     counter [" << first << ", " << n-1 << "] = "
				"[" << first << ", " << tprime-1 << "] + " << n-1);
	 }
	 size_t i= 0;
	 size_t carry= x[n-1];
	 size_t s= 0;
	 while (i < d2.size()) {
		add_full_adder(cnf, d1[i], d2[i], carry, s, carry);
		d.push_back(s);
		TRACE("        out[" << i << "]= " << s);
		++i;
	 }
	 while (i < d1.size()) {
		add_half_adder(cnf, d1[i], carry, s, carry);
		d.push_back(s);
		TRACE("        out[" << i << "]= " << s);
		++i;
	 }
	 d.push_back(carry);
	 TRACE("        out[" << i << "]= " << carry);
  }
  TRACE("      #outputs = " << d.size());
};

static void
generate_comparator(pedcnf_t& cnf,
						  const std::vector<var_t>& d,
						  const size_t k,
						  my_logger& logger) {
  TRACE("Generate comparator with " << d.size() << " variables (<="
		  << k << ")...");
  std::vector<pedcnf_t::clause_t> clauses;
  size_t quotient= k;
  size_t i= 0;
  while (i < d.size()) {
	 const size_t remainder= quotient % 2;
	 TRACE("  k[" << i << "]= " << remainder);
	 if (remainder == 0) {
		pedcnf_t::clause_t c;
		c.insert(-d[i]);
		clauses.push_back(c);
	 } else { // remainder == 1
		BOOST_FOREACH( pedcnf_t::clause_t& c, clauses ) {
		  c.insert(-d[i]);
		}
	 }
	 quotient= quotient >> 1;
	 ++i;
  }
  BOOST_FOREACH( pedcnf_t::clause_t& c, clauses ) {
	 cnf.add_clause(pedcnf_t::clause_t(c.begin(), c.end()));
  }
}

static inline void
divide_vector(const std::vector<var_t>& V,
				  std::vector<var_t>& odd,
				  std::vector<var_t>& even) {
  for(std::vector<var_t>::const_iterator it= V.begin();
		it != V.end();) {
	 odd.push_back(*it);
	 ++it;
	 even.push_back(*it);
	 ++it;
  }
};


static void
generate_hmerge(pedcnf_t& cnf,
					 const std::vector<var_t>& A,
					 const std::vector<var_t>& B,
					 std::vector<var_t>& vars,
					 my_logger& logger) {
  MY_ASSERT_DBG(A.size() == B.size());
  const size_t n= A.size();
  MY_ASSERT_DBG(n > 0);
  MY_ASSERT_DBG(n == pow2_of_floor_log2(n));
  vars.clear();
  TRACE(" --> HMerge on two sequences of " << n << " variables.");
  if (n == 1) {
	 const lit_t c1= cnf.generate_dummy();
	 const lit_t c2= cnf.generate_dummy();
	 vars.push_back(c1);
	 vars.push_back(c2);
	 cnf.add_clause<3>((lit_t[]){ -A[0], -B[0], c2 });
	 cnf.add_clause<2>((lit_t[]){ -A[0], c1 });
	 cnf.add_clause<2>((lit_t[]){ -B[0], c1 });
  } else {
	 std::vector<var_t> A_odd, A_even;
	 divide_vector(A, A_odd, A_even);
	 std::vector<var_t> B_odd, B_even;
	 divide_vector(B, B_odd, B_even);
	 std::vector<var_t> V_odd, V_even;
	 generate_hmerge(cnf, A_odd,  B_odd,  V_odd,  logger);
	 generate_hmerge(cnf, A_even, B_even, V_even, logger);
	 std::vector<var_t> C;
	 for (size_t i= 2; i<2*n; ++i) {
		C.push_back(cnf.generate_dummy());
	 }
	 for (size_t i= 1; i<n; ++i) {
		cnf.add_clause<3>((lit_t[]){ -V_odd[i], -V_even[i-1], C[2*i-1] });
		cnf.add_clause<2>((lit_t[]){ -V_odd[i], C[2*i-2] });
		cnf.add_clause<2>((lit_t[]){ -V_even[i-1], C[2*i-2] });
	 }
	 vars.push_back(V_odd[0]);
	 vars.insert(vars.end(), C.begin(), C.end());
	 vars.push_back(V_even.back());
  }

}

static void
generate_hsort(pedcnf_t& cnf,
					const std::vector<var_t>& S,
					std::vector<var_t>& vars,
					my_logger& logger) {
  const size_t n= S.size();
  MY_ASSERT_DBG(n > 1);
  MY_ASSERT_DBG(n == pow2_of_floor_log2(n));
  vars.clear();
  TRACE(" --> HSort a sequence of " << n << " variables.");
  if (n == 2) {
	 std::vector<var_t> A, B;
	 A.push_back(S[0]);
	 B.push_back(S[1]);
	 generate_hmerge(cnf, A, B, vars, logger);
  } else {
	 const size_t halfn= n>>1;
	 std::vector<var_t> S1, S2;
	 S1.insert(S1.end(), S.begin(), S.begin()+halfn);
	 S2.insert(S2.end(), S.begin()+halfn, S.end());
	 std::vector<var_t> Var1, Var2;
	 generate_hsort(cnf, S1, Var1, logger);
	 generate_hsort(cnf, S2, Var2, logger);
	 generate_hmerge(cnf, Var1, Var2, vars, logger);
  }
}

static void
generate_smerge(pedcnf_t& cnf,
					 const std::vector<var_t>& A,
					 const std::vector<var_t>& B,
					 std::vector<var_t>& vars,
					 my_logger& logger) {
  MY_ASSERT_DBG(A.size() == B.size());
  const size_t n= A.size();
  MY_ASSERT_DBG(n > 0);
  MY_ASSERT_DBG(n == pow2_of_floor_log2(n));
  vars.clear();
  TRACE(" --> SMerge on two sequences of " << n << " variables.");
  if (n == 1) {
	 const lit_t c1= cnf.generate_dummy();
	 const lit_t c2= cnf.generate_dummy();
	 vars.push_back(c1);
	 vars.push_back(c2);
	 cnf.add_clause<3>((lit_t[]){ -A[0], -B[0], c2 });
	 cnf.add_clause<2>((lit_t[]){ -A[0], c1 });
	 cnf.add_clause<2>((lit_t[]){ -B[0], c1 });
  } else {
	 std::vector<var_t> A_odd, A_even;
	 divide_vector(A, A_odd, A_even);
	 std::vector<var_t> B_odd, B_even;
	 divide_vector(B, B_odd, B_even);
	 std::vector<var_t> V_odd, V_even;
	 generate_smerge(cnf, A_odd,  B_odd,  V_odd,  logger);
	 generate_smerge(cnf, A_even, B_even, V_even, logger);
	 std::vector<var_t> C;
	 for (size_t i= 2; i<n+2; ++i) {
		C.push_back(cnf.generate_dummy());
	 }
	 for (size_t i= 1; i<=(n>>1); ++i) {
		cnf.add_clause<3>((lit_t[]){ -V_odd[i], -V_even[i-1], C[2*i-1] });
		cnf.add_clause<2>((lit_t[]){ -V_odd[i], C[2*i-2] });
		cnf.add_clause<2>((lit_t[]){ -V_even[i-1], C[2*i-2] });
	 }
	 vars.push_back(V_odd[0]);
	 vars.insert(vars.end(), C.begin(), C.end());
  }

}



static void
le_parallel_counter(pedcnf_t& cnf,
						  const std::vector<var_t>& in_vars,
						  const size_t k) {
  my_logger logger(get_my_logger("card_constraints"));
  DEBUG("Generating a parallel counter for encoding the cardinality constraint.");
  std::vector<var_t> out_vars;
  generate_counter(cnf, in_vars, 0, in_vars.size(), out_vars, logger);
  generate_comparator(cnf, out_vars, k, logger);
};


static void
le_half_sorting_network(pedcnf_t& cnf,
								const std::vector<var_t>& in_vars,
								const size_t k) {
  my_logger logger(get_my_logger("card_constraints"));
  DEBUG("Generating an half sorting network for encoding the cardinality constraint.");
  std::vector<var_t> out_vars;
  std::vector<var_t> all_vars(in_vars.begin(), in_vars.end());
  size_t n= all_vars.size();
  while (n != pow2_of_floor_log2(n)) {
	 const lit_t dummy= cnf.generate_dummy();
	 all_vars.push_back(dummy);
	 cnf.add_clause<1>((lit_t[]){ -dummy });
	 n= all_vars.size();
  }
  generate_hsort(cnf, all_vars, out_vars, logger);
  cnf.add_clause<1>((lit_t[]){ -out_vars[k] });
};


static void
le_card_network(pedcnf_t& cnf,
					 const std::vector<var_t>& in_vars,
					 const size_t k) {
  my_logger logger(get_my_logger("card_constraints"));
  DEBUG("Generating a cardinality network for encoding the cardinality constraint.");
  const size_t p= pow2_of_ceiling_log2(k+1);
  DEBUG("Dividing " << in_vars.size() << " variables in blocks of size "
		  << p << " since the upper bound is " << k << ".");
// Ensure that the variables are multiples of p
  std::vector<var_t> all_vars(in_vars.begin(), in_vars.end());
  while (all_vars.size() % p != 0) {
	 const lit_t dummy= cnf.generate_dummy();
	 all_vars.push_back(dummy);
	 cnf.add_clause<1>((lit_t[]){ -dummy });
  }
  std::vector<var_t> prev_ris;
  size_t next_in_var_counter= p;
  std::vector<var_t>::iterator next_in_var=
	 all_vars.begin() + next_in_var_counter;
  std::vector<var_t> input(all_vars.begin(), next_in_var);
  generate_hsort(cnf, input, prev_ris, logger);
  while (next_in_var != all_vars.end()) {
	 input.clear();
	 input.insert(input.end(), next_in_var, next_in_var+p);
	 std::vector<var_t> hsort_out;
	 generate_hsort(cnf, input, hsort_out, logger);
	 std::vector<var_t> smerge_out;
	 generate_smerge(cnf, prev_ris, hsort_out, smerge_out, logger);
	 prev_ris.clear();
	 prev_ris.insert(prev_ris.end(), smerge_out.begin(), smerge_out.begin()+p);
	 next_in_var += p;
	 next_in_var_counter += p;
	 for (size_t limit= k; limit<=p; ++limit) {
		cnf.add_clause<1>((lit_t[]){ -smerge_out[limit] });
	 }
  }
};


static void
eq_zero(pedcnf_t& cnf,
		  const std::vector<var_t>& in_vars) {
  BOOST_FOREACH(const lit_t& var, in_vars) {
	 cnf.add_clause<1>((lit_t[]){ -var });
  }
};

void
add_card_constraint_less_or_equal_than(pedcnf_t& cnf,
													const std::vector<var_t>& in_vars,
													const size_t k) {
  my_logger logger(get_my_logger("card_constraints"));
  if (k==0) {
	 eq_zero(cnf, in_vars);
  } else if (k >= in_vars.size()) {
	 INFO("The number of variables (" << in_vars.size() << ") "
			"is not greater than the numeric upper bound (" << k << "). "
			"No clauses have been added.");
  } else {
// Use Cardinality Networks (preferred)
	 le_card_network(cnf, in_vars, k);

// Use Parallel Counter
//	 le_parallel_counter(cnf, in_vars, k);

//	Use Half Sorting Network
//	 le_half_sorting_network(cnf, in_vars, k);
  }
};
