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
 * pedcnf_card.cpp
 *
 * Functions implementing the encoding of cardinality constraints in SAT instances.
 *
 **/

#include "pedcnf.hpp"

#include <boost/foreach.hpp>
#include "utility.hpp"

inline static void
add_half_adder(pedcnf_t& cnf,
					const var_t x, const var_t y,
					var_t& s, var_t& cout) {
  s= cnf.generate_dummy();
  cout= cnf.generate_dummy();
  cnf.add_clause(pedcnf_t::clause_t{ x, -y,  s});
  cnf.add_clause(pedcnf_t::clause_t{-x,  y,  s});
  cnf.add_clause(pedcnf_t::clause_t{-x, -y,  cout});
};

inline static void
add_full_adder(pedcnf_t& cnf,
					const var_t x, const var_t y, const var_t cin,
					var_t& s, var_t& cout) {
  s= cnf.generate_dummy();
  cout= cnf.generate_dummy();
  cnf.add_clause(pedcnf_t::clause_t{ x,  y, -cin,  s});
  cnf.add_clause(pedcnf_t::clause_t{ x, -y,  cin,  s});
  cnf.add_clause(pedcnf_t::clause_t{-x,  y,  cin,  s});
  cnf.add_clause(pedcnf_t::clause_t{-x, -y, -cin,  s});
  cnf.add_clause(pedcnf_t::clause_t{-x, -y,   cout});
  cnf.add_clause(pedcnf_t::clause_t{-x, -cin, cout});
  cnf.add_clause(pedcnf_t::clause_t{-y, -cin, cout});
};


static void
generate_counter(pedcnf_t& cnf,
					  const std::vector<var_t>& x,
					  const size_t first, const size_t n,
					  std::vector<var_t>& d,
					  my_logger& logger) {
  const size_t len= n-first;
  d.clear();
  if (len == 1) {
	 d.push_back(x[first]);
  } else if (len == 2) {
	 var_t d1, d2;
	 add_half_adder(cnf, x[first], x[first+1], d1, d2);
	 d.push_back(d1);
	 d.push_back(d2);
  } else {
	 const size_t tprime= first + pow2_of_floor_log2(len) - 1;
	 std::vector<var_t> d1, d2;
	 generate_counter(cnf, x, first, tprime, d1, logger);
	 if (tprime < n-1) {
		generate_counter(cnf, x, tprime, n-1, d2, logger);
	 } else {
	 }
	 size_t i= 0;
	 var_t carry= x[n-1];
	 var_t s= 0;
	 while (i < d2.size()) {
		add_full_adder(cnf, d1[i], d2[i], carry, s, carry);
		d.push_back(s);
		++i;
	 }
	 while (i < d1.size()) {
		add_half_adder(cnf, d1[i], carry, s, carry);
		d.push_back(s);
		++i;
	 }
	 d.push_back(carry);
  }
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
  if (n == 1) {
	 const lit_t c1= cnf.generate_dummy();
	 const lit_t c2= cnf.generate_dummy();
	 vars.push_back(c1);
	 vars.push_back(c2);
	 cnf.add_clause(pedcnf_t::clause_t{ -A[0], -B[0], c2 });
	 cnf.add_clause(pedcnf_t::clause_t{ -A[0], c1 });
	 cnf.add_clause(pedcnf_t::clause_t{ -B[0], c1 });
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
		cnf.add_clause(pedcnf_t::clause_t{ -V_odd[i], -V_even[i-1], C[2*i-1] });
		cnf.add_clause(pedcnf_t::clause_t{ -V_odd[i], C[2*i-2] });
		cnf.add_clause(pedcnf_t::clause_t{ -V_even[i-1], C[2*i-2] });
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
  if (n == 1) {
	 const lit_t c1= cnf.generate_dummy();
	 const lit_t c2= cnf.generate_dummy();
	 vars.push_back(c1);
	 vars.push_back(c2);
	 cnf.add_clause(pedcnf_t::clause_t{ -A[0], -B[0], c2 });
	 cnf.add_clause(pedcnf_t::clause_t{ -A[0], c1 });
	 cnf.add_clause(pedcnf_t::clause_t{ -B[0], c1 });
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
		cnf.add_clause(pedcnf_t::clause_t{ -V_odd[i], -V_even[i-1], C[2*i-1] });
		cnf.add_clause(pedcnf_t::clause_t{ -V_odd[i], C[2*i-2] });
		cnf.add_clause(pedcnf_t::clause_t{ -V_even[i-1], C[2*i-2] });
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
	 cnf.add_clause(pedcnf_t::clause_t{ -dummy });
	 n= all_vars.size();
  }
  generate_hsort(cnf, all_vars, out_vars, logger);
  cnf.add_clause(pedcnf_t::clause_t{ -out_vars[k] });
};


static void
le_card_network_int(pedcnf_t& cnf,
						  const std::vector<var_t>& in_vars,
						  const size_t k1, const size_t k2) {
  my_logger logger(get_my_logger("card_constraints"));
  DEBUG("Generating a cardinality network for encoding the cardinality constraint.");
  const size_t p= pow2_of_ceiling_log2(k2+1);
  DEBUG("Dividing " << in_vars.size() << " variables in blocks of size "
		  << p << " since the upper bound is " << k2 << ".");
// Ensure that the variables are multiples of p
  std::vector<var_t> all_vars(in_vars.begin(), in_vars.end());
  while (all_vars.size() % p != 0) {
	 const lit_t dummy= cnf.generate_dummy();
	 all_vars.push_back(dummy);
	 cnf.add_clause(pedcnf_t::clause_t{ -dummy });
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
	 cnf.add_clause(pedcnf_t::clause_t{ -hsort_out[k2] });
	 // for (size_t limit= k; limit<p; ++limit) {
	 // 	cnf.add_clause(pedcnf_t::clause_t{ -hsort_out[limit] });
	 // }
	 std::vector<var_t> smerge_out;
	 generate_smerge(cnf, prev_ris, hsort_out, smerge_out, logger);
	 prev_ris.clear();
	 prev_ris.insert(prev_ris.end(), smerge_out.begin(), smerge_out.begin()+p);
	 next_in_var += p;
	 next_in_var_counter += p;
	 cnf.add_clause(pedcnf_t::clause_t{ -smerge_out[k2] });
	 // for (size_t limit= k; limit<=p; ++limit) {
	 // 	cnf.add_clause(pedcnf_t::clause_t{ -smerge_out[limit] });
	 // }
  }
  if (k1 > 0) {
	 DEBUG("Setting the lower bound >= " << k1);
	 for (size_t i= 0; i<k1; ++i)
		cnf.add_clause(pedcnf_t::clause_t{ prev_ris[i] });
  }
};

static void
le_card_network(pedcnf_t& cnf,
					 const std::vector<var_t>& in_vars,
					 const size_t k) {
  le_card_network_int(cnf, in_vars, 0, k);
}

static void
le_card_network_tree_rec(pedcnf_t& cnf,
								 std::vector<var_t>& out_vars,
								 const std::vector<var_t>& in_vars,
								 const size_t k, const size_t p,
								 my_logger& logger) {
  const size_t n= in_vars.size();
  std::vector<var_t> firsthalf_in(in_vars.begin(), in_vars.begin()+(n>>1));
  std::vector<var_t> secondhalf_in(in_vars.begin()+(n>>1), in_vars.end());
  std::vector<var_t> out_vars_int;
  if (n == 2*p) {
	 generate_smerge(cnf, firsthalf_in, secondhalf_in, out_vars_int, logger);
  } else {
	 std::vector<var_t> firsthalf_out;
	 std::vector<var_t> secondhalf_out;
	 le_card_network_tree_rec(cnf, firsthalf_out, firsthalf_in, k, p, logger);
	 le_card_network_tree_rec(cnf, secondhalf_out, secondhalf_in, k, p, logger);
	 generate_smerge(cnf, firsthalf_out, secondhalf_out, out_vars_int, logger);
  }
  for (size_t limit= k; limit<=p; ++limit) {
	 cnf.add_clause(pedcnf_t::clause_t{ -out_vars_int[limit] });
  }
  out_vars.insert(out_vars.end(), out_vars_int.begin(), out_vars_int.end()-1);
};

static void
le_card_network_tree(pedcnf_t& cnf,
							std::vector<var_t>& out_vars,
							const std::vector<var_t>& in_vars,
							const size_t k, const size_t p,
							my_logger& logger) {
  std::vector<var_t>::const_iterator next_in_var=
	 in_vars.begin();
  std::vector<var_t> all_hsort_out;
  std::vector<var_t> input;
  while (next_in_var != in_vars.end()) {
	 input.clear();
	 input.insert(input.end(), next_in_var, next_in_var+p);
	 std::vector<var_t> hsort_out;
	 generate_hsort(cnf, input, hsort_out, logger);
	 all_hsort_out.insert(all_hsort_out.end(), hsort_out.begin(), hsort_out.end());
	 next_in_var += p;
	 for (size_t limit= k; limit<p; ++limit) {
		cnf.add_clause(pedcnf_t::clause_t{ -hsort_out[limit] });
	 }
  }
  if (in_vars.size() > p) {
	 le_card_network_tree_rec(cnf, out_vars, all_hsort_out, k, p, logger);
  } else {
	 out_vars.insert(out_vars.end(), all_hsort_out.begin(), all_hsort_out.end());
  }
};


// Generate a cardinality network by using a complete tree of sorters
// and comparators
static void
le_card_network_new(pedcnf_t& cnf,
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
	 cnf.add_clause(pedcnf_t::clause_t{ -dummy });
  }
// The first 2^i * p <= n variables are arranged as a binary tree
  const size_t n= all_vars.size();
  std::vector<var_t> prev_ris;
  size_t next_in_var_counter= pow2_of_floor_log2(n/p)*p;
  DEBUG("The first " << next_in_var_counter << " variables are in a tree over " << all_vars.size() << ".");
  std::vector<var_t>::iterator next_in_var=
	 all_vars.begin() + next_in_var_counter;
  std::vector<var_t> input(all_vars.begin(), next_in_var);
  le_card_network_tree(cnf, prev_ris, input, k, p, logger);

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
		cnf.add_clause(pedcnf_t::clause_t{ -smerge_out[limit] });
	 }
  }
};


static void
eq_zero(pedcnf_t& cnf,
		  const std::vector<var_t>& in_vars) {
  BOOST_FOREACH(const lit_t& var, in_vars) {
	 cnf.add_clause(pedcnf_t::clause_t{ -var });
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

void
add_card_constraint_between(pedcnf_t& cnf,
									 const std::vector<var_t>& in_vars,
									 const size_t k1, const size_t k2) {
  my_logger logger(get_my_logger("card_constraints"));
  if (k2==0) {
	 eq_zero(cnf, in_vars);
  } else if (k2 >= in_vars.size()) {
	 INFO("The number of variables (" << in_vars.size() << ") "
			"is not greater than the numeric upper bound (" << k2 << "). "
			"Setting to " << in_vars.size() << ".");
	 le_card_network_int(cnf, in_vars,
								std::min(k1, in_vars.size()), in_vars.size());
  } else {
// Use Cardinality Networks
	 le_card_network_int(cnf, in_vars, k1, k2);
  }
};


static void
add_uniform_card_constraint_less_or_equal_than_int(pedcnf_t& cnf,
																	const std::vector<var_t>& in_vars,
																	const size_t window_size,
																	const size_t k) {
  my_logger logger(get_my_logger("card_constraints"));
  const size_t p= window_size>>1;
  std::vector<var_t>::const_iterator next_in_var=
	 in_vars.begin();
  std::vector<var_t> all_hsort_out;
  std::vector<var_t> input;
  while (next_in_var != in_vars.end()) {
	 input.clear();
	 input.insert(input.end(), next_in_var, next_in_var+p);
	 std::vector<var_t> hsort_out;
	 generate_hsort(cnf, input, hsort_out, logger);
	 all_hsort_out.insert(all_hsort_out.end(), hsort_out.begin(), hsort_out.end());
	 next_in_var += p;
  }

  std::vector<var_t>::iterator next_hsort_out=
	 all_hsort_out.begin()+p;
  std::vector<var_t> prev_input(all_hsort_out.begin(), next_hsort_out);
  while (next_hsort_out != all_hsort_out.end()) {
	 std::vector<var_t> cur_input(next_hsort_out, next_hsort_out+p);
	 std::vector<var_t> smerge_out;
	 generate_smerge(cnf, prev_input, cur_input, smerge_out, logger);
	 next_hsort_out += p;
	 prev_input= cur_input;
	 cnf.add_clause(pedcnf_t::clause_t{ -smerge_out[k] });
  }

};

void
add_uniform_card_constraint_less_or_equal_than(pedcnf_t& cnf,
															  const std::vector<var_t>& in_vars,
															  const size_t window_size,
															  const size_t k) {
  my_logger logger(get_my_logger("card_constraints"));
  if (k==0) {
	 eq_zero(cnf, in_vars);
  } else if (k >= in_vars.size()) {
	 INFO("The number of variables (" << in_vars.size() << ") "
			"is not greater than the numeric upper bound (" << k << "). "
			"No clauses have been added.");
  } else {
	 MY_ASSERT( window_size == pow2_of_floor_log2(window_size) );
	 MY_ASSERT( k < window_size );
	 std::vector<var_t> all_vars= in_vars;
	 while (all_vars.size() % window_size != 0) {
		var_t v= cnf.generate_dummy();
		all_vars.push_back(v);
		cnf.add_clause(pedcnf_t::clause_t{ -v });
	 }
	 add_uniform_card_constraint_less_or_equal_than_int(cnf, all_vars,
																		 window_size, k);
  }

}
