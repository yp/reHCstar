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
 * ped2cnf.hpp
 *
 * Functions to convert a pedigree into a SAT instance
 *
 **/

#ifndef __PED2CNF_HPP__
#define __PED2CNF_HPP__


#include "data.hpp"
#include "log.hpp"
#include "utility.hpp"
#include "pedcnf.hpp"

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include <boost/foreach.hpp>




template <>
inline const char* logger_name<pedcnf_t>(void) {
  return "pedcnf_t";
}

template <
  typename T_GENOTYPE,
  typename T_HAPLOTYPE,
  typename T_PHENOTYPE,
  typename T_ID
  >
class ped2cnf_conv_t:
  public log_able_t< pedcnf_t >
{
// Types
private:

  typedef basic_pedigree_t<T_GENOTYPE,
									T_HAPLOTYPE,
									T_PHENOTYPE,
									T_ID> pedigree_t;
  typedef typename pedigree_t::individual_t individual_t;
  typedef typename pedigree_t::g g;



public:

// Data members
private:

public:

// Methods
private:

  void prepare_individual(pedcnf_t& cnf, const pedigree_t& ped,
								  const individual_t& ind) {
// Create the variables for the individual
//   -> p-variables
	 for (size_t l= 0; l < ped.genotype_length(); ++l) {
		cnf.get_p(ind.progr_id(), l);
	 }
//   -> m-variables
	 for (size_t l= 0; l < ped.genotype_length(); ++l) {
		cnf.get_m(ind.progr_id(), l);
	 }
//   -> r-variables
// (start at 1, no recombinations in the first locus)
	 if (ind.has_father()) {
		for (size_t l= 1; l < ped.genotype_length(); ++l) {
		  cnf.get_rp(ind.progr_id(), l);
		}
	 }
	 if (ind.has_mother()) {
		for (size_t l= 1; l < ped.genotype_length(); ++l) {
		  cnf.get_rm(ind.progr_id(), l);
		}
	 }
//   -> e-variables
	 for (size_t l= 0; l < ped.genotype_length(); ++l) {
//    errors only for genotyped loci
		if (is_genotyped(ind.obs_g(l))) {
		  cnf.get_e(ind.progr_id(), l);
		}
	 }
  }

  void add_individual_constraint_bi(pedcnf_t& cnf,
												const individual_t& ind,
												const g& gen,
												const size_t l,
												const size_t i) {
/************
 * Clauses for errors variables
 ************/
	 if (is_genotyped(gen)) {
		lit_t e= cnf.get_e(i, l);
		lit_t p= cnf.get_p(i, l);
		lit_t m= cnf.get_m(i, l);
// Alleles are sorted ! (i.e. gen.allele1() <= gen.allele2())
		if (gen.allele1() == 1 && gen.allele2() == 1) {
// gen == HOMO1
// -e p m,  e -p,  e -m
		  cnf.add_clause(pedcnf_t::clause_t{-e,  p,  m});
		  cnf.add_clause(pedcnf_t::clause_t{ e, -p});
		  cnf.add_clause(pedcnf_t::clause_t{ e,     -m});
		} else if (gen.allele1() == 2 && gen.allele2() == 2) {
// gen == HOMO2
// -e -p -m,  e p,  e m
		  cnf.add_clause(pedcnf_t::clause_t{-e, -p, -m});
		  cnf.add_clause(pedcnf_t::clause_t{ e,  p});
		  cnf.add_clause(pedcnf_t::clause_t{ e,      m});
		} else if (gen.allele1() == 1 && gen.allele2() == 2) {
// gen == HETERO12
#ifdef AVOID_XOR_CLAUSES
		  cnf.add_clause(pedcnf_t::clause_t{ e,  p,  m});
		  cnf.add_clause(pedcnf_t::clause_t{ e, -p, -m});
		  cnf.add_clause(pedcnf_t::clause_t{-e,  p, -m});
		  cnf.add_clause(pedcnf_t::clause_t{-e, -p,  m});
#else
		  cnf.add_xor_clause(pedcnf_t::xor_clause_t{e, p, m});
#endif
		} else if (gen == g::MISS) {
// Do nothing
		} else {
		  MY_FAIL;
		}
	 }
/************
 * End clauses for errors variables
 ************/
/************
 * Clauses for recombination variables
 ************/
// (only for loci after the first)
	 if (l > 0) {
// a recombination event occurs only on heterozygous loci, thus:
//    ( p_{i,l} == m_{i,l} )  implies  not r_{i,l}
		if (ind.has_father()) {
		  lit_t p= cnf.get_p(ind.father().progr_id(), l);
		  lit_t m= cnf.get_m(ind.father().progr_id(), l);
		  lit_t rp= cnf.get_rp(i, l);
		  cnf.add_clause(pedcnf_t::clause_t{-rp,  p,  m});
		  cnf.add_clause(pedcnf_t::clause_t{-rp, -p, -m});
		}
		if (ind.has_mother()) {
		  lit_t p= cnf.get_p(ind.mother().progr_id(), l);
		  lit_t m= cnf.get_m(ind.mother().progr_id(), l);
		  lit_t rm= cnf.get_rm(i, l);
		  cnf.add_clause(pedcnf_t::clause_t{-rm,  p,  m});
		  cnf.add_clause(pedcnf_t::clause_t{-rm, -p, -m});
		}
	 }
/************
 * End clauses for recombination variables
 ************/
  };

/************
 * MULTIALLELIC VERSION
 ************/
  void add_individual_constraint_multi(pedcnf_t& cnf,
													const individual_t& ind,
													const g& gen,
													const size_t l,
													const size_t i,
													const size_t no_of_alleles) {
/************
 * Extract pm and mm variables
 ************/
	 lit_t* p= new lit_t[no_of_alleles];
	 lit_t* m= new lit_t[no_of_alleles];
	 for (size_t j= 0; j<no_of_alleles; ++j) {
		p[j]= cnf.get_pm(i, l, j);
		m[j]= cnf.get_mm(i, l, j);
	 }
/************
 * Consistency of haplotypes
 * (exactly only one TRUE among {p,m}_i,j[l])
 ************/
// at least one
	 cnf.add_clause(p, no_of_alleles);
	 cnf.add_clause(m, no_of_alleles);
// at most one (quadratic version)
	 for (size_t j1= 0; j1<no_of_alleles-1; ++j1) {
		for (size_t j2= j1+1; j2<no_of_alleles; ++j2) {
		  cnf.add_clause(pedcnf_t::clause_t{ -p[j1], -p[j2]});
		  cnf.add_clause(pedcnf_t::clause_t{ -m[j1], -m[j2]});
		}
	 }
/************
 * Clauses for errors variables
 ************/
	 if (is_genotyped(gen)) {
		const lit_t e= cnf.get_e(i, l);
		const size_t a= ((size_t)gen.allele1())-1;
		const size_t b= ((size_t)gen.allele2())-1;
		if (is_homozygous(gen)) {
		  cnf.add_clause(pedcnf_t::clause_t{  e, p[a] });
		  cnf.add_clause(pedcnf_t::clause_t{  e, m[a] });
		  cnf.add_clause(pedcnf_t::clause_t{ -e, -p[a], -m[a] });
		} else if (is_heterozygous(gen)) {
		  cnf.add_clause(pedcnf_t::clause_t{  e,  p[a],  m[a] });
		  cnf.add_clause(pedcnf_t::clause_t{  e, -p[a], -m[a] });
		  cnf.add_clause(pedcnf_t::clause_t{  e,  p[b],  m[b] });
		  cnf.add_clause(pedcnf_t::clause_t{  e, -p[b], -m[b] });

		  cnf.add_clause(pedcnf_t::clause_t{  e, -p[a], -p[b] });
		  cnf.add_clause(pedcnf_t::clause_t{  e, -m[a], -m[b] });

		  cnf.add_clause(pedcnf_t::clause_t{ -e, -p[a],  m[a],  p[b], -m[b] });
		  cnf.add_clause(pedcnf_t::clause_t{ -e,  p[a], -m[a], -p[b],  m[b] });
		} else {
		  MY_FAIL;
		}
	 }
/************
 * End clauses for errors variables
 ************/
/************
 * Clauses for recombination variables
 ************/
// (only for loci after the first)
	 if (l > 0) {
// a recombination event occurs only on heterozygous loci (on parent), thus:
//    ( p_{i,l} == m_{i,l} )  implies  not r_{i,l}
		if (ind.has_father()) {
		  lit_t* clause= new lit_t[1+(2*no_of_alleles)];
		  clause[0]= -cnf.get_rp(i, l);
		  for (size_t j= 0; j<no_of_alleles; ++j) {
			 clause[1+(2*j)]= cnf.get_pm(ind.father().progr_id(), l, j);
			 clause[2+(2*j)]= cnf.get_mm(ind.father().progr_id(), l, j);
		  }
		  for (size_t j= 0; j<no_of_alleles; ++j) {
			 clause[1+(2*j)]= -clause[1+(2*j)];
			 clause[2+(2*j)]= -clause[2+(2*j)];
			 cnf.add_clause(clause, 1+(2*no_of_alleles));
			 clause[1+(2*j)]= -clause[1+(2*j)];
			 clause[2+(2*j)]= -clause[2+(2*j)];
		  }
		  delete [] clause;
		}
		if (ind.has_mother()) {
		  lit_t* clause= new lit_t[1+(2*no_of_alleles)];
		  clause[0]= -cnf.get_rm(i, l);
		  for (size_t j= 0; j<no_of_alleles; ++j) {
			 clause[1+(2*j)]= cnf.get_pm(ind.mother().progr_id(), l, j);
			 clause[2+(2*j)]= cnf.get_mm(ind.mother().progr_id(), l, j);
		  }
		  for (size_t j= 0; j<no_of_alleles; ++j) {
			 clause[1+(2*j)]= -clause[1+(2*j)];
			 clause[2+(2*j)]= -clause[2+(2*j)];
			 cnf.add_clause(clause, 1+(2*no_of_alleles));
			 clause[1+(2*j)]= -clause[1+(2*j)];
			 clause[2+(2*j)]= -clause[2+(2*j)];
		  }
		  delete [] clause;
		}
	 }
/************
 * End clauses for recombination variables
 ************/
	 delete [] p;
	 delete [] m;
  };

  void add_parental_constraint(pedcnf_t& cnf,
										 const g& parent_g,
										 const g& individual_g,
										 const size_t locus,
										 const size_t progr_id_parent,
										 const size_t progr_id_ind,
										 const bool is_mother,
										 const size_t no_of_alleles) {
	 if (no_of_alleles<=2) {
		add_parental_constraint_bi(cnf, parent_g, individual_g,
											locus, progr_id_parent, progr_id_ind,
											is_mother);
	 } else {
		add_parental_constraint_multi(cnf, parent_g, individual_g,
												locus, progr_id_parent, progr_id_ind,
												is_mother, no_of_alleles);
	 }
  };

// maximum 2 alleles !!
  void add_parental_constraint_bi(pedcnf_t& cnf,
													  const g& parent_g,
													  const g& individual_g,
													  const size_t locus,
													  const size_t progr_id_parent,
													  const size_t progr_id_ind,
													  const bool is_mother) {
/************
 * Clauses for Mendelian consistency (s-variables)
 ************/
	 lit_t s= (!is_mother) ?
		cnf.get_sp(progr_id_ind, locus) : cnf.get_sm(progr_id_ind, locus);
	 lit_t p= cnf.get_p(progr_id_parent, locus);
	 lit_t m= cnf.get_m(progr_id_parent, locus);
	 lit_t c= (!is_mother) ?
		cnf.get_p(progr_id_ind, locus) : cnf.get_m(progr_id_ind, locus);
	 cnf.add_clause(pedcnf_t::clause_t{ s,  p, -c});
	 cnf.add_clause(pedcnf_t::clause_t{ s, -p,  c});
	 cnf.add_clause(pedcnf_t::clause_t{-s,  m, -c});
	 cnf.add_clause(pedcnf_t::clause_t{-s, -m,  c});
	 cnf.add_clause(pedcnf_t::clause_t{-p, -m,  c});
	 cnf.add_clause(pedcnf_t::clause_t{ p,  m, -c});
/************
 * End clauses for Mendelian consistency (s-variables)
 ************/
/************
 * Clauses for Recombination events (r-variables)
 ************/
	 if (locus>0) {
		lit_t prevs= (!is_mother) ?
		  cnf.get_sp(progr_id_ind, locus-1) : cnf.get_sm(progr_id_ind, locus-1);
		lit_t r= (!is_mother) ?
		  cnf.get_rp(progr_id_ind, locus) : cnf.get_rm(progr_id_ind, locus);
// s == prevs + r
#ifdef AVOID_XOR_CLAUSES
		  cnf.add_clause(pedcnf_t::clause_t{-s,  r,  prevs});
		  cnf.add_clause(pedcnf_t::clause_t{-s, -r, -prevs});
		  cnf.add_clause(pedcnf_t::clause_t{ s,  r, -prevs});
		  cnf.add_clause(pedcnf_t::clause_t{ s, -r,  prevs});
#else
		  cnf.add_xor_clause(pedcnf_t::xor_clause_t{-s, r, prevs});
#endif
	 }
/************
 * End clauses for Recombination events (r-variables)
 ************/
  };

// minimum 3 alleles !!
  void add_parental_constraint_multi(pedcnf_t& cnf,
												 const g& parent_g,
												 const g& individual_g,
												 const size_t locus,
												 const size_t progr_id_parent,
												 const size_t progr_id_ind,
												 const bool is_mother,
												 const size_t no_of_alleles) {
/************
 * Clauses for Mendelian consistency (s-variables)
 ************/
	 lit_t s= (!is_mother) ?
		cnf.get_sp(progr_id_ind, locus) : cnf.get_sm(progr_id_ind, locus);
	 for (size_t j= 0; j<no_of_alleles; ++j) {
		lit_t p= cnf.get_pm(progr_id_parent, locus, j);
		lit_t m= cnf.get_mm(progr_id_parent, locus, j);
		lit_t c= (!is_mother) ?
		  cnf.get_pm(progr_id_ind, locus, j) : cnf.get_mm(progr_id_ind, locus, j);
		cnf.add_clause(pedcnf_t::clause_t{ s,  p, -c});
		cnf.add_clause(pedcnf_t::clause_t{ s, -p,  c});
		cnf.add_clause(pedcnf_t::clause_t{-s,  m, -c});
		cnf.add_clause(pedcnf_t::clause_t{-s, -m,  c});
		cnf.add_clause(pedcnf_t::clause_t{-p, -m,  c});
		cnf.add_clause(pedcnf_t::clause_t{ p,  m, -c});
	 }
/************
 * End clauses for Mendelian consistency (s-variables)
 ************/
/************
 * Clauses for Recombination events (r-variables)
 ************/
	 if (locus>0) {
		lit_t prevs= (!is_mother) ?
		  cnf.get_sp(progr_id_ind, locus-1) : cnf.get_sm(progr_id_ind, locus-1);
		lit_t r= (!is_mother) ?
		  cnf.get_rp(progr_id_ind, locus) : cnf.get_rm(progr_id_ind, locus);
// s == prevs + r
#ifdef AVOID_XOR_CLAUSES
		  cnf.add_clause(pedcnf_t::clause_t{-s,  r,  prevs});
		  cnf.add_clause(pedcnf_t::clause_t{-s, -r, -prevs});
		  cnf.add_clause(pedcnf_t::clause_t{ s,  r, -prevs});
		  cnf.add_clause(pedcnf_t::clause_t{ s, -r,  prevs});
#else
		  cnf.add_xor_clause(pedcnf_t::xor_clause_t{-s, r, prevs});
#endif
	 }
/************
 * End clauses for Recombination events (r-variables)
 ************/
  };

public:

  void convert(const pedigree_t& ped, pedcnf_t& cnf) {
// Compute the number of alleles of each locus
	 const size_t* const max_alleles= ped.no_of_alleles();
	 BOOST_FOREACH( const individual_t& ind,
						 ped.individuals() ) {
		L_TRACE("Considering individual " << ind.progr_id());
		prepare_individual(cnf, ped, ind);
		for (size_t l= 0; l < ped.genotype_length(); ++l) {
		  L_TRACE("  locus = " << l <<
					 ", g_i = " << ind.obs_g(l));
		  if (max_alleles[l] <= 2) {
			 add_individual_constraint_bi(cnf, ind, ind.obs_g(l), l, ind.progr_id());
		  } else {
			 add_individual_constraint_multi(cnf, ind, ind.obs_g(l), l, ind.progr_id(), max_alleles[l]);
		  }
		}
		if (ind.has_father()) {
		  L_TRACE(" --> father " << ind.father().progr_id());
		  const individual_t& parent= ind.father();
		  for (size_t l= 0; l < ped.genotype_length(); ++l) {
			 L_TRACE("  locus = " << l <<
						", g_f = " << parent.obs_g(l) <<
						", g_i = " << ind.obs_g(l));
			 add_parental_constraint(cnf,
											 parent.obs_g(l), ind.obs_g(l),
											 l,
											 parent.progr_id(), ind.progr_id(),
											 false, max_alleles[l]);
		  }
		}
		if (ind.has_mother()) {
		  L_TRACE(" --> mother " << ind.mother().progr_id());
		  const individual_t& parent= ind.mother();
		  for (size_t l= 0; l < ped.genotype_length(); ++l) {
			 L_TRACE("  locus = " << l <<
						", g_f = " << parent.obs_g(l) <<
						", g_i = " << ind.obs_g(l));
			 add_parental_constraint(cnf,
											 parent.obs_g(l), ind.obs_g(l),
											 l,
											 parent.progr_id(), ind.progr_id(),
											 true, max_alleles[l]);
		  }
		}
	 }
	 delete [] max_alleles;
  };
};

template <
  typename T_GENOTYPE,
  typename T_HAPLOTYPE,
  typename T_PHENOTYPE,
  typename T_ID
>
void
ped2cnf(const basic_pedigree_t<T_GENOTYPE,
										 T_HAPLOTYPE,
										 T_PHENOTYPE,
										 T_ID>& ped,
		  pedcnf_t& cnf) {
  ped2cnf_conv_t<T_GENOTYPE, T_HAPLOTYPE, T_PHENOTYPE, T_ID> conv;
  conv.convert(ped, cnf);
}


#endif // __PED2CNF_HPP__
