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

  void add_individual_constraint(pedcnf_t& cnf,
											const individual_t& ind,
											const g& gen,
											const size_t l,
											const size_t i) {
/************
 * Clauses for errors variables
 ************/
	 if (gen != g::MISS) {
		lit_t e= cnf.get_e(i, l);
		lit_t p= cnf.get_p(i, l);
		lit_t m= cnf.get_m(i, l);
		if (gen == g::HOMO1) {
// -e p m,  e -p,  e -m
		  cnf.add_clause<3>((lit_t[]){-e,  p,  m});
		  cnf.add_clause<2>((lit_t[]){ e, -p});
		  cnf.add_clause<2>((lit_t[]){ e,     -m});
		} else if (gen == g::HOMO2) {
// -e -p -m,  e p,  e m
		  cnf.add_clause<3>((lit_t[]){-e, -p, -m});
		  cnf.add_clause<2>((lit_t[]){ e,  p});
		  cnf.add_clause<2>((lit_t[]){ e,      m});
		} else if (gen == g::HETER12) {
#ifdef AVOID_XOR_CLAUSES
		  cnf.add_clause<3>((lit_t[]){ e,  p,  m});
		  cnf.add_clause<3>((lit_t[]){ e, -p, -m});
		  cnf.add_clause<3>((lit_t[]){-e,  p, -m});
		  cnf.add_clause<3>((lit_t[]){-e, -p,  m});
#else
		  cnf.add_xor_clause<3>((lit_t[]){e, p, m});
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
		  cnf.add_clause<3>((lit_t[]){-rp,  p,  m});
		  cnf.add_clause<3>((lit_t[]){-rp, -p, -m});
		}
		if (ind.has_mother()) {
		  lit_t p= cnf.get_p(ind.mother().progr_id(), l);
		  lit_t m= cnf.get_m(ind.mother().progr_id(), l);
		  lit_t rm= cnf.get_rm(i, l);
		  cnf.add_clause<3>((lit_t[]){-rm,  p,  m});
		  cnf.add_clause<3>((lit_t[]){-rm, -p, -m});
		}
	 }
/************
 * End clauses for recombination variables
 ************/
  };

  void add_parental_constraint(pedcnf_t& cnf,
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
	 cnf.add_clause<3>((lit_t[]){ s,  p, -c});
	 cnf.add_clause<3>((lit_t[]){ s, -p,  c});
	 cnf.add_clause<3>((lit_t[]){-s,  m, -c});
	 cnf.add_clause<3>((lit_t[]){-s, -m,  c});
	 cnf.add_clause<3>((lit_t[]){-p, -m,  c});
	 cnf.add_clause<3>((lit_t[]){ p,  m, -c});
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
		  cnf.add_clause<3>((lit_t[]){-s,  r,  prevs});
		  cnf.add_clause<3>((lit_t[]){-s, -r, -prevs});
		  cnf.add_clause<3>((lit_t[]){ s,  r, -prevs});
		  cnf.add_clause<3>((lit_t[]){ s, -r,  prevs});
#else
		  cnf.add_xor_clause<3>((lit_t[]){-s, r, prevs});
#endif
	 }
/************
 * End clauses for Recombination events (r-variables)
 ************/
  };

public:

  void convert(const pedigree_t& ped, pedcnf_t& cnf) {
	 BOOST_FOREACH( const individual_t& ind,
						 ped.individuals() ) {
		L_TRACE("Considering individual " << ind.progr_id());
		prepare_individual(cnf, ped, ind);
		for (size_t l= 0; l < ped.genotype_length(); ++l) {
		  L_TRACE("  locus = " << l <<
					 ", g_i = " << ind.obs_g(l));
		  add_individual_constraint(cnf, ind, ind.obs_g(l), l, ind.progr_id());
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
											 false);
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
											 true);
		  }
		}
	 }
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
