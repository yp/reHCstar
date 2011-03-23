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
								  const individual_t& ind,
								  const double error_rate) {
// Create the variables for the individual
//   -> p-variables
	 for (size_t l= 0; l < ped.genotype_length(); ++l) {
		cnf.get_p(ind.progr_id(), l);
	 }
//   -> m-variables
	 for (size_t l= 0; l < ped.genotype_length(); ++l) {
		cnf.get_m(ind.progr_id(), l);
	 }
//   -> e-variables
	 std::vector<var_t> evars;
	 for (size_t l= 0; l < ped.genotype_length(); ++l) {
//    errors only for genotyped loci
		if (is_genotyped(ind.obs_g(l))) {
		  const var_t& e= cnf.get_e(ind.progr_id(), l);
		  evars.push_back(e);
		}
	 }
/************
 * Clauses for limiting the number of errors
 ************/
	 const size_t k= std::ceil(error_rate*evars.size());
	 L_DEBUG("Generating cardinality constraints for " << evars.size() <<
				" error variables (<= " << k << ")...");
	 add_card_constraint_less_or_equal_than(cnf, evars, k);
/************
 * End clauses for limiting the number of errors
 ************/
  }

  void add_individual_constraint(pedcnf_t& cnf,
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
		} else if (gen == g::HETER) {
		  cnf.add_xor_clause<3>((lit_t[]){e, p, m});
		  // cnf.add_clause<3>((lit_t[]){ e,  p,  m});
		  // cnf.add_clause<3>((lit_t[]){ e, -p, -m});
		  // cnf.add_clause<3>((lit_t[]){-e,  p, -m});
		  // cnf.add_clause<3>((lit_t[]){-e, -p,  m});
		} else if (gen == g::MISS) {
// Do nothing
		} else {
		  MY_FAIL;
		}
	 }
/************
 * End clauses for errors variables
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
	 lit_t s= cnf.get_s(progr_id_parent, progr_id_ind);
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
  };

public:

  pedcnf_t* convert(const pedigree_t& ped, const double error_rate) {
	 pedcnf_t* pcnf= new pedcnf_t;
	 pedcnf_t& cnf= *pcnf;
	 BOOST_FOREACH( const individual_t& ind,
						 ped.individuals() ) {
		L_TRACE("Considering individual " << ind.progr_id());
		prepare_individual(cnf, ped, ind, error_rate);
		for (size_t l= 0; l < ped.genotype_length(); ++l) {
		  L_TRACE("  locus = " << l <<
					 ", g_i = " << ind.obs_g(l));
		  add_individual_constraint(cnf, ind.obs_g(l), l, ind.progr_id());
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
	 L_INFO("The SAT instance is composed by " <<
			  std::setw(8) << cnf.vars().size() << " variables and " <<
			  std::setw(8) << cnf.no_of_clauses() << " clauses");
	 return pcnf;
  };
};

template <
  typename T_GENOTYPE,
  typename T_HAPLOTYPE,
  typename T_PHENOTYPE,
  typename T_ID
>
pedcnf_t*
ped2cnf(const basic_pedigree_t<T_GENOTYPE,
										 T_HAPLOTYPE,
										 T_PHENOTYPE,
										 T_ID>& ped,
		  const double error_rate) {
  ped2cnf_conv_t<T_GENOTYPE, T_HAPLOTYPE, T_PHENOTYPE, T_ID> conv;
  return conv.convert(ped, error_rate);
}


#endif // __PED2CNF_HPP__
