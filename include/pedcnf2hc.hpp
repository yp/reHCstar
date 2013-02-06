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
 * pedcnf2hc.hpp
 *
 * Functions to convert a satisfying assignment to a (r,e)-haplotype
 * configuration.
 *
 **/

#ifndef __PEDCNF2HC_HPP__
#define __PEDCNF2HC_HPP__

#include "pedigree.hpp"
#include "pedcnf.hpp"
#include "log.hpp"

#include <log4cxx/logger.h>

#include <boost/foreach.hpp>


template <typename T_GENOTYPE,
			 typename T_HAPLOTYPE,
			 typename T_PHENOTYPE,
			 typename T_ID>
void compute_reHC_from_SAT(basic_pedigree_t<T_GENOTYPE, T_HAPLOTYPE, T_PHENOTYPE, T_ID>& ped,
									const pedcnf_t& cnf) {
  log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pedcnf2hc"));
  typedef basic_pedigree_t<T_GENOTYPE, T_HAPLOTYPE, T_PHENOTYPE, T_ID> family_t;
  INFO("Computing the (r,e)-haplotype configuration...");
  const size_t* const no_of_alleles= ped.no_of_alleles();
  size_t no_of_errors= 0;
  size_t no_of_imputation= 0;
// For each locus in each individual:
//   (1) the locus is genotyped and it is homozygous (thus the haplotype
//       is 'fixed'), or
//   (2) there exists a H_i[j] variable in the SAT instance.
  BOOST_FOREACH( typename family_t::individual_t& ind,
					  ped.individuals() ) {
	 TRACE("Considering individual " << ind.progr_id());
	 for (size_t locus= 0; locus < ped.genotype_length(); ++locus) {
		allele_t gp, gm;
		if (no_of_alleles[locus] <= 2) {
		  gp = !cnf.p(ind.progr_id(), locus) ? 1 : 2;
		  gm = !cnf.m(ind.progr_id(), locus) ? 1 : 2;
		} else {
		  for (gp= 0; gp<no_of_alleles[locus] && !cnf.pm(ind.progr_id(), locus, gp); ++gp);
		  MY_ASSERT(gp < no_of_alleles[locus]);
		  gp += 1;
		  for (gm= 0; gm<no_of_alleles[locus] && !cnf.mm(ind.progr_id(), locus, gm); ++gm);
		  MY_ASSERT(gm < no_of_alleles[locus]);
		  gm += 1;
		}
		if ( ! is_genotyped(ind.obs_g(locus)) ) {
//        Individual not genotyped ->
//          -> imputing genotype based on variables p_i_l and m_i_l
		  TRACE("Individual " << ind.progr_id() << " at locus " << locus
				  << " is not genotyped.");
		  TRACE("pil " << gp << "   mil " << gm);
		  DEBUG("Not-genotyped individual " << ind.progr_id() <<
				  " at locus " << locus << " is imputed as (" << gp << ", " << gm << ").");
		  ind.real_g(locus).set_alleles(gp, gm);
		  ++no_of_imputation;
		} else {
		  ind.real_g(locus)= ind.obs_g(locus);
		}
		if ( !is_genotyped(ind.real_g(locus)) ) {
// Impossible: we performed genotype imputation in the previous step
		  MY_FAIL;
		} else {
		  ind.hp(locus)= family_t::h::ALLELE(gp);
		  ind.hm(locus)= family_t::h::ALLELE(gm);
		  typename family_t::g real_g= ind.real_g(locus);
		  real_g.set_alleles(gp, gm);
		  if (real_g != ind.real_g(locus)) {
//   Not-genotyped loci cannot have errors
			 MY_ASSERT( is_genotyped(ind.obs_g(locus)) );
			 DEBUG("Genotyped individual " << ind.progr_id() <<
					 " at locus " << locus << " has an error: changing " <<
					 ind.obs_g(locus) << " into " << real_g << ".");
			 MY_ASSERT( cnf.e(ind.progr_id(), locus) );
			 ind.real_g(locus)= real_g;
			 ++no_of_errors;
		  }
		}
	 }
  }
  delete [] no_of_alleles;
  INFO("Number of imputed genotypes: " << no_of_imputation);
  INFO("Number of corrected errors:  " << no_of_errors);
  INFO("(r,e)-haplotype configuration successfully computed.");
};



#endif // __PEDCNF2HC_HPP__
