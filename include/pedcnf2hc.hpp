/**
 *
 *                              ZRHC-*
 * Zero-Recombinant Haplotype Configuration with missing genotypes
 *
 * Copyright (C) 2010,2011  Yuri Pirola <yuri.pirola(-at-)gmail.com>
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
 * pedcnf2hc.hpp
 *
 * Functions to convert a satisfying assignment to a zero-recombinant haplotype
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
void compute_ZRHC_from_SAT(basic_pedigree_t<T_GENOTYPE, T_HAPLOTYPE, T_PHENOTYPE, T_ID>& ped,
									const pedcnf_t& cnf) {
  log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pedcnf2hc"));
  typedef basic_pedigree_t<T_GENOTYPE, T_HAPLOTYPE, T_PHENOTYPE, T_ID> family_t;
  INFO("Computing the zero-recombinant haplotype configuration...");
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
		const bool pil= cnf.p(ind.progr_id(), locus);
		const bool mil= cnf.m(ind.progr_id(), locus);
		if ( ! is_genotyped(ind.obs_g(locus)) ) {
//        Individual not genotyped ->
//          -> imputing genotype based on variables p_i_l and m_i_l
		  TRACE("Individual " << ind.progr_id() << " at locus " << locus
				  << " is not genotyped.");
		  TRACE("pil " << pil << "   mil " << mil);
		  if ( pil == mil ) {
			 DEBUG("Not-genotyped individual " << ind.progr_id() <<
					 " at locus " << locus << " is imputed as homozygous.");
			 if ( !pil ) {
				ind.real_g(locus)= family_t::g::HOMO1;
			 } else {
				ind.real_g(locus)= family_t::g::HOMO2;
			 }
		  } else {
			 DEBUG("Not-genotyped individual " << ind.progr_id() <<
					 " at locus " << locus << " is imputed as heterozygous.");
			 ind.real_g(locus)= family_t::g::HETER;
		  }
		  ++no_of_imputation;
		} else {
		  ind.real_g(locus)= ind.obs_g(locus);
		}
		if ( !is_genotyped(ind.real_g(locus)) ) {
// Impossible: we performed genotype imputation in the previous step
		  MY_FAIL;
		} else {
		  ind.hp(locus)= (!pil) ? (family_t::h::ALLELE1) : (family_t::h::ALLELE2);
		  ind.hm(locus)= (!mil) ? (family_t::h::ALLELE1) : (family_t::h::ALLELE2);
		  typename family_t::g real_g= ind.real_g(locus);
		  if ( pil == mil ) {
			 if ( !pil ) {
				real_g= family_t::g::HOMO1;
			 } else {
				real_g= family_t::g::HOMO2;
			 }
		  } else {
			 real_g= family_t::g::HETER;
		  }
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
  INFO("Number of imputed genotypes: " << no_of_imputation);
  INFO("Number of corrected errors:  " << no_of_errors);
  INFO("Zero-recombinant haplotype configuration successfully computed.");
};



#endif // __PEDCNF2HC_HPP__
