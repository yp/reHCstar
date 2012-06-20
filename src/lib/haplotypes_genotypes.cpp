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
 * haplotypes_genotypes.cpp
 *
 * Classes to represent haplotypes and genotypes.
 *
 **/

#include "haplotypes_genotypes.hpp"



const single_multiallelic_genotype_t
single_multiallelic_genotype_t::MISS(0, 0);

bool operator==(const single_multiallelic_genotype_t& g1,
					 const single_multiallelic_genotype_t& g2) {
  return (g1.allele1()==g2.allele1()) &&
	 (g1.allele2()==g2.allele2());
};

bool operator!=(const single_multiallelic_genotype_t& g1,
					 const single_multiallelic_genotype_t& g2) {
  return (g1.allele1()!=g2.allele1()) ||
	 (g1.allele2()!=g2.allele2());
};

bool operator<=(const single_multiallelic_genotype_t& g1,
					 const single_multiallelic_genotype_t& g2) {
  return (g1.allele1()<g2.allele1()) ||
	 (g1.allele1()==g2.allele1() && g1.allele2()<=g2.allele2());
};

bool operator<(const single_multiallelic_genotype_t& g1,
					const single_multiallelic_genotype_t& g2) {
  return (g1.allele1()<g2.allele1()) ||
	 (g1.allele1()==g2.allele1() && g1.allele2()<g2.allele2());
};

std::ostream& operator<<(std::ostream& out,
								 const single_multiallelic_genotype_t& g) {
  return (out << g.allele1() << " " << g.allele2());
};

std::istream& operator>>(std::istream& in,
								 single_multiallelic_genotype_t& g) {
  allele_t a1(0), a2(0);
  if (in >> a1 >> a2) {
	 g.set_alleles(a1, a2);
  }
  return in;
};


bool
is_genotyped(const single_multiallelic_genotype_t& g) {
  return (g.allele1() != 0)&&(g.allele2()!=0);
};

bool
is_homozygous(const single_multiallelic_genotype_t& g) {
  return is_genotyped(g) &&
	 (g.allele1() == g.allele2());
};

bool
is_heterozygous(const single_multiallelic_genotype_t& g) {
  return is_genotyped(g) &&
	 (g.allele1() != g.allele2());
};


const single_multiallelic_haplotype_t
single_multiallelic_haplotype_t::MISS(0);


bool operator==(const single_multiallelic_haplotype_t& h1,
					 const single_multiallelic_haplotype_t& h2) {
  return (h1.allele()==h2.allele());
};

bool operator!=(const single_multiallelic_haplotype_t& h1,
					 const single_multiallelic_haplotype_t& h2) {
  return (h1.allele()!=h2.allele());
};

bool operator<=(const single_multiallelic_haplotype_t& h1,
					 const single_multiallelic_haplotype_t& h2) {
  return (h1.allele()<=h2.allele());
};

bool operator<(const single_multiallelic_haplotype_t& h1,
					const single_multiallelic_haplotype_t& h2) {
  return (h1.allele()<h2.allele());
};

std::ostream& operator<<(std::ostream& out,
								 const single_multiallelic_haplotype_t& h) {
  return (out << h.allele());
};

std::istream& operator>>(std::istream& in,
								 single_multiallelic_haplotype_t& h) {
  allele_t a;
  if (in >> a) {
	 h= single_multiallelic_haplotype_t(a);
  }
  return in;
};

bool
is_missing(const single_multiallelic_haplotype_t& h) {
  return h.allele() == 0;
};


const single_multiallelic_haplotype_t
homozygous_to_haplotype(const single_multiallelic_genotype_t& g) {
  MY_ASSERT_DBG( is_homozygous(g) );
  return single_multiallelic_haplotype_t(g.allele1());
};

bool
haplotype_genotype_consistent(const single_multiallelic_haplotype_t& h1,
										const single_multiallelic_haplotype_t& h2,
										const single_multiallelic_genotype_t& g) {
  if (is_missing(h1) && is_missing(h2))
	 return true;
  if (! is_genotyped(g))
	 return true;
  if (h1.allele() > h2.allele())
	 return haplotype_genotype_consistent(h2, h1, g);

  MY_ASSERT( !is_missing(h2) );
  MY_ASSERT( is_genotyped(g) );

  if (is_missing(h1)) {
	 return
		(h2.allele() == g.allele1()) ||
		(h2.allele() == g.allele2());
  } else {
	 return
		(h1.allele() == g.allele1()) &&
		(h2.allele() == g.allele2());
  }
};

bool
strict_haplotype_genotype_consistent(const single_multiallelic_haplotype_t& h1,
												 const single_multiallelic_haplotype_t& h2,
												 const single_multiallelic_genotype_t& g) {
  if (is_missing(h1) || is_missing(h2)) {
	 return false;
  } else {
	 return haplotype_genotype_consistent(h1, h2, g);
  }
};


