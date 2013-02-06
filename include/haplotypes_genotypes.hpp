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
 * haplotypes_genotypes.hpp
 *
 * Classes to represent haplotypes and genotypes.
 *
 **/

#ifndef __HAPLOTYPES_GENOTYPES_HPP__
#define __HAPLOTYPES_GENOTYPES_HPP__

#include <cstdlib>
#include <vector>
#include <iostream>
#include <algorithm>
#include <boost/static_assert.hpp>

#include "log.hpp"
#include "assertion.hpp"
#include "utility.hpp"


/**
 *
 *
 *
 * Representation of single-locus genotypes
 *
 *
 *
 **/

typedef unsigned int allele_t;


class single_multiallelic_genotype_t {

private:
  allele_t _allele1;
  allele_t _allele2;

public:

  static const single_multiallelic_genotype_t MISS;

  single_multiallelic_genotype_t()
		:_allele1(0), _allele2(0)
  {};

  explicit single_multiallelic_genotype_t(const allele_t allele1,
														const allele_t allele2)
		:_allele1(std::min(allele1, allele2)),
		 _allele2(std::max(allele1, allele2))
  {};

  single_multiallelic_genotype_t(const single_multiallelic_genotype_t& g)
		:_allele1(g._allele1), _allele2(g._allele2)
  {};

  allele_t allele1() const {
	 return _allele1;
  };

  allele_t allele2() const {
	 return _allele2;
  };

  void set_alleles(const allele_t allele1,
						 const allele_t allele2) {
	 if (((allele1 == 0)&&(allele2 != 0))||
		  ((allele1 != 0)&&(allele2 == 0))) {
		MY_FAIL;
	 }
	 _allele1= std::min(allele1, allele2);
	 _allele2= std::max(allele1, allele2);
  };

};

bool
operator==(const single_multiallelic_genotype_t& g1,
			  const single_multiallelic_genotype_t& g2);

bool
operator!=(const single_multiallelic_genotype_t& g1,
			  const single_multiallelic_genotype_t& g2);

bool
operator<=(const single_multiallelic_genotype_t& g1,
			  const single_multiallelic_genotype_t& g2);

bool
operator<(const single_multiallelic_genotype_t& g1,
			 const single_multiallelic_genotype_t& g2);

std::ostream&
operator<<(std::ostream& out,
			  const single_multiallelic_genotype_t& g);

std::istream&
operator>>(std::istream& in,
			  single_multiallelic_genotype_t& g);

bool
is_genotyped(const single_multiallelic_genotype_t& g);

bool
is_homozygous(const single_multiallelic_genotype_t& g);

bool
is_heterozygous(const single_multiallelic_genotype_t& g);


/**
 *
 *
 *
 * Representation of single-locus haplotypes
 *
 *
 *
 **/



class single_multiallelic_haplotype_t {
private:
  allele_t _allele;
public:

  static const single_multiallelic_haplotype_t MISS;

  single_multiallelic_haplotype_t()
		:_allele(0)
  {};

  explicit single_multiallelic_haplotype_t(const allele_t allele)
		:_allele(allele)
  {};

  single_multiallelic_haplotype_t(const single_multiallelic_haplotype_t& h)
		:_allele(h._allele)
  {};

  allele_t allele() const {
	 return _allele;
  };

  static single_multiallelic_haplotype_t ALLELE(const allele_t allele) {
	 MY_ASSERT_DBG(allele != single_multiallelic_haplotype_t::MISS.allele());
	 return single_multiallelic_haplotype_t(allele);
  };

};

bool
operator==(const single_multiallelic_haplotype_t& h1,
			  const single_multiallelic_haplotype_t& h2);

bool
operator!=(const single_multiallelic_haplotype_t& h1,
			  const single_multiallelic_haplotype_t& h2);

bool
operator<=(const single_multiallelic_haplotype_t& h1,
			  const single_multiallelic_haplotype_t& h2);

bool
operator<(const single_multiallelic_haplotype_t& h1,
			 const single_multiallelic_haplotype_t& h2);

std::ostream&
operator<<(std::ostream& out,
			  const single_multiallelic_haplotype_t& h);

std::istream&
operator>>(std::istream& in,
			  single_multiallelic_haplotype_t& h);

bool
is_missing(const single_multiallelic_haplotype_t& h);


const single_multiallelic_haplotype_t
homozygous_to_haplotype(const single_multiallelic_genotype_t& g);

bool
haplotype_genotype_consistent(const single_multiallelic_haplotype_t& h1,
										const single_multiallelic_haplotype_t& h2,
										const single_multiallelic_genotype_t& g);

bool
strict_haplotype_genotype_consistent(const single_multiallelic_haplotype_t& h1,
												 const single_multiallelic_haplotype_t& h2,
												 const single_multiallelic_genotype_t& g);



/**
 *
 *
 *
 * Representation of multi-locus genotypes and haplotype
 *
 *
 *
 **/

template <class base_t>
class default_base_writer_t {
private:
  std::ostream& _out;
public:
  default_base_writer_t(std::ostream& out)
		: _out(out)
  {}

  std::ostream& operator()(const base_t& val) {
	 return (_out << val);
  }
};

template <class base_t>
class default_base_reader_t {
private:
  std::istream& _in;
public:
  default_base_reader_t(std::istream& in)
		: _in(in)
  {}

  std::istream& operator()(base_t& val) {
	 return (_in >> val);
  }
};


template <class _base_t,
			 class _reader=default_base_reader_t<_base_t>,
			 class _writer=default_base_writer_t<_base_t> >
class generic_fixlen_vector_t {
public:

  typedef _base_t base;
  typedef _reader reader;
  typedef _writer writer;

  typedef base* iterator;
  typedef const base* const_iterator;

private:
  const size_t _len;
  base* _v;

public:

  generic_fixlen_vector_t(const size_t len)
		:_len(len), _v(new base[len])
  {}

  generic_fixlen_vector_t(const std::vector<base> v)
		:_len(v.size()), _v(new base[v.size()])
  {
	 std::copy(v.begin(), v.end(), _v);
  }

  ~generic_fixlen_vector_t() {
	 delete [] _v;
  }

  const base& operator[](const size_t pos) const {
	 MY_ASSERT_DBG(pos < _len);
	 return _v[pos];
  }

  base& operator[](const size_t pos) {
	 MY_ASSERT_DBG(pos < _len);
	 return _v[pos];
  }

  base const * begin() const {
	 return _v;
  }

  base* begin() {
	 return _v;
  }

  base const * end() const {
	 return _v+_len;
  }

  base * end() {
	 return _v+_len;
  }

  friend std::ostream& operator<<(std::ostream& out,
											 const generic_fixlen_vector_t& val) {
	 std::for_each(val._v, val._v+val._len, writer(out));
	 return out;
  }

  friend std::istream& operator>>(std::istream& in, generic_fixlen_vector_t& val) {
	 std::for_each(val._v, val._v+val._len, reader(in));
	 return in;
  }

  size_t size() const {
	 return _len;
  };

  bool is_compatible_with(const generic_fixlen_vector_t<_base_t,_reader,_writer>& v) const {
	 typedef generic_fixlen_vector_t<_base_t,_reader,_writer> v_t;
	 if (size() != v.size())
		return false;
	 const typename v_t::base* v1it= begin();
	 const typename v_t::base* v2it= v.begin();
	 for (; v1it != end(); ++v1it, ++v2it) {
		if ( (is_genotyped(*v1it)) &&
			  (is_genotyped(*v2it)) &&
			  (*v1it != *v2it) ) {
		  return false;
		}
	 }
	 return true;
  };

};

template <typename h_t, typename g_t>
bool
multilocus_haplotype_genotype_consistent(const h_t& h1,
													  const h_t& h2,
													  const g_t& g) {
  MY_ASSERT_DBG(g.size() == h1.size());
  MY_ASSERT_DBG(g.size() == h2.size());
  const typename g_t::base* git= g.begin();
  const typename h_t::base* h1it= h1.begin();
  const typename h_t::base* h2it= h2.begin();
  for (; git != g.end(); ++git, ++h1it, ++h2it) {
	 if (! haplotype_genotype_consistent(*h1it, *h2it, *git)) {
		return false;
	 }
  }
  return true;
}

template <typename h_t, typename g_t>
bool
strict_multilocus_haplotype_genotype_consistent(const h_t& h1,
																const h_t& h2,
																const g_t& g) {
  MY_ASSERT_DBG(g.size() == h1.size());
  MY_ASSERT_DBG(g.size() == h2.size());
  typename g_t::const_iterator git= g.begin();
  typename h_t::const_iterator h1it= h1.begin();
  typename h_t::const_iterator h2it= h2.begin();
  for (; git != g.end(); ++git, ++h1it, ++h2it) {
	 if (! strict_haplotype_genotype_consistent(*h1it, *h2it, *git)) {
		return false;
	 }
  }
  return true;
}

template <typename h_t>
int
strict_multilocus_mendelian_consistent(const h_t& hpf,
													const h_t& hpm,
													const h_t& h) {
  MY_ASSERT_DBG(hpf.size() == hpm.size());
  MY_ASSERT_DBG(hpf.size() == h.size());
  typename h_t::const_iterator hpfit= hpf.begin();
  typename h_t::const_iterator hpmit= hpm.begin();
  typename h_t::const_iterator* its[]= { &hpfit, &hpmit };
  typename h_t::const_iterator hit= h.begin();
  int phase1= 0; int recomb1= 0;
  int phase2= 1; int recomb2= 0;
  int locus= 0;
  for (; hit != h.end(); ++hpfit, ++hpmit, ++hit, ++locus) {
	 if (*hit != **its[phase1]) {
		phase1= (phase1+1)%2;
		if (*hit != **its[phase1]) {
		  return -locus;
		} else {
		  ++recomb1;
		}
	 }
	 if (*hit != **its[phase2]) {
		phase2= (phase2+1)%2;
		if (*hit != **its[phase2]) {
		  return -locus;
		} else {
		  ++recomb2;
		}
	 }
  }
  return std::min(recomb1, recomb2);
}


typedef generic_fixlen_vector_t<single_multiallelic_genotype_t> genotype_t;
typedef generic_fixlen_vector_t<single_multiallelic_haplotype_t> haplotype_t;

template <class _base_t,
			 class _reader,
			 class _writer>
bool operator==(const generic_fixlen_vector_t<_base_t,_reader,_writer>& v1,
					 const generic_fixlen_vector_t<_base_t,_reader,_writer>& v2) {
  typedef generic_fixlen_vector_t<_base_t,_reader,_writer> v_t;
  if (v1.size() != v2.size())
	 return false;
  const typename v_t::base* v1it= v1.begin();
  const typename v_t::base* v2it= v2.begin();
  for (; v1it != v1.end(); ++v1it, ++v2it) {
	 if ( (*v1it) != (*v2it) ) {
		return false;
	 }
  }
  return true;
};

template <class _base_t,
			 class _reader,
			 class _writer>
bool operator!=(const generic_fixlen_vector_t<_base_t,_reader,_writer>& v1,
					 const generic_fixlen_vector_t<_base_t,_reader,_writer>& v2) {
  typedef generic_fixlen_vector_t<_base_t,_reader,_writer> v_t;
  if (v1.size() != v2.size())
	 return true;
  const typename v_t::base* v1it= v1.begin();
  const typename v_t::base* v2it= v2.begin();
  for (; v1it != v1.end(); ++v1it, ++v2it) {
	 if ( (*v1it) != (*v2it) ) {
		return true;
	 }
  }
  return false;
};

#endif // __HAPLOTYPES_GENOTYPES_HPP__
