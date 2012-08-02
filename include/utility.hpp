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
 * utility.hpp
 *
 * Various utility functions and classes.
 *
 **/

#ifndef __UTILITY_HPP__
#define __UTILITY_HPP__

#include <functional>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/version.hpp>

#include "assertion.hpp"

#ifndef EXIT_NO_reHC
#define EXIT_NO_reHC (2)
#endif

#ifndef EXIT_reHC_ERROR
#define EXIT_reHC_ERROR (4)
#endif


template <typename T>
std::string tostr(const std::vector<T>& x) {
  std::ostringstream out;
  bool first= true;
  for (typename std::vector<T>::const_iterator it= x.begin();
		 it != x.end();
		 ++it) {
	 if (!first) {
		out << " ";
	 }
	 out << *it;
	 first= false;
  }
  return out.str();
};

template <typename T>
std::string tostr(const T& x) {
  std::ostringstream ostr;
  ostr << x;
  return ostr.str();
};

template <class T>
class deleter_t
  : public std::unary_function<T*, void> {
public:
  void operator()(T* pt) {
	 delete pt;
  }
};

class generic_raise_error_class_t
  :public log_able_t<generic_raise_error_class_t>
{
public:
  void raise_error(const std::string& reason="no detail specified") {
	 L_FATAL("Generic error. Reason: " << reason << ". Exiting...");
	 MY_FAIL;
  };
};

template <int i, int default_index, class derived_enum, class raise_error_class_t>
class static_str_switch {
public:
  static inline derived_enum EXEC(const std::string& s, raise_error_class_t& err) {
	 if (s == derived_enum::str_values[i]) {
		return derived_enum::enum_values[i];
	 } else {
		return static_str_switch< i-1, default_index, derived_enum, raise_error_class_t>::EXEC(s, err);
	 }
  }
};

template <int default_index, class derived_enum, class raise_error_class_t>
class static_str_switch<0, default_index, derived_enum, raise_error_class_t> {
public:
  static inline derived_enum EXEC(const std::string& s, raise_error_class_t& err) {
	 if (s == derived_enum::str_values[0]) {
		return derived_enum::enum_values[0];
	 } else {
		err.raise_error("string >" + s + "< not recognized");
		return derived_enum::enum_values[default_index];
	 }
  }
};

template <int i, class derived_enum>
class static_int_switch {
public:
  static inline int EXEC(const int v) {
	 if (v == derived_enum::int_values[i]) {
		return i;
	 } else {
		return static_int_switch< i-1, derived_enum>::EXEC(v);
	 }
  }
};

template <class derived_enum>
class static_int_switch<0, derived_enum> {
public:
  static inline int EXEC(const int v) {
	 if (v == derived_enum::int_values[0]) {
		return 0;
	 } else {
		MY_FAIL;
		return -1;
	 }
  }
};



template <class derived_enum, size_t n_values, int default_index>
class enum_like_t {
private:
  int _d;

protected:

  class raise_input_error_t {
  private:
	 std::istream& in;
  public:
	 raise_input_error_t(std::istream& _in)
		  :in(_in)
	 {};
	 void raise_error(const std::string& reason="no deatail provided") {
		in.setstate(std::ios::failbit);
	 };
  };

  enum_like_t(const size_t d)
		:_d(d)
  {
	 MY_ASSERT_DBG((0 <= d) && (d<n_values));
  }

  size_t get_ordinal_data() const {
	 return _d;
  }

public:

  enum_like_t(const derived_enum& copy)
		: _d(copy.get_ordinal_data())
  {
	 MY_ASSERT_DBG( copy.get_ordinal_data() <n_values );
  };

  const derived_enum&
  operator=(const derived_enum& copy) {
	 _d= copy.get_ordinal_data();
	 return *this;
  };

  const static size_t N_VALUES= n_values;

  friend bool operator==(const derived_enum& e1, const derived_enum& e2) {
	 return e1.get_ordinal_data() == e2.get_ordinal_data();
  }

  friend bool operator!=(const derived_enum& e1, const derived_enum& e2) {
	 return e1.get_ordinal_data() != e2.get_ordinal_data();
  }

  friend bool operator<=(const derived_enum& e1, const derived_enum& e2) {
	 return e1.get_ordinal_data() <= e2.get_ordinal_data();
  }

  friend bool operator<(const derived_enum& e1, const derived_enum& e2) {
	 return e1.get_ordinal_data() < e2.get_ordinal_data();
  }

  friend std::ostream& operator<<(std::ostream& out, const derived_enum& val) {
	 return (out << derived_enum::str_values[val.get_ordinal_data()]);
  }

  friend std::istream& operator>>(std::istream& in, derived_enum& val) {
	 std::string h;
	 typename derived_enum::raise_input_error_t err(in);
	 if (in >> h) {
		val= static_str_switch< n_values-1, default_index, derived_enum, raise_input_error_t>::EXEC(h, err);
	 } else {
		err.raise_error("generic error");
		val= derived_enum::enum_values[default_index];
	 }
	 return in;
  }
};




class file_utility:
  public log_able_t<file_utility> {

private:

  file_utility() {};

public:

  typedef boost::shared_ptr<std::istream> pistream;
  typedef boost::shared_ptr<std::ostream> postream;

  static const file_utility&
  get_file_utility() {
	 static const file_utility fu;
	 return fu;
  };

  pistream get_ifstream(const std::string& file_name,
								const bool compressed=false) const {
	 L_DEBUG("Opening file '" << file_name << "' for reading...");
	 boost::iostreams::file_source source(file_name);
	 if (!source.is_open()) {
		L_ERROR("Impossible to open file '" << file_name << "' for reading.");
		throw std::logic_error(std::string("Impossible to open file '")
									  + file_name + "' for reading.");
	 }
	 L_DEBUG("File '" << file_name << "' successfully opened.");
	 std::istream* is;
	 if (compressed) {
		L_DEBUG("...setting file '" << file_name << "' as gzip-ed.");
		boost::iostreams::filtering_istream* gis= new boost::iostreams::filtering_istream;
		gis->push(boost::iostreams::gzip_decompressor());
		gis->push(source);
		is= gis;
	 } else {
		is= new boost::iostreams::stream<boost::iostreams::file_source>(source);
	 }
	 return pistream(is);
  };

  postream get_ofstream(const std::string& file_name,
								const bool compressed=false) const {
	 L_DEBUG("Opening file '" << file_name << "' for writing...");
	 boost::iostreams::file_sink sink(file_name);
	 if (!sink.is_open()) {
		L_ERROR("Impossible to open file '" << file_name << "' for writing.");
		throw std::logic_error(std::string("Impossible to open file '")
									  + file_name + "' for writing.");
	 }
	 L_DEBUG("File '" << file_name << "' successfully opened.");
	 std::ostream* os;
	 if (compressed) {
		L_DEBUG("...setting file '" << file_name << "' as gzip-ed.");
		boost::iostreams::filtering_ostream* gos= new boost::iostreams::filtering_ostream;
		gos->push(boost::iostreams::gzip_compressor());
		gos->push(sink);
		os= gos;
	 } else {
		os= new boost::iostreams::stream<boost::iostreams::file_sink>(sink);
	 }
	 return postream(os);
  };

  postream get_tmp_ostream(const std::string file_template,
									std::string& real_name,
									const bool compressed=false) const {
	 L_DEBUG("Opening a temporary file with hint '" << file_template <<
				"' for writing...");
	 MY_ASSERT(boost::ends_with(file_template, "XXXXXX"));
	 const size_t name_len= strlen(file_template.c_str());
	 char* name= (char*)calloc(name_len+1, sizeof(char));
	 strncpy(name, file_template.c_str(), name_len+1);
	 const int fd= mkstemp(name);
	 real_name= std::string(name);
	 free(name);
	 if (fd == -1) {
		L_ERROR("Impossible to open a temporary file with hint '" << file_template
				  << "' for writing.");
		throw std::logic_error(std::string("Impossible to open a temporary file with hint '")
									  + file_template + "' for writing.");
	 }
#if (BOOST_VERSION >= 104400)
	 boost::iostreams::file_descriptor_sink sink(fd, boost::iostreams::close_handle);
#else
	 boost::iostreams::file_descriptor_sink sink(fd);
#endif
	 if (!sink.is_open()) {
		L_ERROR("Impossible to open file '" << real_name << "' for writing.");
		throw std::logic_error(std::string("Impossible to open file '")
									  + real_name + "' for writing.");
	 }
	 L_INFO("File '" << real_name << "' successfully opened.");
	 std::ostream* os;
	 if (compressed) {
		L_DEBUG("...setting file '" << real_name << "' as gzip-ed.");
		boost::iostreams::filtering_ostream* gos= new boost::iostreams::filtering_ostream;
		gos->push(boost::iostreams::gzip_compressor());
		gos->push(sink);
		os= gos;
	 } else {
		os= new boost::iostreams::stream<boost::iostreams::file_descriptor_sink>(sink);
	 }
	 postream pos(os);
	 return pos;
  };

  postream get_ostream_from_file(FILE* file) const {
	 MY_ASSERT(file != NULL);
	 const int fd = fileno(file);
	 if (fd == -1) {
		L_ERROR("Impossible to open a stream associated to the given C-file.");
		throw std::logic_error(std::string("Impossible to open a stream associated to the given C-file."));
	 }
#if (BOOST_VERSION >= 104400)
	 boost::iostreams::file_descriptor_sink sink(fd, boost::iostreams::close_handle);
#else
	 boost::iostreams::file_descriptor_sink sink(fd);
#endif
	 if (!sink.is_open()) {
		L_ERROR("Impossible to open the given C-file.");
		throw std::logic_error(std::string("Impossible to open the given C-file."));
	 }
	 std::ostream* os= new boost::iostreams::stream<boost::iostreams::file_descriptor_sink>(sink);
	 postream pos(os);
	 return pos;
  };

};


class percent_t {

public:
  const double val;
  const double tot;

  percent_t(const double val_,
				const double tot_)
		:val(val_), tot(tot_)
  {};
};

std::ostream&
operator<<(std::ostream& out, const percent_t& val);

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


#endif // __UTILITY_HPP__
