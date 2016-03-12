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
#include "assumptions.hpp"

#include "log.hpp"
#include "assertion.hpp"
#include "data.hpp"
#include "utility.hpp"


#include <boost/algorithm/string/trim.hpp>

std::istream&
operator>>(std::istream& in, pedcnf_t::pedvar_t& var) {
  ped_var_kind var_kind= ped_var_kind::DUMMY;
  size_t i1, i2, i3;
  in >> var_kind;
  in >> i1;
  in >> i2;
  in >> i3;
  var= boost::make_tuple(var_kind, i1, i2, i3);
  return in;
};


class assumption_manager_t
  :public log_able_t<assumption_manager_t> {
public:
  typedef plink_reader_t<>::multifamily_pedigree_t::pedigree_t pedigree_t;

private:
  std::istream& in;
  const pedigree_t& ped;
  pedcnf_t& cnf;

  class invalid_line_t {
  public:
	 const std::string msg;
	 invalid_line_t(const std::string& _msg= "unspecified error")
		  :msg(_msg)
	 {}
  };

public:

  assumption_manager_t(std::istream& _in,
							  const pedigree_t& _ped,
							  pedcnf_t& _cnf)
		:in(_in), ped(_ped), cnf(_cnf)
  {
  };

  void add_assumptions() {
	 L_INFO("Adding assumptions...");
	 size_t n_assumptions= 0;
	 std::string buff;
	 std::string prev_family_id= "";
	 while (std::getline(in, buff)) {
		L_TRACE("Read line starting with >" << buff.substr(0, 30) << "<");
		boost::trim(buff);
		if ((buff.length()==0) || boost::starts_with(buff, "#")) {
		  L_TRACE("Comment or empty line ignored.");
		} else {
// Process line
		  try {
			 std::istringstream is(buff);
			 ped_var_kind vk= ped_var_kind::DUMMY;
			 if (!(is >> vk)) throw invalid_line_t("unrecognized variable kind");
			 size_t ii, ij, ik;
			 if (!(is >> ii)) throw invalid_line_t("unrecognized variable first index");
			 if (!(is >> ij)) throw invalid_line_t("unrecognized variable second index");
			 if (vk==ped_var_kind::PM || vk==ped_var_kind::MM) {
				if (!(is >> ik)) throw invalid_line_t("unrecognized variable third index");
			 } else {
				ik= 0;
			 }
			 pedcnf_t::pedvar_t v= boost::make_tuple(vk, ii, ij, ik);
			 bool value;
			 if (!(is >> value)) throw invalid_line_t("unrecognized boolean value");
			 L_DEBUG("Read assumption '" << buff << "' ==> ("
						<< v <<  " == " << (value?"True":"False") << ").");
			 v.get<1>()= ped.get_by_id(v.get<1>()).progr_id();
			 L_DEBUG("...and transformed into " << v <<  " == " << (value?"True":"False") << ".");
			 if (v.get<0>() == ped_var_kind::SP) {
				if (!cnf.has_sp(v.get<1>(), v.get<2>())) {
				  L_WARN("Variable '" << v << "' does not exist in the SAT instance."
							" Adding anyway...");
				}
				lit_t iv= cnf.get_sp(v.get<1>(), v.get<2>());
				cnf.add_clause(pedcnf_t::clause_t{ (value? +1 : -1)  *  iv } );
				++n_assumptions;
			 } else if (v.get<0>() == ped_var_kind::SM) {
				if (!cnf.has_sm(v.get<1>(), v.get<2>())) {
				  L_WARN("Variable '" << v << "' does not exist in the SAT instance."
							" Adding anyway...");
				}
				lit_t iv= cnf.get_sm(v.get<1>(), v.get<2>());
				cnf.add_clause(pedcnf_t::clause_t{ (value? +1 : -1)  *  iv } );
				++n_assumptions;
			 } else if (v.get<0>() == ped_var_kind::P) {
				if (!cnf.has_p(v.get<1>(), v.get<2>())) {
				  L_WARN("Variable '" << v << "' does not exist in the SAT instance."
							" Adding anyway...");
				}
				lit_t iv= cnf.get_p(v.get<1>(), v.get<2>());
				cnf.add_clause(pedcnf_t::clause_t{ (value? +1 : -1)  *  iv } );
				++n_assumptions;
			 } else if (v.get<0>() == ped_var_kind::M) {
				if (!cnf.has_m(v.get<1>(), v.get<2>())) {
				  L_WARN("Variable '" << v << "' does not exist in the SAT instance."
							" Adding anyway...");
				}
				lit_t iv= cnf.get_m(v.get<1>(), v.get<2>());
				cnf.add_clause(pedcnf_t::clause_t{ (value? +1 : -1)  *  iv } );
				++n_assumptions;
			 } else if (v.get<0>() == ped_var_kind::RP) {
				if (!cnf.has_rp(v.get<1>(), v.get<2>())) {
				  L_WARN("Variable '" << v << "' does not exist in the SAT instance."
							" Adding anyway...");
				}
				lit_t iv= cnf.get_rp(v.get<1>(), v.get<2>());
				cnf.add_clause(pedcnf_t::clause_t{ (value? +1 : -1)  *  iv } );
				++n_assumptions;
			 } else if (v.get<0>() == ped_var_kind::RM) {
				if (!cnf.has_rm(v.get<1>(), v.get<2>())) {
				  L_WARN("Variable '" << v << "' does not exist in the SAT instance."
							" Adding anyway...");
				}
				lit_t iv= cnf.get_rm(v.get<1>(), v.get<2>());
				cnf.add_clause(pedcnf_t::clause_t{ (value? +1 : -1)  *  iv } );
				++n_assumptions;
			 } else if (v.get<0>() == ped_var_kind::PM) {
				if (!cnf.has_pm(v.get<1>(), v.get<2>(), v.get<3>())) {
				  L_WARN("Variable '" << v << "' does not exist in the SAT instance."
							" Adding anyway...");
				}
				lit_t iv= cnf.get_pm(v.get<1>(), v.get<2>(), v.get<3>());
				cnf.add_clause(pedcnf_t::clause_t{ (value? +1 : -1)  *  iv } );
				++n_assumptions;
			 } else if (v.get<0>() == ped_var_kind::MM) {
				if (!cnf.has_mm(v.get<1>(), v.get<2>(), v.get<3>())) {
				  L_WARN("Variable '" << v << "' does not exist in the SAT instance."
							" Adding anyway...");
				}
				lit_t iv= cnf.get_mm(v.get<1>(), v.get<2>(), v.get<3>());
				cnf.add_clause(pedcnf_t::clause_t{ (value? +1 : -1)  *  iv } );
				++n_assumptions;
			 } else if (v.get<0>() == ped_var_kind::E) {
				if (!cnf.has_e(v.get<1>(), v.get<2>())) {
				  L_WARN("Variable '" << v << "' does not exist in the SAT instance."
							" Adding anyway...");
				}
				lit_t iv= cnf.get_e(v.get<1>(), v.get<2>());
				cnf.add_clause(pedcnf_t::clause_t{ (value? +1 : -1)  *  iv } );
				++n_assumptions;
			 } else if (v.get<0>() == ped_var_kind::DUMMY) {
				L_FATAL("Error while processing assumption '" << buff << "' ==> ("
						  << v <<  " == " << (value?"True":"False") << ").");
				L_FATAL("Assumptions over dummy variables are NOT allowed. Exiting...");
				MY_FAIL;
			 } else {
				L_FATAL("Unrecognized assumption '" << buff << "'. Exiting...");
				MY_FAIL;
			 };
		  } catch (invalid_line_t& e) {
			 L_FATAL("!! Discarded invalid assumption >" << buff << "<. Reason: " << e.msg << ".");
			 MY_FAIL;
		  }
		}
	 }
	 L_INFO("Added " << n_assumptions << " assumptions.");
  };
};

void
add_assumptions(std::istream& is,
					 const plink_reader_t<>::multifamily_pedigree_t::pedigree_t& ped,
					 pedcnf_t& cnf) {
  assumption_manager_t mgr(is, ped, cnf);
  mgr.add_assumptions();
};
