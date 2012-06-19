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
#include "configuration.h"

#include "application.hpp"
#include "log.hpp"
#include "assertion.hpp"
#include "data.hpp"
#include "utility.hpp"

#include "io-pedigree.hpp"
#include "ped2cnf.hpp"
#include "ped2cnf-constraints.hpp"
#include "pedcnf2hc.hpp"

#include "assumptions.hpp"

#include <iostream>
#include <boost/program_options.hpp>

using namespace std;


class rehcstar_t:
  public log_able_t< rehcstar_t >
{

public:

  typedef plink_reader_t<>::multifamily_pedigree_t pedigree_t;
  typedef pedigree_t::pedigree_t family_t;

private:

  composite_constraints_t err_constraints;
  composite_constraints_t global_recomb_constraints;
  composite_constraints_t separated_recomb_constraints;
  bool has_errors;
  bool has_global_recombinations;
  bool has_separated_recombinations;
  bool has_recombinations;

  bool has_assumptions;
  std::string assumption_file;

public:

  void prepare_program_options(const boost::program_options::variables_map& vm) {
// Analyze error-related program options
	 has_errors= false;
	 if (vm["global-error"].as<bool>()) {
// Use default 'global-error-rate' if 'global-error-number' nor 'global-error-rate'
// have been specified
		if (vm.count("global-error-number") && !vm["global-error-number"].defaulted()) {
		  const unsigned int error_number= vm["global-error-number"].as<unsigned int>();
		  L_INFO("Enabling *GLOBAL* error handling ("
					"error-number=" << error_number << ")");
		  err_constraints.add(new at_most_global_constraints_abs_t(error_number, "errors"));
		} else {
		  const double err_rate= vm["global-error-rate"].as<double>();
		  L_INFO("Enabling *GLOBAL* error handling ("
					"error-rate=" << err_rate << ")");
		  err_constraints.add(new at_most_global_constraints_t(err_rate, "errors"));
		}
		has_errors= true;
	 }
	 if (vm["individual-error"].as<bool>()) {
		const double err_rate= vm["individual-error-rate"].as<double>();
		L_INFO("Enabling *INDIVIDUAL* error handling ("
				 "error-rate=" << err_rate << ")");
		err_constraints.add(new at_most_individual_constraints_t(err_rate, "errors"));
		has_errors= true;
	 }
	 if (vm["uniform-error"].as<bool>()) {
		const unsigned int winerr= vm["max-errors-in-window"].as<unsigned int>();
		const unsigned int winlen= vm["error-window-length"].as<unsigned int>();
		L_INFO("Enabling *WINDOWED* error handling ("
				 "max-errors-in-windows=" << winerr << ", "
				 "error-window-length=" << winlen << ")");
		err_constraints.add(new at_most_windowed_constraints_t(winerr, winlen, "errors"));
		has_errors= true;
	 }
	 if (!has_errors) {
		L_INFO("*DISABLING* genotyping errors");
		err_constraints.add(new all_false_constraints_t("errors"));
	 }
// Analyze recombination-related program options
	 has_global_recombinations=
		has_separated_recombinations=
		has_recombinations= false;
	 if (vm["global-recomb"].as<bool>()) {
// Use default 'global-recomb-rate' if 'global-recomb-number' nor 'global-recomb-rate'
// have been specified
		if (vm.count("global-recomb-number") && !vm["global-recomb-number"].defaulted()) {
		  const unsigned int recomb_number= vm["global-recomb-number"].as<unsigned int>();
		  if (vm.count("global-recomb-min-number") && !vm["global-recomb-min-number"].defaulted()) {
			 const unsigned int recomb_min_number= vm["global-recomb-min-number"].as<unsigned int>();
			 L_INFO("Enabling *GLOBAL* recombination handling ("
					  "recomb-number=" << recomb_number << ", "
					  "recomb-min-number=" << recomb_min_number << ")");
			 global_recomb_constraints.add(new interval_global_constraints_abs_t(recomb_min_number, recomb_number, "recombinations"));
		  } else {
			 L_INFO("Enabling *GLOBAL* recombination handling ("
					  "recomb-number=" << recomb_number << ")");
			 global_recomb_constraints.add(new at_most_global_constraints_abs_t(recomb_number, "recombinations"));
		  }
		} else {
		  const double recomb_rate= vm["global-recomb-rate"].as<double>();
		  L_INFO("Enabling *GLOBAL* recombination handling ("
					"recomb-rate=" << recomb_rate << ")");
		  global_recomb_constraints.add(new at_most_global_constraints_t(recomb_rate, "recombinations"));
		}
		has_global_recombinations= true;
	 }
	 if (vm["individual-recomb"].as<bool>()) {
		const double recomb_rate= vm["individual-recomb-rate"].as<double>();
		L_INFO("Enabling *INDIVIDUAL* recombination handling ("
				 "recomb-rate=" << recomb_rate << ")");
		global_recomb_constraints.add(new at_most_individual_constraints_t(recomb_rate, "recombinations"));
		has_global_recombinations= true;
	 }
	 if (vm["uniform-recomb"].as<bool>()) {
		const unsigned int winrecomb= vm["max-recombs-in-window"].as<unsigned int>();
		const unsigned int winlen= vm["recomb-window-length"].as<unsigned int>();
		L_INFO("Enabling *WINDOWED* recombination handling ("
				 "max-recombs-in-windows=" << winrecomb << ", "
				 "recomb-window-length=" << winlen << ")");
		separated_recomb_constraints.add(new at_most_windowed_constraints_t(winrecomb, winlen, "recombinations"));
		has_separated_recombinations= true;
	 }
	 has_recombinations= has_global_recombinations || has_separated_recombinations;
	 if (!has_recombinations) {
		L_INFO("*DISABLING* recombinations");
		global_recomb_constraints.add(new all_false_constraints_t("recombinations"));
	 }

	 if (vm.count("assumptions")>0 && !vm["assumptions"].defaulted()) {
		has_assumptions= true;
		assumption_file= vm["assumptions"].as<std::string>();
	 } else {
		has_assumptions= false;
		assumption_file= ".";
	 }
  };

  void prepare_pedigree_and_sat(std::istream& ped_is,
										  pedigree_t& mped,
										  pedcnf_t& cnf) const {
	 L_INFO("Reading pedigree...");
	 multiallelic_genotype_reader_t gr;
	 plink_reader_t<> reader(gr);
	 reader.read(ped_is, mped);
	 L_INFO("Pedigree successfully read.");

	 if (mped.families().empty()) {
		throw std::logic_error(std::string("No family has been read."));
	 }
	 if (mped.families().size() > 1) {
		throw std::logic_error(std::string("The pedigree has more than one family."));
	 }

	 mped.print_stats();

// Prepare the SAT instance
	 L_INFO("Preparing SAT instance from pedigree...");
	 ped2cnf(mped.families().front(), cnf);
	 L_DEBUG("So far the SAT instance is composed by " <<
				std::setw(8) << cnf.vars().size() << " variables and " <<
				std::setw(8) << cnf.no_of_clauses() << " clauses");
	 L_INFO("Adding clauses for managing genotyping errors...");
	 error_handler_t err_handler(err_constraints);
	 err_handler.handle_errors(cnf, mped.families().front().size(),
										mped.families().front().genotype_length());
	 L_DEBUG("So far the SAT instance is composed by " <<
				std::setw(8) << cnf.vars().size() << " variables and " <<
				std::setw(8) << cnf.no_of_clauses() << " clauses");
	 L_INFO("Adding clauses for managing recombination events...");
	 global_recombination_handler_t global_recomb_handler(global_recomb_constraints);
	 global_recomb_handler.handle_recombinations(cnf, mped.families().front().size(),
																mped.families().front().genotype_length());
	 separated_recombination_handler_t separated_recomb_handler(separated_recomb_constraints);
	 separated_recomb_handler.handle_recombinations(cnf, mped.families().front().size(),
																	mped.families().front().genotype_length());
	 if (has_assumptions) {
		std::ifstream assumptions(assumption_file.c_str());
		if (assumptions) {
		  add_assumptions(assumptions, mped.families().front(), cnf);
		  assumptions.close();
		} else {
		  L_FATAL("Impossible to open assumption file '" << assumption_file << "'. Aborting...");
		  MY_FAIL;
		}
	 }

	 L_INFO("SAT instance successfully prepared.");
	 L_INFO("The SAT instance is composed by " <<
			  std::setw(8) << cnf.vars().size() << " variables and " <<
			  std::setw(8) << cnf.no_of_clauses() << " clauses");
  }



public:

  explicit rehcstar_t()
		:has_errors(false), has_recombinations(false)
  {
  };

  void save_reHC(pedigree_t& ped,
					  std::ostream& hap_os) const {
	 L_INFO("Saving haplotype configuration...");
// FIXME: Improve template instantiation
	 multiallelic_haplotype_pair_writer_t hpw;
	 plink_haplotype_writer_t<> writer(hpw, "\t", "|");
	 writer.write(hap_os, ped);
	 L_INFO("Haplotype configuration successfully saved.");
  };

  boost::tribool read_SAT_results(pedcnf_t& cnf,
								std::istream& res_is) const {
	 L_INFO("Reading SAT results...");
	 const boost::tribool is_sat= cnf.assignment_from_minisat_format(res_is);
	 L_INFO("SAT results successfully read.");
	 return is_sat;
  };


#ifndef ONLY_INTERNAL_SAT_SOLVER
  void create_SAT_instance_from_pedigree(std::istream& ped_is,
													  std::ostream& sat_os,
													  const std::vector<std::string>& headers) const {
	 pedigree_t ped;
	 pedcnf_t cnf;
	 create_SAT_instance_from_pedigree(ped_is, sat_os, headers,
												  ped, cnf);
  }

  void create_SAT_instance_from_pedigree(std::istream& ped_is,
													  std::ostream& sat_os,
													  const std::vector<std::string>& headers,
													  pedigree_t& ped,
													  pedcnf_t& cnf) const {
	 prepare_pedigree_and_sat(ped_is, ped, cnf);
// Output the instance
	 L_INFO("Saving SAT instance...");
	 cnf.clauses_to_dimacs_format(sat_os, headers);
	 L_INFO("SAT instance successfully saved.");
  };



  boost::tribool
  compute_HC_from_SAT_results(std::istream& ped_is,
										std::istream& res_is,
										pedigree_t& ped,
										pedcnf_t& cnf) const {
	 prepare_pedigree_and_sat(ped_is, ped, cnf);
	 return compute_HC_from_SAT_results(res_is, ped, cnf);
  };

  boost::tribool
  compute_HC_from_SAT_results(std::istream& res_is,
										pedigree_t& ped,
										pedcnf_t& cnf) const {
// Open the result file and read the assignment
	 const boost::tribool is_sat= read_SAT_results(cnf, res_is);
	 if (is_sat) {
		L_INFO("The pedigree can be realized by a (r,e)-haplotype "
				 "configuration.");
		const bool ok= compute_HC_from_model(ped, cnf);
		if (!ok) {
		  MY_FAIL;
		}
	 } else if (!is_sat) {
		L_INFO("The pedigree CANNOT be realized by a (r,e)-haplotype "
			  "configuration.");
// Do nothing
	 } else {
		L_WARN("We do NOT know if the pedigree can be realized by a "
				 "(r,e)-haplotype configuration. "
				 "The SAT solver did not give a valid result.")
	 }
	 return is_sat;
  }
#endif // not defined ONLY_INTERNAL_SAT_SOLVER

  bool compute_HC_from_model(pedigree_t& ped,
									  pedcnf_t& cnf) const {
	 pedigree_t::pedigree_t& family= ped.families().front();
// Compute the actual haplotype configuration
	 compute_reHC_from_SAT(family, cnf);
// Check the haplotype configuration
	 const bool ok=
		family.is_completely_haplotyped() &&
		family.is_consistent(false);
	 const int n_recomb= family.is_mendelian_consistent();
	 L_INFO("The computed haplotype configuration has " << n_recomb << " recombinations.");
	 if (ok) {
		L_INFO("The computed haplotype configuration is valid.");
	 } else {
		L_ERROR("The computed haplotype configuration is not valid.");
	 }
	 return ok;
  };

  bool compute_HC_from_model_and_save(pedigree_t& ped,
												  pedcnf_t& cnf,
												  std::ostream& hap_os) const {
	 const bool ok= compute_HC_from_model(ped, cnf);
	 if (ok) {
// Output the haplotype configuration
		save_reHC(ped, hap_os);
	 }
	 return ok;
  };


#ifndef ONLY_INTERNAL_SAT_SOLVER

  boost::tribool
  compute_HC_from_SAT_results(std::istream& ped_is,
										std::istream& res_is,
										std::ostream& hap_os) const {
	 pedigree_t ped;
	 pedcnf_t cnf;
	 const boost::tribool is_sat=
		compute_HC_from_SAT_results(ped_is, res_is,
											 ped, cnf);
	 if (is_sat) {
// Output the haplotype configuration
		save_reHC(ped, hap_os);
	 }
	 return is_sat;
  }

  boost::tribool
  compute_HC_from_SAT_results(pedigree_t& ped,
										pedcnf_t& cnf,
										std::istream& res_is,
										std::ostream& hap_os) const {
	 const boost::tribool is_sat=
		compute_HC_from_SAT_results(res_is,
											 ped, cnf);
	 if (is_sat) {
// Output the haplotype configuration
		save_reHC(ped, hap_os);
	 }
	 return is_sat;
  }

#endif // not defined ONLY_INTERNAL_SAT_SOLVER

};
