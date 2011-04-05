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
#include "zrhc_app.hpp"

#include "configuration.h"

#include "application.hpp"
#include "log.hpp"
#include "assertion.hpp"
#include "utility.hpp"

#include <iostream>
#include <fstream>

#include <boost/static_assert.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>

#ifndef EXIT_NO_ZRHC
#define EXIT_NO_ZRHC (2)
#endif

BOOST_STATIC_ASSERT(EXIT_FAILURE != EXIT_NO_ZRHC);

using namespace std;


class zrhcstar_application_t: public application_t {

public:

  zrhcstar_application_t()
		:application_t(APPLICATION_CODENAME " " APPLICATION_VERSION_STRING
							" " BOOST_PP_STRINGIZE(SAT_SOLVER))
  {}

private:


protected:

  virtual po::options_description
  get_named_options() const {
	 po::options_description desc("Program options");
	 desc.add_options()
#ifndef ONLY_INTERNAL_SAT_SOLVER
		("create,1", po::bool_switch(),
		 "Create the SAT instance from the pedigree file.")
		("read,2", po::bool_switch(),
		 "Read the results produced by the SAT solver.")
		("create-read,3", po::bool_switch(),
		 "Create the SAT instance from the pedigree file, execute the SAT solver, and read the results.")
#endif // ONLY_INTERNAL_SAT_SOLVER
#ifdef INTERNAL_SAT_SOLVER
		("solve-internal,4", po::bool_switch(),
		 "Execute the integrated SAT solver.")
#endif // INTERNAL_SAT_SOLVER
		("pedigree,p",
		 po::value< std::string >()->default_value("pedigree.ped"),
		 "File storing the genotyped pedigree.")
		("sat,s",
		 po::value< std::string >()->default_value("instance.cnf"),
		 "File storing the SAT instance.")
		("result,r",
		 po::value< std::string >()->default_value("sat-result.txt"),
		 "File storing the results produced by the SAT solver.")
		("haplotypes,h",
		 po::value< std::string >()->default_value("haplotypes.txt"),
		 "File storing the computed haplotype configuration.")
		("sat-cmdline,c",
		 po::value< std::string >(),
		 "The command line used to execute the SAT solver: "
		 "%%INPUT%% and %%OUTPUT%% are markers used to represent the input and the output file "
		 "of the solver and they are automatically substituted by the program into the "
		 "corresponding filenames.")
		("compress,z", po::bool_switch()->default_value(false),
		 "Use compressed input and output files.")
		("compress-input", po::bool_switch()->default_value(false),
		 "Use compressed input files.")
		("compress-output", po::bool_switch()->default_value(false),
		 "Use compressed output files.")
		("keep,k", po::bool_switch()->default_value(false),
		 "Keep temporary files (such as 'cnf-instance-*' and 'res-cnf-instance-*' "
		 "files for '--create-read'/'-3' mode) after the execution.")
		("global-error,g", po::bool_switch()->default_value(false),
		 "Enable GLOBAL error handling (i.e., the error rate in each genotype is less or "
		 "equal than the specified error rate, computed over genotyped loci).")
		("global-error-rate,e", po::value< double >()->default_value(0.05),
		 "Maximum error rate in each genotype, computed only over genotyped loci "
		 "(used only if '--global-error' is specified).")
		("uniform-error,u", po::bool_switch()->default_value(false),
		 "Enable UNIFORM error handling (i.e., the number of errors in each window is less "
		 "or equal than the maximum number of errors).")
		("max-errors-in-window,m", po::value< unsigned int >()->default_value(4),
		 "Maximum number of errors in each window "
		 "(used only if '--uniform-error' is specified).\n"
		 "*MUST* be less or equal than half window size.")
		("window-length,w", po::value< unsigned int >()->default_value(16),
		 "number of typed loci that compose a window "
		 "(used only if '--uniform-error' is specified).\n"
		 "*MUST* be a power of 2 and *MUST* be greater than 2.\n"
		 "Windows overlap each other by half their length.")
		;
	 return desc;
  };

  virtual bool check_options(const po::variables_map& vm) {
#if defined(NO_INTERNAL_SAT_SOLVER)
	 mode_options(vm, "create", "read", "create-read");
#endif
#if defined(INTERNAL_SAT_SOLVER) && !defined(ONLY_INTERNAL_SAT_SOLVER)
	 mode_options(vm, "create", "read", "create-read", "solve-internal");
#endif
#if defined(ONLY_INTERNAL_SAT_SOLVER)
	 mode_options(vm, "solve-internal");
#endif
#ifndef ONLY_INTERNAL_SAT_SOLVER
	 conflicting_options(vm, "create", "read");
	 conflicting_options(vm, "create", "create-read");
	 conflicting_options(vm, "create-read", "read");
	 option_dependency(vm, "create", "pedigree");
	 option_dependency(vm, "create", "sat");
	 option_dependency(vm, "read", "pedigree");
	 option_dependency(vm, "read", "result");
	 option_dependency(vm, "read", "haplotypes");
	 option_dependency(vm, "create-read", "pedigree");
	 option_dependency(vm, "create-read", "haplotypes");
	 option_dependency(vm, "create-read", "sat-cmdline");
#endif
#ifdef INTERNAL_SAT_SOLVER
#ifndef ONLY_INTERNAL_SAT_SOLVER
	 conflicting_options(vm, "solve-internal", "create");
	 conflicting_options(vm, "solve-internal", "read");
	 conflicting_options(vm, "solve-internal", "create-read");
#endif
	 option_dependency(vm, "solve-internal", "pedigree");
	 option_dependency(vm, "solve-internal", "haplotypes");
#endif
	 option_dependency(vm, "global-error-rate", "global-error");
	 option_dependency(vm, "max-errors-in-window", "uniform-error");
	 option_dependency(vm, "window-length", "uniform-error");
	 if (vm["uniform-error"].as<bool>()) {
		const unsigned int wlen= vm["window-length"].as<unsigned int>();
		const unsigned int merr= vm["max-errors-in-window"].as<unsigned int>();
		if (wlen != pow2_of_floor_log2(wlen)) {
		  throw std::logic_error(std::string("The window length must be a power of 2."));
		}
		if (wlen < 4) {
		  throw std::logic_error(std::string("The window length must be greater than 2."));
		}
		if (merr > (wlen>>1)) {
		  throw std::logic_error(std::string("The maximum number of errors in a single window must be less or equal then the half the window length."));
		}
	 }
	 return true;
  }

  virtual int execution(int argc, char** argv,
								const po::variables_map& vm) {

	 int main_ris= EXIT_SUCCESS;


	 const bool in_compress=
		vm["compress"].as<bool>() ||
		vm["compress-input"].as<bool>();
	 const bool out_compress=
		vm["compress"].as<bool>() ||
		vm["compress-output"].as<bool>();

	 zrhcstar_t zrhcstar;
	 zrhcstar.prepare_program_options(vm);

// Dispatch the work depending on the program parameters
#ifndef ONLY_INTERNAL_SAT_SOLVER
	 if (vm["create"].as<bool>()) {
//    Creation of the SAT instance
		INFO("Creation of the SAT instance from the pedigree of file '"
			  << vm["pedigree"].as<string>() << "'...");

		file_utility::pistream ped_is=
		  file_utility::get_file_utility().
		  get_ifstream(vm["pedigree"].as<string>(), in_compress);
		file_utility::postream sat_os=
		  file_utility::get_file_utility().
		  get_ofstream(vm["sat"].as<string>(), out_compress);
		const std::string headers[] = {
		  "SAT instance",
		  std::string("pedigree: ") + vm["pedigree"].as<string>(),
		  std::string("sat: ") + vm["sat"].as<string>(),
		  std::string("source version: ") + APPLICATION_SOURCE_VERSION
		};
		zrhcstar.create_SAT_instance_from_pedigree(*ped_is,
																 *sat_os,
																 vector<string>(headers,
																					 headers+4));

		INFO("SAT instance successfully created and saved.");

		main_ris= EXIT_SUCCESS;

	 } else if (vm["read"].as<bool>()) {
//    Computation of the haplotype configuration
		INFO("Computation of the haplotype configuration from the "
			  "pedigree of file '"
			  << vm["pedigree"].as<string>()
			  << "' and the results of the SAT solver stored in file '"
			  << vm["result"].as<string>() << "'...");
		file_utility::pistream ped_is=
		  file_utility::get_file_utility().
		  get_ifstream(vm["pedigree"].as<string>(), in_compress);
		file_utility::pistream res_is=
		  file_utility::get_file_utility().
		  get_ifstream(vm["result"].as<string>(), false); // The SAT result is always not compressed
		file_utility::postream hap_os=
		  file_utility::get_file_utility().
		  get_ofstream(vm["haplotypes"].as<string>(), out_compress);
		bool is_zrhc=
		  zrhcstar.compute_HC_from_SAT_results(*ped_is, *res_is, *hap_os);

		if (is_zrhc) {
		  INFO("Zero-Recombinant Haplotype Configuration successfully "
				 "computed and saved.");
		  main_ris= EXIT_SUCCESS;
		} else {
		  WARN("No Zero-Recombinant Haplotype Configuration can exist. "
				 "Exiting without haplotype configuration.");
		  main_ris= EXIT_NO_ZRHC;
		}

	 } else if (vm["create-read"].as<bool>()) {
		INFO("Computation of the haplotype configuration from the "
			  "pedigree of file '"
			  << vm["pedigree"].as<string>() << "' by direct invocation of the SAT solver...");
		string sat_name;
		zrhcstar_t::pedigree_t ped;
		pedcnf_t cnf;
// Block for reading the pedigree and writing the SAT instance
// The block is needed to close the SAT instance stream before executing
// the solver
		{
		  file_utility::pistream ped_is=
			 file_utility::get_file_utility().
			 get_ifstream(vm["pedigree"].as<string>(), in_compress);
		  file_utility::postream sat_os=
			 file_utility::get_file_utility().
			 get_tmp_ostream("cnf-instance-XXXXXX", sat_name, out_compress);
		  const std::string headers[] = {
			 "SAT instance",
			 std::string("pedigree: ") + vm["pedigree"].as<string>(),
			 std::string("sat: ") + sat_name,
			 std::string("source version: ") + APPLICATION_SOURCE_VERSION
		  };
		  zrhcstar.create_SAT_instance_from_pedigree(*ped_is, *sat_os,
																	vector<string>(headers,
																						headers+4),
																	ped, cnf);
		}

// Execute the SAT solver
		INFO("Execution of the SAT solver...");
		INFO("Given command line: '"<< vm["sat-cmdline"].as<string>() << "'");
		const string res_name= "res-" + sat_name;
		string cmdline= vm["sat-cmdline"].as<string>();
		boost::replace_all(cmdline, "%%INPUT%%", "\"" + sat_name + "\"");
		boost::replace_all(cmdline, "%%OUTPUT%%", "\"" + res_name + "\"");
		INFO("Real command line: '"<< cmdline <<"'");

		boost::filesystem::remove(res_name);
		const int ret_value= system(cmdline.c_str());
// We cannot trust the return value
		DEBUG("The SAT solver returned: '" << ret_value << "'.");

// Check if the result file exists
		if (!boost::filesystem::exists(res_name)) {
		  WARN("Impossible to find SAT result file '" << res_name << "'. "
				 "Assuming that the SAT instance is not satisfiable.");
		  main_ris= EXIT_NO_ZRHC;
		} else {
// Read the results and compute the haplotype configuration
		  file_utility::pistream res_is=
			 file_utility::get_file_utility().
			 get_ifstream(res_name, false); // The SAT result is always not compressed
		  file_utility::postream hap_os=
			 file_utility::get_file_utility().
			 get_ofstream(vm["haplotypes"].as<string>(), out_compress);
		  bool is_zrhc= zrhcstar.compute_HC_from_SAT_results(ped, cnf,
																			  *res_is, *hap_os);

		  if (!vm["keep"].as<bool>()) {
			 DEBUG("Removing temporary files...");
			 bool remove_res= boost::filesystem::remove(sat_name);
			 remove_res= boost::filesystem::remove(res_name) && remove_res;
			 if (!remove_res) {
				WARN("File '" << sat_name << "' or '" << res_name <<
					  "' (or both) have NOT been removed successfully. Continuing anyway...");
			 } else {
				DEBUG("Temporary files removed.");
			 }
		  }

		  if (is_zrhc) {
			 INFO("Zero-Recombinant Haplotype Configuration successfully "
					"computed and saved.");
			 main_ris= EXIT_SUCCESS;
		  } else {
			 WARN("No Zero-Recombinant Haplotype Configuration can exist. "
					"Exiting without haplotype configuration.");
			 main_ris= EXIT_NO_ZRHC;
		  }
		}

	 } else
#endif
#ifdef INTERNAL_SAT_SOLVER
		if (vm["solve-internal"].as<bool>()) {
		INFO("Computation of the haplotype configuration from the "
			  "pedigree of file '"
			  << vm["pedigree"].as<string>() << "' by using the internal SAT solver...");
		zrhcstar_t::pedigree_t ped;
		pedcnf_t cnf;
// Block for reading the pedigree and writing the SAT instance
// The block is needed to close the SAT instance stream before executing
// the solver
		{
		  file_utility::pistream ped_is=
			 file_utility::get_file_utility().
			 get_ifstream(vm["pedigree"].as<string>(), in_compress);
		  zrhcstar.prepare_pedigree_and_sat(*ped_is,
														ped, cnf);
		}

// Execute the SAT solver
		INFO("Execution of the internal SAT solver...");
		const bool ret_value= cnf.solve();
// We have to trust the return value
		DEBUG("The SAT solver returned: '" << ret_value << "'.");

		if (ret_value) {
		  file_utility::postream hap_os=
			 file_utility::get_file_utility().
			 get_ofstream(vm["haplotypes"].as<string>(), out_compress);
		  bool is_zrhc= zrhcstar.compute_HC_from_model_and_save(ped, cnf,
																				  *hap_os);
		  if (is_zrhc) {
			 INFO("Zero-Recombinant Haplotype Configuration successfully "
					"computed and saved.");
			 main_ris= EXIT_SUCCESS;
		  } else {
			 WARN("A Haplotype Configuration has been computed but it is not valid!!");
			 main_ris= EXIT_NO_ZRHC;
		  }
		} else {
		  WARN("No Zero-Recombinant Haplotype Configuration can exist. "
				 "Exiting without haplotype configuration.");
		  main_ris= EXIT_NO_ZRHC;
		}

	 } else
#endif // INTERNAL_SAT_SOLVER
	 {
// We should not arrive here
		MY_FAIL;
		throw logic_error(string("Modes have not been recognized."));
	 }

	 return main_ris;
  }

};



int main(int argc, char** argv) {
  zrhcstar_application_t app;
  return app.execute(argc, argv);
}
