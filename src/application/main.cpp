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
#include "rehc_app.hpp"

#include "configuration.h"

#include "application.hpp"
#include "log.hpp"
#include "assertion.hpp"
#include "utility.hpp"
#include "timelimit.hpp"

#include <iostream>
#include <fstream>

#include <boost/static_assert.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>

BOOST_STATIC_ASSERT(EXIT_FAILURE != EXIT_NO_reHC);
BOOST_STATIC_ASSERT(EXIT_FAILURE != EXIT_reHC_ERROR);

using namespace std;


class rehcstar_application_t: public application_t {

public:

  rehcstar_application_t()
		:application_t(APPLICATION_CODENAME " " APPLICATION_VERSION_STRING
							" " BOOST_PP_STRINGIZE(SAT_SOLVER))
  {}

private:


protected:

  virtual po::options_description
  get_named_options() const {
	 po::options_description modes("Program Modes",
											 po::options_description::m_default_line_length+12,
											 po::options_description::m_default_line_length-12);
#ifndef ONLY_INTERNAL_SAT_SOLVER
	 modes.add_options()
		("create,1", po::bool_switch(),
		 "Create the SAT instance from the pedigree file.")
		("read,2", po::bool_switch(),
		 "Read the results produced by the SAT solver.")
		("create-read,3", po::bool_switch(),
		 "Create the SAT instance from the pedigree file, execute the SAT solver, and read the results.");
#endif // ONLY_INTERNAL_SAT_SOLVER
#ifdef INTERNAL_SAT_SOLVER
	 modes.add_options()
		("solve-internal,4", po::bool_switch(),
		 "Execute the integrated SAT solver.");
#endif // INTERNAL_SAT_SOLVER
	 po::options_description files("Input/Output",
											 po::options_description::m_default_line_length+12,
											 po::options_description::m_default_line_length-12);
	 files.add_options()
		("pedigree,p",
		 po::value< std::string >()->default_value("pedigree.ped"),
		 "File storing the genotyped pedigree.");
#ifndef ONLY_INTERNAL_SAT_SOLVER
	 files.add_options()
		("sat,s",
		 po::value< std::string >()->default_value("instance.cnf"),
		 "File storing the SAT instance.")
		("result,r",
		 po::value< std::string >()->default_value("sat-result.txt"),
		 "File storing the results produced by the SAT solver.");
#endif // ONLY_INTERNAL_SAT_SOLVER
	 files.add_options()
		("haplotypes,h",
		 po::value< std::string >()->default_value("haplotypes.txt"),
		 "File storing the computed haplotype configuration.");
	 files.add_options()
		("assumptions,a",
		 po::value< std::string >(),
		 "File storing additional assumptions/constraints that must hold in the "
		 "reconstructed haplotype configuration.\n"
		 "Each assumption is in a single row composed by 4 white-spaced fields:\n"
		 "\t<kind of variable> <individual index> <locus index> <bool value (0/1)>");
#ifndef ONLY_INTERNAL_SAT_SOLVER
	 files.add_options()
		("sat-cmdline,c",
		 po::value< std::string >(),
		 "The command line used to execute the SAT solver: "
		 "%%INPUT%% and %%OUTPUT%% are markers used to represent the input and the output file "
		 "of the solver and they are automatically substituted by the program into the "
		 "corresponding filenames.");
#endif // ONLY_INTERNAL_SAT_SOLVER
	 files.add_options()
		("compress,z", po::bool_switch()->default_value(false),
		 "Use compressed input and output files.")
		("compress-input", po::bool_switch()->default_value(false),
		 "Use compressed input files.")
		("compress-output", po::bool_switch()->default_value(false),
		 "Use compressed output files.");
#ifndef ONLY_INTERNAL_SAT_SOLVER
	 files.add_options()
		("compress-sat", po::bool_switch()->default_value(false),
		 "Compress the file containing the SAT instance.")
		("keep,k", po::bool_switch()->default_value(false),
		 "Keep temporary files (such as 'cnf-instance-*' and 'res-cnf-instance-*' "
		 "files for '--create-read'/'-3' mode) after the execution.")
		("pipe", po::bool_switch()->default_value(false),
		 "Pipe the SAT instance to the external solver instead of using an intermediate file.");
#endif // ONLY_INTERNAL_SAT_SOLVER
	 po::options_description errors("Error Management Options",
											  po::options_description::m_default_line_length+12,
											  po::options_description::m_default_line_length-12);
	 errors.add_options()
		("global-error", po::bool_switch()->default_value(false),
		 "Enable GLOBAL error handling (i.e., the global error rate in the whole pedigree is "
		 "less than or equal to the specified error rate, computed over genotyped loci).")
		("global-error-rate", po::value< double >()->default_value(0.005),
		 "Maximum error rate in all the genotypes, computed only over genotyped loci "
		 "(used only if '--global-error' is specified, cannot be used with '--global-error-number').")
		("global-error-number", po::value< unsigned int >()->default_value(1),
		 "Maximum number of errors in the genotyped loci"
		 "(used only if '--global-error' is specified, cannot be used with '--global-error-rate').")
		("individual-error", po::bool_switch()->default_value(false),
		 "Enable INDIVIDUAL error handling (i.e., the error rate in each genotype is less than "
		 "or equal to the specified error rate, computed over genotyped loci).")
		("individual-error-rate", po::value< double >()->default_value(0.01),
		 "Maximum error rate in each genotype, computed only over genotyped loci "
		 "(used only if '--individual-error' is specified).")
		("uniform-error", po::bool_switch()->default_value(false),
		 "Enable UNIFORM error handling (i.e., the number of errors in each window is less than "
		 "or equal to the maximum number of errors).")
		("max-errors-in-window", po::value< unsigned int >()->default_value(4),
		 "Maximum number of errors in each window "
		 "(used only if '--uniform-error' is specified).\n"
		 "NOTE: It *must* be less than or equal to half window size.")
		("error-window-length", po::value< unsigned int >()->default_value(16),
		 "Number of typed loci that compose a window "
		 "(used only if '--uniform-error' is specified).\n"
		 "NOTE: It *must* be a power of 2 and *must* be greater than 2.\n"
		 "Windows overlap each other by half their length.");
	 po::options_description recombs("Recombination Management Options",
												po::options_description::m_default_line_length+12,
												po::options_description::m_default_line_length-12);
	 recombs.add_options()
		("global-recomb", po::bool_switch()->default_value(false),
		 "Enable GLOBAL recombination handling (i.e., the global recombination rate in the whole "
		 "pedigree is less than or equal to the specified recombination rate, computed over *ALL* loci).")
		("global-recomb-rate", po::value< double >()->default_value(0.01),
		 "Maximum recombination rate in all the genotypes "
		 "(used only if '--global-recomb' is specified, cannot be used with '--global-recomb-number').")
		("global-recomb-number", po::value< unsigned int >()->default_value(1),
		 "Maximum number of recombinations in all the genotypes "
		 "(used only if '--global-recomb' is specified, cannot be used with '--global-recomb-rate').")
		("global-recomb-min-number", po::value< unsigned int >()->default_value(0),
		 "Minimum number of recombinations in all the genotypes "
		 "(used only if '--global-recomb' is specified, cannot be used with '--global-recomb-rate').\n"
		 "NOTE: This lower bound is not strictly enforced, since the SAT solver could compute a "
		 "solution with unnecessary recombinations. This option should be used to improve the "
		 "efficiency of the search of the optimum and should be set to a value such that "
		 "a solution with this number of recombinations does not exist.")
		("individual-recomb", po::bool_switch()->default_value(false),
		 "Enable INDIVIDUAL recombination handling (i.e., the recombination rate in each genotype is "
		 "less than or equal to the specified recombination rate, computed over *ALL* loci).")
		("individual-recomb-rate", po::value< double >()->default_value(0.01),
		 "Maximum recombination rate in each genotype "
		 "(used only if '--individual-recomb' is specified).")
		("uniform-recomb", po::bool_switch()->default_value(false),
		 "Enable UNIFORM recombination handling (i.e., the number of recombinations in each window is "
		 "less than or equal to the specified maximum number of recombinations).")
		("max-recombs-in-window", po::value< unsigned int >()->default_value(4),
		 "Maximum number of recombinations in each window "
		 "(used only if '--uniform-recomb' is specified).\n"
		 "NOTE: It *must* be less than or equal to half window size.")
		("recomb-window-length", po::value< unsigned int >()->default_value(16),
		 "Number of loci that compose a window "
		 "(used only if '--uniform-recomb' is specified).\n"
		 "NOTE: It *must* be a power of 2 and *must* be greater than 2.\n"
		 "Windows overlap each other by half their length.")
		;

	 po::options_description exec_opt("Execution Management Options",
												 po::options_description::m_default_line_length+12,
												 po::options_description::m_default_line_length-12);
	 exec_opt.add_options()
		("time-limit", po::value< unsigned int >()->default_value(0),
		 "Maximum (approximated) execution time in seconds (0=no limit).")
		;

	 modes.add(files);
	 modes.add(errors);
	 modes.add(recombs);
	 modes.add(exec_opt);
	 return modes;
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
	 option_dependency(vm, "pipe", "create-read");
#endif
#ifdef INTERNAL_SAT_SOLVER
#ifndef ONLY_INTERNAL_SAT_SOLVER
	 conflicting_options(vm, "solve-internal", "create");
	 conflicting_options(vm, "solve-internal", "read");
	 conflicting_options(vm, "solve-internal", "create-read");
#endif
	 conflicting_options(vm, "compress-sat", "pipe");
	 option_dependency(vm, "solve-internal", "pedigree");
	 option_dependency(vm, "solve-internal", "haplotypes");
#endif
	 option_dependency(vm, "global-error-rate", "global-error");
	 option_dependency(vm, "global-error-number", "global-error");
	 option_dependency(vm, "individual-error-rate", "individual-error");
	 option_dependency(vm, "max-errors-in-window", "uniform-error");
	 option_dependency(vm, "error-window-length", "uniform-error");
	 option_dependency(vm, "global-recomb-rate", "global-recomb");
	 option_dependency(vm, "global-recomb-number", "global-recomb");
	 option_dependency(vm, "global-recomb-min-number", "global-recomb");
	 conflicting_options(vm, "global-recomb-rate", "global-recomb-number");
	 conflicting_options(vm, "global-recomb-rate", "global-recomb-min-number");
	 option_dependency(vm, "individual-recomb-rate", "individual-recomb");
	 option_dependency(vm, "max-recombs-in-window", "uniform-recomb");
	 option_dependency(vm, "recomb-window-length", "uniform-recomb");
	 if (vm["global-recomb"].as<bool>()) {
		const unsigned int rmin= vm["global-recomb-min-number"].as<unsigned int>();
		const unsigned int rmax= vm["global-recomb-number"].as<unsigned int>();
		if (rmin > rmax) {
		  throw std::logic_error(std::string("The minimum number of recombinations must be not greater than the maximum number of recombinations."));
		}
	 }
	 if (vm["uniform-error"].as<bool>()) {
		const unsigned int wlen= vm["error-window-length"].as<unsigned int>();
		const unsigned int merr= vm["max-errors-in-window"].as<unsigned int>();
		if (wlen != pow2_of_floor_log2(wlen)) {
		  throw std::logic_error(std::string("The window length must be a power of 2."));
		}
		if (wlen < 4) {
		  throw std::logic_error(std::string("The window length must be greater than 2."));
		}
		if (merr >= wlen) {
		  throw std::logic_error(std::string("The maximum number of errors in a single window must be less than the window length."));
		}
	 }
	 if (vm["uniform-recomb"].as<bool>()) {
		const unsigned int wlen= vm["recomb-window-length"].as<unsigned int>();
		const unsigned int merr= vm["max-recombs-in-window"].as<unsigned int>();
		if (wlen != pow2_of_floor_log2(wlen)) {
		  throw std::logic_error(std::string("The window length must be a power of 2."));
		}
		if (wlen < 4) {
		  throw std::logic_error(std::string("The window length must be greater than 2."));
		}
		if (merr >= wlen) {
		  throw std::logic_error(std::string("The maximum number of recombinations in a single window must be less than the window length."));
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
#ifndef ONLY_INTERNAL_SAT_SOLVER
	 const bool sat_compress=
		out_compress ||
		vm["compress-sat"].as<bool>();
#endif

	 rehcstar_t rehcstar;
	 rehcstar.prepare_program_options(vm);

// Check if a time-limit is specified
	 if (vm["time-limit"].as<unsigned int>()>0) {
		register_time_limit(vm["time-limit"].as<unsigned int>());
	 }

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
		  get_ofstream(vm["sat"].as<string>(), sat_compress);
		const std::string headers[] = {
		  "SAT instance",
		  std::string("pedigree: ") + vm["pedigree"].as<string>(),
		  std::string("sat: ") + vm["sat"].as<string>(),
		  std::string("source version: ") + APPLICATION_SOURCE_VERSION
		};
		rehcstar.create_SAT_instance_from_pedigree(*ped_is,
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
		boost::tribool is_rehc=
		  rehcstar.compute_HC_from_SAT_results(*ped_is, *res_is, *hap_os);

		if (is_rehc) {
		  INFO("(r,e)-Haplotype Configuration successfully "
				 "computed and saved.");
		  main_ris= EXIT_SUCCESS;
		} else if (!is_rehc) {
		  WARN("No (r,e)-Haplotype Configuration can exist. "
				 "Exiting without haplotype configuration.");
		  main_ris= EXIT_NO_reHC;
		} else {
		  WARN("We do NOT know if a (r,e)-Haplotype Configuration can exist. "
				 "The SAT solver did not give a valid result.");
		  main_ris= EXIT_reHC_ERROR;
		}

	 } else if (vm["create-read"].as<bool>()) {
		INFO("Computation of the haplotype configuration from the "
			  "pedigree of file '"
			  << vm["pedigree"].as<string>() << "' by direct invocation of the SAT solver...");
		string sat_name;
		string res_name= ".";
		rehcstar_t::pedigree_t ped;
		pedcnf_t cnf;
		int ret_value;
//
// ** USE "NORMAL" FILES ** BEGIN
		if (! vm["pipe"].as<bool>()) {
// Block for reading the pedigree and writing the SAT instance
// The block is needed to close the SAT instance stream before executing
// the solver
		  {
			 file_utility::pistream ped_is=
				file_utility::get_file_utility().
				get_ifstream(vm["pedigree"].as<string>(), in_compress);
			 file_utility::postream sat_os=
				file_utility::get_file_utility().
				get_tmp_ostream("cnf-instance-XXXXXX", sat_name, sat_compress);
			 const std::string headers[] = {
				"SAT instance",
				std::string("pedigree: ") + vm["pedigree"].as<string>(),
				std::string("sat: ") + sat_name,
				std::string("source version: ") + APPLICATION_SOURCE_VERSION
			 };
			 rehcstar.create_SAT_instance_from_pedigree(*ped_is, *sat_os,
																	  vector<string>(headers,
																						  headers+4),
																	  ped, cnf);
		  }

// Execute the SAT solver
		  INFO("Execution of the SAT solver...");
		  INFO("Given command line: '"<< vm["sat-cmdline"].as<string>() << "'");
		  res_name= "res-" + sat_name;
		  string cmdline= vm["sat-cmdline"].as<string>();
		  boost::replace_all(cmdline, "%%INPUT%%", "\"" + sat_name + "\"");
		  boost::replace_all(cmdline, "%%OUTPUT%%", "\"" + res_name + "\"");
		  INFO("Real command line: '"<< cmdline <<"'");

		  boost::filesystem::remove(res_name);
		  ret_value= system(cmdline.c_str());
// ** USE "NORMAL" FILES ** END
//
		} else {
//
// ** USE "PIPE" ** BEGIN
		  { // used only to get a temporary filename
			 file_utility::postream sat_os=
				file_utility::get_file_utility().
				get_tmp_ostream("cnf-instance-XXXXXX", sat_name, sat_compress);
		  }
		  boost::filesystem::remove(sat_name);
		  INFO("Starting the SAT solver (piping the instance)...");
		  INFO("Given command line: '"<< vm["sat-cmdline"].as<string>() << "'");
		  res_name= "res-" + sat_name;
		  string cmdline= vm["sat-cmdline"].as<string>();
		  boost::replace_all(cmdline, "%%INPUT%%", "");
		  boost::replace_all(cmdline, "%%OUTPUT%%", "\"" + res_name + "\"");
		  INFO("Real command line: '"<< cmdline <<"'");

		  boost::filesystem::remove(res_name);
		  FILE* sat_pipe= popen(cmdline.c_str(), "w");

		  if (sat_pipe == NULL) {
			 WARN("Impossible to execute the command (or another error has occurred).");
			 ret_value= -1;
		  } else {
			 INFO("Creating the SAT instance and writing to SAT solver's standard input...");
			 file_utility::pistream ped_is=
				file_utility::get_file_utility().
				get_ifstream(vm["pedigree"].as<string>(), in_compress);
			 file_utility::postream sat_os=
				file_utility::get_file_utility().
				get_ostream_from_file(sat_pipe);
			 rehcstar.create_SAT_instance_from_pedigree(*ped_is, *sat_os,
																	  vector<string>(),
																	  ped, cnf);
			 INFO("Finished creating and writing the SAT instance. Waiting for the solver...");
			 ret_value= pclose(sat_pipe);
		  }
// ** USE "PIPE" ** END
//
		}
// We cannot trust the return value
		DEBUG("The SAT solver returned: '" << ret_value << "'.");

// Check if the result file exists
		if (!boost::filesystem::exists(res_name)) {
		  WARN("Impossible to find SAT result file '" << res_name << "'. "
				 "Assuming that the SAT instance is not satisfiable.");
		  main_ris= EXIT_NO_reHC;
		} else {
// Read the results and compute the haplotype configuration
		  file_utility::pistream res_is=
			 file_utility::get_file_utility().
			 get_ifstream(res_name, false); // The SAT result is always not compressed
		  file_utility::postream hap_os=
			 file_utility::get_file_utility().
			 get_ofstream(vm["haplotypes"].as<string>(), out_compress);
		  boost::tribool
			 is_rehc= rehcstar.compute_HC_from_SAT_results(ped, cnf,
																		  *res_is, *hap_os);

		  if (!vm["keep"].as<bool>()) {
			 DEBUG("Removing temporary files...");
			 bool remove_res= true;
			 if (!vm["pipe"].as<bool>()) {
				remove_res= boost::filesystem::remove(sat_name);
			 }
			 remove_res= boost::filesystem::remove(res_name) && remove_res;
			 if (!remove_res) {
				WARN("File '" << sat_name << "' or '" << res_name <<
					  "' (or both) have NOT been removed successfully. Continuing anyway...");
			 } else {
				DEBUG("Temporary files removed.");
			 }
		  }

		  if (is_rehc) {
			 INFO("(r,e)-Haplotype Configuration successfully "
					"computed and saved.");
			 main_ris= EXIT_SUCCESS;
		  } else if (!is_rehc) {
			 WARN("No (r,e)-Haplotype Configuration can exist. "
					"Exiting without haplotype configuration.");
			 main_ris= EXIT_NO_reHC;
		  } else {
			 WARN("We do NOT know if a (r,e)-Haplotype Configuration can exist. "
					"The SAT solver did not give a valid result.");
			 main_ris= EXIT_reHC_ERROR;
		  }
		}

	 } else
#endif
#ifdef INTERNAL_SAT_SOLVER
		if (vm["solve-internal"].as<bool>()) {
		INFO("Computation of the haplotype configuration from the "
			  "pedigree of file '"
			  << vm["pedigree"].as<string>() << "' by using the internal SAT solver...");
		rehcstar_t::pedigree_t ped;
		pedcnf_t cnf;
// Block for reading the pedigree and writing the SAT instance
// The block is needed to close the SAT instance stream before executing
// the solver
		{
		  file_utility::pistream ped_is=
			 file_utility::get_file_utility().
			 get_ifstream(vm["pedigree"].as<string>(), in_compress);
		  rehcstar.prepare_pedigree_and_sat(*ped_is,
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
		  bool is_rehc= rehcstar.compute_HC_from_model_and_save(ped, cnf,
																				  *hap_os);
		  if (is_rehc) {
			 INFO("(r,e)-Haplotype Configuration successfully "
					"computed and saved.");
			 main_ris= EXIT_SUCCESS;
		  } else {
			 WARN("A Haplotype Configuration has been computed but it is not valid!!");
			 main_ris= EXIT_reHC_ERROR;
		  }
		} else {
		  WARN("No (r,e)-Haplotype Configuration can exist. "
				 "Exiting without haplotype configuration.");
		  main_ris= EXIT_NO_reHC;
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
  rehcstar_application_t app;
  return app.execute(argc, argv);
}
