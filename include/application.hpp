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
#ifndef __APPLICATION_HPP__
#define __APPLICATION_HPP__

#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>

#include "log.hpp"
#include "assertion.hpp"

namespace po = boost::program_options;

class application_t {

protected:

  const std::string _name;
  log4cxx::LoggerPtr logger;

  virtual po::options_description
  get_named_options() const {
	 po::options_description desc;
	 return desc;
  };

  virtual po::positional_options_description
  get_positional_options() const {
	 po::positional_options_description p;
	 return p;
  };

// Function used to check that at least one of 'opt1' or 'opt2' are
// specified.
  void mode_options(const po::variables_map& vm,
						  const char* opt) const {
	 TRACE("Checking required option '" << opt <<
			 "'.");
	 if (!vm[opt].as<bool>())
		throw std::logic_error(std::string("Option '")
									  + opt + "' must be specified.");
  }
  void mode_options(const po::variables_map& vm,
						  const char* opt1,
						  const char* opt2) const {
	 TRACE("Checking alternative options '" << opt1 <<
			 "' and '" << opt2 << "'.");
	 if (!vm[opt1].as<bool>() && !vm[opt2].as<bool>())
		throw std::logic_error(std::string("At least one of the options '")
									  + opt1 + "' or '" + opt2 + "' must be specified.");
  }
  void mode_options(const po::variables_map& vm,
						  const char* opt1,
						  const char* opt2,
						  const char* opt3) const {
	 TRACE("Checking alternative options '" << opt1 <<
			 "', '" << opt2 << "', and '" << opt3 << "'.");
	 if (!vm[opt1].as<bool>() && !vm[opt2].as<bool>() && !vm[opt3].as<bool>())
		throw std::logic_error(std::string("At least one of the options '")
									  + opt1 + "', '" + opt2 + "', or '" + opt3 +
									  "' must be specified.");
  }
  void mode_options(const po::variables_map& vm,
						  const char* opt1,
						  const char* opt2,
						  const char* opt3,
						  const char* opt4) const {
	 TRACE("Checking alternative options '" << opt1 <<
			 "', '" << opt2 << "', '" << opt3 << "', and'" <<
			 opt4 << "'.");
	 if (!vm[opt1].as<bool>() && !vm[opt2].as<bool>()
		  && !vm[opt3].as<bool>() && !vm[opt4].as<bool>())
		throw std::logic_error(std::string("At least one of the options '")
									  + opt1 + "', '" + opt2 + "', '" + opt3 +
									  ", or '" + opt4 +
									  "' must be specified.");
  }
// Function used to check that 'opt1' and 'opt2' are not specified
// at the same time.
  void conflicting_options(const po::variables_map& vm,
									const char* opt1,
									const char* opt2) const {
	 TRACE("Checking conflicting options '" << opt1 <<
			 "' and '" << opt2 << "'.");
	 if (vm.count(opt1) && !vm[opt1].defaulted()
		  && vm.count(opt2) && !vm[opt2].defaulted())
		throw std::logic_error(std::string("Conflicting options '")
									  + opt1 + "' and '" + opt2 + "'.");
  }

// Function used to check that of 'for_what' is specified, then
//	'required_option' is specified too.
  void option_dependency(const po::variables_map& vm,
								 const char* for_what,
								 const char* required_option) const {
	 TRACE("Checking dependency between options '" << for_what <<
			 "' and '" << required_option << "'.");
	 if (vm.count(for_what) && !vm[for_what].defaulted())
		if (vm.count(required_option) == 0 || vm[required_option].defaulted())
		  throw std::logic_error(std::string("Option '") + for_what
										 + "' requires option '" + required_option + "'.");
  }

  virtual bool check_options(const po::variables_map& vm) {
	 return true;
  }

  virtual int execution(int argc,
								char** argv,
								const po::variables_map& vm) = 0;


public:
  explicit application_t(const std::string& name)
		:_name(name),
		 logger(log4cxx::Logger::getLogger("application"))
  {
	 initialize_logger();
  }

  virtual ~application_t() {
  }

  int execute(int argc, char** argv) {
	 boost::posix_time::ptime time_start(boost::posix_time::second_clock::local_time());
	 INFO(_name << " - started at " << time_start);
	 DEBUG("Source code version: " APPLICATION_SOURCE_VERSION);
	 DEBUG("Compiled on:         " __DATE__ "  -  " __TIME__);

	 int result= EXIT_SUCCESS;
	 try {

		DEBUG("Parsing program parameters...");
// Parse options
		po::variables_map vm;
		po::options_description desc("Help Options",
											  po::options_description::m_default_line_length,
											  po::options_description::m_default_line_length-16);
		desc.add_options()
		  ("help,?", po::bool_switch(),
			"Produce (this) help message.");
		desc.add(get_named_options());
		po::positional_options_description p= get_positional_options();

		po::command_line_parser parser(argc, argv);
		parser.style(po::command_line_style::default_style ^ po::command_line_style::allow_guessing);
		po::store(parser.options(desc).positional(p).run(), vm);

		po::notify(vm);
		DEBUG("Program parameters successfully parsed.");

		if (vm["help"].as<bool>()) {
// Generate the help message and exit
		  std::cout << _name << std::endl;
		  std::cout << desc << std::endl;
		  result= EXIT_SUCCESS;
		} else {
// Check parameter values
		  DEBUG("Checking program parameters...");
		  bool po_check= true;
		  std::string po_failing_desc= "";
		  try {
			 po_check= check_options(vm);
		  } catch (std::logic_error& e) {
			 FATAL("EXCEPTION OCCURRED: " << e.what());
			 po_failing_desc= e.what();
			 po_check= false;
		  }
		  if (!po_check) {
			 DEBUG("Check failed.");
			 std::cout << _name << std::endl;
			 std::cout << std::endl <<
				"*** Invalid parameters!" << std::endl;
			 if (po_failing_desc != "") {
				std::cout << "Reason: " << po_failing_desc << std::endl;
			 }
			 std::cout << std::endl;
			 std::cout << desc << std::endl;
			 result= EXIT_FAILURE;
		  } else {
// Execute
			 DEBUG("Check successfully completed.");
			 DEBUG("Beginning execution...");
			 result= execution(argc, argv, vm);
		  }
		  DEBUG("Execution successfully completed.");
		}


	 } catch (std::exception& e) {
		FATAL("EXCEPTION OCCURRED: " << e.what());
		result= EXIT_FAILURE;
	 } catch (assertion_failed_exception& e) {
		FATAL("ASSERTION FAILED: " << e.what());
		result= EXIT_FAILURE;
		throw e;
	 } catch (...) {
		FATAL("GENERIC EXCEPTION OCCURRED.");
		result= EXIT_FAILURE;
	 }

	 boost::posix_time::ptime time_end(boost::posix_time::second_clock::local_time());
	 INFO("Running time:        " << (time_end-time_start));
	 INFO(_name << " - terminated at " << time_end);

	 return result;
  }

};



#endif //__APPLICATION_HPP__
