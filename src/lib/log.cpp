/**
 *
 *                              reHC-*
 * Haplotyping with Recombinations, Errors, and Missing Genotypes
 *
 * Copyright (C) 2010  Yuri Pirola <yuri.pirola(-at-)gmail.com>
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
#include "log.hpp"

#include <boost/filesystem.hpp>
#include <log4cxx/propertyconfigurator.h>
#include <log4cxx/helpers/exception.h>
#include <log4cxx/logstring.h>
#include <log4cxx/patternlayout.h>
#include <log4cxx/consoleappender.h>
#include <log4cxx/logmanager.h>
#include <log4cxx/logger.h>

#include "configuration.h"

using namespace std;
using namespace log4cxx;
using namespace log4cxx::helpers;
using namespace boost::filesystem;

int initialize_logger() {
  static bool configured= false;
  if (configured)
	 return EXIT_SUCCESS;
  int result = EXIT_SUCCESS;
  try {
	 bool file_found= false;
	 path path(current_path());
	 while (!(file_found= exists(path / log_configuration_default_file))
				 && (path != path.root_path())) {
		path= path.parent_path();
	 }
	 if (file_found) {
		try {
		  PropertyConfigurator::configure((path / log_configuration_default_file).string());
		  configured= true;
		} catch (Exception& e) {
		  configured= false;
		}
	 }
	 if (!configured) {
           LogManager::getLoggerRepository()->setConfigured(true);
           LoggerPtr root = Logger::getRootLogger();
           static const LogString TTCC_CONVERSION_PATTERN(LOG4CXX_STR("%5p [%15r] (%30.30F:%-4L) - %m%n"));
           LayoutPtr layout(new PatternLayout(TTCC_CONVERSION_PATTERN));
           AppenderPtr appender(new ConsoleAppender(layout, LOG4CXX_STR("System.err")));
           root->addAppender(appender);
           root->setLevel(Level::getInfo());
           configured= true;
	 }
	 if (!configured) {
		result = EXIT_FAILURE;
	 }
  } catch(Exception&) {
	 result = EXIT_FAILURE;
  }
  return result;
}

