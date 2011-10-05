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

#include "timelimit.hpp"

#include "log.hpp"
#include "utility.hpp"

#include <csignal>
#include <cstdlib>
#include <sys/resource.h>

void exceeded_time_limit(int sig, siginfo_t * info, void * context) {
  ROOT_FATAL("Time limit reached. Exiting now...");
  exit(EXIT_reHC_ERROR);
};

void register_time_limit(unsigned long int secs) {
  if (secs>0){
	 if (secs<=MIN_TIME_LIMIT) secs= MIN_TIME_LIMIT;
	 ROOT_INFO("Setting time-limit to " << secs << " seconds.");
// Setting time-limit
	 struct sigaction act;
	 act.sa_sigaction= exceeded_time_limit;
	 sigemptyset(&act.sa_mask);
	 act.sa_flags= SA_SIGINFO;
	 sigaction(SIGXCPU, &act, NULL);
	 struct rlimit cpu_limit;
	 getrlimit(RLIMIT_CPU, &cpu_limit);
	 cpu_limit.rlim_cur= secs;
	 const int lim_res= setrlimit(RLIMIT_CPU, &cpu_limit);
	 if (lim_res!=0) {
		ROOT_ERROR("Failed to set the \"hard\" time-limit! "
					  "Error: " << lim_res << ". Continuing anyway...");
	 }

  }
};
