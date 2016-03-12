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
 * sat_solver_interface.hpp
 *
 * Structures to represent a light-weight interface to an internal SAT solver.
 * It only wants to narrow down the interface to the SAT solver.
 *
 **/

#ifndef __SAT_SOLVER_INTERFACE_HPP__
#define __SAT_SOLVER_INTERFACE_HPP__

#include <boost/cstdint.hpp>

typedef boost::int_fast32_t lit_t;
typedef boost::int_fast32_t var_t;


// Check that the 'right' preprocessor symbols have been defined:
// NO_INTERNAL_SAT_SOLVER
// INTERNAL_SAT_SOLVER
// ONLY_INTERNAL_SAT_SOLVER (implies INTERNAL_SAT_SOLVER)
//
// 1- check at least one is defined
#if !defined(ONLY_INTERNAL_SAT_SOLVER) && !defined(INTERNAL_SAT_SOLVER) && !defined(NO_INTERNAL_SAT_SOLVER)
// default= ONLY_INTERNAL_SAT_SOLVER
#define ONLY_INTERNAL_SAT_SOLVER
#endif

// 2- ONLY_INTERNAL_SAT_SOLVER => INTERNAL_SAT_SOLVER
#if defined(ONLY_INTERNAL_SAT_SOLVER) && !defined(INTERNAL_SAT_SOLVER)
#define INTERNAL_SAT_SOLVER
#endif

// 3- NO_INTERNAL_SAT_SOLVER => not INTERNAL_SAT_SOLVER
#if defined(NO_INTERNAL_SAT_SOLVER) && defined(INTERNAL_SAT_SOLVER)
#error "We have to choose only one alternative: NO_INTERNAL_SAT_SOLVER or INTERNAL_SAT_SOLVER"
#endif


#ifndef INTERNAL_SAT_SOLVER
#ifdef SAT_SOLVER
#undef SAT_SOLVER
#endif
#define SAT_SOLVER NO internal SAT solver
#endif

// Provide the interface only if asked to do so
#ifdef INTERNAL_SAT_SOLVER

// Enable the "right" interface depending on the macro definitions

#if defined(USE_CRYPTOMINISAT) && defined(USE_MINISAT)
#error "Only one SAT solver can be integrated: please choose among CryptoMiniSat and MiniSat"
#endif

#if !defined(USE_CRYPTOMINISAT) && !defined(USE_MINISAT)
#message "No SAT solver specified: enabling 'CryptoMiniSat' by default"
#define USE_CRYPTOMINISAT
#endif


#ifdef USE_CRYPTOMINISAT

#include "cms_sat_solver_interface.hpp"

#endif


#ifdef USE_MINISAT

#include "ms_sat_solver_interface.hpp"

#endif


#endif

#endif // __SAT_SOLVER_INTERFACE_HPP__
