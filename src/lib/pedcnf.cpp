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
 * pedcnf.cpp
 *
 * Structures to represent SAT instances derived from pedigrees.
 *
 **/

#include "pedcnf.hpp"
#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/trim.hpp>

const ped_var_kind_reHCstar ped_var_kind_reHCstar::SP(0);
const ped_var_kind_reHCstar ped_var_kind_reHCstar::SM(1);
const ped_var_kind_reHCstar ped_var_kind_reHCstar::P(2);
const ped_var_kind_reHCstar ped_var_kind_reHCstar::M(3);
const ped_var_kind_reHCstar ped_var_kind_reHCstar::RP(4);
const ped_var_kind_reHCstar ped_var_kind_reHCstar::RM(5);
const ped_var_kind_reHCstar ped_var_kind_reHCstar::E(6);
const ped_var_kind_reHCstar ped_var_kind_reHCstar::DUMMY(7);
const int ped_var_kind_reHCstar::int_values[]={0, 1, 2, 3, 4, 5, 6, 7};
const std::string ped_var_kind_reHCstar::str_values[]={"sp", "sm", "p", "m", "rp", "rm", "e", "dummy"};
const ped_var_kind_reHCstar ped_var_kind_reHCstar::enum_values[]={SP, SM, P, M, RP, RM, E, DUMMY};

/*
const ped_var_kind_r0HCstar ped_var_kind_r0HCstar::SP(0);
const ped_var_kind_r0HCstar ped_var_kind_r0HCstar::SM(1);
const ped_var_kind_r0HCstar ped_var_kind_r0HCstar::H(2);
const ped_var_kind_r0HCstar ped_var_kind_r0HCstar::W(3);
const ped_var_kind_r0HCstar ped_var_kind_r0HCstar::RP(4);
const ped_var_kind_r0HCstar ped_var_kind_r0HCstar::RM(5);
const ped_var_kind_r0HCstar ped_var_kind_r0HCstar::DUMMY(6);
const int ped_var_kind_r0HCstar::int_values[]={0, 1, 2, 3, 4, 5, 6};
const std::string ped_var_kind_r0HCstar::str_values[]={"sp", "sm", "h", "w", "rp", "rm", "dummy"};
const ped_var_kind_r0HCstar ped_var_kind_reHCstar::enum_values[]={SP, SM, H, W, RP, RM, DUMMY};

*/
