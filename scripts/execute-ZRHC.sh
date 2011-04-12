#!/bin/bash

##########
#
#                               ZRHC-*
#  Zero-Recombinant Haplotype Configuration with missing genotypes
#
#  Copyright (C) 2010  Yuri Pirola <yuri.pirola(-at-)gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#
#  This file is part of ZRHC-* (ZRHCstar).
#
#  ZRHC-* is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ZRHC-* is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with ZRHC-*.  If not, see <http://www.gnu.org/licenses/>.
#
##########

##########
#
#  execute-ZRHC.sh
#
#  A simple script that execute ZRHC in the 'create-read' mode on
#  every genotyped pedigree of the current directory.
#
#  Useful for 'serial' experimentations.
#
##########

LANG=C

# Check existence of required files
if [ ! -f "./config-execution.ini"  ];
then
    echo "Configuration file 'config-execution.ini' not found. Aborting..."
    exit 1
fi

# Set parameters to default values
ZRHC_exe=""
ZRHC_sat_mode=""
error_handler=""
recomb_handler=""
pedigrees=gen-ped-*.txt
dest_dir="./zrhcstar-out/"

# Read configuration
source ./config-execution.ini

# Check required parameters...
if [ -z "${ZRHC_exe}" ]; then echo "Parameter 'ZRHC_exe' not defined. Aborting..."; exit 1; fi
if [ -z "${ZRHC_sat_mode}" ]; then echo "Parameter 'ZRHC_sat_mode' not defined. Aborting..."; exit 1; fi
if [ -z "${pedigrees}" ]; then echo "Parameter 'pedigrees' not defined. Aborting..."; exit 1; fi


if [ ! -x "${ZRHC_exe}" ]; then echo "Main program '${ZRHC_exe}' not found or not executable. Aborting..."; exit 1; fi

for pedigree in ${pedigrees}; do
    if [ ! -f ${pedigree} ]; then
        echo "Pedigree '${pedigree}' refers to a non-existent file. Aborting..."
        exit 1
    fi
done

if [ ! -e "${dest_dir}" ]; then
    mkdir -p "${dest_dir}"
fi
if [ ! -d "${dest_dir}" ]; then
    echo "Impossible to create output directory '${dest_dir}'. Aborting..."
    exit 1
fi


ZRHC_partial="${ZRHC_exe} ${ZRHC_sat_mode} ${error_handler} ${recomb_handler}"

echo "`date`  --  Starting experimentation" > execution.log
for full_pedigree in ${pedigrees}; do
    pedigree=`basename ${full_pedigree}`
    echo "`date`  --  Executing on ${pedigree}"
    echo "`date`  --  Executing on ${pedigree}" >> execution.log
    nice time -f "%U %S %E %x %M %C" -o ${dest_dir}/time-${pedigree} \
        ${ZRHC_partial} \
        -p ${full_pedigree} -h ${dest_dir}/hap-${pedigree} > ${dest_dir}/log-${pedigree}
done
