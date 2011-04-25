#!/bin/bash

##########
#
#                               ZRHC-*
#  Zero-Recombinant Haplotype Configuration with missing genotypes
#
#  Copyright (C) 2010,2011  Yuri Pirola <yuri.pirola(-at-)gmail.com>
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
#  generate-instances-JPSGCS.sh
#
#  A simple script that generate various instances using 'JPSGCS' and
#  'hc-generator.py' with various parameters.
#
#  Useful for 'serial' experimentations.
#
##########


LANG=C

if [ "${0:0:1}" == "/" ]; then
    script_dir=${1:-`dirname $0`}
else
    script_dir=${1:-`dirname ${PWD}/$0`}
fi

# Check existence of required files
if [ ! -f "./config-generation-JPSGCS.ini"  ];
then
    echo "Configuration file 'config-generation-JPSGCS.ini' not found. Aborting..."
    exit 1
fi
if [ ! -x ${script_dir}/hc-generator.py ]; then
    echo "Script 'hc-generator.py' not found in directory '${script_dir}'.  Aborting..."
    exit 1
fi

jpsgcs_dir=${script_dir}
if [ ! -f ${jpsgcs_dir}/JPSGCS.jar ]; then
    jpsgcs_dir=${script_dir}/../thirdparty
    echo "Program 'JPSGCS.jar' not found in directory '${script_dir}'.  Trying '${jpsgcs_dir}'..."
fi
if [ ! -f ${jpsgcs_dir}/JPSGCS.jar ]; then
    echo "Program 'JPSGCS.jar' not found in directory '${jpsgcs_dir}'.  Aborting..."
    exit 1
fi

# Set parameters to default values
#   Required
pedigree_sizes=-1
avg_family_sizes=-1
genotype_lengths=-1
missing_rates=-1
error_rates=-1
recomb_rates=-1
n_ped=-1
n_hc=-1

#   Optional
additional_JPSGCS_params=""
additional_hg_params=""
fmt_dl="%.5d"
fmt_ds="%.2d"
fmt_fl="%5.4f"
fmt_fs="%5.3f"
dest_dir="./"

# Read configuration
source ./config-generation-JPSGCS.ini

# Check required parameters...
if [ "${pedigree_sizes}" = "-1" ]; then echo "Parameter 'pedigree_sizes' not defined. Aborting..."; exit 1; fi
if [ "${avg_family_sizes}" = "-1" ]; then echo "Parameter 'avg_family_sizes' not defined. Aborting..."; exit 1; fi
if [ "${n_ped}" = "-1" ]; then echo "Parameter 'n_ped' not defined. Aborting..."; exit 1; fi
if [ "${genotype_lengths}" = "-1" ]; then echo "Parameter 'genotype_lengths' not defined. Aborting..."; exit 1; fi
if [ "${missing_rates}" = "-1" ]; then echo "Parameter 'missing_rates' not defined. Aborting..."; exit 1; fi
if [ "${recomb_rates}" = "-1" ]; then echo "Parameter 'recomb_rates' not defined. Aborting..."; exit 1; fi
if [ "${error_rates}" = "-1" ]; then echo "Parameter 'error_rates' not defined. Aborting..."; exit 1; fi
if [ "${n_hc}" = "-1" ]; then echo "Parameter 'n_hc' not defined. Aborting..."; exit 1; fi

# Check semantic of parameter values
if [ ! -e "${dest_dir}" ]; then
    mkdir -p "${dest_dir}"
fi
if [ ! -d "${dest_dir}" ]; then
    echo "Impossible to create output directory '${dest_dir}'. Aborting..."
    exit 1
fi

seed=2095
for size in ${pedigree_sizes}; do
    size_str=`printf "${fmt_dl}" ${size}`
    for fam in ${avg_family_sizes}; do
        fam_str=`printf "${fmt_ds}" ${fam}`
        mat_nod=$(((size-1)/(fam-1)))
        for ped in `seq 0 $((n_ped - 1))`; do
            ped_str=`printf "${fmt_ds}" ${ped}`
            prefix_id="size${size_str}-fam${fam_str}-conf${ped_str}"
            echo "Generating pedigree structure: pedigree-${prefix_id}.txt"
            echo "Generating pedigree structure: pedigree-${prefix_id}.txt" >&2
            java -cp ${jpsgcs_dir}/JPSGCS.jar ZeroLoopPed \
                ${mat_nod} ${size} ${additional_JPSGCS_params} \
                | awk '{print "1", $0}'  >  ${dest_dir}/pedigree-${prefix_id}.txt
            for length in ${genotype_lengths}; do
                length_str=`printf "${fmt_dl}" ${length}`
                for miss in ${missing_rates}; do
                    miss_str=`printf "${fmt_fs}" ${miss}`
                    for rec in ${recomb_rates}; do
                        rec_str=`printf "${fmt_fs}" ${rec}`
                        for err in ${error_rates}; do
                            err_str=`printf "${fmt_fl}" ${err}`
                            for hc in `seq 0 $((n_hc - 1))`; do
                                hc_str=`printf "${fmt_ds}" ${hc}`
                                id="${prefix_id}-length${length_str}-miss${miss_str/./_}-rec${rec_str/./_}-err${err_str/./_}-hc${hc_str}"
                                echo "Generating instance: gen-ped-${id}.txt"
                                echo "Generating instance: gen-ped-${id}.txt" >&2
                                ${script_dir}/hc-generator.py \
                                    -l ${length} \
                                    -m ${miss} \
                                    -e ${err} \
                                    -r ${rec} \
                                    -s ${seed} \
                                    ${additional_hg_params} \
                                    < ${dest_dir}/pedigree-${prefix_id}.txt \
                                    > ${dest_dir}/gen-ped-${id}.txt
                                seed=$((seed + 122))
                            done
                        done
                    done
                done
            done
        done
    done
done 2> generation.log
