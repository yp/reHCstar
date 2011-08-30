#!/usr/bin/python3

##########
#
#                               reHC-*
#  Haplotyping with Recombinations, Errors, and Missing Genotypes
#
#  Copyright (C) 2010,2011  Yuri Pirola <yuri.pirola(-at-)gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#
#  This file is part of reHC-* (reHCstar).
#
#  reHC-* is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  reHC-* is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with reHC-*.  If not, see <http://www.gnu.org/licenses/>.
#
##########

##########
#
#  reHCstar-big.py
#
#  A program that decomposes an instance with long genotypes is several fixed-length
#  blocks that are solved by reHCstar and then recombined.
#
##########

import array
import logging
import math
import optparse
import os
import random
import shlex
import subprocess
import sys
import time

gen2code= { '0 0': 0, '1 1': 1, '2 2': 2, '1 2': 3 }
code2gen= { 0: '0 0', 1: '1 1', 2: '2 2', 3: '1 2' }

def parse_command_line():
    usage= "usage: %prog [options]"
    parser= optparse.OptionParser(usage=usage)
    parser.add_option("-v", "--verbose",
                      action="count", dest="verbose",
                      default=0,
                      help="print additional log messages "
                      "(specified more than once increase verbosity)")

    io_group= optparse.OptionGroup(parser, "Input/Output Options",
                                   "Options that specify the input and "
                                   "the output files.")
    io_group.add_option("-p", "--pedigree",
                        action="store", dest="pedigree",
                        type="string", default=None,
                        help="the file containing the genotyped pedigree",
                        metavar="FILE")
    io_group.add_option("-r", "--results",
                        action="store", dest="haplotypes",
                        type="string", default=None,
                        help="the file that will contain the haplotype "
                        "configuration computed by reHC-*",
                        metavar="FILE")
    parser.add_option_group(io_group)

    search_group= optparse.OptionGroup(parser, "Block Options",
                                       "Options that regulates the division of "
                                       "the full genotypes into smaller blocks.")
    search_group.add_option("-l", "--block-length",
                            action="store", dest="length",
                            type="int", default=50,
                            help="the length of the blocks which the instance "
                            "is divided into  "
                            "(default: %default)",
                            metavar="LENGTH")
    search_group.add_option("-a", "--lookahead-length",
                            action="store", dest="lookahead",
                            type="int", default=0,
                            help="the length of the additional overlapping "
                            "blocks  (default: %default)",
                            metavar="LENGTH")
    parser.add_option_group(search_group)

    cmd_group= optparse.OptionGroup(parser, "Interoperability Options",
                                    "Options that specify how the 'reHCstar' "
                                    "is invoked.")
    cmd_group.add_option("--cmd",
                         action="store", dest="cmd", type="string",
                         default="./reHCstar -4 -p \"{pedigree}\" -h \"{haplotypes}\" -a \"{assumptions}\"",
                         help="the command-line used to invoke the 'reHCstar' program  "
                         "(default: %default)",
                         metavar="CMD-LINE")
    cmd_group.add_option("--cmd-rec",
                         action="store", dest="cmdrec", type="string",
                         default="--global-recomb --global-recomb-number \"{number}\"",
                         help="the options of the 'reHCstar' program used to specify "
                         "the maximum number of RECOMBINATIONS  (default: %default)",
                         metavar="PROGR-OPT")
    parser.add_option_group(cmd_group)

    (options, args) = parser.parse_args()

    return options

def write_haplotypes(filename, order, ped, haplotypes, stop, verbose= False):
    source_vector_enc= { 0 : ".",   1 : "O" }
    with open(filename, 'w') as haplotypes_file:
        if verbose:
            for ind in order:
                haplotypes_file.write("## INDIVIDUAL IND[" + ind + "] ")
                haplotypes_file.write("\t".join([str(x) for x in ped[ind]]))
                haplotypes_file.write("\n## PAT_HAPLOT IND[" + ind + "] ")
                haplotypes_file.write("".join([str(x) for x in haplotypes[ind][0][0:stop]]))
                haplotypes_file.write("\n## MAT_HAPLOT IND[" + ind + "] ")
                haplotypes_file.write("".join([str(x) for x in haplotypes[ind][1][0:stop]]))
                haplotypes_file.write("\n## PAT_SOUR_V IND[" + ind + "] ")
                haplotypes_file.write("".join([source_vector_enc[x] for x in haplotypes[ind][2][0:stop]]))
                haplotypes_file.write("\n## MAT_SOUR_V IND[" + ind + "] ")
                haplotypes_file.write("".join([source_vector_enc[x] for x in haplotypes[ind][3][0:stop]]))
                haplotypes_file.write("\n")

        for ind in order:
            haplotypes_file.write("\t".join([str(x) for x in ped[ind]]))
            haplotypes_file.write("\t")
            haplotypes_file.write("\t".join([ "|".join((str(x), str(y))) for x, y in zip(haplotypes[ind][0][0:stop], haplotypes[ind][1][0:stop]) ]))
            haplotypes_file.write("\n")


def compute_recombinations_and_gps(h, ph, known_gps=None):
    lh= len(h)
    (n_rec1, n_rec2)= (0, 0)
    (phase1, phase2)= (0, 1)
    (gps1, gps2)= ([0]*lh, [0]*lh)
    for (l, a, a1, a2) in zip(range(lh), h, ph[0], ph[1]):
        av= (a1, a2)
        if a != av[phase1]:
            phase1= 1-phase1
            assert(a == av[phase1])
            n_rec1 = n_rec1 + 1
        if a != av[phase2]:
            phase2= 1-phase2
            assert(a == av[phase2])
            n_rec2 = n_rec2 + 1
        gps1[l]= phase1
        gps2[l]= phase2
    if known_gps == None:
        if n_rec1 <= n_rec2:
            known_gps= 0
        else:
            known_gps= 1
    if known_gps==0:
        return (n_rec1, gps1)
    else:
        return (n_rec2, gps2)

def read_and_process_hc(haplotypes_filename, pedigree, complete_haplotypes, good_chunk_length):
    logging.debug("Reading and parsing the haplotype configuration computed for "
                  "the current chunk...")
    current_haplotypes= {}
    with open(haplotypes_filename, 'r') as haplotypes_file:
        for r in haplotypes_file:
            if not r.startswith('#'):
                row= r.strip().split("\t", 6)
                comb_hap= [ [int(x) for x in single.split("|")] for single in row[6].split("\t") ]
                (pat_hap, mat_hap)= [list(x) for x in zip(*comb_hap)]
                current_haplotypes[row[1]]= (pat_hap, mat_hap)
                if not ( row[1] in complete_haplotypes ):
                    complete_haplotypes[row[1]]= [ pat_hap[0:good_chunk_length],
                                                   mat_hap[0:good_chunk_length],
                                                   [], []]
                else:
                    complete_haplotypes[row[1]][0].extend( pat_hap[1:good_chunk_length] )
                    complete_haplotypes[row[1]][1].extend( mat_hap[1:good_chunk_length] )

    tot_rec= 0
    for ind in pedigree:
        row= pedigree[ind]
        if int(row[2]) > 0:
            # father is present
            (n_rec, gps)= compute_recombinations_and_gps(complete_haplotypes[row[1]][0],
                                                         complete_haplotypes[row[2]],
                                                         None)
            complete_haplotypes[row[1]][2]= gps
            tot_rec= tot_rec + n_rec

        if int(row[3]) > 0:
            # mother is present
            (n_rec, gps)= compute_recombinations_and_gps(complete_haplotypes[row[1]][1],
                                                         complete_haplotypes[row[3]],
                                                         None)
            complete_haplotypes[row[1]][3]= gps
            tot_rec= tot_rec + n_rec
    return tot_rec

def read_and_process_partial_hc(haplotypes_filename, pedigree, complete_haplotypes, fixed_index):
    logging.debug("Reading the haplotype configuration computed for "
                  "the current chunk...")
    hap_str= []
    with open(haplotypes_filename, 'r') as haplotypes_file:
        hap_str= [ line.strip()
                   for line in haplotypes_file
                   if not line.startswith('#')  ]

    logging.debug("Parsing the haplotype configuration...")
    current_haplotypes= {}
    for r in hap_str:
        split_r= r.split("\t", 6)
        comb_hap= split_r[6].split("\t")
        pat_hap= [ int(single.split('|')[0]) for single in comb_hap ]
        mat_hap= [ int(single.split('|')[1]) for single in comb_hap ]
        current_haplotypes[split_r[1]]= (pat_hap, mat_hap)

    tot_rec= 0
    for ind in pedigree:
        row= pedigree[ind]
        if int(row[2]) != 0:
            # father is present
            known_gps= None
            if row[1] in complete_haplotypes and len(complete_haplotypes[row[1]][2]) > 0:
                known_gps= complete_haplotypes[row[1]][2][fixed_index]
            (n_rec, gps)= compute_recombinations_and_gps(current_haplotypes[row[1]][0],
                                                         current_haplotypes[row[2]],
                                                         known_gps)
            tot_rec= tot_rec + n_rec

        if int(row[3]) != 0:
            # mother is present
            known_gps= None
            if row[1] in complete_haplotypes and len(complete_haplotypes[row[1]][3]) > 0:
                known_gps= complete_haplotypes[row[1]][3][fixed_index]
            (n_rec, gps)= compute_recombinations_and_gps(current_haplotypes[row[1]][1],
                                                         current_haplotypes[row[3]],
                                                         known_gps)
            tot_rec= tot_rec + n_rec
    return tot_rec

def basic_exec_reHCstar(pedigree_filename,
                        haplotypes_filename,
                        assumptions_filename,
                        max_recombs,
                        cmd_templ):
    rehcstar_success= False
    terminate= False
    full_cmd_str= cmd_templ['cmd'].format(pedigree=    pedigree_filename,
                                          haplotypes=  haplotypes_filename,
                                          assumptions= assumptions_filename)
    cmd= []
    try:
        if max_recombs>0:
            full_cmd_str+= " " + cmd_templ['rec'].format(number= str(max_recombs))
        cmd= shlex.split(full_cmd_str)
        logging.info("Invoking >%s<...", " ".join(cmd))
        retcode = subprocess.call(cmd, shell=False)
        rehcstar_success= (retcode==0)
    except OSError as e:
        logging.warn("reHCstar execution failed. "
                     "Command-line: >%s<. "
                     "Additional information: '%s'",
                     " ".join(cmd), e)

    return rehcstar_success


def exec_reHCstar(pedigree_filename, haplotypes_filename, assumptions_filename,
                  pedigree, complete_haplotypes, fixed_index):
    rehcstar_success= False
    terminate= False
    max_recombs= 0
    cmd_bin="./reHCstar"
    logging.info("Trying to solve the instance in file '%s'.", pedigree_filename)
    logging.info("Step 1. Searching an upper bound on the number of recombinations...")
    new_haplotypes_filename= haplotypes_filename
    while not terminate and not rehcstar_success:
        logging.debug("Step 1. Trying with at most %d recombinations.", max_recombs)
        rehcstar_success= basic_exec_reHCstar(pedigree_filename,
                                              haplotypes_filename,
                                              assumptions_filename,
                                              max_recombs, cmd)
        if rehcstar_success:
            logging.info("Step 1. Found a solution with at most %d recombinations.", max_recombs)
            new_haplotypes_filename="{}-success-{}".format(haplotypes_filename, str(os.getpid()))
            os.rename(haplotypes_filename, new_haplotypes_filename)
        else:
            if max_recombs > 2*len(pedigree):
                logging.info("Step 1. Solution NOT found. The instance requires "
                             "too much recombinations (more than %d). Aborting...", max_recombs)
                terminate= True
            else:
                logging.debug("Step 1. Solution not found yet. Increasing the maximum "
                              "number of recombinations.")
                max_recombs=2*(max_recombs+1)-1

    if rehcstar_success:
        if max_recombs == 0:
            logging.info("Step 2. A solution without recombinations has been found. It is optimal, thus bisection is skipped...")
        else:
            # Read the haplotype configuration and compute the real upper-bound
            lb= ((max_recombs+1)/2)-1
            ub= read_and_process_partial_hc(new_haplotypes_filename, pedigree, complete_haplotypes, fixed_index)
            logging.info("Step 2. Performing a bisection to find an optimal solution...")
            while lb+1 < ub:
                mid= math.floor((lb+ub)/2)
                logging.debug("Step 2. Bisecting interval (%d-%d]", lb, ub)
                logging.debug("Step 2. Trying with at most %d recombinations.", mid)
                rehcstar_success= basic_exec_reHCstar(pedigree_filename,
                                                      haplotypes_filename,
                                                      assumptions_filename,
                                                      mid, cmd)
                if rehcstar_success:
                    logging.debug("Step 2. Found a solution with at most %d recombinations.", mid)
                    new_haplotypes_filename="{}-success-{}".format(haplotypes_filename, str(os.getpid()))
                    os.rename(haplotypes_filename, new_haplotypes_filename)
                    ub= read_and_process_partial_hc(new_haplotypes_filename, pedigree, complete_haplotypes, fixed_index)
                    max_recombs= ub
                else:
                    logging.debug("Step 2. Solution NOT found in this half interval. Using the other...")
                    lb= mid
            logging.info("Step 2. Found an optimal solution with at most %d recombinations...", max_recombs)
            rehcstar_success= True

    return rehcstar_success




options = parse_command_line()

log_level= logging.DEBUG if options.verbose>0 else logging.INFO
verbose=  options.verbose>0
verbose1= options.verbose>0
verbose2= options.verbose>1
verbose3= options.verbose>2

logging.basicConfig(level=log_level,
                    format='%(levelname)5s [%(relativeCreated)15d] (%(filename)30s:%(lineno)-4d) - %(message)s',
                    stream=sys.stderr)

logging.info("reHCstar manager - started at %s", time.asctime())

# Check program options
#  input genotyped pedigree
if ( not options.pedigree or
     not os.path.isfile(options.pedigree) or
     not os.access(options.pedigree, os.R_OK) ):
    logging.fatal("Pedigree file not specified or not accessible. Given '%s'",
                  options.pedigree)
    sys.exit(2)
#  output haplotype configuration
if ( not options.haplotypes ):
    logging.fatal("Result file not specified. Given '%s'",
                  options.haplotypes)
    sys.exit(2)

#  block lengths
if options.length <= 0:
    logging.fatal("Maximum block length must be a positive integer. Given '%d'",
                  options.length)
    sys.exit(2)
if options.lookahead < 0:
    logging.fatal("Overlapping block length must be a non-negative integer. Given '%d'",
                  options.lookahead)
    sys.exit(2)

#  command-line related options
cmd= {}
if not ( '{pedigree}' in options.cmd and
         '{haplotypes}' in options.cmd and
         '{assumptions}' in options.cmd ):
    logging.fatal("The command-line *MUST* include the following placeholders: "
                  "'{pedigree}' '{haplotypes}' '{assumptions}'. "
                  "Given: '%s'",
                  options.cmd)
    sys.exit(2)
if not '{number}' in options.cmdrec:
    logging.fatal("The options for handling recombinations *MUST* include "
                  "the following placeholder: '{number}'. "
                  "Given: '%s'",
                  options.cmdrec)
    sys.exit(2)
cmd['cmd']= options.cmd
cmd['rec']= options.cmdrec



logging.info("CONFIGURATION:")
logging.info("Input pedigree:       '%s'", options.pedigree)
logging.info("Resulting haplotypes: '%s'", options.haplotypes)
logging.info("Max. block length:     %d", options.length)
logging.info("Lookahead length:      %d", options.lookahead)


# Read the original pedigree
logging.info("Reading pedigree file '%s'...", options.pedigree)
ped_str= []
with open(options.pedigree, 'r') as ped_file:
    ped_str= [ line.strip()
               for line in ped_file
               if not line.startswith('#')  ]
logging.debug("Parsing the pedigree...")
genotypes= {}
ped= {}
order= []
backorder= {}
gen_len= 0
i= 0
for r in ped_str:
    split_r= r.split("\t", 6)
    genotypes[ split_r[1] ]= array.array('b', [gen2code[x] for x in split_r[6].split("\t")])
    gen_len= len(genotypes[ split_r[1] ])
    order.append(split_r[1])
    backorder[split_r[1]]= i
    ped[split_r[1]]= split_r[0:6]
    i= i+1
del ped_str
logging.info("Genotyped pedigree '%s' successfully read.", options.pedigree)

logging.info("Starting the computation of the haplotype configuration...")
assumptions= []
complete_haplotypes= {}
tot_rec= 0
suffix="{}-{}".format(os.path.basename(options.pedigree), os.path.basename(options.haplotypes))
for start in range(0, gen_len, options.length):
    good_stop= min( gen_len, start + options.length + 1 )
    stop=      min( gen_len, good_stop + options.lookahead )
    good_chunk_length= good_stop - start
    logging.info("Considering pedigree chunk between the markers [%d-%d)...", start, stop)
    # Write genotype block
    genfile= "tmp-pedigree-{}-{}-{}".format(start, stop, suffix)
    with open(genfile, 'w') as out_file:
        for ind in order:
            out_file.write("\t".join([str(x) for x in ped[ind]]) + "\t" +
                           "\t".join([code2gen[x] for x in genotypes[ind][start:stop]]) + "\n")

    # Write assumptions
    assumptions_filename= "tmp-assumptions-{}-{}-{}".format(start, stop, suffix)
    with open(assumptions_filename, 'w') as assumptions_file:
        assumptions_file.write("\n".join(assumptions))
        assumptions_file.write("\n")

    # Execute reHCstar
    haplotypes_filename= "tmp-haplotypes-{}-{}-{}".format(start, stop, suffix)
    rehcstar_success= exec_reHCstar(genfile, haplotypes_filename, assumptions_filename,
                                    ped, complete_haplotypes, start)

    if rehcstar_success:
        # Reading the haplotype configuration
        logging.info("reHCstar has successfully computed a solution on the chunk [%d-%d).", start, stop)
        logging.debug("Reading the haplotype configuration computed for "
                      "the current chunk [%d-%d)...", start, stop)
        tot_rec= read_and_process_hc(haplotypes_filename+"-success-"+str(os.getpid()), ped, complete_haplotypes, good_chunk_length)
        logging.debug("The haplotype configuration has %d recombinations so far", tot_rec)

        if verbose2:
            incremental_haplotypes_filename= "tmp-incremental-{}-{}-{}".format(start,
                                                                               stop,
                                                                               suffix)
            write_haplotypes(incremental_haplotypes_filename, order, ped,
                             complete_haplotypes, good_stop, verbose3)

        assumptions= []
        if start + options.length < gen_len:
            for ind in order:
                assumptions.append(" ".join(("p", str(ind), "0", str(complete_haplotypes[ind][0][good_stop-1]-1))))
                assumptions.append(" ".join(("m", str(ind), "0", str(complete_haplotypes[ind][1][good_stop-1]-1))))
                if int(ped[ind][2]) > 0:
                    assumptions.append(" ".join(("sp", str(ind), "0", str(complete_haplotypes[ind][2][good_stop-1]))))
                if int(ped[ind][3]) > 0:
                    assumptions.append(" ".join(("sm", str(ind), "0", str(complete_haplotypes[ind][3][good_stop-1]))))

    else:
        # No solution of reHCstar
        logging.info("reHCstar did NOT computed a solution on the chunk [%d-%d]. Aborting...", start, stop-1)
        sys.exit(1)

logging.info("Found a complete haplotype configuration with (at most) %d recombinations.", tot_rec)

logging.info("Writing the haplotype configuration to '%s'...", options.haplotypes)
write_haplotypes(options.haplotypes, order, ped, complete_haplotypes, gen_len, verbose2)

logging.info("reHCstar manager - completed at %s", time.asctime())

