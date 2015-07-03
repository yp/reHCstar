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
#  reHCstar-mgr.py
#
#  A program that decomposes an instance with long genotypes is several fixed-length
#  blocks that are solved by reHCstar and then recombined.
#
##########

import logging
import math
import optparse
import os
import random
import shlex
import subprocess
import sys
import time
import resource

MISSING_G = '0 0'
source_vector_enc= { 0 : ".",   1 : "*" }

def is_genotyped(g):
    return g != MISSING_G

def is_same_genotype(g, a1, a2):
    return (g == "{} {}".format(a1, a2)) if (a1 <= a2) else (g == "{} {}".format(a2, a1))

rehcstar_license_msgs= [
    "reHCstar-mgr -- reHC-* manager",
    "Haplotyping on Pedigrees with Recombinations, Errors, and Missing Genotypes.",
    "Copyright (C) 2010, 2011  Yuri Pirola <yuri.pirola(-at-)gmail.com>.",
    "This program is distributed under the terms of the GNU General Public License (GPL), either version 3 of the License, or (at your option) any later version.",
    "This program comes with ABSOLUTELY NO WARRANTY. See the GNU General Public License for more details.",
    "This is free software, and you are welcome to redistribute it under the conditions specified by the license."
]

class MyIndentedHelpFormatter(optparse.IndentedHelpFormatter):
    def format_usage(self, usage):
        return "%s\n" % usage

class REHCstarMgrError(Exception):
    """Base class for exceptions of reHCstar-mgr."""
    pass


class REHCstarMgrOutOfTimeError(REHCstarMgrError):
    """Exception raised if the execution exceeds the time limit.

    Attributes:
        msg    -- explanation of the error
    """

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return '{0}'.format(self.msg)

class TimeLimit:
    """A class representing a CPU running time limit."""

    def __init__(self, time_limit= None):
        # Compute the minimum system limit (None= unlimited)
        (soft_sys_limit, hard_sys_limit)= resource.getrlimit(resource.RLIMIT_CPU)
        assert soft_sys_limit != -1 or hard_sys_limit == -1
        if hard_sys_limit == -1:
            hard_sys_limit= soft_sys_limit
        assert soft_sys_limit <= hard_sys_limit
        sys_limit= None if soft_sys_limit <= 0 else soft_sys_limit
        # Check the given time limit (None= unlimited)
        if time_limit is not None and time_limit <= 0:
            time_limit= None
        assert sys_limit is None or sys_limit > 0
        assert time_limit is None or time_limit > 0
        if sys_limit is None:
            self.limit= time_limit
        elif time_limit is None:
            self.limit= sys_limit
        else:
            self.limit= min(time_limit, sys_limit)



    def get_remaining_time(self):
        time_used= ( resource.getrusage(resource.RUSAGE_CHILDREN)[0] +
                     resource.getrusage(resource.RUSAGE_SELF)[0] )
        logging.debug("Time limit: %s",
                      "unlimited" if self.limit is None
                      else "{0:.3f} s".format(self.limit))
        logging.debug("Time used:  %.3f s", time_used)
        if self.limit is None:
            return (None, self.limit, time_used)
        else:
            return (self.limit - time_used, self.limit, time_used)

    def check_remaining_time(self):
        (rem_time, time_limit, time_used)= self.get_remaining_time()

        if rem_time is None:
            return None
        if rem_time <= 0:
            msg= "Maximum CPU time exceeded. (used {:.2f} s, given {:d} s)".format( time_used, time_limit )
            logging.warn(msg)
            raise REHCstarMgrOutOfTimeError(msg)

        return max(time_limit - time_used, 1)




class Individual:
    '''A class representing an individual and its parents.'''
    family= "0"
    id=     "0"
    father= "0"
    mother= "0"
    gender= "0"

    def __init__(self, indiv_id):
        self.id= indiv_id

    def parse_list(self, l):
        assert len(l)==6, "The initializer list must have 6 elements. Given {}".format(len(l))
        assert l[1] == self.id, "Individual ID cannot change. Was {}. Given {}.".format(self.id,
                                                                                        l[1])
        (self.family, self.id, self.father,
         self.mother, self.gender, self.phenotype)= l

    def __str__(self):
        return "\t".join((self.family, self.id, self.father,
                          self.mother, self.gender, self.phenotype))

    def has_father(self):
        return (int(self.father)!=0)

    def has_mother(self):
        return (int(self.mother)!=0)

def parse_genotype(genstr):
    genv = genstr.split()
    assert len(genv)%2==0
    return [ " ".join(sorted((genv[2*i], genv[2*i+1]))) for i in range(int(len(genv)/2)) ]

class Solution:
    '''A class representing both the input and the output.'''
    pedigree=    {}
    genotypes=   {}
    individuals= []
    haplotypes=  {}
    gps=         {}
    recombs=     []
    errors=      []
    genotype_length= 0
    chunk_info=  {}
    notes=       []

    def __init__(self):
        pass

    def read_genotyped_pedigree(self, pedigree_filename):
        '''Read a genotyped pedigree from a file.'''

        logging.debug("Reading pedigree file '%s'...", pedigree_filename)

        with open(pedigree_filename, 'r') as ped_file:
            for row in ped_file:
                if row.startswith('#'):
                    continue
                split_row= row.strip().split(None, 6)
                indiv_id= split_row[1]
                indiv= Individual(indiv_id)
                indiv.parse_list(split_row[0:6])
                self.pedigree[ indiv_id ]= indiv
                self.genotypes[ indiv_id ]= parse_genotype(split_row[6])
                self.genotype_length= len(self.genotypes[ indiv_id ])
                self.individuals.append(indiv)

        to_add = []
        for indiv in self.pedigree.values():
            if indiv.has_father() and indiv.father not in self.pedigree:
                father = Individual(indiv.father)
                father.parse_list([ indiv.family, indiv.father, "0", "0", "1", indiv.phenotype ])
                to_add.append((indiv.father, father))
                self.genotypes[ indiv.father ] = [ MISSING_G ] * len(self.genotypes[ indiv.id ])
                self.individuals.append(father)
            if indiv.has_mother() and indiv.mother not in self.pedigree:
                mother = Individual(indiv.mother)
                mother.parse_list([ indiv.family, indiv.mother, "0", "0", "1", indiv.phenotype ])
                to_add.append((indiv.mother, mother))
                self.genotypes[ indiv.mother ] = [ MISSING_G ] * len(self.genotypes[ indiv.id ])
                self.individuals.append(mother)


        self.pedigree.update(to_add)
        del to_add


        # Initialize other information
        haplotypes= {}
        gps=        {}
        recombs=    []
        errors=     []
        chunk_info= {}

        logging.info("Genotyped pedigree successfully read from file '%s' "
                     "(individuals=%d, genotype length=%d).",
                     pedigree_filename, len(self.individuals), self.genotype_length)

    def write_genotype_block(self, filename, start=None, stop=None):
        if start is None:
            start= 0
        if stop is None:
            stop= self.genotype_length
        with open(filename, 'w') as genotype_file:
            for ind in self.individuals:
                genotype_file.write(str(ind))
                genotype_file.write("\t")
                genotype_file.write("\t".join(self.genotypes[ind.id][start:stop]))
                genotype_file.write("\n")


    def write_haplotypes(self, filename, stop=None, verbose= False):

        logging.debug("Writing haplotypes to file '%s'%s...", filename, " (verbose)" if verbose else "")

        if stop is None:
            stop= self.genotype_length ## Genotype length should always be greater than haplotype length

        with open(filename, 'w') as haplotypes_file:
            haplotypes_file.write("#### BEGIN NOTES\n")
            for note in self.notes:
                haplotypes_file.write("#### {}\n".format(note))
            haplotypes_file.write("#### END NOTES\n")

            if verbose:
                # Chunk information
                haplotypes_file.write("### CHUNK INFORMATION\n")
                for chunk in self.chunk_info:
                    chunk_str="## CHUNK: '{}'".format(chunk)
                    for key in self.chunk_info[chunk]:
                        haplotypes_file.write("{}@'{}': '{}'\n".format(chunk_str, key,
                                                                       self.chunk_info[chunk][key]))
                # Recombinations
                haplotypes_file.write("### RECOMBINATIONS (number={})\n".format(len(self.recombs)))
                for recomb in self.recombs:
                    haplotypes_file.write("## REC {}\n".format(" ".join([str(x) for x in recomb])))
                # Errors
                haplotypes_file.write("### ERRORS (number={})\n".format(len(self.errors)))
                for error in self.errors:
                    haplotypes_file.write("## ERR {}\n".format(" ".join([str(x) for x in error])))
                # Haplotypes and grand-parental sources
                for ind in self.individuals:
                    ind_id= "{:10s}".format(ind.id)
                    haplotypes_file.write("## INDIVIDUAL IND[" + ind_id + "] ")
                    haplotypes_file.write( str(ind) )
                    haplotypes_file.write("\n## GENOTYPE   IND[" + ind_id + "] ")
                    haplotypes_file.write(" ".join( [ x.replace(" ", "")
                                                     for x in self.genotypes[ind.id][:stop] ] ) )
                    haplotypes_file.write("\n## MASK_ERR   IND[" + ind_id + "] ")
                    haplotypes_file.write("".join( [ "X" if is_genotyped(g) and not is_same_genotype(g, a1, a2)
                                                     else " "
                                                     for g, a1, a2 in zip(self.genotypes[ind.id][:stop],
                                                                          *self.haplotypes[ind.id]) ] ) )
                    haplotypes_file.write("\n## PAT_HAPLOT IND[" + ind_id + "] ")
                    haplotypes_file.write("".join( [ str(x)
                                                     for x in self.haplotypes[ind.id][0][:stop] ] ) )
                    haplotypes_file.write("\n## MAT_HAPLOT IND[" + ind_id + "] ")
                    haplotypes_file.write("".join( [ str(x)
                                                     for x in self.haplotypes[ind.id][1][:stop] ] ) )
                    haplotypes_file.write("\n## PAT_SOUR_V IND[" + ind_id + "] ")
                    haplotypes_file.write("".join( [ source_vector_enc[x]
                                                     for x in self.gps[ind.id][0][:stop] ] ) )
                    haplotypes_file.write("\n## MAT_SOUR_V IND[" + ind_id + "] ")
                    haplotypes_file.write("".join( [ source_vector_enc[x]
                                                     for x in self.gps[ind.id][1][:stop] ] ) )
                    haplotypes_file.write("\n")

            # Haplotypes
            for ind in self.individuals:
                haplotypes_file.write( str(ind) )
                haplotypes_file.write("\t")
                haplotypes_file.write("\t".join( [ "|".join((str(x), str(y)))
                                                   for x, y in zip(self.haplotypes[ind.id][0][:stop],
                                                                   self.haplotypes[ind.id][1][:stop]) ] ) )
                haplotypes_file.write("\n")

        logging.debug("Haplotypes successfully written to file '%s'%s...",
                      filename, " (verbose)" if verbose else "")

    def compute_recombinations_and_gps(self, indiv_id, parent_id):
        gender_idx= int(self.pedigree[parent_id].gender) - 1
        h= self.haplotypes[indiv_id][gender_idx]
        ph= self.haplotypes[parent_id]
        lh= len(h)
        (n_rec1, n_rec2)= (0, 0)
        (phase1, phase2)= (0, 1)
        gps= (gps1, gps2)= ([0]*lh, [0]*lh)
        for (l, a, a1, a2) in zip(range(lh), h, ph[0], ph[1]):
            av= (a1, a2)
            if a != av[phase1]:
                phase1= 1-phase1
                assert a == av[phase1]
                n_rec1= n_rec1 + 1
            if a != av[phase2]:
                phase2= 1-phase2
                assert a == av[phase2]
                n_rec2= n_rec2 + 1
            gps1[l]= phase1
            gps2[l]= phase2
        best_gps= 0 if n_rec1 <= n_rec2 else 1
        rec= [ (indiv_id, locus, gender_idx)
               for locus in range(1, lh)
               if gps[best_gps][locus] != gps[best_gps][locus-1] ]
        return (rec, gps[best_gps])

    def compute_best_solution(self):
        logging.debug("Recomputing the assignment of grand parental sources which minimizes recombinations...")
        self.recombs= []
        self.errors= []
        for ind in self.individuals:
            self.gps[ind.id]= [ [], [] ]
            if ind.has_father():
                # father is present
                (rec, gps)= self.compute_recombinations_and_gps(ind.id, ind.father)
                self.recombs.extend(rec)
                self.gps[ind.id][0]= gps

            if ind.has_mother():
                # father is present
                (rec, gps)= self.compute_recombinations_and_gps(ind.id, ind.mother)
                self.recombs.extend(rec)
                self.gps[ind.id][1]= gps

            (h1, h2)= self.haplotypes[ind.id]
            gv= self.genotypes[ind.id]
            assert len(h1) == len(h2)
            assert len(h1) <= len(gv)
            lh= len(h1)
            self.errors.extend([ (ind.id, locus)
                                 for locus, g, a1, a2 in zip(range(lh), gv, h1, h2)
                                 if is_genotyped(g) and not is_same_genotype(g, a1, a2) ])

    def analyze_optimality(self):
        valid_chunk_info= [ chunk for chunk, info in self.chunk_info.items() if info["chunk included in solution"] ]
        self.notes.append("VALID CHUNKS {}".format(valid_chunk_info))
        is_single_chunk= len(valid_chunk_info)==1
        is_all_optimal= all([ self.chunk_info[valid_chunk]["optimum found"] for valid_chunk in valid_chunk_info ])
        is_optimal= is_single_chunk and is_all_optimal
        solution_start= min((chunk[0] for chunk in valid_chunk_info))
        solution_end= max((chunk[1] for chunk in valid_chunk_info))
        if (solution_start, solution_end) == (0, self.genotype_length):
            self.notes.append("COMPLETE SOLUTION")
            logging.info("The solution is COMPLETE.")
        else:
            self.notes.append("PARTIAL SOLUTION")
            logging.info("The solution is PARTIAL.")
        if is_optimal:
            self.notes.append("OPTIMAL SOLUTION")
            logging.info("The solution is OPTIMAL.")
        else:
            self.notes.append("SUBOPTIMAL SOLUTION")
            if is_all_optimal:
                logging.info("The solution is *NOT* OPTIMAL since it is composed by several (optimal) chunks.")
            else:
                logging.info("The partial solution is *NOT* OPTIMAL since it at least one chunk is not optimal.")




def parse_command_line():
    usage= "\n".join( rehcstar_license_msgs +
                      [ "", "Usage:", "  %prog --pedigree 'PEDIGREE FILE' --results 'RESULT FILE' [options]"] )
    parser= optparse.OptionParser(usage= usage,
                                  formatter= MyIndentedHelpFormatter())
    parser.add_option("-v", "--verbose",
                      action="count", dest="verbose",
                      default=0,
                      help="print additional log messages "
                      "(specified more than once increase verbosity)")

    io_group= optparse.OptionGroup(parser, "Input/Output Options (MANDATORY)",
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
                         default="reHCstar -4 -p \"{pedigree}\" -h \"{haplotypes}\" -a \"{assumptions}\"",
                         help="the command-line used to invoke the 'reHCstar' program  "
                         "(default: %default)",
                         metavar="CMD-LINE")
    cmd_group.add_option("--cmd-rec",
                         action="store", dest="cmdrec", type="string",
                         default="--global-recomb --global-recomb-number \"{number}\" --global-recomb-min-number \"{min_number}\"",
                         help="the options of the 'reHCstar' program used to specify "
                         "the maximum (and optionally minimum) number of "
                         "RECOMBINATIONS  (default: %default)",
                         metavar="PROGR-OPT")
    cmd_group.add_option("--cmd-time",
                         action="store", dest="cmdtime", type="string",
                         default="--time-limit {time}",
                         help="the option of the 'reHCstar' program used to specify "
                         "the maximum CPU time  (default: %default)",
                         metavar="PROGR-OPT")
    parser.add_option_group(cmd_group)

    oth_group= optparse.OptionGroup(parser, "Search Options",
                                    "Other options that specify how 'reHCstar-mgr' behaves.")
    oth_group.add_option("--recomb-upper-limit",
                         action="store", dest="ulrecomb",
                         type="int", default=-1,
                         help="haplotype configurations with more that this number of "
                         "recombinations are not searched for  (default: %default)",
                         metavar="int")
    oth_group.add_option("--initial-recomb-lb",
                         action="store", dest="recomb_lb",
                         type="int", default=-1,
                         help="an initial pre-determined lower bound of the number of "
                         "recombinations (i.e. a haplotype configuration with that number "
                         "of recombinations does not exist)  (default: %default)",
                         metavar="int")
    oth_group.add_option("--initial-recomb-ub",
                         action="store", dest="recomb_ub",
                         type="int", default=-1,
                         help="an initial pre-determined upper bound of the number of "
                         "recombinations (i.e. a haplotype configuration with that number "
                         "of recombinations certainly exists) or '-1' to disable  "
                         "(default: %default)",
                         metavar="int")
    oth_group.add_option("--initial-haplotype-configuration",
                         action="store", dest="init_hc",
                         type="string", default=None,
                         help="a file containing an initial haplotype configuration "
                         "that possibly induces more recombinations than the optimum. "
                         "It conflicts with '--initial-recomb-ub' and cannot be used "
                         "with the automatic genotype segmentation (i.e. block length "
                         "must be greater than genotype length)  (default: %default)",
                         metavar="FILE")
    oth_group.add_option("--bootstrap",
                         action="store_true", dest="bootstrap",
                         default=False,
                         help="enable the bootstrap phase if an initial lower bound is not provided. "
                         "It conflicts with '--initial-recomb-lb')  "
                         "(default: %default)")
    oth_group.add_option("--bootstrap-time-limit",
                         action="store", dest="bootstrap_time_limit",
                         type="int", default=120,
                         help="the maximum CPU time (in seconds) for the bootstrap phase  "
                         "(default: %default)",
                         metavar="SECONDS")
    oth_group.add_option("--time-limit",
                         action="store", dest="time_limit",
                         type="int", default=0,
                         help="the maximum CPU time (in seconds) for 'reHCstar-mgr' execution "
                         "or '0' to disable  (default: %default)",
                         metavar="SECONDS")
    parser.add_option_group(oth_group)

    (options, args) = parser.parse_args()

    options.cmdline = " ".join(sys.argv)

    return options

def check_program_options(options):
    #  input genotyped pedigree
    if ( not options.pedigree or
         not os.path.isfile(options.pedigree) or
         not os.access(options.pedigree, os.R_OK) ):
        logging.fatal("Pedigree file not specified or not accessible. Given '%s'. Aborting...",
                      options.pedigree)
        sys.exit(2)
    #  output haplotype configuration
    if not options.haplotypes:
        logging.fatal("Result file not specified. Aborting...")
        sys.exit(2)

    #  block lengths
    if options.length <= 0:
        logging.fatal("Maximum block length must be a positive integer. Given '%d'. Aborting...",
                      options.length)
        sys.exit(2)
    if options.lookahead < 0:
        logging.fatal("Overlapping block length must be a non-negative integer. Given '%d'. Aborting...",
                      options.lookahead)
        sys.exit(2)

    #  command-line related options
    if not ( '{pedigree}' in options.cmd and
             '{haplotypes}' in options.cmd and
             '{assumptions}' in options.cmd ):
        logging.fatal("The command-line *MUST* include the following placeholders: "
                      "'{pedigree}' '{haplotypes}' '{assumptions}'. "
                      "Given: '%s'. Aborting...",
                      options.cmd)
        sys.exit(2)
    if not '{number}' in options.cmdrec:
        logging.fatal("The options for handling recombinations *MUST* include "
                      "the following placeholder: '{number}'. "
                      "Given: '%s'. Aborting...",
                      options.cmdrec)
        sys.exit(2)
    if not '{min_number}' in options.cmdrec:
        logging.warn("The options do NOT include the placeholder '{min_number}' "
                     "for lower bounding the number of recombinations. "
                     "Given: '%s'",
                     options.cmdrec)

    if not '{time}' in options.cmdtime:
        logging.warn("The options do NOT include the placeholder '{time}'. "
                     "Time limit is *DISABLED*. "
                     "Given: '%s'",
                     options.cmdtime)

    #  other options
    if options.recomb_lb < -1 or options.recomb_ub < -1:
        logging.fatal("The initial bounds must be non-negative or equal to '-1' to disable. "
                      "Given: (%d, %d]", options.recomb_lb, options.recomb_ub)
        sys.exit(2)
    if options.recomb_ub != -1 and options.recomb_ub <= options.recomb_lb:
        logging.fatal("The initial upper bound must be greater than the initial lower bound. "
                      "Given: (%d, %d]", options.recomb_lb, options.recomb_ub)
        sys.exit(2)
    if options.time_limit < 0:
        logging.fatal("Running-time limit CANNOT be negative. Set to '0' to disable. "
                      "Given '%d'. Aborting...", options.time_limit)
        sys.exit(2)

    if options.recomb_ub != -1 and options.init_hc is not None:
        logging.fatal("Parameters '--initial-recomb-ub' and "
                      "'--initial-haplotype-configuration' cannot be used togheter.")
        sys.exit(2)
    if options.init_hc is not None and (
        not os.path.isfile(options.init_hc) or
        not os.access(options.init_hc, os.R_OK) ):
        logging.fatal("Initial haplotype configuration file not accessible. "
                      "Given '%s'. Aborting...",
                      options.init_hc)
        sys.exit(2)

    if options.bootstrap_time_limit <= 0:
        logging.fatal("Bootstrap time limit CANNOT be negative or '0'. "
                      "Given '%d'. Aborting...", options.bootstrap_time_limit)
        sys.exit(2)

    if options.bootstrap and options.recomb_lb != -1:
        logging.fatal("Parameters '--initial-recomb-lb' and "
                      "'--bootstrap' cannot be used togheter.")
        sys.exit(2)



def read_haplotypes(haplotypes_filename):
    logging.debug("Reading the haplotype configuration from file '%s'...",
                  haplotypes_filename)
    haplotypes= {}
    with open(haplotypes_filename, 'r') as haplotypes_file:
        for r in haplotypes_file:
            if r.startswith('#'):
                continue
            split_row= r.strip().split('\t', 6)
            comb_hap= split_row[6].split('\t')
            pat_hap= [ int(single.split('|')[0]) for single in comb_hap ]
            mat_hap= [ int(single.split('|')[1]) for single in comb_hap ]
            haplotypes[split_row[1]]= (pat_hap, mat_hap)
    return haplotypes



def integrate_hc(current_haplotypes, solution, good_chunk_length):
    logging.debug("Processing the haplotype configuration computed for "
                  "the current chunk...")
    for ind in current_haplotypes:
        (pat_hap, mat_hap)= current_haplotypes[ind]
        if ind not in solution.haplotypes:
            solution.haplotypes[ind]= [ pat_hap[0:good_chunk_length],
                                        mat_hap[0:good_chunk_length] ]
        else:
            solution.haplotypes[ind][0].extend( pat_hap[1:good_chunk_length] )
            solution.haplotypes[ind][1].extend( mat_hap[1:good_chunk_length] )



def compute_number_of_recombinations(h, ph, gps= None):
    if gps is None:
        nr1= compute_number_of_recombinations(h, ph, 0)
        nr2= compute_number_of_recombinations(h, ph, 1)
        return min(nr1, nr2)
    else:
        lh= len(h)
        nr= 0
        phase= gps
        for (l, a, a1, a2) in zip(range(lh), h, ph[0], ph[1]):
            av= (a1, a2)
            if a != av[phase]:
                phase= 1-phase
                assert a == av[phase]
                nr= nr+1
        return nr


def process_partial_hc(current_haplotypes, solution, fixed_index):
    tot_rec= 0
    for ind in solution.individuals:
        assert ind.id in current_haplotypes
        if ind.has_father():
            # father is present
            known_gps= None
            if ( fixed_index is not None and
                 solution.gps[ind.id][0] ):
                known_gps= solution.gps[ind.id][0][fixed_index]
            nr= compute_number_of_recombinations(current_haplotypes[ind.id][0],
                                                 current_haplotypes[ind.father],
                                                 known_gps)
            tot_rec= tot_rec + nr

        if ind.has_mother():
            # mother is present
            known_gps= None
            if ( fixed_index is not None and
                 solution.gps[ind.id][1] ):
                known_gps= solution.gps[ind.id][1][fixed_index]
            nr= compute_number_of_recombinations(current_haplotypes[ind.id][1],
                                                 current_haplotypes[ind.mother],
                                                 known_gps)
            tot_rec= tot_rec + nr

    logging.info("A partial haplotype configuration with %d recombinations has been found.",
                 tot_rec)
    return tot_rec

def basic_exec_reHCstar(filenames,
                        min_recombs, max_recombs,
                        cmd_templ, time_limit):
    rehcstar_success= False
    full_cmd_str= cmd_templ['cmd'].format(pedigree=    filenames['genotypes'],
                                          haplotypes=  filenames['haplotypes'],
                                          assumptions= filenames['assumptions'])
    cmd= []
    try:
        if max_recombs>0:
            full_cmd_str+= " " + cmd_templ['rec'].format(number= str(max_recombs), min_number= str(min_recombs+1))
        rem_time= time_limit.check_remaining_time()
        if rem_time is not None:
            full_cmd_str+= " " + cmd_templ['time'].format(time= str(int(rem_time)))

        cmd= shlex.split(full_cmd_str)
        logging.debug("Invoking >%s<...", " ".join(cmd))
        retcode = subprocess.call(cmd, shell=False)

        if retcode==4:
            (t_remain, t_limit, t_used)= time_limit.get_remaining_time()
            msg= "'reHCstar' has exceeded the maximum CPU time. (used {:.2f} s, given {:s} s)".format( t_used, str(t_limit) )
            logging.critical(msg)
            raise REHCstarMgrOutOfTimeError(msg)

        rehcstar_success= (retcode==0)

    except OSError as e:
        logging.warn("reHCstar execution failed. "
                     "Command-line: >%s<. "
                     "Additional information: '%s'",
                     " ".join(cmd), e)

    return rehcstar_success

def limited_exec_reHCstar(filenames,
                          min_recombs, max_recombs,
                          cmd_templ, time_limit):
    try:
        rehcstar_success= basic_exec_reHCstar(filenames, min_recombs, max_recombs,
                                              cmd_templ, time_limit)
        return "success" if rehcstar_success else "failure"
    except REHCstarMgrOutOfTimeError as e:
        return "out-of-time"

def step1_bootstrap(filenames, solution, init_n_rec, cmd, time_limit):
    logging.info("Step 1(bootstrap). An upper bound is given but a lower bound is not. "
                 "Trying to quickly find a lower bound.")
    bootstrap_ub= 0
    bootstrap_lb= -1
    bootstrap_end= bootstrap_ub >= init_n_rec
    successful_haplotypes_filename= None
    while not bootstrap_end:
        logging.info("Step 1(bootstrap). Trying with at most %d recombinations.",
                     bootstrap_ub)
        exec_rehcstar_status= limited_exec_reHCstar(filenames,
                                                    bootstrap_lb, bootstrap_ub,
                                                    cmd, time_limit)
        if exec_rehcstar_status == "success":
            logging.info("Step 1(bootstrap). Found a solution with at most %d recombinations.",
                         bootstrap_ub)
            successful_haplotypes_filename= filenames['haplotypes'] + "-success"
            os.rename(filenames['haplotypes'], successful_haplotypes_filename)
            bootstrap_end= True
        elif exec_rehcstar_status == "failure":
            logging.info("Step 1(bootstrap). Solution not found. Increasing the maximum "
                         "number of recombinations...")
            bootstrap_lb= bootstrap_ub
            bootstrap_ub= max(bootstrap_lb+1, int(3*(bootstrap_ub+1)/2-1))
            if bootstrap_ub >= init_n_rec:
                logging.info("Step 1(bootstrap). "
                             "The new maximum is greater than the given upper bound. "
                             "New maximum %d, given upper bound %d. Skipping bootstrap...",
                             bootstrap_ub, init_n_rec)
                bootstrap_end= True
        else: # exec_rehcstar_status == "out-of-time"
            bootstrap_end= True

    if successful_haplotypes_filename is not None:
        current_haplotypes= read_haplotypes(successful_haplotypes_filename)
        n_rec= process_partial_hc(current_haplotypes, solution, None)
        assert n_rec < init_n_rec
        logging.info("Step 1(bootstrap). "
                     "Found a haplotype configuration which improves the initial one. "
                     "New recombination upper bound: %d.", n_rec)
        return (True, bootstrap_lb, n_rec, current_haplotypes)
    elif bootstrap_lb >= 0:
        logging.info("Step 1(bootstrap). "
                     "No solution has been found but a lower bound has been computed. "
                     "New recombination lower bound: %d.", bootstrap_lb)
        return (False, bootstrap_lb, init_n_rec, None)
    else:
        logging.info("Step 1(bootstrap). "
                     "Bounds have not been improved. Continuing...")
        return (False, bootstrap_lb, init_n_rec, None)

def step1_regular(filenames, solution, min_recombs, max_recombs, out_upper_limit_recombinations, fixed_index, cmd, time_limit):
    logging.info("Step 1. Searching an upper bound on the number of recombinations...")
    step1_end= False
    out_of_time= False
    successful_haplotypes_filename= None
    upper_limit_recombinations = out_upper_limit_recombinations
    if not upper_limit_recombinations:
        upper_limit_recombinations = 2*len(solution.genotypes)*solution.genotype_length

    while not step1_end:
        logging.info("Step 1. Trying with at most %d recombinations.", max_recombs)
        exec_rehcstar_status= limited_exec_reHCstar(filenames, min_recombs, max_recombs,
                                                    cmd, time_limit)
        if exec_rehcstar_status == "success":
            logging.info("Step 1. Found a solution with at most %d recombinations.", max_recombs)
            successful_haplotypes_filename= filenames['haplotypes'] + "-success"
            os.rename(filenames['haplotypes'], successful_haplotypes_filename)
            step1_end= True
        elif exec_rehcstar_status == "failure":
            if max_recombs >= upper_limit_recombinations:
                logging.fatal("Step 1. Solution NOT found. The instance requires "
                              "too much recombinations (more than %d). Aborting...", max_recombs)
                step1_end= True
            else:
                logging.info("Step 1. Solution not found. Increasing the maximum "
                             "number of recombinations...")
                min_recombs= max_recombs
                max_recombs= min(2*(max_recombs+1)-1, upper_limit_recombinations)
        else:
            assert exec_rehcstar_status == "out-of-time"
            logging.fatal("Step 1. Time-limit exceeded before computing an "
                          "upper-bound to the number of recombinations. Aborting...")
            step1_end= True
            out_of_time= True

    if successful_haplotypes_filename is not None:
        current_haplotypes= read_haplotypes(successful_haplotypes_filename)
        n_rec= process_partial_hc(current_haplotypes, solution, fixed_index)
        assert n_rec <= max_recombs
        logging.info("Step 1. "
                     "Found a haplotype configuration with %d recombinations.", n_rec)
        return (True, min_recombs, n_rec, current_haplotypes)
    else:
        logging.info("Step 1. No solution has been found.")
        return (False, min_recombs, max_recombs, None)


def exec_reHCstar(filenames, solution, chunk, bounds, upper_limit_recombinations, cmd, time_limit, initial_haplotypes, bootstrap):
    out_of_time= False
    min_recombs= -1
    max_recombs= 0
    (c_start, c_stop, c_stop_la)= chunk
    fixed_index= c_start if c_start > 0 else None
    logging.info("Trying to solve the instance in file '%s'.", filenames['genotypes'])
    ## Search and/or set the initial bounds
    if bounds[0] != -1:
        logging.info("Lower bound on the number of recombinations= %d (given).", bounds[0])
        min_recombs= bounds[0]
        max_recombs= min_recombs+1
    if bounds[1] != -1:
        logging.info("Upper bound on the number of recombinations= %d (given).", bounds[1])
        max_recombs= bounds[1]
    else:
        max_recombs= int(math.ldexp(1, math.frexp(max_recombs)[1]))-1
    assert min_recombs < max_recombs

    successful_haplotypes_filename= None
    chunk_info= {}
    solution.chunk_info[chunk]= chunk_info

    current_haplotypes= None
    if bounds[0] == -1 and (bounds[1] != -1 or initial_haplotypes) and bootstrap['enabled']:
        # A lower bound has not been provided.
        # Try to quickly compute a lower bound...
        (bs_status, bs_lb, bs_ub, bs_haplotypes)= \
            step1_bootstrap(filenames, solution, max_recombs, cmd, TimeLimit(bootstrap['time-limit']))
        assert min_recombs <= bs_lb and bs_ub <= max_recombs
        if bs_status: # New better solution
            initial_haplotypes= None
            current_haplotypes= bs_haplotypes
        min_recombs= bs_lb
        max_recombs= bs_ub

    if initial_haplotypes is not None:
        assert not current_haplotypes
        current_haplotypes= initial_haplotypes

    if current_haplotypes is None: # No initial haplotypes and not successful bootstrap phase
        (s1_status, s1_lb, s1_ub, s1_haplotypes)= \
            step1_regular(filenames, solution, min_recombs, max_recombs, upper_limit_recombinations,
                          fixed_index, cmd, time_limit)
        min_recombs= s1_lb
        max_recombs= s1_ub
        if not s1_status:
            chunk_info["optimum found"]= False
            chunk_info["status"]= "Failed Step 1 "\
                "(time limit or recombination limit has exceeded): no upper bound found"
            chunk_info["last tried upper bound"]= s1_ub
            chunk_info["known lb"]= s1_lb
        else:
            current_haplotypes= s1_haplotypes


    rehcstar_success= current_haplotypes is not None

    if current_haplotypes is not None:
        if max_recombs == 0:
            logging.info("Step 2. A solution without recombinations has been found. "
                         "It is optimal, thus bisection is skipped...")
            chunk_info["optimum found"]= True
            chunk_info["status"]= "Zero-recombinant chunk"
            chunk_info["known lb"]= -1
            chunk_info["known ub"]=  0
        else:
            # Read the haplotype configuration and compute the real upper-bound
            lb= min_recombs
            ub= max_recombs
            mid= ub
            assert lb < ub
            try:

                logging.info("Step 2. Performing a bisection of interval (%d, %d] to find an optimal solution...",
                             lb, ub)
                while lb+1 < ub:
                    assert lb < ub
                    mid= math.floor((lb+ub)/2)
                    logging.info("Step 2. Bisecting interval (%d-%d]", lb, ub)
                    logging.info("Step 2. Trying with at most %d recombinations.", mid)
                    rehcstar_success= basic_exec_reHCstar(filenames, lb, mid, cmd, time_limit)
                    if rehcstar_success:
                        logging.debug("Step 2. Found a solution with at most %d recombinations.", mid)
                        successful_haplotypes_filename= filenames['haplotypes'] + "-success"
                        os.rename(filenames['haplotypes'], successful_haplotypes_filename)
                        current_haplotypes= read_haplotypes(successful_haplotypes_filename)
                        ub= process_partial_hc(current_haplotypes, solution, fixed_index)
                        max_recombs= ub
                    else:
                        logging.debug("Step 2. Solution NOT found in this half interval. Using the other...")
                        lb= mid

                logging.info("Step 2. Found an optimal solution with at most %d recombinations...", max_recombs)
                chunk_info["optimum found"]= True
                chunk_info["status"]= "Bisection correctly terminated"
                chunk_info["known lb"]= lb
                chunk_info["known ub"]= ub
                rehcstar_success= True

            except REHCstarMgrOutOfTimeError as e:
                logging.critical("Step 2. Time-limit exceeded before computing the minimum "
                                 "number of recombinations. Continuing anyway with a suboptimal solution...")
                chunk_info["optimum found"]= False
                chunk_info["status"]= "Time limit exceeded during bisection step (Step 2)"
                chunk_info["known lb"]= lb
                chunk_info["known ub"]= ub
                chunk_info["last tried upper bound"]= mid
                out_of_time= True

    if current_haplotypes is not None:
        chunk_info["chunk included in solution"]= True
        chunk_info["chunk is initial"]= current_haplotypes is initial_haplotypes
        integrate_hc(current_haplotypes, solution, c_stop - c_start)
        solution.compute_best_solution()
        rehcstar_success= True
    else:
        chunk_info["chunk included in solution"]= False

    return (rehcstar_success, out_of_time)



def main(options):
    log_level= logging.DEBUG if options.verbose>0 else logging.INFO
    verbose=  options.verbose>0
    verbose1= options.verbose>0
    verbose2= options.verbose>1
    verbose3= options.verbose>2

    logging.basicConfig(level=log_level,
                        format='%(levelname)5.5s [%(relativeCreated)15d] (%(filename)30s:%(lineno)-4d) - %(message)s',
                        stream=sys.stderr)

    for lic_lin in rehcstar_license_msgs:
        logging.info(lic_lin)
    logging.info("Started at %s", time.asctime())

    # Check program options
    check_program_options(options)

    cmd= {}
    cmd['cmd']= options.cmd
    cmd['rec']= options.cmdrec
    cmd['time']= options.cmdtime

    bootstrap= {
        'enabled': options.bootstrap,
        'time-limit': options.bootstrap_time_limit
       }

    logging.info("CONFIGURATION:")
    logging.info("Input pedigree:       '%s'", options.pedigree)
    logging.info("Resulting haplotypes: '%s'", options.haplotypes)
    logging.info("Max. block length:     %4d", options.length)
    logging.info("Lookahead length:      %4d", options.lookahead)
    logging.info("Recombination bounds: (%4d, %4d]", options.recomb_lb, options.recomb_ub)
    time_limit= TimeLimit(options.time_limit)
    logging.info("Time limit (secs):     %s",
                 "{:4d}".format(time_limit.limit)
                 if time_limit.limit is not None
                 else "unlimited")

    upper_limit_recombinations = options.ulrecomb if options.ulrecomb >= 0 else None

    # Read the original pedigree
    complete_sol= Solution()
    complete_sol.read_genotyped_pedigree(options.pedigree)

    if options.init_hc is not None and complete_sol.genotype_length > options.length:
        logging.fatal("Option '--initial-haplotype-configuration' is not "
                      "compatible with automatic genotype block partitioning. "
                      "Specify a block length greater than the genotype length. "
                      "(Given block length: %d, read genotype length %d).  "
                      "Aborting...", options.length, complete_sol.genotype_length)
        sys.exit(2)

    if len(complete_sol.individuals) == 0 or complete_sol.genotype_length == 0:
        logging.warn("The input genotyped pedigree is empty (individuals=%d, genotype length=%d). "
                     "A haplotype configuration cannot be computed. Exiting...",
                     len(complete_sol.individuals), complete_sol.genotype_length)
        sys.exit(1)

    logging.info("Starting the computation of the haplotype configuration...")
    assumptions= []
    suffix="{}-{}-pid{}".format(os.path.basename(options.pedigree),
                                os.path.basename(options.haplotypes),
                                os.getpid())

    initial_haplotypes= None
    bounds=[options.recomb_lb, options.recomb_ub]
    if options.init_hc is not None:
        initial_haplotypes= read_haplotypes(options.init_hc)
        # Read the haplotype configuration and compute the real upper-bound
        bounds[1]= process_partial_hc(initial_haplotypes, complete_sol, None)

    assert bounds[0] < bounds[1] or (bounds[0] == bounds[1] and bounds[0] == -1)

    complete_sol.notes.extend([
            "COMMAND LINE: {}".format(options.cmdline if options.cmdline else ""),
            "START TIME: {}".format(time.asctime()),
            "INPUT PEDIGREE: '{}'".format(options.pedigree),
            "OUTPUT HAPLOTYPES: '{}'".format(options.haplotypes),
            "MAX BLOCK LENGTH: '{}'".format(options.length),
            "LOOKAHEAD LENGTH: '{}'".format(options.lookahead),
            "INITIAL RECOMBINATION BOUNDS: '({}, {}]'".format(*bounds),
            "INITIAL HAPLOTYPE CONFIGURATION: '{}'".format(options.init_hc),
            "BOOTSTRAP ENABLED: '{}'".format(options.bootstrap),
            "BOOTSTRAP TIME LIMIT: '{}s'".format(options.bootstrap_time_limit),
            "TIME LIMIT: '{}'".format("{}s".format(time_limit.limit)
                                      if time_limit.limit is not None
                                      else "unlimited") ])

    partial_solution_found= initial_haplotypes is not None
    for start in range(0, complete_sol.genotype_length, options.length):
        good_stop= min( complete_sol.genotype_length, start + options.length + 1 )
        stop=      min( complete_sol.genotype_length, good_stop + options.lookahead )

        logging.info("Considering pedigree chunk between the markers [%d-%d)...", start, stop)
        filenames= {}
        # Write genotype block if necessary
        if (start != 0) or (stop!=complete_sol.genotype_length):
            filenames['genotypes']= "tmp-pedigree-{}-{}-{}".format(start, stop, suffix)
            complete_sol.write_genotype_block(filenames['genotypes'], start, stop)
        else:
            filenames['genotypes']= options.pedigree

        # Write assumptions
        filenames['assumptions']= "tmp-assumptions-{}-{}-{}".format(start, stop, suffix)
        with open(filenames['assumptions'], 'w') as assumptions_file:
            assumptions_file.write("\n".join(assumptions))
            assumptions_file.write("\n")

        # Execute reHCstar
        chunk= (start, good_stop, stop)
        filenames['haplotypes']= "tmp-haplotypes-{}-{}-{}".format(start, stop, suffix)
        (rehcstar_success, out_of_time)= exec_reHCstar(filenames, complete_sol,
                                                       chunk, bounds,
                                                       upper_limit_recombinations,
                                                       cmd, time_limit,
                                                       initial_haplotypes, bootstrap)
        if rehcstar_success:
            logging.info("reHCstar has successfully computed %s solution on the chunk [%d-%d).",
                         "an optimal" if complete_sol.chunk_info[chunk]["optimum found"] else "a suboptimal",
                         start, stop)
            logging.info("The haplotype configuration has %d recombinations and %d errors so far.",
                         len(complete_sol.recombs), len(complete_sol.errors))
            partial_solution_found= True

        if out_of_time:
            logging.fatal("Time-limit exceeded. Finalizing...")
            if partial_solution_found:
                partial_sol_filename= options.haplotypes + ".part"
                assert len(complete_sol.chunk_info) > 0
                complete_sol.analyze_optimality()
                logging.info("Saving the partial haplotype configuration to file '%s'", partial_sol_filename)
                complete_sol.notes.append("END TIME: {}".format(time.asctime()))
                complete_sol.write_haplotypes(partial_sol_filename, verbose=verbose1)
            else:
                logging.warn("No partial solution found.")
            logging.fatal("reHCstar manager - aborted at %s", time.asctime())
            sys.exit(1)

        elif not rehcstar_success:
            # No solution of reHCstar
            logging.fatal("reHCstar did NOT computed a solution on the chunk [%d-%d]. Aborting...", start, stop-1)
            logging.fatal("reHCstar manager - aborted at %s", time.asctime())
            sys.exit(1)
        else:
            # Prepare for next chunk
            if verbose2:
                incremental_haplotypes_filename= "tmp-incremental-{}-{}-{}".format(start,
                                                                                   stop,
                                                                                   suffix)
                complete_sol.write_haplotypes(incremental_haplotypes_filename, good_stop, verbose3)
            assumptions= []
            if start + options.length < complete_sol.genotype_length:
                assumption_locus_alleles = set()
                for ind in complete_sol.individuals:
                    curr_g = complete_sol.genotypes[ind.id][good_stop-1]
                    if is_genotyped(curr_g):
                        assumption_locus_alleles.update(set(curr_g.split(' ')))
                assumption_locus_n_alleles = len(assumption_locus_alleles)

                for ind in complete_sol.individuals:
                    if assumption_locus_n_alleles<=2:
                        assumptions.append(" ".join( [ "p", ind.id, "0",
                                                       str(complete_sol.haplotypes[ind.id][0][good_stop-1]-1) ] ))
                        assumptions.append(" ".join( [ "m", ind.id, "0",
                                                       str(complete_sol.haplotypes[ind.id][1][good_stop-1]-1) ] ))
                    else:
                        pat_allele = complete_sol.haplotypes[ind.id][0][good_stop-1]-1
                        mat_allele = complete_sol.haplotypes[ind.id][1][good_stop-1]-1
                        for a in range(assumption_locus_n_alleles):
                            print(a==complete_sol.haplotypes[ind.id][0][good_stop-1]-1)
                            assumptions.append(" ".join( [ "pm", ind.id, "0", str(a),
                                                           "1" if (a == pat_allele) else "0" ] ))
                            assumptions.append(" ".join( [ "mm", ind.id, "0", str(a),
                                                           "1" if (a == mat_allele) else "0" ] ))
                    if ind.has_father():
                        assumptions.append(" ".join( [ "sp", ind.id, "0",
                                                       str(complete_sol.gps[ind.id][0][good_stop-1]) ] ))
                    if ind.has_mother():
                        assumptions.append(" ".join( [ "sm", ind.id, "0",
                                                       str(complete_sol.gps[ind.id][1][good_stop-1]) ] ))


    logging.info("Found a complete haplotype configuration with %d recombinations and %d errors.",
                 len(complete_sol.recombs), len(complete_sol.errors))
    complete_sol.analyze_optimality()

    logging.info("Writing the haplotype configuration to '%s'...", options.haplotypes)
    complete_sol.notes.append("END TIME: {}".format(time.asctime()))
    complete_sol.write_haplotypes(options.haplotypes, verbose=verbose1)

    logging.info("reHCstar manager - completed at %s", time.asctime())

def maincmd():
    options = parse_command_line()
    main(options)

if __name__ == "__main__":
    maincmd()
