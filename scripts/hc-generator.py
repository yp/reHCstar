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
#  hc-generator.py
#
#  A program that generates a random set of genotypes with
#  recombinations, errors and missing genotypes given a pedigree
#  structure.
#
##########

import sys
import random
import logging
import functools
from optparse import OptionParser

def parse_command_line():
    usage= "usage: %prog [options]"
    parser= OptionParser(usage=usage)
    parser.add_option("-l", "--genotype_length",
                      action="store", dest="length",
                      type="int", default=100,
                      help="the desired genotype length",
                      metavar="LENGTH")
    parser.add_option("-m", "--missing_genotype_probability",
                      action="store", dest="missing",
                      type="float", default=0.1,
                      help="the probability that a single-locus genotype has not been called",
                      metavar="PROBABILITY [0,1]")
    parser.add_option("-e", "--error_probability",
                      action="store", dest="error",
                      type="float", default=0.1,
                      help="the probability that a single-locus genotype has been mis-called",
                      metavar="PROBABILITY [0,1]")
    parser.add_option("-r", "--recombination_probability",
                      action="store", dest="recomb",
                      type="float", default=0.1,
                      help="the probability that a recombination event has occurred in each haplotype locus",
                      metavar="PROBABILITY [0,1]")
    parser.add_option("-s", "--seed",
                      action="store", dest="seed",
                      type="int", default=122295,
                      help="the seed of the random generator",
                      metavar="INT")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      default=False,
                      help="print additional log messages")
    (options, args) = parser.parse_args()
    return(options)

options = parse_command_line()

length= options.length
missing_genotype_prob= options.missing
error_prob= options.error
recomb_prob= options.recomb
seed= options.seed
log_level= logging.DEBUG if options.verbose else logging.INFO
allele1= "1"
allele2= "2"
homo1= "1 1"
homo2= "2 2"
heter= "1 2"
missing= "0 0"
pheno= "phenotype"

random.seed(seed)
logging.basicConfig(level=log_level,
                    format='%(levelname)-6s [%(asctime)s]  %(message)s')

logging.info("GENOTYPED PEDIGREE GENERATION WITH RECOMBINATIONS, ERRORS, AND MISSING GENOTYPES")

logging.info("Haplotype encoding: "
             "allele1= '%s', allele2= '%s'",
             allele1, allele2)
logging.info("Genotype encoding:  "
             "homo1= '%s', homo2= '%s', heter= '%s', missing= '%s'",
             homo1, homo2, heter, missing)
logging.info("Genotype length: %d", length)
logging.info("Missing genotype probability: %f", missing_genotype_prob)
logging.info("Genotyping error probability: %f", error_prob)
logging.info("Recombination probability:    %f", recomb_prob)
logging.info("Seed: %d", seed)



# Compute the mapping between haplotypes and genotypes
mapping={}
mapping[(allele1, allele1)]= homo1
mapping[(allele1, allele2)]= heter
mapping[(allele2, allele1)]= heter
mapping[(allele2, allele2)]= homo2

others={}
others[homo1]= (heter, homo2)
others[homo2]= (heter, homo1)
others[heter]= (homo1, homo2)


pedigree= set()
fathers= {}
mothers= {}
genders= {}
for f in sys.stdin:
    f= f.rstrip()
    (family, ind, father, mother, gender)= f.split()
    pedigree.add(ind)
    fathers[ind]= father
    mothers[ind]= mother
    genders[ind]= gender
    logging.debug("Read individual %4s.  Father=%4s  Mother=%4s  Gender= %s",
                  ind, father, mother, "male" if gender == 1 else "female")

logging.info("Pedigree size: %6d", len(pedigree))

hp={}
hm={}
gpp={} # Grand-parental sources
gpm={}
rp={} # Recombinations
rm={}
recombinations=[]
errors=[]

def acc_recomb(l, r):
    if len(l)==0:
        return [r]
    else:
        return l + [ (r+l[-1])% 2]

# Compute haplotypes of founders
founders=[ ind for ind in pedigree if fathers[ind] == "0" and mothers[ind] == "0" ]
logging.info("No. of founders: %4d", len(founders))
for ind in founders:
    logging.debug("Computing the haplotypes of founder %7s", ind)
    hp[ind]= [random.choice((allele1, allele2)) for i in range(length)]
    hm[ind]= [random.choice((allele1, allele2)) for i in range(length)]
    gpp[ind]= 0
    gpm[ind]= 0

genotyped= set(founders)
remaining= pedigree - genotyped

# Compute the haplotypes of the remaining individuals
logging.info("Not-haplotyped individuals: %4d", len(remaining))
while len(remaining)>0:
    for ind in remaining:
        if fathers[ind] in genotyped and mothers[ind] in genotyped:
            logging.debug("Computing the haplotypes of individual %4s", ind)
            genotyped.add(ind)
            gpp[ind]= random.randint(0,1)
            gpm[ind]= random.randint(0,1)
            rp[ind]= [gpp[ind]] + [ 1 if ((random.random() <= recomb_prob) &
                                          (hp[fathers[ind]][l] != hm[fathers[ind]][l]))
                                    else 0
                                    for l in range(length-1) ]
            rm[ind]= [gpm[ind]] + [ 1 if ((random.random() <= recomb_prob) &
                                          (hp[mothers[ind]][l] != hm[mothers[ind]][l]))
                                    else 0
                                    for l in range(length-1) ]
            recombinations= recombinations + [
                "REC {} {} 0".format(ind, l) for l in range(1,length) if rp[ind][l]==1 ] + [
                "REC {} {} 1".format(ind, l) for l in range(1,length) if rm[ind][l]==1 ]
            rp[ind]= functools.reduce(acc_recomb, rp[ind], [])
            rm[ind]= functools.reduce(acc_recomb, rm[ind], [])
            hp[ind]= [ hp[fathers[ind]][l] if rp[ind][l] == 0 else hm[fathers[ind]][l]
                       for l in range(length) ]
            hm[ind]= [ hp[mothers[ind]][l] if rm[ind][l] == 0 else hm[mothers[ind]][l]
                       for l in range(length) ]
    remaining -= genotyped
    logging.info("Not-haplotyped individuals: %4d", len(remaining))

# Computes the corresponding genotypes
logging.info("Computing the corresponding genotypes...")
tot_no_of_err = 0
max_err_x_ind = 0
genotypes={}
for ind in pedigree:
    genotypes[ind]= [mapping[(a,b)] for a,b in zip(hp[ind],hm[ind])]
    # Errors
    generr = [g if (random.random() >= error_prob)
              else random.choice(others[g])
              for g in genotypes[ind] ]
    curr_err = sum( [ a != b for a,b in zip(genotypes[ind],generr) ] )
    errors = errors + [ "ERR {}\t{}\t({}) -> ({})".format(ind, l, a, b)
                        for l,a,b in zip(range(length), genotypes[ind], generr)
                        if a != b ]
    tot_no_of_err = tot_no_of_err + curr_err
    max_err_x_ind = max(max_err_x_ind, curr_err)
    # Missing
    genmiss = [g if (random.random() >= missing_genotype_prob) else missing
               for g in generr ]

    genotypes[ind] = genmiss

# Print haplotype configuration as comment
logging.info("Saving the generated haplotype configuration...")
print("##ERRORS_IN_A_INDIVIDUAL\t{}\t{}".format(max_err_x_ind,max_err_x_ind/length))
print("##TOTAL_ERRORS\t{}\t{}".format(tot_no_of_err, tot_no_of_err/(length*len(genotypes))))
if len(errors)>0 :
    print("## ERR individual locus (original_genotype) -> (miscalled genotype)")
    print("\n".join(["# {}".format(f) for f in errors]))
if len(recombinations)>0 :
    print("## REC individual locus 0/1 (0=on paternal haplotype, 1=on maternal haplotype)")
    print("\n".join(["# {}".format(f) for f in recombinations]))
for ind in pedigree:
    print("# GENERATED_HAPLOTYPES", 0, ind, fathers[ind], mothers[ind], genders[ind], pheno,
          "\t".join(["{}|{}".format(a,b)
                      for a,b in zip(hp[ind], hm[ind]) ] ),
          sep="\t")

# Print the genotyped pedigree
logging.info("Saving the genotyped pedigree...")
for ind in pedigree:
    print(0, ind, fathers[ind], mothers[ind], genders[ind], pheno,
          "\t".join(genotypes[ind]), sep="\t")

logging.info("Terminated!")
