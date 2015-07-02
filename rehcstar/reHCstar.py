#!/usr/bin/python

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

import sys
import pkg_resources
import subprocess
import rehcstar

## Scripts that allows to simulate reHCstar as it was in the PATH by forwarding
## the invocation to the real binary (packaged as data_file)

def main():
    cmdline = sys.argv
    cmdline[0] = pkg_resources.resource_filename('rehcstar', '../bin/reHCstar')
    subprocess.call(cmdline)

if __name__ == "__main__":
    main()
