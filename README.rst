by `Yuri Pirola <http://algolab.eu/pirola>`_

Started: September 27, 2010

Current release: **2.1.5** (March 21, 2016)

--------------

Introduction
------------

This program is based on a reduction of the *Haplotype Configuration
with Recombinations and Errors* problem to *Boolean Satisfiability*,
which is then solved by a SAT solver. A haplotype configuration is
finally recovered from the satisfying assignment.

The algorithm is described in the following papers:

Yuri Pirola, Gianluca Della Vedova, Stefano Biffani, Alessandra Stella,
and Paola Bonizzoni. *A fast and practical approach to genotype phasing
and imputation on a pedigree with erroneous and incomplete information*.
IEEE/ACM Transactions on Computational Biology and Bioinformatics
(2012). `Link <http://dx.doi.org/10.1109/TCBB.2012.100>`__

Yuri Pirola, Gianluca Della Vedova, Stefano Biffani, Alessandra Stella,
and Paola Bonizzoni. *A fast and practical approach to genotype phasing
and imputation on a pedigree with erroneous and incomplete information*.
In: Proc. of IEEE 2nd International Conference on Computational Advances
in Bio and Medical Sciences, ICCABS 2012.
`Link <http://dx.doi.org/10.1109/ICCABS.2012.6182643>`__

Download and Installation
-------------------------

reHC-\* is currently distributed only on source form. It has been
developed on Ubuntu Linux machines (10.04 and later) and has been tested
on both 32 and 64 bit. The program should work on (or should be easily
ported to) on MacOS X but has not been tested and it is not supported on
this operating system.

Dependencies
~~~~~~~~~~~~

-  Python (>= 2.7)
-  CMake (>= 2.8)
-  GNU make
-  Boost FileSystem, System, DateTime, ProgramOptions, IOStreams, and
   other include-only libraries (tested with 1.42)
-  Apache Log4cxx (tested with 0.10.0)

Automatic Download and Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way for having reHC-\* correctly installed on your machine
is through `PyPI <https://pypi.python.org/pypi>`_ with the command:

::

    $ pip install -v reHCstar

(Please be patient, because it can take some time to build the package.)

If ``pip`` is not available on your system (and you cannot install it
following `these
instructions <https://pip.pypa.io/en/latest/installing.html>`_, you can
manually download the reHC-\* package from
https://github.com/yp/reHCstar/tarball/master, unpack it in a directory
of your choice, and then build it with the command:

::

    $ python setup.py install

Alternatively, you can proceed with the manual download, compilation,
and installation as detailed below.

Manual Download
~~~~~~~~~~~~~~~

reHC-\* is developed on the ``yp/reHCstar`` Git repository hosted by
GitHub. The repository can be explored using the GitHub web interface at
https://github.com/yp/reHCstar.

The latest stable version of reHC-\* can be downloaded in either
`.tar.gz <https://github.com/yp/reHCstar/tarball/master>`_ or in
`.zip <https://github.com/yp/reHCstar/zipball/master>`_ format. Previous
stable releases can be downloaded from
https://github.com/yp/reHCstar/archives/master.

It is also possible to clone the entire repository using the following
command:

::

    $ git clone git://github.com/yp/reHCstar.git

Or, if you have a GitHub account, you can fork the project from the
`repository web page <https://github.com/yp/reHCstar>`_.

Manual Compilation
~~~~~~~~~~~~~~~~~~

The program can be compiled by issuing the command at the command
prompt:

::

    $ make STATUS=Release

The program can be compiled in three different variants:

1. **with** an integrated SAT solver and **without** the ability to
   invoke an external SAT solver. This is the default variant, and
   should be the most efficient on both time and memory.
2. **with** an integrated SAT solver and **with** the ability to invoke
   an external SAT solver. This is the most flexible variant, but it is
   also the variant that uses more time and memory than the others.
3. **without** an integrated SAT solver but **with** the ability to
   invoke an external SAT solver. This variant should be preferred when
   one wants to use an external SAT solver, since it is more memory- and
   time-efficient than the previous one.

Currently, the SAT solvers that can be directly integrated into reHC-\*
are `CryptoMiniSat v2.9.1 <http://gitorious.org/cryptominisat>`_ by Mate
Soos and `MiniSat v2.2.0 <http://www.minisat.se/MiniSat.html>`_ by
Niklas Een and Niklas Sorensson. reHC-\* can interact (in the 2nd and
3rd variants) with any SAT solver that follow the standard
`requirements <http://www.satcompetition.org/2004/format-solvers2004.html>`_
of the last SAT competitions. The variant can be specified by modifying
the file ``CMakeOptions.txt`` in the root directory. In particular, two
options have to be used to specify the variant:

-  ``INTEGRATE_SAT_SOLVER``, which specifies if the SAT solver should be
   integrated into reHC-\*;
-  ``DISABLE_EXTERNAL_SAT_SOLVERS``, which specifies if the invocation
   of external SAT solvers is allowed.

The following combinations are allowed:

-  ``INTEGRATE_SAT_SOLVER=ON`` and ``DISABLE_EXTERNAL_SAT_SOLVERS=ON``,
   (**default**), which corresponds to the **first** variant;
-  ``INTEGRATE_SAT_SOLVER=ON`` and ``DISABLE_EXTERNAL_SAT_SOLVERS=OFF``,
   which corresponds to the **second** variant;
-  ``INTEGRATE_SAT_SOLVER=OFF`` and
   ``DISABLE_EXTERNAL_SAT_SOLVERS=OFF``, which corresponds to the
   **third** variant.

The SAT solver that will be integrated (if ``INTEGRATE_SAT_SOLVER`` is
``ON``) can be specified by setting ``USE_CRYPTOMINISAT`` or
``USE_MINISAT`` to ``ON`` in the file ``CMakeOptions.txt``.

By default, reHC-\* generates SAT instances in a augmented CNF format
where the "extended clauses" can also be XORs of literals. The SAT
instance is generally smaller if XORs are allowed, represents better the
"internal structure" of the Boolean formula, and *should be used* if it
is supported by the SAT solver. The SAT solver embedded by default,
CryptoMiniSat, supports this functionality, while MiniSat does not. If
MiniSat is chosen instead of CryptoMiniSat by editing the file
``CMakeOptions.txt``, then XOR-clauses are automatically disabled. To
disable the augmented CNF format, reHC-\* *must be rebuilt* specifying
the preprocessor symbol ``AVOID_XOR_CLAUSES`` in the ``CXXFLAGS``. For
example, if the program is compiled in a bash shell in Linux, it
suffices the following command to enable the "pure" CNF format:

::

    $ CXXFLAGS="-DAVOID_XOR_CLAUSES" make STATUS=Release

Usage
-----

The program takes as input a genotyped pedigree (with missing genotypes)
and returns (if possible) a complete haplotype configuration with at
most *r* recombinations and *e* errors. (The file formats are described
below.) Depending on the variant that has been compiled, the program
works in four different modes that have to be specified on the command
line as program parameter:

1. ``--create`` (short form ``-1``), that, given a genotyped pedigree,
   creates the associated SAT instance. (Available only on variants *2*
   and *3*.)
2. ``--read`` (short form ``-2``), that, given a genotyped pedigree,
   reads a satisfying model of the associated SAT instance (if such a
   model exists) and computes the associated haplotype configuration.
   (Available only on variants *2* and *3*.)
3. ``--create-read`` (short form ``-3``), that, given a genotyped
   pedigree, creates the associated SAT instance, invokes the external
   SAT solver, reads a satisfying model of the SAT instance (if such a
   model exists), and computes the associated haplotype configuration.
   This mode essentially combines the previous two modes by
   automatically invoking the external SAT solver. (Available only on
   variants *2* and *3*.)
4. ``--solve-internal`` (short form ``-4``), that, given a genotyped
   pedigree, creates the associated SAT instance, uses the integrated
   SAT solver for solving the instance, and, if the SAT instance is
   satisfiable, computes the associated haplotype configuration.
   (Available only on variant *1*.)

The following options are used to specify the input/output files:

-  ``--pedigree`` (short form ``-p``), that specifies the file
   containing the genotyped pedigree (input file);
-  ``--sat`` (short form ``-s``), that specifies the file containing the
   SAT instance associated with the genotyped pedigree (output file);
-  ``--result`` (short form ``-r``), that specifies the file containing
   the results computed by the external SAT solver for the SAT instance
   associated with the genotyped pedigree (input file);
-  ``--haplotypes`` (short form ``-h``), that specifies the file that
   will contain the haplotype configuration of the genotyped pedigree
   computed by reHC-\* (output file);
-  ``--assumptions`` (short form ``-a``), that specifies an *optional*
   file that contains additional assumptions that *must* be satisfied by
   the resulting haplotype configuration. Assumptions are specified one
   for each row with the following syntax:

   ::

       <variable kind> <individual id> <locus> <value>

Where ``<variable kind>`` is one of ``sp`` (paternal source), ``sm``
(maternal source), ``p`` (paternal allele), ``m`` (maternal allele),
``rp`` (paternal recombination), ``rm`` (maternal recombination), and
``e`` (genotyping error), ``<individual id>`` is the numerical
identifier of the individual (1-based), ``<locus>`` is the genotype
locus, and ``<value>`` is the boolean value (0/1) that the variable must
have. Please note that biallelic and multi-allelic loci are treated
differently, thus they have different set of variables.

For the ``--create-read`` mode, the command-line that has to be used to
invoke the external SAT must be specified by using the ``--sat-cmdline``
(short form ``-c``) program option. The strings ``%%INPUT%%`` and
``%%OUTPUT%%`` are placeholders for, respectively, the input and the
output files of the SAT solver. If the SAT solver can read the SAT
instance from its standard input, then it is possible to write the SAT
instance to the solver's standard input by specifying the option
``--pipe``. In this case, the placeholder ``%%INPUT%%`` will *not* be
used.

Options for Recombinations and Errors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, reHC-\* search for a haplotype configuration with zero
recombinations and zero errors. To enable recombinations in the
haplotyping process, the program options ``--global-recomb`` and either
``--global-recomb-rate=XX`` or ``--global-recomb-number=YY`` *must be
specified*. Here ``XX`` is a number between ``0.0`` and ``1.0`` that
represents the maximum number of recombinations *r* as a fraction of the
total number of possible recombination loci, while ``YY`` is (directly)
the maximum number of recombinations *r*. Moreover, if option
``--global-recomb`` is enabled and ``--global-recomb-number`` is used,
it is also possible to search for a haplotype configuration with a given
minimum number of recombinations by specifying the option
``--global-recomb-min-number=ZZ``, where ``ZZ`` is the sought lower
bound. This option should only be used to specify a lower bound that has
been already proved since the resulting haplotype configuration could
induce unnecessary recombination in order to satisfy the given lower
bound.

Similarly, to enable genotyping errors in the computed haplotype
configuration, the program options ``--global-error`` and either
``--global-error-rate=XX`` or ``--global-error-number=YY`` *must be
specified*. As before, ``XX`` is a number between ``0.0`` and ``1.0``
that represents the maximum number of errors *e* as a fraction of the
number of non-missing genotypes, while ``YY`` is (directly) the maximum
number of errors *e*.

Other program options allow a finer control over the distribution of
recombinations and errors. Please refer to the help of the program (that
can be obtained by specifying the ``--help`` program option) for their
presentation and explanation.

Other Options
~~~~~~~~~~~~~

reHC-\* can also read and write files compressed by GZip. The GZip
compression allows to save some space and, especially for large
instances and when an external SAT solver is used, it could reduce the
running time, since it greatly reduces to time spent for I/O operations.
It is disabled by default since not all the SAT solvers support it.
Three options regulates the GZip compression:

-  ``--compress-input``, which enables the GZip compression of some
   files that are read by reHC-\* (currently only the ``--pedigree``
   file);
-  ``--compress-output``, which enables the GZip compression of some
   files that are written by reHC-\* (currently the ``--sat`` and
   ``--haplotypes`` files);
-  ``--compress`` (short form ``-z``), which is equivalent to specify
   both ``--compress-input`` and ``--compress-output``;
-  ``--compress-sat``, which enables the GZip compression only for the
   file that contains the computed SAT instance.

Temporary files of the ``--create-read`` mode are automatically removed
by default. To keep them (for example, for manual inspection), the
program option ``--keep`` (short form ``-k``) has to be specified.

A summary of the available program options can be printed by invoking
reHC-\* with the ``--help`` (short form ``-?``) option.

Example
~~~~~~~

For example, if the genotyped pedigree is described in file
``genotyped-pedigree.txt``, the following commands perform the complete
haplotype inference process (saving the resulting haplotype
configuration in file ``haplotype-configuration.txt``).

Using the integrated SAT solver (variant *1* or *2*):

::

    $ ./bin/reHCstar -4  \
          -p genotyped-pedigree.txt  \
          -h haplotype-configuration.txt

Using an external SAT solver (variant *2* or *3*) with *manual*
invocation of the SAT solver:

::

    $ ./bin/reHCstar -1  \
          -p genotyped-pedigree.txt  \
          -s instance.cnf
    # ...execution of the external SAT solver, assuming that
    #    it writes the results in file sat-result.txt
    $ ./bin/reHCstar -2  \
          -p genotyped-pedigree.txt  \
          -r sat-result.txt  \
          -h haplotype-configuration.txt

Using an external SAT solver (variant *2* or *3*) with *automatic*
invocation of the SAT solver:

::

    $ ./bin/reHCstar -3  \
          -p genotyped-pedigree.txt  \
          -h haplotype-configuration.txt  \
          -c "./external-sat-solver %%INPUT%% %%OUTPUT%%"

Or, if the SAT solver reads the SAT instance from its standard input:

::

    $ ./bin/reHCstar -3  \
          -p genotyped-pedigree.txt  \
          -h haplotype-configuration.txt  \
          --pipe  \
          -c "./external-sat-solver %%OUTPUT%%"

Optimization Version
--------------------

reHC-\* also includes a program that uses the basic ``reHCstar``
executable in order to achieve two different aims:

-  finding (by a bisect-like search) the haplotype configuration that
   induces the minimum number of recombinations;
-  splitting long input genotypes into smaller overlapping blocks on
   which a partial haplotype configuration is computed independently and
   then used to reconstruct the complete haplotype configuration.

Please notice that the optimality of the solution (in term of number of
recombinations) is guaranteed if the genotypes are *not* split into
smaller blocks.

These functionalities are provided by the program ``reHCstar-mgr``
written in `Python <http://www.python.org>`_ version 3 and later.

``reHCstar-mgr`` requires two parameters, ``-p`` and ``-r``, that
specify, respectively, the file containing the input genotyped pedigree
and the file on which the computed haplotype configuration will be
saved.

By default, ``reHCstar-mgr`` invokes the ``reHCstar`` executable in the
current directory using the internal SAT solver mode (option
``--solve-internal`` described above). To change the default, the
complete command line must be provided as argument of the program option
``--cmd`` and must contain the following three placeholders
``{pedigree}``, ``{haplotypes}``, and ``{assumptions}`` that will be
replaced, respectively, with the input pedigree file, the output
haplotype configuration file, and the input additional assumption file.

For example, the default value of the ``--cmd`` option (i.e. the default
command line) is:

::

    ./reHCstar -4 -p "{pedigree}" -h "{haplotypes}" -a "{assumptions}"

The command line used to invoke the ``reHCstar`` executable is composed
by concatenating the argument of the previous option with the arguments
of two other options: ``--cmd-rec`` and ``--cmd-time``. The first one,
``--cmd-rec``, specifies the options (of ``reHCstar``) that regulates
the maximum (and, possibly, minimum) number of recombinations. In
particular, the argument must include the placeholder ``{number}`` which
will be replaced before invocation with the actual maximum number of
recombinations. Moreover, the argument may include the placeholder
``{min_number}`` which will be replaced before invocation with the
largest lower bound on the number of recombinations computed so far.

For example, the default value of the ``--cmd-rec`` option is:

::

    --global-recomb --global-recomb-number "{number}" --global-recomb-min-number "{min_number}"

The last option that regulates the final command line of ``reHCstar`` is
``--cmd-time`` and, if specified, must include the placeholder
``{time}`` which will be replaced before invocation with the maximum CPU
time of the ``reHCstar`` execution (in seconds). An empty argument
disables the running time limit control (albeit it could be enforced
anyway via OS services).

For example, the default value of the ``--cmd-time`` option is:

::

    --time-limit {time}

The following sections present the other main features of
``reHCstar-mgr`` while the full list of its options is available in the
integrated help (option ``-h``).

Automatic Genotype Partition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The subdivision of the input genotypes in (smaller) overlapping blocks
is regulated by the following two options: ``--block-length`` (short
form ``-l``, default ``50``) and ``--lookahead-length`` (short form
``-a``, default ``0``). The first option specifies the non-overlapping
(maximum) length of each block which the genotypes are divided into,
while the second option specifies the number of loci (in addition to a
single fixed locus) which two consecutive blocks overlap on. In other
words, a single block can be considered as composed by three parts: the
first part spans ``block-length`` loci, the second is composed by a
single locus, and the third (optional) part spans ``lookahead-length``
loci. (Hence, the total length of a block is ``block-length`` + ``1`` +
``lookahead-length``.) The second part of a block always overlaps with
the first locus of the first part of the next genotype block. Moreover
the haplotype configuration computed on this locus during the solution
of the "current" block is used as assumptions during the solution of the
next block (thus coincide). The third part of a block, the "look-ahead"
part, if it is present overlaps with the next block starting from its
second locus. This part is used to compute a haplotype configuration of
the "current" block, but the solution is then discarded when the next
block is considered (thus it may not coincide). Its purpose is to
provide a hint of the structure of the next block and it should be
particularly useful when the proportion of missing genotypes is
relevant, since when the overlapping locus has many missing genotypes,
the solution of the current block could impute the genotypes in a way
that is locally optimal, but globally sub-optimal.

Please notice that ``reHCstar-mgr`` finds a solution that requires the
minimum number of recombinations only if the genotypes are *not* divided
into blocks.

Initial Bounds on the Number of Recombinations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Initial lower and upper bounds on the number of recombinations may be
specified with the options ``--initial-recomb-lb=XX`` and
``--initial-recomb-ub=YY``, respectively. The options' arguments, ``XX``
and ``YY``, are non-negative numbers such that a haplotype configuration
with ``XX`` recombinations does not exist and a haplotype configuration
with ``YY`` recombinations certainly exists. The default value of both
of them is ``-1`` which means that no bound is known/provided. Moreover
it is possible to specify a file containing an initial haplotype
configuration that ``reHCstar-mgr`` tries to improve (in terms of number
of recombinations). In this case, the initial haplotype configuration is
read and the number of recombinations that it induces is used as initial
upper bound. If not better solution is found (for example, due to time
limits), then ``reHCstar-mgr`` outputs the initial haplotype
configuration. The file containing the initial haplotype configuration
is specified as argument of the ``--initial-haplotype-configuration``
program option. Please notice that options
``--initial-haplotype-configuration`` and ``--initial-recomb-ub`` cannot
be used together. These options could help to speed-up the process of
searching the solution with the minimum number of recombinations since
they provide the initial interval which the bisect-like search is
performed on.

If an initial upper bound is known but an initial lower bound is not, it
is possible to enable a *bootstrap* phase that attempts to quickly
identify an initial lower bound and then the execution continues by
bisecting the interval so determined. The bootstrap phase can be
activated by specifying the ``--bootstrap`` switch, while the maximum
CPU time spent in the bootstrap phase can be specified with the
``--bootstrap-time-limit=XX`` parameter, where ``XX`` is the time limit
expressed in seconds.

Running Time Management
~~~~~~~~~~~~~~~~~~~~~~~

``reHCstar-mgr`` provides basic tools for limiting its total running
time (CPU time). In particular, option ``--time-limit=SS`` specifies the
maximum running time of the program (``SS`` seconds). For the proper
functioning of this feature, the option ``--cmd-time`` must be valid. If
the program execution exceeds the given time limit, then
``reHCstar-mgr`` tries to save the solution computed so far in a file
whose name is the name specified by the option ``--results``
concatenated with the (fixed) extension ``.part``. The saved solution
could be *partial* (if the original instance has been partitioned in
blocks) and/or *suboptimal* (if the minimum number of recombinations has
not been computed within the time limit). The status of the solution is
saved as a comment line in the same file of the solution. We suggest to
enable the verbose mode (with ``-v`` or ``-vv``) for getting additional
information.

File Formats
------------

Input: Genotyped Pedigrees
~~~~~~~~~~~~~~~~~~~~~~~~~~

Genotyped pedigrees are described by a single file with the standard PED
format used in
`plink <http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped>`_.
In particular, each line of the pedigree file fully describes a single
individual and it is composed by at least *six* whitespace-separated
fields. The first (mandatory) six fields are:

-  ``Family ID`` (numeric only)
-  ``Individual ID`` (numeric only, greater than ``0``)
-  ``Paternal ID`` (the ID of the father, ``0`` if unknown/not present)
-  ``Maternal ID`` (the ID of the mother, ``0`` if unknown/not present)
-  ``Sex`` (``1`` = male, ``2`` = female)
-  ``Phenotype`` (ignored, could be any string not containing a
   whitespace)

**Remark:** reHC-\* currently works only on single-family pedigrees,
thus the ``Family ID`` *must be* the same for all the individuals.

The remaining fields (field 7 onwards) represent the genotype of the
individual, where each field represents a single allele of a single SNP
**biallelic** locus. Both the alleles of each locus *must be* specified
(they can be missing alleles), thus the total number of fields of each
row *must be* even. Major and minor alleles are encoded by the
characters ``1`` and ``2``. Missing genotypes are encoded by the pair
``0 0`` (i.e. by two fields containing the missing allele ``0``). The
pairs composed by a valid allele (``1`` or ``2``) and a missing allele
(``0``) *are not valid*. Since reHC-\* 2.0.0, there could also be
**multi-allelic** loci. Alleles are encoded by a number greater than
``0`` (which is always considered the missing allele code).

Rows starting with the character ``#`` are considered as comments and
ignored.

**Remark:** The order of the two alleles on each locus is meaningless
(i.e., the pair ``2 1`` is considered the same as the pair ``1 2``).

A simple single-family pedigree composed by 5 individuals genotyped over
5 biallelic loci is as follows.

::

    0 1 0 0 1 phenotype 1 1 2 2 2 2 2 2 1 1
    0 2 0 0 2 phenotype 2 2 1 1 1 1 1 1 1 1
    0 3 1 2 2 phenotype 1 2 0 0 1 2 1 2 1 1
    0 4 0 0 1 phenotype 1 2 1 2 1 1 1 1 0 0
    0 5 4 3 1 phenotype 1 2 1 2 0 0 1 1 1 2

Output: Haplotype Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The haplotype configuration computed by reHC-\* is represented in a
PED-like format. In particular, the first six fields are equal to the
PED format. The remaining fields represent the computed haplotype pair
of the individual, where *each* field represents the two alleles on a
single locus (separated by the character ``|``). In this case, the order
of the two alleles is important and represents the *phase* of each
locus. The first allele in each pair is the paternal allele, while the
second one is the maternal allele.

For the example, a zero-recombinant haplotype configuration for the
previous genotyped pedigree is as follows.

::

    0 1 0 0 1 phenotype 1|1 2|2 2|2 2|2 1|1
    0 2 0 0 2 phenotype 2|2 1|1 1|1 1|1 1|1
    0 3 1 2 2 phenotype 1|2 2|1 2|1 2|1 1|1
    0 4 0 0 1 phenotype 2|1 1|2 1|1 1|1 2|2
    0 5 4 3 1 phenotype 1|2 2|1 1|1 1|1 2|1

where the two (multi-locus) haplotypes of individual ``5`` are ``12112``
(paternal haplotype) and ``21111`` (maternal haplotype).

License
-------

reHC-\* is released under the terms of the GNU General Public License
(GPL) as published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

reHC-\* is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details.

Please refer to file ``COPYING`` or to the `GNU
website <http://www.gnu.org/licenses/>`_ for a copy of the GNU General
Public License.

Acknowledgments
---------------

The template of reHC-\* is based on the
`cpp-project-template <http://code.google.com/p/cpp-project-template/>`_
by Michael Aaron Safyan.

reHC-\* incorporates the following SAT solvers:

-  `CryptoMiniSat <http://gitorious.org/cryptominisat>`_ version 2.9.1
   (commit e819ab3236e, date 26/May/2011) by Mate Soos, which is
   distributed under the GNU General Public License version 3;
-  `MiniSat <http://www.minisat.se/MiniSat.html>`_ version 2.2.0 by
   Niklas Een and Niklas Sorensson, which is distributed under the MIT
   license.

For extracting source version information from git repository tags,
reHC-\* uses
`autorevision <https://github.com/Autorevision/autorevision>`_ by dak180
and others, which is distributed under the MIT license.

We would like to thank Gianluca Della Vedova for useful discussions.

Contacts
--------

Please contact *Yuri Pirola* for additional information.

E-mail: yuri.pirola@gmail.com

Web page: http://algolab.eu/pirola
