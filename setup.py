"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup
# To use a consistent encoding
import codecs
import os
import glob
import subprocess

here = os.path.abspath(os.path.dirname(__file__))

print("Generating version information...")
version = {}
with open("version.py", "w") as vf:
    subprocess.call("sh thirdparty/autorevision.sh -t py -o VERSION".split(" "), stdout=vf)
with open("version.py", "r") as vf:
    exec(vf.read(), version)
VERSION="{VCS_TAG}".format(**version).lstrip("v")
print("Extracted version: {0}".format(VERSION))

print("Building reHCstar...")
subprocess.call(["make"])


setup(
    name='reHCstar',
    version= VERSION,
    description='A SAT-based program to compute a haplotype configuration on pedigrees with recombinations, genotyping errors, and missing genotypes over biallelic and multi-allelic loci.',
    long_description=codecs.open(os.path.join(here, 'README.rst'), encoding='utf-8').read(),

    url='http://rehcstar.algolab.eu/',

    author='Yuri Pirola',
    author_email='yuri.pirola@gmail.com',
    license=codecs.open(os.path.join(here, 'COPYING'), encoding='utf-8').read(),

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2'
    ],

    keywords='haplotypes SNPs genotypes',

    zip_safe=False,

    packages=['rehcstar'],
    entry_points = {
        'console_scripts': [
            'reHCstar-mgr = rehcstar.reHCstar_mgr:maincmd'
        ]
    },
    data_files=[ (datadir, [f for f in glob.glob(os.path.join(here, datadir, '*'))])
                 for datadir in ('bin', 'docs') ]
)
