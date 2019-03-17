#!/bin/env python3


'''

usage: predixcan magma [options] --gwas=FILE --geneannot=FILE --ref=FILE

options:
 --gwas=FILE       gwas summary statistics
 --geneannot=FILE  geneannot file
 --ref=FILE        genotype reference file
 --out=PREFIX      outname prefix [default: magma]
 --nojob           run in front end
 --dry-run         just show the codes
 --int             submit jobs from front end
 --njobs=NUMBER    number of parallel jobs;applicable only when running in front end
 --cluster=NAME    cluster name  [default: open]

'''
from docopt import docopt
import sys
sys.path.insert(1, sys.path[0] + '/../../../library')
import md

arguments = docopt(__doc__)
if __name__ == '__main__':
    md.main(arguments,['magma'])
