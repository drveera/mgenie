#!/bin/python

'''

usage:
 sumstats impute [options] --sumstats=FILE --ld=PLINK 

options:
 --sumstats=FILE     munged summary statistics file 
 --out=PREFIX        output name prefix [default: clump_out]
 --ld=PLINKFILE      plink file without extension
 --bed=BEDFILE       bed file  [default: |resources/sumstats/impute/fizi/locations.bed]
 --N=NUMBER          sample size 
 --nojob             if should run in front end
 --int
 --njobs=NUMBER      Number of parallel jobs when running in front end
 --dry-run           dry run snakemake
 --cluster=NAME      cluster name  [default: open]

'''

from docopt import docopt
import sys
sys.path.insert(1, sys.path[0] + '/../../../library')
import md
from md import process_list

arguments = docopt(__doc__)
if __name__ == '__main__':
    md.main(arguments,['impute'])
