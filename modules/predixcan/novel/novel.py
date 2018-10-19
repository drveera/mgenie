#!/bin/env python3


'''

usage:
 predixcan novel [options] --gwas=FILE

options:
 --gwas=FILE       gwas summary statistics
 --annot=FILE      gencode annot  [default: |resources/predixcan/gencode.v27.build37.txt]
 --out=PREFIX      outname prefix [default: predixcan]
 --nojob           run in front end
 --dry-run         just show the codes
 --int             submit jobs from front end
 --njobs=NUMBER    number of parallel jobs; applicable only when running 
                    in front end

'''
from docopt import docopt
import sys
sys.path.insert(1, sys.path[0] + '/../../../library')
import md

arguments = docopt(__doc__)
if __name__ == '__main__':
    md.main(arguments,['novel'])
