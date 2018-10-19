#!/bin/env python3


'''

usage:
 predixcan gcta [options] --grm=FILE --pheno=NAME

options:
 --grm=BASENAME    basename of grm file 
 --pheno=FILE      file or list of plink format phenotype files
 --out=PREFIX      outname prefix
 --nojob           run in front end
 --dry-run         just show the codes
 --njobs=NUMBER    number of parallel jobs; applicable only when running 
                    in front end
 --int 

'''
from docopt import docopt
import sys
sys.path.insert(1, sys.path[0] + '/../../../library')
import md

arguments = docopt(__doc__)
if __name__ == '__main__':
    md.main(arguments,['gcta'])
