#!/bin/env python3


'''

usage:
 predixcan covmat [options] --gds=FILE --db=FILE --genes=FILE

options:
 --gds=FILE        genotype gds file
 --db=FILE         prediction model db file
 --genes=FILE      genes files
 --out=PREFIX      outname prefix [default: out]
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
    md.main(arguments,['covmat'])
