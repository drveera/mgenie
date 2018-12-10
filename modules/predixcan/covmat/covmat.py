#!/bin/env python3


'''

usage:
 predixcan covmat [options]

options:
 --gds=FILE        genotype gds file  [default: 1kg.gds]
 --db=FILE         prediction model db file  [default: module.db]
 --genes=FILE      genes files
 --out=PREFIX      outname prefix [default: out]
 --module          if computing for modules
 --cluster=NAME    minerva or genomedk  [default: minerva]
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
