#!/bin/env python3


'''

usage:
 predixcan predict [options] --gds=FILE --db=FILE

options:
 --gds=FILE        genotype gds file
 --db=FILE         prediction model db file
 --genes=FILE      a file with list of genes with co-ordinates
 --out=PREFIX      outname prefix [default: predixcan]
 --nojob           run in front end
 --cluster=NAME    cluster name  [default: minerva]
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
    md.main(arguments,['predict'])
