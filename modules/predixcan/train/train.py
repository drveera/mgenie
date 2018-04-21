#!/bin/env python3


'''

usage:
 predixcan train [options] --gds=FILE --expr=FILE

options:
 --gds=FILE        genotype gds file
 --expr=FILE       expression file
 --snpannot=FILE   snp annotation RDS file
 --genes=FILE      a file with list of genes with co-ordinates
 --out=PREFIX      outname prefix [default: predixcan]
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
    md.main(arguments,['train'])
