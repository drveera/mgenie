#!/bin/env python3


'''

usage:
 predixcan train_module [options] --gds=FILE --expr=FILE --module=FILE --bfile=NAME

options:
 --gds=FILE        genotype gds file
 --bfile=NAME      plink base name
 --expr=FILE       expression file
 --module=FILE     file with list of module genes
 --out=PREFIX      outname prefix [default: predixcan_module]
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
    md.main(arguments,['train_module'])
