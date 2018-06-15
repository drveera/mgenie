#!/bin/env python3


'''

usage:
 predixcan train_module1 [options] --expr=FILE --out=NAME

options:
 --gds=FILE        genotype gds file  [default: merged.gds]
 --expr=FILE       expression file
 --bfile=NAME      plink base name  [default: merged]
 --module=FILE     file with list of module genes  [default: module.names.list]
 --pcutoff=NUMBER  p value threshold  [default: 0.0001]
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
    md.main(arguments,['train_module1'])
