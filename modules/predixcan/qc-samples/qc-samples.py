#!/bin/env python3


'''

usage:
 predixcan qc-samples [options] --bfile=FILE --ref=FILE

options:
 --bfile=FILE      plink file basename or list 
 --ref=FILE        referrence plink file (only basename)
 --out=PREFIX      outname prefix [default: predixcan]
 --skip-sex        skip sex checks
 --overlap         if the ref is in bfile 
 --cluster=NAME    cluster name  [default: open]
 --nojob           run in front end
 --dry-run         just show the codes
 --int             control job submission from frontend
 --njobs=NUMBER    number of parallel jobs; applicable only when running 
                    in front end

'''
from docopt import docopt
import sys
sys.path.insert(1, sys.path[0] + '/../../../library')
import md

arguments = docopt(__doc__)
if __name__ == '__main__':
    md.main(arguments,['qc-samples'])
