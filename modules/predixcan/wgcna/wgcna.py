#!/bin/env python3


'''

usage:
 predixcan wgcna [options] --expr=RDS

options:
 --expr=RDS        expression rds file
 --out=PREFIX      outname prefix [default: WGCNA]
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
    md.main(arguments,['wgcna'])
