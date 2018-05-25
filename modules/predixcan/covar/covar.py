#!/bin/env python3


'''

usage:
 predixcan covar [options] --expr=FILE 

options:
 --expr=FILE       expression matrix file
 --thold=NUMBER    correlation value to filter  [default: 0.8]
 --cutoff=NUMBER   filtering cutoff 20 or 50 percent [default: 0.2]
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
    md.main(arguments,['covar'])
