#!/bin/env python3


'''

usage:
 mvp gdsforpca [options] --pgen=FILE

options:
 --pgen=PREFIX        pgen file or .list of pgen files
 --nvariants=NUMBER   final number of variants  [default: 40000]
 --out=PREFIX         outname prefix [default: predixcan]
 --cluster=NAME       minerva or genomedk  [default: minerva]
 --nojob              run in front end
 --dry-run            just show the codes
 --int                submit jobs from front end
 --njobs=NUMBER       number of parallel jobs; applicable only when running in front end

'''
from docopt import docopt
import sys
sys.path.insert(1, sys.path[0] + '/../../../library')
import md

arguments = docopt(__doc__)
if __name__ == '__main__':
    md.main(arguments,['gdsforpca'])
