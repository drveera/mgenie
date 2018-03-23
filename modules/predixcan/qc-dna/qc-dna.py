#!/bin/env python3


'''

usage:
 predixcan qc-dna [options] --vcf=FILE

options:
 --vcf=FILE        bgzipped and indexed vcf file or list
 --ref=FILE        bgzipped and indexed referrence vcf file
 --out=PREFIX      outname prefix [default: predixcan]
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
    md.main(arguments,['qc-dna'])
