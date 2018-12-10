#!/bin/env python3


'''

usage:
 predixcan metax [options] --gwas=FILE 

options:
 --gwas=FILE       gwas sumstats or file list
 --db=FILE         prediction model db file or file list  [default: |resources/metaxcan/alldb.list]
 --genes=FILE      a file or a list of files with list of genes
 --brain           If you want only brain tissues to be analysed ignored if u provide own db list
 --out=PREFIX      outname prefix [default: predixcan]
 --cluster=NAME    minerva or genomedk  [default: minerva]
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
    md.main(arguments,['metax'])
