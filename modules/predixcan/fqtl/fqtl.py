#!/bin/env python3


'''

usage:
 predixcan fqtl [options] --vcf=FILE --expr=RDS

options:
 --vcf=FILE        vcf file
 --expr=RDS        expression rds file or list
 --oargs=ARGS      other arguments to pass within quotes
 --npeer=NUMBER    number of peer factors  [default: 10]
 --bed=annotation  bed file annotation for genes
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
    md.main(arguments,['fqtl'])
