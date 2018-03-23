#!/bin/env python

'''

usage:
 genie predixcan <method> [<args>...]

The following prs methods are available:

COMMAND  DESCRIPTION
train    train prediction models
predict  predict gene expressions
qc-dna   qc vcf files
qc-rna   qc gene expression matrix/covariates

type 'genie predixcan <method> --help' for more information on specific command
For example, type 'genie predixcan train --help'

'''

########################################################################################################
#IMPORTS
from subprocess import call
from docopt import docopt
import sys
########################################################################################################

########################################################################################################
#DOCOPT
if __name__ == '__main__':
    arguments = docopt(__doc__, options_first=True)
argv = [arguments['<method>']] + arguments['<args>']
########################################################################################################

########################################################################################################
# CHECK IF METHOD IS VALID
method = arguments['<method>']
gwasmethods = ['train','qc-dna','predict','-h','--help','qc-samples','metax','covar']
if method not in gwasmethods:
    exit(method + " is not valid")
########################################################################################################


########################################################################################################
# HELP MESSAGE
if method == '-h' or method == '--help':
    exit(print(__doc__))
########################################################################################################


########################################################################################################
# CALL THE METHOD
exit(call(['python', sys.path[0] + "/" + method + "/" + method + '.py'] + argv))
########################################################################################################
