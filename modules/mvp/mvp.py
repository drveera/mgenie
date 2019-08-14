#!/bin/env python

'''

usage:
 genie mvp <method> [<args>...]

The following prs methods are available:

COMMAND             DESCRIPTION
chunk-pgen-samples  chunk the pgen files sample wise
chunk-pgen-variants chunk the pgen files variant wise

type 'genie mvp <method> --help' for more information on specific command
For example, type 'genie mvp pgen --help'

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
gwasmethods = ['chunk-pgen-samples','chunk-pgen-variants',
               '-h','--help']

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
