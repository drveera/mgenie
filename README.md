
## Installation
1. install anaconda 
2. create a new enviornmnent with python 3.7 and docopt

`conda create -n genie python=3.7 docopt=0.6.2`

3. activate the enviornment

`conda activate genie`

4. install snakemake

`conda install -c bioconda snakemake`

5. clone this full repository

`git clone https://github.com/drveera/mgenie.git`

6. add alias to the script genie.py in the home folder of the repo to your bashfile, then when you call it

`genie --help`

you should see the help message, then its done. 
