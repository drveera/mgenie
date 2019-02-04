
## Installation
1. install anaconda 
2. create a new enviornment with python 3.7 and docopt

`conda create -n genie python=3.7 docopt=0.6.2`

3. activate the enviornment

`conda activate genie`

4. install snakemake

`conda install -c bioconda snakemake`

5. clone this full repository

`git clone https://github.com/drveera/mgenie.git`

6. the main script you need to call is `genie.py`. Add an alias to it. In my case, i added a line,

`alias genie="/home/users/veera/softwares/mgenie/genie.py`

to my ~/.bashrc file. 

7. After adding alias, source the bashrc file and activate the genie environment again (remember: everytime you source the bashrc file, you need to activate the genie environment, unless you have a line, `source activate genie`, in your bashrc file)

8. then type, `genie --help` and you should see the help message, then its done. If you didn't get the help message, repeat the steps 1-8 again. 
