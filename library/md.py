import sys
import os
from os.path import basename
import json
import re

methoddir = sys.path[0]
#here we assume that the method is always 2 levels down mainfolder, ex, ./module/method/.
#maindir = methoddir + "/../../../"
#so we grep search genie instead
maindir = re.search(r'.*genie/',methoddir).group()

def main(args, methods):
    method = next(arg for arg in methods if arg in args and args[arg])

    args = process_arguments(args)
    outfolder = f"genie_{method}/{args['--out']}"
    pfix = f"{outfolder}/{args['--out']}"

    args[f"pfix_{method}"] = pfix
    args['--outfolder'] = outfolder

    debugdir = f"_debug/{outfolder}"

    os.makedirs(debugdir, exist_ok=True)
    os.makedirs(debugdir + '/out', exist_ok=True)
    os.makedirs(debugdir + '/error', exist_ok=True)
    os.makedirs(outfolder, exist_ok=True)

    write_config(args, debugdir)
    try:
        if args['--nojob']:
            if args['--njobs']:
                njobs = args['--njobs']
            else:
                njobs = 1
        else:
            njobs = 100000
    except LookupError:
        if args['--nojob']:
            njobs = 1
        else:
            njobs = 1


    cmds = f"snakemake -j {njobs} --use-conda --keep-going " \
          f"--cluster-config {maindir}/library/cluster.json " \
          f"--configfile {debugdir}/config.json " \
          f"--nolock " \
          f"-s {methoddir}/{method}.snake"

    if args['--dry-run']:
        os.system(cmds + " --dryrun")
    elif args['--nojob']:
        os.system(cmds)
    elif args['--int']:
        run_job(debugdir, args['--out'], cmds, interactive=True)
    else:
        run_job(debugdir, args['--out'], cmds, interactive=False)


def process_arguments(args):
    for arg, val in args.items():
        if isinstance(val, str) and '|' in val:
            args[arg] = maindir + '/' + val[1:]
    return args


def write_config(args, debugdir):
    try:
        if args['--cluster']=='minerva':
            scriptcalls = f"{maindir}/library/minervascripts.json"
        else:
            scriptcalls = f"{maindir}/library/openscripts.json"
    except LookupError:
        scriptcalls = f"{maindir}/library/openscripts.json"
    with open(scriptcalls,'r') as s:
        script_files = json.load(s)
    script_files = script_files[0]
    args.update(script_files)
    with open(debugdir + '/config.json', 'w') as outfile:
        json.dump(args, outfile, indent=4)


def run_job(debugdir, outname, scmds, interactive=False):

    cmds = f"#!/bin/sh \n {scmds} " \
           f"--jobname {outname}.{{rulename}}.{{jobid}} " \
           f"--cluster 'bsub -P {{cluster.P}} -q {{cluster.q}} -W {{cluster.W}} -M {{cluster.M}} -o {debugdir}/out/{{cluster.output}} -n {{cluster.n}}' "
    jobscript = debugdir + '/jobscript.sh'
    with open(jobscript, 'w') as outfile:
        outfile.write(cmds)
    if interactive:
        os.system(f"sh {jobscript}")
    else:
        os.system(f"bsub -n 1 -P acc_epigenAD -q premium -W 12:00 -M 4000 -o logfile < {jobscript}")


def process_list(argument):
    if ".list" in argument:
        d = {}
        argument_list = [x.strip() for x in list(open(argument))]
        argument_list = [x.split(" ") for x in argument_list]
        for i in argument_list:
            d[basename(i[0])] = i
        return(d)
    else:
        argument_list = argument.split(" ")
        d = {}
        d[basename(argument_list[0])] = argument_list
        return(d)

def flen(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
