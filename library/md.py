import sys
import os
from os.path import basename
import json
import re

methoddir = sys.path[0]
#here we assume that the method is always 2 levels down mainfolder, ex, ./module/method/.
##maindir = methoddir + "/../../../"
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


    try:
        clusterjson=f"{args['--cluster]}.cluster.json"
    except LookupError:
        clusterjson="minerva.cluster.json"
    cmds = f"snakemake -j {njobs} --use-conda --keep-going " \
          f"--cluster-config {maindir}/library/{clusterjson} " \
          f"--configfile {debugdir}/config.json " \
          f"--nolock " \
          f"-s {methoddir}/{method}.snake"

    if args['--dry-run']:
        os.system(cmds + " --dryrun")
    elif args['--nojob']:
        os.system(cmds)
    elif args['--int']:
        run_job(debugdir, args, cmds, interactive=True)
    else:
        run_job(debugdir, args, cmds, interactive=False)


def process_arguments(args):
    for arg, val in args.items():
        if isinstance(val, str) and '|' in val:
            args[arg] = maindir + '/' + val[1:]
    return args


def write_config(args, debugdir):
    try:
        clustername = args['--cluster']
        scriptcalls = f"{maindir}/library/scripts.{clustername}.json"
        ##if args['--cluster']=='minerva':
        ##    scriptcalls = f"{maindir}/library/minerva.json"
        ##else:
        ##    scriptcalls = f"{maindir}/library/openscripts.json"
    except LookupError:
        scriptcalls = f"{maindir}/library/scripts.open.json"
    try:
        with open(scriptcalls,'r') as s:
            script_files = json.load(s)
            script_files = script_files[0]
            args.update(script_files)
        with open(debugdir + '/config.json', 'w') as outfile:
            json.dump(args, outfile, indent=4)
    except FileNotFoundError:
        scriptcalls = f"{maindir}/library/scripts.open.json"
        with open(scriptcalls,'r') as s:
            script_files = json.load(s)
            script_files = script_files[0]
            args.update(script_files)
        with open(debugdir + '/config.json', 'w') as outfile:
            json.dump(args, outfile, indent=4)


def run_job(debugdir, args, scmds, interactive=False):
    outname=args['--out']
    try:
        ##edit added on 10th Dec 2018 to make it work both in minerva and genome.dk
        if args['--cluster']=='minerva':
            cmds = f"#!/bin/sh \n {scmds} " \
                   f"--jobname {outname}.{{rulename}}.{{jobid}} " \
                   f"--cluster 'bsub -P {{cluster.P}} -q {{cluster.q}} -W {{cluster.W}} -M {{cluster.M}} -o {debugdir}/out/{{cluster.output}} -n {{cluster.n}}' "
        elif args['--cluster']=="genomedk":
            cmds = f"#!/bin/sh \n {scmds} " \
                   f"--jobname {outname}.{{rulename}}.{{jobid}} " \
                   f"--cluster 'sbatch -e {debugdir}/error/slurm-%A_%a.out.error -o {debugdir}/out/slurm-%A_%a.out " \
                   f"--mem={{cluster.mem}} --time={{cluster.time}} -c {{cluster.cores}}' "
    except LookupError:
            cmds = f"#!/bin/sh \n {scmds} " \
                   f"--jobname {outname}.{{rulename}}.{{jobid}} " \
                   f"--cluster 'bsub -P {{cluster.P}} -q {{cluster.q}} -W {{cluster.W}} -M {{cluster.M}} -o {debugdir}/out/{{cluster.output}} -n {{cluster.n}}' "
    jobscript = debugdir + '/jobscript.sh'
    with open(jobscript, 'w') as outfile:
        outfile.write(cmds)
    if interactive:
        os.system(f"sh {jobscript}")
    else:
        try:
            if args['--cluster']=='minerva':
                os.system(f"bsub -n 1 -P acc_epigenAD -q premium -W 12:00 -M 4000 -o logfile < {jobscript}")
            elif args['--cluster']=='genomedk':
                os.system(f"sbatch --time=12:00:00 -e {debugdir}/master.error -o {debugdir}/master.out {jobscript}")
        except LookupError:
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
