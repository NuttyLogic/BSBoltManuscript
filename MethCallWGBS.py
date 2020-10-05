#!/usr/bin/env python



from collections import defaultdict, namedtuple
import os
import pickle
import random
import subprocess
import sys
import time
from typing import Dict, List, NamedTuple, Tuple

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pysam
import seaborn as sns
from tqdm.notebook import tqdm

sample_name = sys.argv[1]
tool = sys.argv[2]



home_dir = '~/BSBoltPaper/'
index_dir = '~/BSBoltPaper/real_data_comps/indices/'
alignment_dir = f'~/BSBoltPaper/real_data_comps/alignments/'
hg38 = '~/BSBoltPaper/real_data_comps/indices/hg38_lambda.fa'

bismark_dir = '~/BSBoltPaper/Tools/Bismark-0.22.3/'
bsseeker_dir = '~/BSBoltPaper/Tools/BSseeker2-BSseeker2-v2.1.8/'
biscuit = '~/BSBoltPaper/Tools/biscuit-release/biscuit'
bwameth_tabulate = '~/BSBoltPaper/Tools/bwa-meth-0.2.2/scripts/tabulate-methylation.py' 


# set alignment file path 
alignment_file = f'{alignment_dir}{sample_name}_{tool}'
output_file = f'{alignment_dir}{sample_name}_{tool}.meth_out'

if tool == 'bismark':
    alignment_file = f'{alignment_file}/{sample_name}_1_{tool}_bt2_pe.bam'
    output = alignment_file
elif tool != 'bsseeker':
    alignment_file = f'{alignment_file}.bam'
    
fixmate_alignment = f'{alignment_file.replace(".bam", "")}.fixmate.bam'
sorted_alignment = f'{alignment_file.replace(".bam", "")}.sorted.bam'
dup_alignment = f'{alignment_file.replace(".bam", "")}.dup.bam'
name_sorted = f'{alignment_file.replace(".bam", "")}.name.bam'


def time_func(func):
    def inner(*args, **kwargs):
        start = time.time()
        output = func(*args, **kwargs)
        run_time = time.time() - start
        return run_time, output
    return inner
        
fixmate_cmd = ['samtools', 'fixmate', '-p', '-m', alignment_file, fixmate_alignment]
sort_cmd = ['samtools', 'sort', '-@', '12', '-o', sorted_alignment, fixmate_alignment]
mark_cmd = ['samtools', 'markdup', sorted_alignment, dup_alignment]
index_cmd = ['samtools', 'index', dup_alignment]


def check_file(file_path):
    if not os.path.exists(file_path):
        return 1
    if os.path.getsize(file_path) < 100:
        return 1
    return 0

if check_file(dup_alignment):
    fix = subprocess.run(fixmate_cmd)
    sort = subprocess.run(sort_cmd)
    rm_fix = subprocess.run(['rm', fixmate_alignment])
    mark = subprocess.run(mark_cmd)
    rm_sort = subprocess.run(['rm', sorted_alignment])
    index = subprocess.run(index_cmd)
    
if tool == 'bismark':
    name = subprocess.run(['samtools', 'sort', '-n', '-@', '12', '-o', name_sorted, dup_alignment])

methylation_calling_defaults = {
    'BSBolt': [dict(output='-O', defaults=['python3', '-m', 'BSBolt', 'CallMethylation', '-DB', f'{index_dir}bsbolt_index', '-t', '12', '-BQ', '10', '-MQ', '20', '-text', '-min', '5', '-CG'],
                   input_file='-I')],
    
    'bsseeker': [dict(output= '--CGmap', defaults= ['~/BSBoltPaper/Tools/Python-2.7.18rc1/python', f'{bsseeker_dir}bs_seeker2-call_methylation.py', 
                 '--txt', '--rm-overlap', '--txt', '--db', f'{index_dir}bsseeker_index/hg38_lambda.fa.gz_bowtie2'], input_file='-i')],
    
    'biscuit': [dict(output='-o', defaults=[biscuit, 'pileup', '-q', '20', '-b', '10', hg38]), 
                dict(defaults=[biscuit, 'vcf2bed', '-k', '5', '-t', 'cg'])],
    
    'bismark': [dict(output='-o', defaults=[f'{bismark_dir}bismark_methylation_extractor', '--no_overlap', '--bedGraph', '--cytosine_report', '--parallel', '3', '--genome_folder', f'{index_dir}bismark_index'])],
    
    'bwa_meth': [dict(defaults=['python3', bwameth_tabulate, '--reference', hg38, '--base-q', '10', '--rf', '0', '--read-length', '100'])]
    }


@time_func
def run_methylation_calling(alignment_file, calling_defaults):
    meth_cmds = calling_defaults[tool]
    input_file = alignment_file
    output = f'{alignment_file}.meth_out'
    for cmd_count, command in enumerate(meth_cmds):
        command_args = list(command['defaults'])
        if 'output' in command:
            command_args.extend([command['output'], output])
        if 'input_file' in command:
            command_args.extend([command['input_file'], input_file])
        else:
            if cmd_count == 1:
                command_args.extend([output])
            else:
                command_args.extend([input_file])
        print(' '.join(command_args))
        with open(f'{output}_{tool}_out', 'w') as std_out:
            with open(f'{output}_{tool}_err', 'w') as std_err:
                meth_call = subprocess.Popen(command_args, stdout=std_out, cwd=alignment_dir, stderr=std_err)
                meth_call.wait()
    return command_args 

if tool == 'bismark':
    run_info = run_methylation_calling(name_sorted, methylation_calling_defaults)
else:
    run_info = run_methylation_calling(dup_alignment, methylation_calling_defaults)
    

with open(f'{home_dir}run_stats/{sample_name}_{tool}_meth_call_stats.pkl', 'wb') as out:
    pickle.dump(run_info, out)
