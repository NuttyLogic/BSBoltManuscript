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

sim = sys.argv[1]

alignment_dir = f'~/BSBolt_Paper/alignments/'



def format_run_info(run_info):
    run_time, run_stats = run_info
    reference, run_number = run_stats['description'].split('_x_')
    output = run_stats['output'].replace('.out', '')
    alignment_file = f'{output}.bam'
    if run_stats['tool'] =='bismark':
        file_base = output.split('/')[-1]
        if 'se' == sim[0:2]:
            alignment_file = f'{output}/sim_1_bismark_bt2.bam'
        else:
            alignment_file = f'{output}/sim_1_bismark_bt2_pe.bam'
    elif run_stats['tool'] == 'bsseeker':
        alignment_file = output
    sim_reference = f'{simulation_directory}{reference}'
    return dict(alignment_time=run_time, run_stats=run_stats, 
                reference=reference, run_number=run_number, 
                output=output, alignment_file=alignment_file, tool=run_stats['tool'],
                description=run_stats['description'], sim_reference=sim_reference)



home_dir = '~/BSBoltPaper/'
wd = '~/BSBoltPaper/indices/'
simulation_directory = '~/BSBoltPaper/simulated_reads/'

simulation_fasta = '~/BSBoltPaper/sim_genome.fa'

bismark_dir = '~/BSBoltPaper/Tools/Bismark-0.22.3/'
bsseeker_dir = '~/BSBoltPaper/Tools/BSseeker2-BSseeker2-v2.1.8/'
biscuit = '~/BSBoltPaper/Tools/biscuit/biscuit'
bwameth_tabulate = '~/BSBoltPaper/Tools/bwa-meth-0.2.2/scripts/tabulate-methylation.py' 


def time_func(func):
    
    def inner(*args, **kwargs):
        start = time.time()
        output = func(*args, **kwargs)
        run_time = time.time() - start
        return run_time, output
    return inner
        

@time_func
def run_alignment(alignment_command, alignment_tool, description, output_path):
    sam_aligners = {'bwa_meth', 'biscuit'}
    alignment = subprocess.Popen(alignment_command, stdout=subprocess.PIPE, universal_newlines=True)
    if alignment_tool in sam_aligners:
        # if output as sam pipe to samtools and output bam
        samtools_convert = ['samtools', 'view', '-b', '-']
        with open(f'{output_path}.bam', 'wb') as bam_out:
            bam_conversion = subprocess.Popen(samtools_convert,
                                              stdin=alignment.stdout,
                                              stdout=bam_out)
            bam_conversion.wait()
    else:
        with open(output_path, 'w') as out:
            for line in iter(alignment.stdout.readline, ''):
                out.write(line)
    alignment.wait()
    if alignment.returncode:
        return False
    return dict(cmd=alignment_command, tool=alignment_tool, description=description, output=output_path)



sim_dir = '/u/nobackup/mcdbscratch/colinpat/BSBolt_Paper/simulated_reads/'

simulation_parameters = {
'pe_directional_50': {'directory': f'{sim_dir}pe_directional_50', 'local':False, 'directional':True, 'read_length':50},
'pe_directional_100': {'directory': f'{sim_dir}pe_directional_100', 'local':False, 'directional':True, 'read_length':100},
'pe_directional_150': {'directory': f'{sim_dir}pe_directional_150', 'local':False, 'directional':True, 'read_length':150},
'pe_undirectional_50': {'directory': f'{sim_dir}pe_undirectional_50', 'local':False, 'directional':False, 'read_length':50},
'pe_undirectional_100': {'directory': f'{sim_dir}pe_undirectional_100', 'local':False, 'directional':False, 'read_length':100},
'pe_undirectional_150': {'directory': f'{sim_dir}pe_undirectional_150', 'local':False, 'directional':False, 'read_length':150},
'se_directional_50': {'directory': f'{sim_dir}se_directional_50', 'local':False, 'directional':True, 'read_length':50},
'se_directional_100': {'directory': f'{sim_dir}se_directional_100', 'local':False, 'directional':True, 'read_length':100},
'se_directional_150': {'directory': f'{sim_dir}se_directional_150', 'local':False, 'directional':True, 'read_length':150},
'se_undirectional_50': {'directory': f'{sim_dir}se_undirectional_50', 'local':False, 'directional':False, 'read_length':50},
'se_undirectional_100': {'directory': f'{sim_dir}se_undirectional_100', 'local':False, 'directional':False, 'read_length':100},
'se_undirectional_150': {'directory': f'{sim_dir}se_undirectional_150', 'local':False, 'directional':False, 'read_length':150},
'pe_low_coverage_directional_50': {'directory': f'{sim_dir}pe_low_coverage_directional_50', 'local':False, 'directional':True, 'read_length':50},
'pe_low_coverage_directional_100': {'directory': f'{sim_dir}pe_low_coverage_directional_100', 'local':False, 'directional':True, 'read_length':100},
'pe_low_coverage_directional_150': {'directory': f'{sim_dir}pe_low_coverage_directional_150', 'local':False, 'directional':True, 'read_length':150},
'se_low_coverage_directional_50': {'directory': f'{sim_dir}se_low_coverage_directional_50', 'local':False, 'directional':True, 'read_length':50},
'se_low_coverage_directional_100': {'directory': f'{sim_dir}se_low_coverage_directional_100', 'local':False, 'directional':True, 'read_length':100},
'se_low_coverage_directional_150': {'directory': f'{sim_dir}se_low_coverage_directional_150', 'local':False, 'directional':True, 'read_length':150},
'pe_error_directional_50': {'directory': f'{sim_dir}pe_error_directional_50', 'local':False, 'directional':True, 'read_length':50},
'pe_error_directional_100': {'directory': f'{sim_dir}pe_error_directional_100', 'local':False, 'directional':True, 'read_length':100},
'pe_error_directional_150': {'directory': f'{sim_dir}pe_error_directional_150', 'local':False, 'directional':True, 'read_length':150}}

with open(f'{home_dir}run_stats/{sim}_alignment_stats.pkl', 'rb') as stats:
    alignment_run_info = pickle.load(stats)

    
for alignment_run in alignment_run_info:
    formatted_info = format_run_info(alignment_run)
    sort_cmd = ['samtools', 'sort', '-@', '12', '-o', f'{formatted_info["output"]}.sorted.bam', 
                formatted_info['alignment_file']]
    index_cmd = ['samtools', 'index', f'{formatted_info["output"]}.sorted.bam']
    sort = subprocess.Popen(sort_cmd)
    sort.wait()
    index = subprocess.Popen(index_cmd)
    index.wait()

methylation_calling_defaults = {
    'BSBolt': [dict(output='-O', defaults=['python3', '-m', 'BSBolt', 'CallMethylation', '-DB', f'{wd}bsbolt_index', '-t', '12', 
                                           '-BQ', '10', '-MQ', '20', '-text', '-min', '5', '-CG'],
                   input_file='-I')],
    
    'bsseeker': [dict(output= '--CGmap', defaults= ['~/BSBoltPaper/Tools/Python-2.7.18rc1/python', f'{bsseeker_dir}bs_seeker2-call_methylation.py', 
                 '--txt', '--rm-overlap', '--txt', '--db', f'{wd}bsseeker_index/sim_genome_vector.fa_bowtie2'], input_file='-i')],
    
    'biscuit': [dict(output='-o', defaults=[biscuit, 'pileup', '-q', '20', '-b', '10', simulation_fasta]), 
                dict(defaults=[biscuit, 'vcf2bed', '-k', '5', '-t', 'cg'])],
    
    'bismark': [dict(output='-o', defaults=[f'{bismark_dir}bismark_methylation_extractor', '--no_overlap', '--bedGraph', 
                                            '--cytosine_report', '--parallel', '3', '--genome_folder', f'{wd}bismark_index'])],
    
    'bwa_meth': [dict(defaults=['python3', bwameth_tabulate, '--reference', simulation_fasta, '--base-q', '10', '--rf', '0', 
                                '--skip-left', '0', '--skip-right', '0'])]
    }


@time_func
def run_methylation_calling(run, calling_defaults):
    run_info = dict(run)
    tool = run_info['tool']
    meth_cmds = calling_defaults[tool]
    input_file = f'{run_info["output"]}.sorted.bam'
    output = f'{run_info["output"]}.meth_out'
    if tool == 'bismark':
        input_file = run_info['alignment_file']
        output = '/'.join(run_info['alignment_file'].split('/')[0:-1])
    for cmd_count, command in enumerate(meth_cmds):
        command_args = list(command['defaults'])
        if tool == 'bwa_meth':
            if run_info['reference'].split('_')[-1] == 'vec':
                command_args.extend(['--read-length', '150'])
            else:
                command_args.extend(['--read-length', run_info['reference'].split('_')[-1]])
        if 'output' in command:
            command_args.extend([command['output'], output])
        if 'input_file' in command:
            command_args.extend([command['input_file'], input_file])
        else:
            if cmd_count == 1:
                command_args.extend([output])
            else:
                command_args.extend([input_file])
        with open(f'{output}_{tool}_out', 'w') as std_out:
            with open(f'{output}_{tool}_err', 'w') as std_err:
                print(' '.join(command_args))
                meth_call = subprocess.Popen(command_args, stdout=std_out, cwd=alignment_dir, stderr=std_err)
                meth_call.wait()
    if tool == 'bismark':
        if 'se' in sim:
            output = f'{output}/sim_meth_1_bismark_bt2.CpG_report.txt'
        else:
            output = f'{output}/sim_meth_1_bismark_bt2_pe.CpG_report.txt'
    elif tool == 'BSBolt':
        output = f'{output}.CGmap'
    elif tool == 'bwa_meth':
        output = input_file.replace('.bam', '')
        output = f'{output}.methylation.txt'
    elif tool == 'biscuit':
        output = f'{output}_{tool}_out'
    run_info.update(dict(meth_output=output))
    return run_info



meth_calling_info = {}

for run in tqdm(alignment_run_info):
    formatted_run = format_run_info(run)
    meth_stats = run_methylation_calling(formatted_run, methylation_calling_defaults)
    meth_calling_info[f'{formatted_run["description"]}_{formatted_run["tool"]}'] = meth_stats


with open(f'{home_dir}run_stats/{sim}_meth_call_stats.pkl', 'wb') as out:
    pickle.dump(meth_calling_info, out)
