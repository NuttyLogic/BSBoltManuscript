#!/usr/bin/env python


import os
import pickle
import random
import subprocess
import sys
import time

sample = sys.argv[1]
tool = sys.argv[2]

def get_node_info():
    node_stats = {}
    commands = [['lscpu'], ['cat', '/proc/meminfo'], ['hostname'], ['cat',  '/etc/system-release']]
    for count, cmd in enumerate(commands):
        info = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
        for line in iter(info.stdout.readline, ''):
            if count < 2:
                cat, cat_stats = line.strip().split(':')
                node_stats[cat.strip()] = cat_stats.strip()
            elif count == 2:
                node_stats['hostname'] = line.strip()
            else:
                node_stats['os'] = line.strip()
    return node_stats 

base_dir = '~/BSBolt_Paper/real_data_comps/'
index_dir = f'{base_dir}indices/'
output_dir = f'~/BSBoltPaper/real_data_comps/alignments/'
fastq_dir = f'{base_dir}fastqs/'

# tool paths
bismark_dir = '~/BSBoltPaper/Tools/Bismark-0.22.3/'
bsseeker_dir = '~/BSBoltPaper/Tools/BSseeker2-BSseeker2-v2.1.8/'
biscuit = '~/BSBoltPaper/Tools/biscuit-release/biscuit'
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
    alignment = subprocess.Popen([x for x in alignment_command if x], stdout=subprocess.PIPE, universal_newlines=True)
    if alignment_tool in sam_aligners:
        # if output as sam pipe to samtools and output bam
        samtools_convert = ['samtools', 'view', '-@', '2', '-b', '-']
        with open(f'{output_path}.bam', 'wb') as bam_out:
            bam_conversion = subprocess.Popen(samtools_convert,
                                              stdin=alignment.stdout,
                                              stdout=bam_out)
            bam_conversion.wait()
    else:
        with open(f'{output_path}.log', 'w') as out:
            for line in iter(alignment.stdout.readline, ''):
                out.write(line)
    alignment.wait()
    if alignment.returncode:
        return False
    return dict(cmd=alignment_command, tool=alignment_tool, description=description, output=output_path)



alignment_defaults = {
    'BSBolt': dict(output='-O', undirectional=['-UN'], fastq1='-F1', fastq2='-F2', index=['-DB', f'{index_dir}bsbolt_index'], 
                  threads='-t', defaults=['python3', '-m', 'BSBolt', 'Align', '-OT', '2']),
    
    'bsseeker': dict(output='-o', undirectional=['-t', 'Y'], fastq1='-1', fastq2='-2', index=['-g', 'hg38_lambda.fa.gz', '-d', f'{index_dir}bsseeker_index'],
                    threads='--bt2-p', defaults=['~/BSBoltPaper/Tools/Python-2.7.18rc1/python', f'{bsseeker_dir}bs_seeker2-align.py', '--aligner', 'bowtie2', 
                                                f'--temp_dir={base_dir}' ]),
    
    'biscuit': dict(output='>', fastq1='', fastq2='', index=[f'{index_dir}biscuit_index/hg38_lambda.fa.gz'], 
                   threads='-t', defaults=[biscuit, 'align', '-b', '1']),
    
    'bismark': dict(output='-o', undirectional=['--non_directional'], fastq1='-1', fastq2='-2', index=[f'{index_dir}bismark_index'], 
                    threads='-p', defaults=[f'{bismark_dir}bismark', '--bam']),
    
    'bwa_meth': dict(output='>', threads='--threads', index=['--reference', f'{index_dir}bwa_index/hg38_lambda.fa.gz'], fastq1='', fastq2='', defaults=['python3', '-m', 'bwameth'])
    }
output_alignment = f'{output_dir}{sample}_{tool}'
alignment_args = alignment_defaults[tool]
alignment_base = list(alignment_args['defaults'])
alignment_base.extend([alignment_args['threads'], '12'])
  
if alignment_args['output'] != '>':
    alignment_base.extend([alignment_args['output'], output_alignment])
alignment_base.extend(alignment_args['index'])
if tool == 'bismark':
    alignment_base.extend(['--temp_dir', base_dir])
alignment_base.extend([alignment_args['fastq1'], f'{fastq_dir}{sample}_1.fastq.gz'])
alignment_base.extend([alignment_args['fastq2'], f'{fastq_dir}{sample}_2.fastq.gz'])

print(' '.join(alignment_base))


with open(f'~/BSBoltPaper/run_stats/{sample}_{tool}_alignment_stats.pkl', 'wb') as out:
    pickle.dump((None, {'cmd':' '.join(alignment_base)}), out)
    
node_stats = get_node_info()

with open(f'~/BSBoltPaper/run_stats/{sample}_{tool}_alignment_node_stats.pkl', 'wb') as out:
    pickle.dump(node_stats, out)

run_info = run_alignment(alignment_base, tool, sample, output_alignment)
    
with open(f'~/BSBoltPaper/run_stats/{sample}_{tool}_alignment_stats.pkl', 'wb') as out:
    pickle.dump(run_info, out)
    
