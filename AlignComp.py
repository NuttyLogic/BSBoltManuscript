#!/usr/bin/env python


import os
import pickle
import random
import subprocess
import sys
import time

sim = sys.argv[1]

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



wd = '~/BSBoltPaper/indices/'


simulation_fasta = os.getcwd() + '/sim_genome_vector.fa'

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

sim_dir = '~/BSBoltPaper/simulated_reads/'

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



alignment_defaults = {
    'BSBolt': dict(output='-O', undirectional=['-UN'], fastq1='-F1', fastq2='-F2', index=['-DB', f'{wd}bsbolt_index'], 
                  threads='-t', defaults=['python3', '-m', 'BSBolt', 'Align', '-M', '-T', '0'], read_length=True),
    
    'bsseeker': dict(output='-o', undirectional=['-t', 'Y'], fastq1='-1', fastq2='-2', index=['-g', 'sim_genome_.fa', '-d', f'{wd}bsseeker_index'], end_to_end=['--bt2--end-to-end'], 
                    threads='--bt2-p', defaults=['~/BSBoltPaper/Tools/Python-2.7.18rc1/python', f'{bsseeker_dir}bs_seeker2-align.py', '--aligner', 'bowtie2']),
    
    'biscuit': dict(output='>', fastq1='', fastq2='', index=[f'{wd}biscuit_index/sim_genome.fa'], 
                   threads='-t', defaults=[biscuit, 'align'], directional=['-b', '1']),
    
    'bismark': dict(output='-o', undirectional=['--non_directional'], fastq1='-1', fastq2='-2', index=[f'{wd}bismark_index'], 
                    threads='-p', defaults=[f'{bismark_dir}bismark', '--bam']),
    
    'bwa_meth': dict(output='>', threads='--threads', index=['--reference', f'{wd}bwa_index/sim_genome.fa'], fastq1='', fastq2='', defaults=['python3', '-m', 'bwameth'])
    }




alignment_dir = f'~/BSBoltPaper/alignments/'


alignment_commands = []

run = 0
for tool, alignment_args in alignment_defaults.items(): 
    sim_params = simulation_parameters[sim]
    sim_label = f'{sim}_x_{run + 1}'
    alignment_base = list(alignment_args['defaults'])
    alignment_base.extend([alignment_args['threads'], '12'])
    if sim_params['directional']:
        if 'directional' in alignment_args:
            alignment_base.extend(alignment_args['directional'])
    else:
        if 'undirectional' in alignment_args:
            alignment_base.extend(alignment_args['undirectional'])
    else:
        if 'end_to_end' in alignment_args:
            alignment_base.extend(alignment_args['end_to_end'])
    output = f'{alignment_dir}{sim_label}_x_{tool}.out'
    if alignment_args['output'] != '>':
        alignment_base.extend([alignment_args['output'], f'{alignment_dir}{sim_label}_x_{tool}'])
    else:
        output = f'{alignment_dir}{sim_label}_x_{tool}'
    alignment_base.extend(alignment_args['index'])
    if tool == 'bismark':
        alignment_base.extend(['--temp_dir', sim_params['directory']])
    fastq1_arg = alignment_args['fastq1']
    if not os.path.exists(f'{sim_params["directory"]}/sim_2.fq'):
        if tool == 'bsseeker':
            fastq1_arg = '-i'
        elif tool == 'bismark':
            fastq1_arg = ''
    alignment_base.extend([fastq1_arg, f'{sim_params["directory"]}/sim_1.fq'])
    if os.path.exists(f'{sim_params["directory"]}/sim_2.fq'):
        alignment_base.extend([alignment_args['fastq2'], f'{sim_params["directory"]}/sim_2.fq'])
    alignment_commands.append((alignment_base, tool, sim_label, output))


random.shuffle(alignment_commands)

print(alignment_commands)

alignment_run_info = []

for alignment in alignment_commands:
    run_info = run_alignment(*alignment)
    if not run_info[1]:
        print(' '.join(alignment[0]))
        print(alignment)
    alignment_run_info.append(run_info)
 
with open(f'~/BSBoltPaper/run_stats/{sim}_alignment_stats.pkl', 'wb') as out:
    pickle.dump(alignment_run_info, out)
    
node_stats = get_node_info()

with open(f'~/BSBoltPaper/run_stats/{sim}_alignment_node_stats.pkl', 'wb') as out:
    pickle.dump(node_stats, out)

