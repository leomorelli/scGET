import gzip
import os
import re
import glob
import importlib.util
import scanpy as sc
import anndata

abs_path=os.path.abspath('')

# PARAMS NEEDED

CELL_NUMBER=config['cell_number']
#BC_PATTERN=config['bc_pattern']
SAMPLE= config['sample']
READS= config['reads']
THREADS=config['threads']
INPUT_PATH=config['input_path']
INPUT_LIST=config['input_list']
OUTPUT_PATH=config['output_path']
BC_PATTERN=config['bc_pattern']


samples_tot=[x.strip() for x in open(f'{INPUT_LIST}').read().split('\n')[:-1]]
SAMPLES=list(set([x.split(' ')[-1] for x in samples_tot if len(x)>0]))

wildcard_constraints:
    sample='.*[a-zA-Z0-9_]'


#1b) umi_tools
rule umi_tools:
    input:
        file='{output}/{sample}/{sample}_READ2.fastq.gz',
    params:
        method='reads',
        extract_method='string',
        bc_pattern=BC_PATTERN, #'CCCCCCCCCCCCCCCC',
        cell_number=CELL_NUMBER, # 5000 by default
        subset_reads='10000000000',
        output_path=OUTPUT_PATH
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 18000
    output:
        '{output}/{sample}/{sample}_whitelist.tsv',
        '{output}/{sample}/{sample}_cell_barcode_knee.png',
        '{output}/{sample}/{sample}_cell_barcode_counts.png'
    shell:
        """
        umi_tools whitelist --method={params.method} --extract-method={params.extract_method} --bc-pattern={params.bc_pattern} -I {input.file} -S {params.output_path}/{wildcards.sample}/{wildcards.sample}_whitelist_full.tsv --plot-prefix={params.output_path}/{wildcards.sample}/{wildcards.sample} --set-cell-number={params.cell_number} --subset-reads={params.subset_reads}
        awk '$1 !~ /N/' {params.output_path}/{wildcards.sample}/{wildcards.sample}_whitelist_full.tsv > {params.output_path}/{wildcards.sample}/{wildcards.sample}_whitelist.tsv
        """
