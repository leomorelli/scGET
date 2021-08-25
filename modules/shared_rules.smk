import gzip
import os
import re
import glob
import importlib.util
import scanpy as sc
import anndata

abs_path=os.path.abspath('')

# [B] CONFIG HANDLING
configfile: f'{abs_path}/config.yaml'

# [b.2] OTHER VARIABLES
CELL_NUMBER=config['cell_number']
SAMPLE= config['sample']
READS= config['reads']
THREADS=config['threads']
INPUT_PATH=config['input_path']
INPUT_LIST=config['input_list']
SAMPLE_NAME=utilities.sample_name(SAMPLE,INPUT_PATH)
OUTPUT_PATH=utilities.output(config['output_path'],utilities.sample_name(SAMPLE,INPUT_PATH))


#1b) umi_tools
rule umi_tools:
    input:
        '{output}/{sample}_READ2.fastq.gz'
    params:
        prefix=SAMPLE_NAME,
        method='reads',
        extract_method='string',
        bc_pattern='CCCCCCCCCCCCCCCC',
        cell_number=CELL_NUMBER, # 5000 by default
        subset_reads='10000000000',
        output_path=OUTPUT_PATH
    resources:
        mem_mb=35000
    output:
        '{output}/{sample}_whitelist.tsv',
        '{output}/{sample}_cell_barcode_knee.png',
        '{output}/{sample}_cell_barcode_counts.png'
    shell:
        'umi_tools whitelist --method={params.method} --extract-method={params.extract_method} --bc-pattern={params.bc_pattern} -I {input} -S {params.output_path}/{params.prefix}_whitelist.tsv --plot-prefix={params.output_path}/{params.prefix} --set-cell-number={params.cell_number} --subset-reads={params.subset_reads}'


