import gzip
import os
import re
import glob
import importlib.util
import scanpy as sc
import anndata

abs_path=os.path.abspath('')

# [A] SCRIPTS IMPORT
spec_u = importlib.util.spec_from_file_location("utility_functions",f"{abs_path}/scripts/utility_functions.py")
utilities = importlib.util.module_from_spec(spec_u)
spec_u.loader.exec_module(utilities)

spec_l = importlib.util.spec_from_file_location("tn_layers",f"{abs_path}/scripts/tn_layers.py")
layers = importlib.util.module_from_spec(spec_l)
spec_l.loader.exec_module(layers)

# [B] CONFIG HANDLING  
# [b.1] OTHER VARIABLES (no barcode in atac analysis)

CELL_NUMBER=config['cell_number']  
SAMPLE= config['sample']
READS= config['reads']
THREADS=config['threads']
INPUT_PATH=config['input_path']
INPUT_LIST=config['input_list']
OUTPUT_PATH=config['output_path']

# [b.2] SELECTION OF OUTPUT MATRICES CHARACTERISTICS
binary_dictionary={True:'-B',False:''}
BINARY=binary_dictionary[config['binary']]

wildcard_constraints:
    sample='.*[a-zA-Z0-9_]'

# [C] WORKFLOW      

#1) alignment
rule bwa:
    input:
        config['genome'], 
        '{output}/{sample}/{sample}_READ1.fastq.gz',
        '{output}/{sample}/{sample}_READ3.fastq.gz'
    params:
        center='COSR',
        platform='Illumina',
        lib='not_specified',
        threads_bwa=THREADS-2,
        threads_samtools=THREADS-6
    resources:
        cpus=8,
        mem_mb=lambda wildcards, attempt: attempt  * 20000
    output:
        '{output}/{sample}/{sample}_aligned.bam'
    shell:
        'bwa mem -R "@RG\\tID:{wildcards.sample}\\tPL:{params.platform}\\tPU:{wildcards.sample}\\tLB:{params.lib}\\tSM:{wildcards.sample}\\tCN:{params.center}" -t {params.threads_bwa} {input} |  samtools sort -T {output}_tmp -@ {params.threads_samtools} -o {output}'

#2) indexing alignment
rule index_alignment:
    input:
        '{output}/{sample}/{sample}_aligned.bam'
    params:
        threads=THREADS
    resources:
        cpus=8,
        mem_mb=lambda wildcards, attempt: attempt  * 15000
    output:
        '{output}/{sample}/{sample}_aligned.bam.bai'
    shell:
        'samtools index {input}'

# 3) deduplication
rule dedup:
    input:
        whitelist='{output}/{sample}/{sample}_whitelist.tsv',
        read2='{output}/{sample}/{sample}_READ2.fastq.gz',
        bamfile='{output}/{sample}/{sample}_aligned.bam',
        indexed_bam='{output}/{sample}/{sample}_aligned.bam.bai'
    resources:
        cpus=8,
        mem_mb=lambda wildcards, attempt: attempt  * 100000
    params:
        tn='atac',
        scatACC_path=config['scatacc_path']
    output:
        '{output}/{sample}/{sample}_bcdedup.bam'
    shell:
        'python {params.scatACC_path}/bc2rg.py -w {input.whitelist} -i {input.read2} -b {input.bamfile} -O -k -G {wildcards.sample}_{params.tn} | python {params.scatACC_path}/cbdedup.py -I -o {output}'
        

# 4) indexing dedup files
rule index_dedup:
    input:
        '{output}/{sample}/{sample}_bcdedup.bam'
    params:
        threads=THREADS
    resources:
        cpus=8,
        mem_mb=lambda wildcards, attempt: attempt  * 20000
    output:
        '{output}/{sample}/{sample}_bcdedup.bam.bai'
    shell:
        'samtools index {input}'


rule bam_metrics:
    input:
        bai='{output}/{sample}/{sample}_bcdedup.bam.bai',
        bam='{output}/{sample}/{sample}_bcdedup.bam',
    params:
        scatACC_path=config['scatacc_path'],
    output:
        stats='{output}/{sample}/{sample}_B_stats.txt',
    shell:
        'python {params.scatACC_path}/bam_metrics.py {input.bam} > {output.stats}'

rule merge_metrics:
    input:
        st_file='{output}/{sample}/{sample}_B_stats.txt',
    output:
        stats='{output}/{sample}/{sample}_stats.txt',
    params:
        sample=lambda wildcards: {wildcards.sample}
    shell:
        'paste {input.st_file} | cut -f1,2,4,6,8,10,12,14,16 > {output.stats}'


# 5) Peak_count
def spl(file):
        name=file
        return name[:-12]

rule peak_count:
    input:
        bai='{output}/{sample}/{sample}_bcdedup.bam.bai',
        bam='{output}/{sample}/{sample}_bcdedup.bam',
        bed=config['bed_file']
    params:
        out=spl('{output}/{sample}/{sample}_bcdedup.bam'),
        threads=THREADS,
        scatACC_path=config['scatacc_path'],
        binary=BINARY
    resources:
        cpus=8,
        mem_mb=lambda wildcards, attempt: attempt  * 20000
    output:
        '{output}/{sample}/{sample}.h5ad'
    shell:
        'python {params.scatACC_path}/para_count.py -p {input.bed} -b {input.bam} -o {params.out} -t {params.threads} {params.binary}'

