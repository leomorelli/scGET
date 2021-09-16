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

# [b.1] DEFINITION OF TRANSPOSASES AND BARCODES FOR SCGET ANALYSIS
if config['tn5']!=True:
    TN_BARCODE=config['barcodes']
    TN_BARCODES={'tnh':TN_BARCODE['tnh']}
    BARCODES=TN_BARCODES['tnh']
elif config['tnh']!=True:
    TN_BARCODE=config['barcodes']
    TN_BARCODES={'tn5':TN_BARCODE['tn5']}
    BARCODES=TN_BARCODES['tn5']
else:
    TN_BARCODES=config['barcodes']
    BARCODES=TN_BARCODES['tn5']+TN_BARCODES['tnh']
    
# [b.2] OTHER VARIABLES
CELL_NUMBER=config['cell_number']  
SAMPLE= config['sample']
READS= config['reads']
THREADS=config['threads']
INPUT_PATH=config['input_path']
INPUT_LIST=config['input_list']
OUTPUT_PATH=config['output_path']

# [b.3] SELECTION OF OUTPUT MATRICES CHARACTERISTICS
binary_dictionary={True:'-B',False:''}
BINARY=binary_dictionary[config['binary']]

wildcard_constraints:
    sample='.*[a-zA-Z0-9_]'

# [C] WORKFLOW

def get_mem_mb_bwa(wildcards, attempt):
    return attempt * 20000

#1a) Classify each read using its barcode
rule tag_dust:
    input:
        expand('{output}/{{sample}}/{{sample}}_READ{read}.fastq.gz',output=OUTPUT_PATH, read=READS)
    params:
        output_path=OUTPUT_PATH,
        list_BCs=','.join(BARCODES),
        threads=THREADS
    resources:
        cpus=8,
        mem_mb=get_mem_mb_bwa
    output:
        expand('{output}/{{sample}}/{{sample}}_logfile.txt',output=OUTPUT_PATH),        
        expand('{output}/{{sample}}/{{sample}}_un_READ{read}.fq',output=OUTPUT_PATH, read=READS),        
        expand('{output}/{{sample}}/{{sample}}_BC_{barcodes}_READ{read}.fq',output=OUTPUT_PATH, barcodes=BARCODES,read=READS)
    shell:
        'tagdust -1 B:{params.list_BCs} -2 S:AGATGTGTATAAGAGACAG -3 R:N {input} -o {params.output_path}/{wildcards.sample}/{wildcards.sample} -t {params.threads}'


#1b) umi_tools (shared_rules)

#2a) compress files
rule compress:
    input:
        '{output}/{sample}/{sample}_BC_{barcode}_READ{read}.fq'
    output:
        '{output}/{sample}/{sample}_BC_{barcode}_READ{read}.fq.gz'
    shell:
        r'gzip {input} {output}'

#3a) alignement

rule bwa:
    input:
        config['genome'], 
        '{output}/{sample}/{sample}_BC_{barcode}_READ1.fq.gz',
        '{output}/{sample}/{sample}_BC_{barcode}_READ3.fq.gz'
    params:
        ids=utilities.exp_bc('merged_sample_BC_{barcode}_READ1.fq.gz'),
        center='COSR',
        platform='Illumina',
        prefix=get_mem_mb_bwa,
        lib='not_specified',
        threads_bwa=THREADS-2,
        threads_samtools=THREADS-6
    resources:
        cpus=8,
        mem_mb=get_mem_mb_bwa
    output:
        '{output}/{sample}/{sample}_BC_{barcode}.bam'
    wildcard_constraints:
        barcode="[A-Z]+"
    shell:
        'bwa mem -R "@RG\\tID:BC_{params.ids}\\tPL:{params.platform}\\tPU:{params.ids[1]}\\tLB:{params.lib}\\tSM:{params.prefix}\\tCN:{params.center}" -t {params.threads_bwa} {input} |  samtools sort -T {output}_tmp -@ {params.threads_samtools} -o {output}'

#4a) indexing alignement
rule index_alignement:
    input:
        '{output}/{sample}/{sample}_BC_{barcode}.bam'
    params:
        threads=THREADS
    resources:
        cpus=8,
        mem_mb=get_mem_mb_bwa
    output:
        '{output}/{sample}/{sample}_BC_{barcode}.bam.bai'
    wildcard_constraints:
        barcode="[A-Z]+"
    shell:
        'samtools index {input}'

# 5) deduplication
def tn_id(bc_in):
    if bc_in in config['barcodes']['tn5']:
        return 'tn5'        
    elif bc_in in config['barcodes']['tnh']:
        return 'tnh' 
    else:
        return 'ATAC'

rule dedup:
    input:
        whitelist='{output}/{sample}/{sample}_whitelist.tsv',
        read2='{output}/{sample}/{sample}_BC_{barcode}_READ2.fq.gz',
        bamfile='{output}/{sample}/{sample}_BC_{barcode}.bam',
        indexed_bam='{output}/{sample}/{sample}_BC_{barcode}.bam.bai'
    resources:
        cpus=8,
        mem_mb=get_mem_mb_bwa
    params:
        tn=lambda wildcards: tn_id(wildcards.barcode),
        prefix=lambda wildcards: {wildcards.sample},
        scatACC_path=config['scatacc_path']
    output:
        '{output}/{sample}/{sample}_BC_{barcode}_bcdedup.bam'
    shell:
        'python {params.scatACC_path}/bc2rg.py -w {input.whitelist} -i {input.read2} -b {input.bamfile} -O -k -G {params.prefix}_{params.tn} | python {params.scatACC_path}/cbdedup.py -I -o {output}'
        

# 6) indexing dedup files
rule index_dedup:
    input:
        '{output}/{sample}/{sample}_BC_{barcode}_bcdedup.bam'
    params:
        threads=THREADS
    resources:
        cpus=8,
        mem_mb=get_mem_mb_bwa
    output:
        '{output}/{sample}/{sample}_BC_{barcode}_bcdedup.bam.bai'
    shell:
        'samtools index {input}'


# 7) Peak_count
def spl(file):
        name=file
        return name[:-12]

rule peak_count:
    input:
        bai='{output}/{sample}/{sample}_BC_{barcode}_bcdedup.bam.bai',
        bam='{output}/{sample}/{sample}_BC_{barcode}_bcdedup.bam',
        bed=config['bed_file']
    params:
        out=spl('{output}/{sample}/{sample}_BC_{barcode}_bcdedup.bam'),
        threads=THREADS,
        scatACC_path=config['scatacc_path'],
        binary=BINARY
    resources:
        cpus=8,
        mem_mb=get_mem_mb_bwa
    output:
        '{output}/{sample}/{sample}_BC_{barcode}.h5ad'
    shell:
        'python {params.scatACC_path}/para_count.py -p {input.bed} -b {input.bam} -o {params.out} -t {params.threads} {params.binary}'

#8) single anndata with tnh and tn5 layers
rule layers:
    input:
        expand('{output}/{{sample}}/{{sample}}_BC_{barcode}.h5ad',output=OUTPUT_PATH, barcode=BARCODES)
    params:
        sample=lambda wildcards: {wildcards.sample}
    resources:
        cpus=8,
        mem_mb=get_mem_mb_bwa
    output:
        '{output}/{sample}/{sample}.h5ad'
    run:
        layers.add_layer(OUTPUT_PATH,{wildcards.sample}.pop(),TN_BARCODES)
