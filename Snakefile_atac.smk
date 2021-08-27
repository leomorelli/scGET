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

spec_i = importlib.util.spec_from_file_location("get_info",f"{abs_path}/scripts/get_info.py")
infos = importlib.util.module_from_spec(spec_i)
spec_i.loader.exec_module(infos)

spec_l = importlib.util.spec_from_file_location("tn_layers",f"{abs_path}/scripts/tn_layers.py")
layers = importlib.util.module_from_spec(spec_l)
spec_l.loader.exec_module(layers)

# [B] CONFIG HANDLING
configfile: f'{abs_path}/config.yaml'

include: f'{abs_path}/modules/init_lane'  #initialization (handling of input files)
include: f'{abs_path}/modules/shared_rules' #umi_tools 
  
# [b.1] OTHER VARIABLES (no barcode in atac analysis)
CELL_NUMBER=config['cell_number']  
SAMPLE= config['sample']
READS= config['reads']
THREADS=config['threads']
INPUT_PATH=config['input_path']
INPUT_LIST=config['input_list']
SAMPLE_NAME=utilities.sample_name(SAMPLE,INPUT_PATH)
OUTPUT_PATH=utilities.output(config['output_path'],utilities.sample_name(SAMPLE,INPUT_PATH))

# [b.2] SELECTION OF OUTPUT MATRICES CHARACTERISTICS
binary_dictionary={True:'-B',False:''}
BINARY=binary_dictionary[config['binary']]

# [C] WORKFLOW
rule all:
    input:
        expand('{output}/{sample}.h5ad', output=OUTPUT_PATH,sample=SAMPLE_NAME)
       

#  create log dir
log_path = OUTPUT_PATH+"/logs_slurm"
try:
    os.mkdir(OUTPUT_PATH+"/logs_slurm")
    os.mkdir(log_path)
except OSError:
    try:
        os.makedirs(log_path)
        print("Successfully created the directory %s " % log_path)
    except OSError:
        print ("Creation of the directory %s failed" % log_path)
else:
    print ("Successfully created the directory %s " % log_path)


# 0) PREMERGE OF DIFFERENT INPUT FILES

# 0a) Create a text file where each file is assigned to its read (rule file_read)

#0b) Rename each file with a tag expressing the number of read (rule merge_reads)



#1) umi_tools (shared_rules)


#2) allignment
rule bwa:
    input:
        config['genome'], 
        '{output}/{sample}_READ1.fastq.gz',
        '{output}/{sample}_READ3.fastq.gz'
    params:
        center='COSR',
        platform='Illumina',
        prefix=SAMPLE_NAME,
        lib='not_specified',
        threads_bwa=THREADS-2,
        threads_samtools=THREADS-6
    resources:
        cpus=8,
        mem_mb=10000
    output:
        '{output}/{sample}_alligned.bam'
    shell:
        'bwa mem -R "@RG\\tID:{params.prefix}\\tPL:{params.platform}\\tPU:{params.prefix}\\tLB:{params.lib}\\tSM:{params.prefix}\\tCN:{params.center}" -t {params.threads_bwa} {input} |  samtools sort -T {output}_tmp -@ {params.threads_samtools} -o {output}'

#3) indexing allignment
rule index_allignment:
    input:
        '{output}/{sample}_alligned.bam'
    params:
        threads=THREADS
    resources:
        cpus=8,
        mem_mb=10000
    output:
        '{output}/{sample}_alligned.bam.bai'
    shell:
        'samtools index {input}'

# 4) deduplication
rule dedup:
    input:
        whitelist='{output}/{sample}_whitelist.tsv',
        read2='{output}/{sample}_READ2.fastq.gz',
        bamfile='{output}/{sample}_alligned.bam',
        indexed_bam='{output}/{sample}_alligned.bam.bai'
    resources:
        cpus=8,
        mem_mb=30000
    params:
        tn='atac',
        prefix=SAMPLE_NAME,
        scatACC_path=config['scatacc_path']
    output:
        '{output}/{sample}_bcdedup.bam'
    shell:
        'python {params.scatACC_path}/bc2rg.py -w {input.whitelist} -i {input.read2} -b {input.bamfile} -O -k -G {params.prefix}_{params.tn} | python {params.scatACC_path}/cbdedup.py -I -o {output}'
        

# 5) indexing dedup files
rule index_dedup:
    input:
        '{output}/{sample}_bcdedup.bam'
    params:
        threads=THREADS
    resources:
        cpus=8,
        mem_mb=24000
    output:
        '{output}/{sample}_bcdedup.bam.bai'
    shell:
        'samtools index {input}'


# 6) Peak_count
def spl(file):
        name=file
        return name[:-12]

rule peak_count:
    input:
        bai='{output}/{sample}_bcdedup.bam.bai',
        bam='{output}/{sample}_bcdedup.bam',
        bed=config['bed_file']
    params:
        out=spl('{output}/{sample}_bcdedup.bam'),
        threads=THREADS,
        scatACC_path=config['scatacc_path'],
        binary=BINARY
    resources:
        cpus=8,
        mem_mb=24000
    output:
        '{output}/{sample}.h5ad'
    shell:
        'python {params.scatACC_path}/para_count.py -p {input.bed} -b {input.bam} -o {params.out} -t {params.threads} {params.binary}'

