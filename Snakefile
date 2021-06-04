import gzip
import os
import glob
import importlib.util

abs_path=os.path.abspath('')


spec_u = importlib.util.spec_from_file_location("utility_functions",f"{abs_path}/scripts/utility_functions.py")
utilities = importlib.util.module_from_spec(spec_u)
spec_u.loader.exec_module(utilities)

spec_i = importlib.util.spec_from_file_location("get_info",f"{abs_path}/scripts/get_info.py")
infos = importlib.util.module_from_spec(spec_i)
spec_i.loader.exec_module(infos)

configfile: f'{abs_path}/config.yaml'

if config['atac']==True:
    TN_BARCODES=config['barcodes']
    BARCODES=TN_BARCODES['tn5']
else:
    TN_BARCODES=config['barcodes']
    BARCODES=TN_BARCODES['tn5']+TN_BARCODES['tnh']
    
CELL_NUMBER=config['cell_number']  
SAMPLE= config['sample']
READS= config['reads']
THREADS=config['threads']
INPUT_PATH=config['input_path']
INPUT_LIST=config['input_list']
SAMPLE_NAME=utilities.sample_name(SAMPLE,INPUT_PATH)
OUTPUT_PATH=utilities.output(config['output_path'],utilities.sample_name(SAMPLE,INPUT_PATH))

binary_dictionary={True:'',False:'-c'}
BINARY=binary_dictionary[config['binary']]

rule all:
    input:
        expand('{output}/{sample}_BC_{barcode}.h5ad',output=OUTPUT_PATH, sample=SAMPLE_NAME, barcode=BARCODES),
       

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
# 0a) Create a text file where each file is assigned to its read

rule file_read:
    input:
        files=utilities.list_of_files(INPUT_PATH,INPUT_LIST),
    output:
        temp(expand('{output}/info_READ{read}.txt',output=OUTPUT_PATH,read=READS)),
    run:
        info_read=infos.df_info(input)
        infos.create_read_files(OUTPUT_PATH,info_read)


#0b) Rename each file with a tag expressing the number of read
rule merge_reads:
    input:
        '{output_path}/info_READ{read}.txt'
    output:
        '{output_path}/{sample}_READ{read}.fastq.gz'
    shell:
        "cat $(cat {input}) > {output}"


#1a) Classify each read using its barcode
rule tag_dust:
    input:
        expand('{output}/{sample}_READ{read}.fastq.gz',output=OUTPUT_PATH, sample=SAMPLE_NAME, read=READS)
    params:
        output_path=OUTPUT_PATH,
        prefix=SAMPLE_NAME,
        list_BCs=','.join(BARCODES),
        threads=THREADS
    resources:
        cpus=8,
        mem_mb=15000
    output:
        expand('{output}/{sample}_logfile.txt',output=OUTPUT_PATH, sample=SAMPLE_NAME),        
        expand('{output}/{sample}_un_READ{read}.fq',output=OUTPUT_PATH, sample=SAMPLE_NAME, read=READS),        
        expand('{output}/{sample}_BC_{barcodes}_READ{read}.fq',output=OUTPUT_PATH, sample=SAMPLE_NAME, barcodes=BARCODES,read=READS)
    shell:
        'tagdust -1 B:{params.list_BCs} -2 S:AGATATATATAAGGAGACAG -3 R:N {input} -o {params.output_path}/{params.prefix} -t {params.threads}'


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


#2a) compress files
rule compress:
    input:
        expand('{output}/{sample}_BC_{barcode}_READ{read}.fq',output=OUTPUT_PATH, sample=SAMPLE_NAME, read=READS,barcode=BARCODES)
    output:
        expand('{output}/{sample}_BC_{barcode}_READ{read}.fq.gz',output=OUTPUT_PATH, sample=SAMPLE_NAME, read=READS,barcode=BARCODES)
    shell:
        r'gzip {input} {output}'

#3a) allignement
rule bwa:
    input:
        config['genome'], 
        '{output}/{sample}_BC_{barcode}_READ1.fq.gz',
        '{output}/{sample}_BC_{barcode}_READ3.fq.gz'
    params:
        ids=utilities.exp_bc('merged_sample_BC_{barcode}_READ1.fq.gz'),
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
        '{output}/{sample}_BC_{barcode}.bam'
    shell:
        'bwa mem -R "@RG\\tID:BC_{params.ids}\\tPL:{params.platform}\\tPU:{params.ids[1]}\\tLB:{params.lib}\\tSM:{params.prefix}\\tCN:{params.center}" -t {params.threads_bwa} {input} |  samtools sort -T {output}_tmp -@ {params.threads_samtools} -o {output}'
#4a) indexing allignement

rule index_allignement:
    input:
        '{output}/{sample}_BC_{barcode}.bam'
    params:
        threads=THREADS
    resources:
        cpus=8,
        mem_mb=10000
    output:
        '{output}/{sample}_BC_{barcode}.bam.bai'
    wildcard_constraints:
        barcode="[A-Z]+"
    shell:
        'samtools index {input}'

# 5) deduplication
def tn_id(file):
    BCs=TN_BARCODES
    bc_in=utilities.exp_bc(file)
    if bc_in in BCs['tn5']:
        return 'tn5'        
    elif bc_in in BCs['tnh']:
        return 'tnh' 
    else:
        return 'nan'

rule dedup:
    input:
        whitelist='{output}/{sample}_whitelist.tsv',
        read2='{output}/{sample}_BC_{barcode}_READ2.fq.gz',
        bamfile='{output}/{sample}_BC_{barcode}.bam',
        indexed_bam='{output}/{sample}_BC_{barcode}.bam.bai'
    resources:
        cpus=8,
        mem_mb=30000
    params:
        tn=tn_id('{sample}_BC_{barcode}_READ2.fq.gz'),
        prefix=SAMPLE_NAME,
        scatACC_path=config['scatacc_path']
    output:
        '{output}/{sample}_BC_{barcode}_bcdedup.bam'
    shell:
        'python {params.scatACC_path}/bc2rg.py -w {input.whitelist} -i {input.read2} -b {input.bamfile} -O -k -G {params.prefix}_{params.tn}_{wildcards.barcode} | python {params.scatACC_path}/cbdedup.py -I -o {output}'
        

# 6) indexing dedup files
rule index_dedup:
    input:
        '{output}/{sample}_BC_{barcode}_bcdedup.bam'
    params:
        threads=THREADS
    resources:
        cpus=8,
        mem_mb=24000
    output:
        '{output}/{sample}_BC_{barcode}_bcdedup.bam.bai'
    shell:
        'samtools index {input}'


# 7) Peak_count
def spl(file):
        name=file
        return name[:-12]

rule peak_count:
    input:
        bai='{output}/{sample}_BC_{barcode}_bcdedup.bam.bai',
        bam='{output}/{sample}_BC_{barcode}_bcdedup.bam',
        bed=config['bed_file']
    params:
        out=spl('{output}/{sample}_BC_{barcode}_bcdedup.bam'),
        threads=THREADS,
        scatACC_path=config['scatacc_path'],
        binary=BINARY
    resources:
        cpus=8,
        mem_mb=24000
    output:
        '{output}/{sample}_BC_{barcode}.h5ad'
    shell:
        'python {params.scatACC_path}/peak_count.py -p {input.bed} -b {input.bam} -o {params.out} -A {params.binary}'


#8) naming files with tn tag

def tn_tag(file):
    bc=file.split('.')[0].split('_')[-1]
    if bc in TN_BARCODES['tn5']:
        return 'tn5'
    elif bc in TN_BARCODES['tnh']:
        return 'tnh'
    else:
        return 'nan'

#rule rename_h5ad:
#    input:
#        '{output}/{sample}_BC_{barcode}.h5ad'
#    params:
#        tn=tn_tag('{output}/{sample}_BC_{barcode}.h5ad')
#    output:
#        '{output}/{sample}_{params.tn}_{barcode}.h5ad'
#    shell:
#        'mv {input} {output}'
