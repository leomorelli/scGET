import gzip
import glob
import importlib.util

spec = importlib.util.spec_from_file_location("utility_functions", "scripts/utility_functions.py")
utilities = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utilities)


configfile: 'config.yaml'

SAMPLE= config['sample']
READS= config['reads']
TN_BARCODES=config['barcodes']
BARCODES=TN_BARCODES['tn5']+TN_BARCODES['tnh']
GENOME=config['genome']
THREADS=config['threads']
CELL_NUMBER=config['cell_number']
INPUT_PATH=config['input_path']
INPUT_LIST=config['input_list']
OUTPUT_PATH=utilities.output(config['output_path'])
SAMPLE_NAME=utilities.sample_name(SAMPLE,INPUT_PATH)


rule all:
    input:              
        expand('{sample}_{tn}_merged.bam',sample=SAMPLE_NAME, tn=TN_BARCODES.keys())


# 0) PREMERGE OF DIFFERENT INPUT FILES
# 0a) Create a text file where each file is assigned to its read

rule file_read:
    input:
        files=utilities.list_of_files(INPUT_PATH,INPUT_LIST)
    output:
        expand('info_READ{read}.txt', read=READS)
    script:
        'scripts/get_info.py'


#0b) Rename each file with a tag expressing the number of read
rule merge_reads:
    input:
        'info_READ{read}.txt'
    output:
        '{sample}_READ{read}.fastq.gz'
    shell:
        "for F in $(cat {input}) ; do cat $F>>{output}; done"


#1a) Classify each read using its barcode
rule tag_dust:
    input:
        expand('{sample}_READ{read}.fastq.gz', sample=SAMPLE_NAME, read=READS)
    params:
        prefix=SAMPLE_NAME,
        list_BCs=','.join(BARCODES),
        threads=THREADS
    output:
        expand('{sample}_logfile.txt',sample=SAMPLE_NAME),        
        expand('{sample}_un_READ{read}.fq', sample=SAMPLE_NAME, read=READS),        
        expand('{sample}_BC_{barcodes}_READ{read}.fq', sample=SAMPLE_NAME, barcodes=BARCODES,read=READS)
    shell:
        'tagdust -1 B:{params.list_BCs} -2 S:AGATATATATAAGGAGACAG -3 R:N {input} -o {params.prefix} -t {params.threads}'


#1b) umi_tools
rule umi_tools:
    input:
        '{sample}_READ2.fastq.gz'
    params:
        prefix=SAMPLE_NAME,
        method='reads',
        extract_method='string',
        bc_pattern='CCCCCCCCCCCCCCCC',
        cell_number=CELL_NUMBER, # 5000 by default
        subset_reads='10000000000'
    output:
        '{sample}_whitelist.tsv',
        '{sample}_cell_barcode_knee.png',
        '{sample}_cell_barcode_counts.png'
    shell:
        'umi_tools whitelist --method={params.method} --extract-method={params.extract_method} --bc-pattern={params.bc_pattern} -I {input} -S {params.prefix}_whitelist.tsv --plot-prefix={params.prefix} --set-cell-number={params.cell_number} --subset-reads={params.subset_reads}'


#2a) compress files
rule compress:
    input:
        expand('{sample}_BC_{barcode}_READ{read}.fq', sample=SAMPLE_NAME, read=READS,barcode=BARCODES)
    output:
        expand('{sample}_BC_{barcode}_READ{read}.fq.gz',sample=SAMPLE_NAME, read=READS,barcode=BARCODES)
    shell:
        r'gzip {input} {output}'

#3a) allignement
rule bwa:
    input:
        GENOME, 
        '{sample}_BC_{barcode}_READ1.fq.gz',
        '{sample}_BC_{barcode}_READ3.fq.gz'
    params:
        ids=utilities.exp_bc('merged_sample_BC_{barcode}_READ1.fq.gz'),
        center='COSR',
        platform='Illumina',
        prefix=SAMPLE_NAME,
        lib='not_specified',
        threads_bwa=THREADS-2,
        threads_samtools=THREADS-6
    output:
        '{sample}_BC_{barcode}.bam'
    shell:
        'bwa mem -R "@RG\\tID:BC_{params.ids}\\tPL:{params.platform}\\tPU:{params.ids[1]}\\tLB:{params.lib}\\tSM:{params.prefix}\\tCN:{params.center}" -t {params.threads_bwa} {input} |  samtools sort -T {output}_tmp -@ {params.threads_samtools} -o {output}'
#4a) indexing allignement

rule index_allignement:
    input:
        '{sample}_BC_{barcode}.bam'
    params:
        threads=THREADS
    output:
        '{sample}_BC_{barcode}.bam.bai'
    shell:
        'samtools index -@ {params.threads} {input}'

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
        whitelist='{sample}_whitelist.tsv',
        read2='{sample}_BC_{barcode}_READ2.fq.gz',
        bamfile='{sample}_BC_{barcode}.bam',
    params:
        tn=tn_id('{sample}_BC_{barcode}_READ2.fq.gz'),
        prefix=SAMPLE_NAME
    output:
        '{sample}_BC_{barcode}_bcdedup.bam'
    shell:
        'python ~/scatACC/bc2rg.py -w {input.whitelist} -i {input.read2} -b {input.bamfile} -O -k -G {params.prefix}_{params.tn}_{wildcards.barcode} | python ~/scatACC/cbdedup.py -I -o {output}'
        

#6) merge
rule merge_bam_tn5:
    input:
        expand('{sample}_BC_{barcode}_bcdedup.bam', sample=SAMPLE_NAME, barcode=TN_BARCODES['tn5'])
    params:
        threads=THREADS
    output:
        expand('{sample}_tn5_merged.bam',sample=SAMPLE_NAME)
    shell:
        'samtools merge -c -p -@ {params.threads} {output} {input}'

rule merge_bam_tnh:
    input:
        expand('{sample}_BC_{barcode}_bcdedup.bam', sample=SAMPLE_NAME, barcode=TN_BARCODES['tnh'])
    params:
        threads=THREADS
    output:
        expand('{sample}_tnh_merged.bam',sample=SAMPLE_NAME)
    shell:
        'samtools merge -c -p -@ {params.threads} {output} {input}'