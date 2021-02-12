import gzip
import glob

configfile: 'config.yaml'

SAMPLE= config['sample']
READS= config['reads']
TN_BARCODES=config['barcodes']
BARCODES=TN_BARCODES['tn5']+TN_BARCODES['tnh']
GENOME=config['genome']
THREADS=config['threads']
CELL_NUMBER=config['cell_number']
SAMPLE_PATH=config['path']


rule all:
    input:              
        expand('{sample}_{tn}_merged.bam',sample=SAMPLE, tn=TN_BARCODES.keys())


# 0) PREMERGE OF DIFFERENT INPUT FILES
# 0a) Create a text file where each file is assigned to its read
def list_of_files(path):
    if len(path)<1:
        return glob.glob('*fastq.gz')
    else:
        if path[-1]=='/':
            return glob.glob(path+'*fastq.gz')
        else:
            return glob.glob(path+'/'+'*fastq.gz')


rule file_read:
    input:
        files=list_of_files(SAMPLE_PATH)
    output:
        expand('merged_READ{read}.txt', read=READS)
    script:
        'scripts/get_info.py'


#0b) Rename each file with a tag expressing the number of read
rule merge_reads:
    input:
        'merged_READ{read}.txt'
    output:
        'merged_sample_READ{read}.fastq.gz'
    shell:
        "for F in $(cat {input}) ; do cat $F>>{output}; done"


#1a) Classify each read using its barcode
rule tag_dust:
    input:
        expand('merged_sample_READ{read}.fastq.gz', read=READS)
    params:
        prefix='merged_sample',
        list_BCs=','.join(BARCODES),
        threads=THREADS
    output:
        'merged_sample_logfile.txt',        
        expand('merged_sample_un_READ{read}.fq', read=READS),        
        expand('merged_sample_BC_{barcodes}_READ{read}.fq', barcodes=BARCODES,read=READS)
    shell:
        'tagdust -1 B:{params.list_BCs} -2 S:AGATATATATAAGGAGACAG -3 R:N {input} -o {params.prefix} -t {params.threads}'


#1b) umi_tools
rule umi_tools:
    input:
        'merged_sample_READ2.fastq.gz'
    params:
        prefix='merged_sample',
        method='reads',
        extract_method='string',
        bc_pattern='CCCCCCCCCCCCCCCC',
        cell_number=CELL_NUMBER, # 5000 by default
        subset_reads='10000000000'
    output:
        'merged_sample_whitelist.tsv',
        'merged_sample_cell_barcode_knee.png',
        'merged_sample_cell_barcode_counts.png'
    shell:
        'umi_tools whitelist --method={params.method} --extract-method={params.extract_method} --bc-pattern={params.bc_pattern} -I {input} -S {params.prefix}_whitelist.tsv --plot-prefix={params.prefix} --set-cell-number={params.cell_number} --subset-reads={params.subset_reads}'


#2a) compress files
rule compress:
    input:
        expand('merged_sample_BC_{barcode}_READ{read}.fq', read=READS,barcode=BARCODES)
    output:
        expand('merged_sample_BC_{barcode}_READ{read}.fq.gz',read=READS,barcode=BARCODES)
    shell:
        r'gzip {input} {output}'

#3a) allignement
def exp_bc(file):
    bc_in=str(file).split('_')[-2]
    return bc_in

rule bwa:
    input:
        GENOME, 
        'merged_sample_BC_{barcode}_READ1.fq.gz',
        'merged_sample_BC_{barcode}_READ3.fq.gz'
    params:
        ids=exp_bc('merged_sample_BC_{barcode}_READ1.fq.gz'),
        center='COSR',
        platform='Illumina',
        prefix='merged_sample',
        lib='not_specified',
        threads_bwa=THREADS-2,
        threads_samtools=THREADS-6
    output:
        'merged_sample_BC_{barcode}.bam'
    shell:
        'bwa mem -R "@RG\\tID:BC_{params.ids}\\tPL:{params.platform}\\tPU:{params.ids[1]}\\tLB:{params.lib}\\tSM:{params.prefix}\\tCN:{params.center}" -t {params.threads_bwa} {input} |  samtools sort -T {output}_tmp -@ {params.threads_samtools} -o {output}'
#4a) indexing allignement

rule index_allignement:
    input:
        'merged_sample_BC_{barcode}.bam'
    params:
        threads=THREADS
    output:
        'merged_sample_BC_{barcode}.bam.bai'
    shell:
        'samtools index -@ {params.threads} {input}'

# 5) deduplication
def tn_id(file):
    BCs=TN_BARCODES
    bc_in=str(file).split('_')[-2]
    if bc_in in BCs['tn5']:
        return 'tn5'        
    elif bc_in in BCs['tnh']:
        return 'tnh' 
    else:
        return 'nan'


rule dedup:
    input:
        whitelist='merged_sample_whitelist.tsv',
        read2='merged_sample_BC_{barcode}_READ2.fq.gz',
        bamfile='merged_sample_BC_{barcode}.bam',
    params:
        tn=tn_id('merged_sample_BC_{barcode}_READ2.fq.gz'),
        prefix='merged_sample'
    output:
        'merged_sample_BC_{barcode}_bcdedup.bam'
    shell:
        'python ~/scatACC/bc2rg.py -w {input.whitelist} -i {input.read2} -b {input.bamfile} -O -k -G {params.prefix}_{params.tn}_{wildcards.barcode} | python ~/scatACC/cbdedup.py -I -o {output}'
        

#6) merge
rule merge_bam_tn5:
    input:
        expand('merged_sample_BC_{barcode}_bcdedup.bam', barcode=TN_BARCODES['tn5'])
    params:
        threads=THREADS
    output:
        expand('{sample}_tn5_merged.bam',sample=SAMPLE)
    shell:
        'samtools merge -c -p -@ {params.threads} {output} {input}'

rule merge_bam_tnh:
    input:
        expand('merged_sample_BC_{barcode}_bcdedup.bam', barcode=TN_BARCODES['tnh'])
    params:
        threads=THREADS
    output:
        expand('{sample}_tnh_merged.bam',sample=SAMPLE)
    shell:
        'samtools merge -c -p -@ {params.threads} {output} {input}'