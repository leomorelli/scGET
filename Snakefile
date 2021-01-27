import gzip
configfile: 'config.yaml'

SAMPLE= config['sample']
READS= config['reads']
TN_BARCODES=config['barcodes']
BARCODES=TN_BARCODES['tn5']+TN_BARCODES['tnh']
GENOME=config['genome']
THREADS=config['threads']
CELL_NUMBER=config['cell_number']


rule all:
    input:       
        expand('{sample}_un_READ{read}.fq', sample=SAMPLE,read=READS),        
        expand('{sample}_logfile.txt', sample=SAMPLE), 
        expand('{sample}_whitelist.tsv',sample=SAMPLE),
        expand('{sample}_cell_barcode_knee.png',sample=SAMPLE),
        expand('{sample}_cell_barcode_counts.png',sample=SAMPLE),
        expand('{sample}_tn5_merged.bam',sample=SAMPLE),
        expand('{sample}_tnh_merged.bam',sample=SAMPLE)

#1a) Classify each read using its barcode
rule tag_dust:
    input:
        expand('{sample}_R{read}.fastq.gz', sample=SAMPLE,read=READS)
    params:
        prefix=SAMPLE,
        list_BCs=','.join(BARCODES),
        threads=THREADS
    output:
        expand('{sample}_logfile.txt', sample=SAMPLE),        
        expand('{sample}_un_READ{read}.fq', sample=SAMPLE,read=READS),        
        expand('{sample}_BC_{barcodes}_READ{read}.fq', sample=SAMPLE, barcodes=BARCODES,read=READS) 
    shell:
        'tagdust -1 B:{params.list_BCs} -2 S:AGATATATATAAGGAGACAG -3 R:N {input} -o {params.prefix} -t {params.threads}'

#1b) umi_tools
rule umi_tools:
    input:
        expand('{sample}_R2.fastq.gz', sample=SAMPLE)
    params:
        prefix=SAMPLE,
        method='reads',
        extract_method='string',
        bc_pattern='CCCCCCCCCCCCCCCC',
        cell_number=CELL_NUMBER, # it must be 5000 i set it at 100
        subset_reads='10000000000'
    output:
        expand('{sample}_whitelist.tsv',sample=SAMPLE) ,
        expand('{sample}_cell_barcode_knee.png',sample=SAMPLE),
        expand('{sample}_cell_barcode_counts.png',sample=SAMPLE) 
    shell:
        'umi_tools whitelist --method={params.method} --extract-method={params.extract_method} --bc-pattern={params.bc_pattern} -I {input} -S {params.prefix}_whitelist.tsv --plot-prefix={params.prefix} --set-cell-number={params.cell_number} --subset-reads={params.subset_reads}'

#2a) compress files
rule compress:
    input:
        expand('{sample}_BC_{barcode}_READ{read}.fq', sample=SAMPLE,read=READS,barcode=BARCODES)
    output:
        expand('{sample}_BC_{barcode}_READ{read}.fq.gz', sample=SAMPLE,read=READS,barcode=BARCODES)
    shell:
        r'gzip {input} {output}'

#3a) allignement
def exp_id(file):
#   f=open(file).read()
    f=gzip.open(file,'r').read().decode() #quando userò file gz
    f=f.split('\n')
    f1=f[0].split(':')
    ids=[f1[4],f'{f1[2]}_{f1[4]}']
    return ids

rule bwa:
    input:
        GENOME, 
        lambda wildcards: expand('{sample}_BC_{barcode}_READ1.fq.gz', sample=SAMPLE,barcode=wildcards.barcode),
        lambda wildcards: expand('{sample}_BC_{barcode}_READ3.fq.gz', sample=SAMPLE,barcode=wildcards.barcode)
    params:
        ids=exp_id(expand('{sample}_R1.fastq.gz',sample=SAMPLE)[0]),
        center='COSR',
        platform='Illumina',
        prefix=SAMPLE,
        lib='not_specified',
        threads_bwa=THREADS-2,
        threads_samtools=THREADS-6
    output:
        '{sample}_BC_{barcode}.bam'
    shell:
        'bwa mem -R "@RG\\tID:{params.ids[0]}_BC\\tPL:{params.platform}\\tPU:{params.ids[1]}\\tLB:{params.lib}\\tSM:{params.prefix}\\tCN:{params.center}" -t {params.threads_bwa} {input} |  samtools sort -T {output}_tmp -@ {params.threads_samtools} -o {output}'
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
    bc_in=str(file).split('_')[-2]
    if bc_in in BCs['tn5']:
        return 'tn5'        
    elif bc_in in BCs['tnh']:
        return 'tnh' 
    else:
        return 'nan'

rule dedup:
    input:
        whitelist='{sample}_whitelist.tsv',
        read2=lambda wildcards: expand('{sample}_BC_{barcode}_READ2.fq.gz', sample=SAMPLE,barcode=wildcards.barcode),
        bamfile=lambda wildcards: expand('{sample}_BC_{barcode}.bam', sample=SAMPLE,barcode=wildcards.barcode)
    params:
        tn=tn_id('{sample}_BC_{barcode}_READ2.fq.gz'),
        prefix=SAMPLE
    output:
        '{sample}_BC_{barcode}_bcdedup.bam'
    shell:
        'python ~/scatACC/bc2rg.py -w {input.whitelist} -i {input.read2} -b {input.bamfile} -O -k -G {params.prefix}_{params.tn}_{wildcards.barcode} | python ~/scatACC/cbdedup.py -I -o {output}'
        

#6) merge
rule merge_bam_tn5:
    input:
        expand('{sample}_BC_{barcode}_bcdedup.bam', sample=SAMPLE,barcode=TN_BARCODES['tn5'])
    params:
        threads=THREADS
    output:
        expand('{sample}_tn5_merged.bam',sample=SAMPLE)
    shell:
        'samtools merge -c -p -@ {params.threads} {output} {input}'

rule merge_bam_tnh:
    input:
        expand('{sample}_BC_{barcode}_bcdedup.bam', sample=SAMPLE,barcode=TN_BARCODES['tnh'])
    params:
        threads=THREADS
    output:
        expand('{sample}_tnh_merged.bam',sample=SAMPLE)
    shell:
        'samtools merge -c -p -@ {params.threads} {output} {input}'