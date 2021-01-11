configfile: "config.yaml"

SAMPLE= config['sample']
READS= config['reads']
BARCODES=config['barcodes']

#bc_tn5=['CGTACTAG','TCCTGAGC','TCATGAGC','CCTGAGAT']
#bc_tnh=['TAAGGCGA','GCTACGCT','AGGCTCCG','CTGCGCAT']


rule all:
    input:
        expand("{sample}_BC_{barcodes}_READ{read}.fq", sample=SAMPLE, barcodes=BARCODES,read=READS),        
        expand("{sample}_un_READ{read}.fq", sample=SAMPLE,read=READS),        
        expand("{sample}_logfile.txt", sample=SAMPLE), 
        expand('{sample}_whitelist.tsv',sample=SAMPLE) ,
        expand('{sample}_cell_barcode_knee.png',sample=SAMPLE),
        expand('{sample}_cell_barcode_counts.png',sample=SAMPLE),
        expand('{sample}_BC_{barcode}_READ{read}.fq.gz', sample=SAMPLE,read=READS,barcode=BARCODES),  
        expand('{sample}_BC_{barcode}.bam', sample=SAMPLE, barcode=BARCODES),
        expand('{sample}_BC_{barcode}.bam.bai', sample=SAMPLE,barcode=BARCODES),
        expand('{sample}_BC_{barcode}_bcdedup.bam', sample=SAMPLE,barcode=BARCODES),
        expand('{sample}_merged.bam',sample=SAMPLE)

#1a) Classify each read using its barcode
rule tag_dust:
    input:
        expand("{sample}_R{read}.fastq", sample=SAMPLE,read=READS)
    params:
        prefix=SAMPLE
    output:
        expand("{sample}_logfile.txt", sample=SAMPLE),        
        expand("{sample}_un_READ{read}.fq", sample=SAMPLE,read=READS),        
        expand("{sample}_BC_{barcodes}_READ{read}.fq", sample=SAMPLE, barcodes=BARCODES,read=READS) 
    shell:
        'tagdust -1 B:TAAGGCGA,GCTACGCT,AGGCTCCG,CTGCGCAT,CGTACTAG,TCCTGAGC,TCATGAGC,CCTGAGAT -2 S:AGATATATATAAGGAGACAG -3 R:N {input} -o {params.prefix}'

#1b) umi_tools
rule umi_tools:
    input:
        expand('{sample}_R2.fastq', sample=SAMPLE)
    params:
        prefix=SAMPLE,
        method='reads',
        extract_method='string',
        bc_pattern='CCCCCCCCCCCCCCCC',
        cell_number='100', #5000
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
	f=open(file).read()
#	f=gzip.open(file,'rb').read() quando userò file gz
	f=f.split('\n')
	f1=f[0].split(':')
	ids=[f1[4],f'{f1[2]}_{f1[4]}']
	return ids

rule bwa:
    input:
        '../genome/hg38.fa', 
        lambda wildcards: expand('{sample}_BC_{barcode}_READ1.fq.gz', sample=SAMPLE,barcode=wildcards.barcode),
        lambda wildcards: expand('{sample}_BC_{barcode}_READ3.fq.gz', sample=SAMPLE,barcode=wildcards.barcode)
    params:
        ids=exp_id(expand('{sample}_R1.fastq',sample=SAMPLE)[0]),
        center='COSR',
        platform='Illumina',
        prefix=SAMPLE,
        lib='not_specified'
    output:
        '{sample}_BC_{barcode}.bam'
    shell:
        'bwa mem -R "@RG\\tID:{params.ids[0]}_BC\\tPL:{params.platform}\\tPU:{params.ids[1]}\\tLB:{params.lib}\\tSM:{params.prefix}\\tCN:{params.center}" -t 4 {input} |  samtools sort -T {output}_tmp -o {output}'

#4a) indexing allignement

rule index_allignement:
    input:
        '{sample}_BC_{barcode}.bam'
    output:
        '{sample}_BC_{barcode}.bam.bai'
    shell:
        'samtools index {input}'

# 5) deduplication
def tn_id(file):
    BCs={'tn5':['CGTACTAG','TCCTGAGC','TCATGAGC','CCTGAGAT'],'tnh':['TAAGGCGA','GCTACGCT','AGGCTCCG','CTGCGCAT']}
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
rule merge_bam:
    input:
        expand('{sample}_BC_{barcode}_bcdedup.bam', sample=SAMPLE,barcode=BARCODES)
    output:
        expand('{sample}_merged.bam',sample=SAMPLE)
    shell:
        'samtools merge -c -p {output} {input}'
