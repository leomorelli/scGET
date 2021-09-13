import os

abs_path=os.path.abspath('')

# PARAMS NEEDED
OUTPUT_PATH=config['output_path']
INPUT_LIST=config['input_list']
READS=config['reads']
INPUT_PATH=config['input_path']
samples_tot=[x.strip() for x in open(f'{INPUT_LIST}').read().split('\n')[:-1]]
SAMPLES=list(set([x.split(' ')[-1] for x in samples_tot if len(x)>0]))



#CREATE DIFFERENT DIR FOR EACH SAMPLE AND GENERATE TEXT FILES DESCRIBING EACH FILE OF EACH SAMPLE
rule info_text:
    input:
        text=INPUT_LIST
    output:
        expand('{out}/{samples}/info_{samples}.txt', out=OUTPUT_PATH, samples=SAMPLES),
        expand('{out}/{samples}/info_READ{read}.txt', out=OUTPUT_PATH, samples=SAMPLES,read=READS)
    run:
        paral.prep_input(INPUT_LIST,OUTPUT_PATH,INPUT_PATH)



#MERGE OF DIFFERENT INPUT FILES OF THE SAME SAMPLE (Rename each file with a tag expressing the number of read)

rule merge_reads:
    input:
        f'{OUTPUT_PATH}/{{sample}}/info_READ{{read}}.txt'    
    output:
        f'{OUTPUT_PATH}/{{sample}}/{{sample}}_READ{{read}}.fastq.gz'
    wildcard_constraints:
        sample='.*[a-zA-Z0-9_]'
    shell:
        "cat $(cat {input}) > {output}"

