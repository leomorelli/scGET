import os

abs_path=os.path.abspath('')

# [A] SCRIPTS IMPORT
# spec_u = importlib.util.spec_from_file_location("utility_functions",f"{abs_path}/scripts/utility_functions.py")
# utilities = importlib.util.module_from_spec(spec_u)
# spec_u.loader.exec_module(utilities)
#
# spec_i = importlib.util.spec_from_file_location("get_info",f"{abs_path}/scripts/get_info.py")
# infos = importlib.util.module_from_spec(spec_i)
# spec_i.loader.exec_module(infos)
#
# spec_l = importlib.util.spec_from_file_location("tn_layers",f"{abs_path}/scripts/tn_layers.py")
# layers = importlib.util.module_from_spec(spec_l)
# spec_l.loader.exec_module(layers)
#
# # [B] CONFIG HANDLING
# configfile: f'{abs_path}/config.yaml'
#
# include: f'{abs_path}/modules/init_lane'
#
# # [b.1] DEFINITION OF TRANSPOSASES AND BARCODES FOR SCGET ANALYSIS
# if config['tn5']!=True:
#     TN_BARCODE=config['barcodes']
#         TN_BARCODES={'tnh':TN_BARCODE['tnh']}
#             BARCODES=TN_BARCODES['tnh']
#             elif config['tnh']!=True:
#                 TN_BARCODE=config['barcodes']
#                     TN_BARCODES={'tn5':TN_BARCODE['tn5']}
#                         BARCODES=TN_BARCODES['tn5']
#                         else:
#                             TN_BARCODES=config['barcodes']
#                                 BARCODES=TN_BARCODES['tn5']+TN_BARCODES['tnh']
#
#                                 # [b.2] OTHER VARIABLES
#                                 CELL_NUMBER=config['cell_number']
SAMPLE= config['sample']
READS= config['reads']
THREADS=config['threads']
INPUT_PATH=config['input_path']
INPUT_LIST=config['input_list']
SAMPLE_NAME=utilities.sample_name(SAMPLE,INPUT_PATH)
OUTPUT_PATH=utilities.output(config['output_path'],utilities.sample_name(SAMPLE,INPUT_PATH))



#PREMERGE OF DIFFERENT INPUT FILES OF THE SAME SAMPLE


#Create a text file where each file is assigned to its read

rule file_read:
    input:
        files=utilities.list_of_files(INPUT_PATH,INPUT_LIST),
    output:
        temp(expand('{output}/info_READ{read}.txt',output=OUTPUT_PATH,read=READS)),
    run:
        info_read=infos.df_info(input)
        infos.create_read_files(OUTPUT_PATH,info_read)



#Rename each file with a tag expressing the number of read

rule merge_reads:
    input:
        '{output_path}/info_READ{read}.txt'
    output:
        '{output_path}/{sample}_READ{read}.fastq.gz'
    shell:
        "cat $(cat {input}) > {output}"

