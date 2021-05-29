# scGET sequencing
f
# Configuration management
`scGET` allows the configuration of different parameters. The `config.yaml` file contains the default setting of each parameter. However, each default setting can be configured differently, directly from terminal. Here you can see a visual representation of configurable parameters needed for a personalyzed analysis.

![conf_params](conf_params.png)


### 1. input_list



The `input_list` must be a text file. Each line has to indicate the name of one of the fastq files. The `input_list` is accepted in two fromats:
1. *name_of_the_file* + *read_number*: read numbers are explicitly expressed for each file:
> sample_S1_L001_R1_001.fastq.gz 1
> 
> sample_S1_L001_R2_001.fastq.gz 2
> 
> sample_S1_L001_R3_001.fastq.gz 3
> 
> sample_S1_L002_R1_001.fastq.gz 1
> 
> sample_S1_L002_R2_001.fastq.gz 2
> 
> sample_S1_L002_R3_001.fastq.gz 3
    

2. *name_of_the_file*: each line of the `input_list` indicates only the name of the file. A simple script will be in charge of searching for the pattern `_R1`,`_R2`, or `_R3` within each file name, in order to assign the number of read:
> sample_S1_L001_R1_001.fastq.gz 
> 
> sample_S1_L001_R2_001.fastq.gz 
> 
> sample_S1_L001_R3_001.fastq.gz 
> 
> sample_S1_L002_R1_001.fastq.gz 
> 
> sample_S1_L002_R2_001.fastq.gz 
> 
> sample_S1_L002_R3_001.fastq.gz 


### 2. input_path

The path for the input files must be expressed.

### 3. cell_number
The `cell_number` parameter accounts for the number of cellular barcodes extracted by `umi_tools`. By default `cell_number`= 5000.
### 4. genome
The path for an indexed genome is necessary for the allignment step. It may be more convenient to modify the `genome` parameter directly from the `config.yaml` file.
### 5. scatacc_path
It is necessary to clarify the path for `scatACC` directory in order to perform the deduplication step. It may be more convenient to modify the `scatacc_path` parameter directly from the `config.yaml` file.
### 6. bed_file
This parameter is mandatory for the peak count step. It may be more convenient to modify the `bed_file` parameter directly from the `config.yaml` file.
### 7. binary

### 8. output_path
### 9. sample
