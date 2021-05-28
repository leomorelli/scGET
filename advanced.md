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
    

2. *name_of_the_file*: each line of the `input_list` indicates the name of the file only. The read is assigned:
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
### 3. cell_number
### 4. genome
### 5. scatacc_path
### 6. bed_file
### 7. binary
### 8. output_path
### 9. sample
