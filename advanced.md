# scGET sequencing
**scGET sequencing** approach enables the analysis of both accessible and compacted chromatin, exploiting the ability of two transposases: **tn5** recognizes open chromatin, while **tnh** recognizes closed chromatin. **scGET** protocol generates 3 different reads each genomic fragment:
- read 1 (R1) represents the **forward read**
- read 2 (R2) accounts for the **cellular barcode**
- read 3 (R3) represents the **reverse read**
<img src="img/scget_workflow.png" alt="img/scget_workflow" width="700"/>

Once the **sequencing** protocol is finished, R1, R2 and R3 reads will be stored separately: each fragment R1 will be stored in the same file, while each fragment R2 will converge in another file, as well as each fragment R3.


From these files starts the analysis with `scGET` library.



# Configuration management
`scGET` allows the configuration of different parameters. The `config.yaml` file contains the default setting of each parameter. However, each default setting can be configured differently, directly from terminal. Here you can see a visual representation of configurable parameters needed for a personalyzed analysis.

![img/scget_params](img/scget_params.png)


### 1. input_list

The `input_list` must be a text file. Each line has to indicate the name of one of the fastq files. Up to now, the `input_list` is accepted only, using the following format: *name_of_the_file* + *read_number* + *sample_name* must be expressed for each file.
> sample_S1_L001_R1_001.fastq.gz 1 S1
> 
> sample_S1_L001_R2_001.fastq.gz 2 S1
> 
> sample_S1_L001_R3_001.fastq.gz 3 S1
> 
> sample_S2_L001_R1_001.fastq.gz 1 S2
> 
> sample_S2_L001_R2_001.fastq.gz 2 S2
> 
> sample_S2_L001_R3_001.fastq.gz 3 S2
    
This is necessary in order to enable the simultaneous processing of different samples, using a single command.

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
`binary`= False by default:
- If `binary`= False, it enables the selection of output counts over intervals instead of binary data. 
- If `binary`= True, matrix counts are binary.
### 8. tn5 & tnh
`tn5`= True & `tnh`= True by default.
If `tn5` or `tnh` = False, the corresponding barcodes are excluded from the analysis.
### 9. atac
`atac`= False by default.
If `atac`= True, the worflow chenges:
- from the `modules/` folder `atac_wf.smk` is imported in stead of `scget_wf.smk`.
- `Tagdust` step is no more necessary, since scATAC-seq does not make use of scGET-seq barcodes
![img/atac_wf](img/atac_wf.png)
### 10. output_path
It indicates where `scGET` results will be stored: resulting files will be stored in a folder named after the sample name.


Example:
- `output_path`= ${HOME}/results
- `sample`= scget_files


Finally, each file can be found at `${HOME}/results/scget_files` path.
