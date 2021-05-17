# 1. scGET
**s**ingle **c**ell **G**enome and **E**pigenome by **T**ransposases sequencing (**scGET**-seq) aims to discriminate between open accessible chormatin regions and closed compacted chromatin regions within each single cell, taking advantage of two different transposases: transposase-5 binds to *accessible chromatin* (**tn5**) and transposase-H, a chimeric form of tn5 (**tnH**), which binds *compacted chromatin*.


The workflow for analysis of scGET-generated data is based on the *Snakemake*: a workflow management system, which guarantees the possibility to **parallelize** independent jobs. **scGET** workflow is described by the image below: it starts from the *FASTQ* files generated by scGET-seq and it terminates with the production of 2 *BAM* files, accounting for open chromatin regions (**tn5_bam**) and closed chromatin regions (**tnH_bam**).


 ![dag](dag.svg)
 <img src="dag.svg" alt="dag" width="4000"/>
# 2. Installing scGET
Before getting your hands dirty with scGET analyses, it is necessary to create a suitable conda environment. However, some of the packages requires different installation procedures. Therefore, we have designed a 4-step process, allowing an easy and quick generation of the *scget environment*.


1. The conda *environment* can be automatically generated, along with the installation of the *majority of required packages* thanks to the **scget.yaml** file:
```
conda env create -f scget.yaml
conda activate scget
```
2. After having activated the *scget* environment, it is necessary to install the **TagDust** package. First, the package must be downloaded and compiled; second, from the *tagdust* directory we can copy the binary tagdust file in our *scget* environment: 
```
wget https://sourceforge.net/projects/tagdust/files/tagdust-2.33.tar.gz
tar -zxvf tagdust-2.33.tar.gz 
cd tagdust-2.33
./configure 
make
make check
cp ./src/tagdust $CONDA_PREFIX/bin
```
3. Similarly to the previous step, also **samtools** must be installed:
   - Firstly, git repositories of **samtools** and **htslib** must be cloned
   - Secondly, **htslib** must be compiled and installed
   - Finally, **samtools** must be compiled and installed
```
git clone https://github.com/samtools/samtools
git clone https://github.com/samtools/htslib

cd htslib
autoreconf -i
git submodule update --init --recursive
./configure --prefix=$CONDA_PREFIX
make 
make install

cd samtools
autoheader
autoconf -Wno-syntax
./configure --prefix=$CONDA_PREFIX --without-curses
make
make install
```
4. The last step accounts for the cloning of **scatACC** repository from github:
```
git clone https://github.com/dawe/scatACC
```
In order to perform the analysis through the calculus cluster, it may be useful to check if `screen` package is already installed:
```
screen --version
```
Output:
>Screen version 4.08.00 (GNU) 05-Feb-20


If `screen` has **not** been **installed** yet, it could be easily installed via `sudo`:
```
sudo apt update
sudo apt install screen
```

 
# 3. Set up
Finally, there are three small tricks left to perform, before you can finally use scGET.


### a. Cluster set up
 - In your home directory check if you have a snakemake folder inside ``${HOME}/.config``. Inside this folder create a slurm folder add in a *config.yaml* file:
```
mkdir -p ${HOME}/.config/snakemake/slurm
vi ${HOME}/.config/snakemake/slurm/config.yaml
```
 - compile the *config.yaml* file with the content specified below (remember to update the queue name specified by the ``-p`` option and  your ``mail-user``):
```
jobs: 34
cluster: "sbatch --mem={resources.mem_mb} -c {resources.cpus} -o logs_slurm/{SAMPLE_NAME}/{rule}_{wildcards} -e logs_slurm/{SAMPLE_NAME}/{rule}_{wildcards} --mail-type=FAIL --mail-user=user.mail@hsr.it"
default-resources: [cpus=1, mem_mb=2000]
resources: [cpus=10, mem_mb=50000]
use-conda: true
```

### b. scatACC path & genome path
The scGET analysis must **start** from the **scGET directory**. The *scatACC library* is used for one the last step of scGET analysis as well as the *genome file* and the *bed file*: therefore, the location of **scatACC directory** together with the path for the **genome** and the **bed_file** must be clarified in order to perform a succesful analysis.


*WHERE SHOULD I CLARIFY SUCH PATHS?*


Let's exemplify! 


Let's assume that the **scGET** directory is located in our **home** directory (*${HOME}/scGET*); while our **scatACC** directory is situated in a directory called *repositories* (*${HOME}/repositories/scatACC*); on the other hand, the **genome**  file, which we have named *hg38.fa*, lays in the *references* directory (*${HOME}/references/hg.38*) as well as the **bed_file** (*${HOME}/references/hg385kbin.bed*):
- First, we should open the *config.yaml* file present in **scGET** directory:
```
cd ${HOME}/scGET
vi config.yaml
```
Output:
>sample: ''
>
>reads: [1,2,3]
>
>barcodes: {'tn5':['CGTACTAG','TCCTGAGC','TCATGAGC','CCTGAGAT'],'tnh':['TAAGGCGA','GCTACGCT','AGGCTCCG','CTGCGCAT']}
>
>genome: ${HOME}/genome.fa
>
>bed_file: ${HOME}/genome.bed
>
>threads: 8
>
>cell_number: 5000
>
>scatacc_path: '${HOME}/scatACC'
>
>input_path: ''
>
>input_list: ''
>
>output_path: ''


- After that, we must modify the field *scatacc_path*, specifying our **actual scatACC path**, the field *genome*, clarifying the **genome path** with the **genome file name** and the field *bed_file* with the path for the **bed file**:


Output:
>sample: ''
>
>reads: [1,2,3]
>
>barcodes: {'tn5':['CGTACTAG','TCCTGAGC','TCATGAGC','CCTGAGAT'],'tnh':['TAAGGCGA','GCTACGCT','AGGCTCCG','CTGCGCAT']}
>
>genome: ${HOME}/**references**/**hg38.fa**
>
>bed_file: ${HOME}/**references**/**hg385kbin.bed**
>
>threads: 8
>
>cell_number: 5000
>
>scatacc_path: '${HOME}/**repositories**/scatACC'
>
>input_path: ''
>
>input_list: ''
>
>output_path: ''


**N.B. the REFERENCE GENOME must be INDEXED before the analysis**


If the genome has not been indexed yet, you can make up for this in three steps:
- Activate the **scget** conda **environment**
- Open the **directory** where the reference genome is stored
- **Index** the genome, using **samtools** library
```
conda activate scget
cd ${HOME}/references
samtools index hg38.fa
```
# 4. How to use
Now it's time to start the analysis! It is important to remember that the scGET analysis must be performed from the **scGET directory**. Therefore, before starting the workflow, you should reach the **scGET directory** and activate the **scGET environment**.
```
cd ${HOME}/scGET
conda activate scget
```
### a. Standard input 
Three different inputs are necessary to start the scGET analysis:
- Of course, **FASTQ files** generated by scGET sequencing step are mandatory:


Example:
>sample_S1_L001_R1_001.fastq.gz
>
>sample_S1_L001_R2_001.fastq.gz
>
>sample_S1_L001_R3_001.fastq.gz
>
>sample_S1_L002_R1_001.fastq.gz
>
>sample_S1_L002_R2_001.fastq.gz
>
>sample_S1_L002_R3_001.fastq.gz


- Since the analysis starts from **scGET** directory, it is necesary to indicate the **path** to **FASTQ files**


Example: *FASTQ files* are stored in ``/home/files/experiment_test`` directory.


- Finally, it must be generated a ``.txt`` file, containing names of input *FASTQ files* and the corresponding number of *read*: first, the **input_file.txt** must be created.
```
vi input_file.txt
```
After that, it must be modified as explained below:
>sample_S1_L001_R1_001.fastq.gz 1
>
>sample_S1_L001_R2_001.fastq.gz 2 
>
>sample_S1_L001_R3_001.fastq.gz 3
> 
>sample_S1_L002_R1_001.fastq.gz 1
>
>sample_S1_L002_R2_001.fastq.gz 2
>
>sample_S1_L002_R3_001.fastq.gz 3


### b. Standard command
In order to start with scGET analysis, you must run the following command, specifying the *input_path* and the *input_list*, generated above:
```
snakemake --cores 8 --config input_path=/home/files/experiment_test input_list=input_file.txt --profile slurm
```
### c. Standard output
Two type of outputs are generated by the **scGET** analysis:
- **Output files**: files generated by scGET analysis. They are stored in a directory named as the directory of input files, located in the scGET directory. Here you can find also the final output of the workflow: the ``tn5.bam`` and the ``tnH.bam`` files.


Example: *output files* are stored in ``${HOME}/scGET/experiment_test``


- **Log files**: log files generated during the analysis, describing each step of the workflow:


Example: *log files* are stored in ``${HOME}/scGET/logs_slurm/experiment_test``






**N.B.**


If you need to dig more into **scGET settings**, you can find more info about scGET usage in the **advanced.md** file.
