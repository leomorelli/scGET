# scGET
**s**ingle **c**ell **G**enome and **E**pigenome by **T**ransposases sequencing (**scGET**-seq) aims to discriminate between open accessible chormatin regions and closed compacted chromatin regions within each single cell, taking advantage of two different transposases: transposase-5 binds to *accessible chromatin* (**tn5**) and transposase-H, a chimeric form of tn5 (**tnH**), which binds *compacted chromatin*.


The workflow for analysis of scGET-generated data is based on the *Snakemake*: a workflow management system, which guarantees the possibility to **parallelize** independent jobs. **scGET** workflow is described by the image below: it starts from the *FASTQ* files generated by scGET-seq and it terminates with the production of 2 *BAM* files, accounting for open chromatin regions (**tn5_bam**) and closed chromatin regions (**tnH_bam**).


 ![dag](dag.svg)
# Installing scGET
Before getting your hands dirty with scGET analyses, it is necessary to create a suitable conda environment. However, some of the packages requires different installation procedures. Therefore, we have designed a 3-step process, allowing an easy and quick generation of the *scget environment*.


1. The conda *environment* can be automatically generated, along with the installation of the *majority of required packages* thanks to the **scget.yaml** file:
```
conda env create -f scget.yaml
conda activate scget
```
2. After having activated the *scget* environment, it is necessary to install the **TagDust** package. First, the package must be downloaded and compiled; second, from the *tagdust* directory we can copy the binary tagdust file in our *scget* environment: 
```
tar -zxvf tagdust-2.33.tar.gz 
cd tagdust
./configure 
make
make check
cp ./src/tagdust $CONDA_PREFIX/bin
```
3. The last step accounts for the cloning of **scatACC** repository from github:
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

 
# Set up
Before starting using scGET


### Cluster set up
scGET workflow is intended to be used with the support of a calculus cluster, therefore

### scatACC path

# How to use
standard input 

standard output

for more read advanced.md
