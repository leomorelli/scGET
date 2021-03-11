# scGET
**s**ingle **c**ell **G**enome and **E**pigenome by **T**ransposases sequencing (**scGET**-seq) aims to discriminate between open accessible chormatin regions and closed compacted chromatin regions within each single cell, taking advantage of two different transposases: transposase-5 binds to *accessible chromatin* (**tn5**) and transposase-H, a chimeric form of tn5 (**tnH**), which binds *compacted chromatin*.


The workflow for analysis of scGET-generated data is based on the *Snakemake*: a workflow management system, which guarantees the possibility to **parallelize** independent jobs. **scGET** workflow is described by the image below: it starts from the *FASTQ* files generated by scGET-seq and it terminates with the production of 2 *BAM* files, accounting for open chromatin regions (**tn5_bam**) and closed chromatin regions (**tnH_bam**).


 ![dag](dag.svg)
# Installation


# Set up


# How to use
