# illuminaCallHPVInt
Workflow for calling HPV integration sites and events in Illumina short-read sequencing data. 

# Installation
This will clone the repository. You can run the IMPALA within this directory.
```
git clone https://github.com/bcgsc/IMPALA.git
```

### Dependencies
> To run this workflow, you must have snakemake (v6.12.3) and singularity (v3.5.2-1.1.el7). You can install snakemake using [this guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and singularity using [this guide](https://docs.sylabs.io/guides/3.5/admin-guide/installation.html). The remaining dependencies will be downloaded automatically within the snakemake workflow.

# Input Files

### **Method 1**: Whole genome short reads <br />
- WGS alignment (bam file)

### **Method 2**: RNA-seq <br />
- RNA alignment (bam file)

# Running Workflow

### **Edit the config files**

#### **Example parameters.yaml:** <br />
Config files to specify parameters and paths needed for the workflow. The main parameter to include is the genome path.

```
genome_path: /path/to/genome/fasta
```

#### **Example samples.yaml:** <br />
Main config file to specify input files.

```
bams:
    sampleName_1: /path/to/bam/file
    sampleName_2: /path/to/bam/file

```

# converting sample paths to yaml file
A text file can be converted to the samples.yaml file using the scripts/sampletsvtoyaml.py script. The tsv file should have the sample name in one column and the path in another and be tab delimited (no header). 

```
scripts/sampletsvtoyaml.py -t samples.txt -o config/samples.yaml

```

### **Run snakemake**
This is the command to run it with conda. The `-c` parameter can be used to specify maximum number of threads. 

```
snakemake -c 30 --use-conda
```
