# Nextflow workflows

This folder contains the workflows for analyzing metagenomic data 
using Nextflow as workflow manager. The workflows are based on those in the 
[bhattlab_workflows Github](https://github.com/bhattlab/bhattlab_workflows)
repository.


## Preparations

### Adjusting the `params.yml` file

The `params.yml` file holds all the important settings for running the 
workflows on your system. You will need to adjust the file paths to the 
databases (see below for downloading those to your file system). 
You can also adjust the run parameters for some tools, such as the parameters
to run mOTUs, for example. The `params.yml` file should hold some explanation
for each of the parameters.

### Indexing the host reference genome

You will need a host reference genome on your file system so that you can 
remove all reads matching to the host genome from your sample, since we are
interested in the bacteria, not the host genome. We map reads to the host 
genome with the `bwa` algorithm, so you will need to index the host 
genome before being able to use the workflow.

You can download for example the human reference genome from the
[UCSC Genome Browser website](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/):
```bash
cd <your-reference-genome-location>
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

Now, you can index the genome with `bwa`. However, if you want to use the 
singularity container that we use for the workflow (because you might not have
`bwa` installed on your system), you can do so like that by using the `--bind`
command when running the singularity container:
```bash
cd <your-reference-genome-location>
singularity shell --bind ./:/mnt docker://ghcr.io/jakob-wirbel/micromamba-focal-preprocessing:latest
cd /mnt
bwa index hg38.fa
```

### Downloading the MetaPhlAn4 database

We also need to download the MetaPhlAn4 database, again using the already
prepared singularity container:

```bash
cd <your-metaphlan-database-location>
singularity shell --bind ./:/mnt docker://ghcr.io/jakob-wirbel/micromamba-focal-classification:latest
cd /mnt
metaphlan --install --bowtie2db ./
```

### Downloading the mOTUs3 database

The mOTU tool needs a custom database of marker genes. Unfortunately, you 
cannot download the database through the tool as of now, but maybe it will be
available as a feature in the future (see 
[this issue](https://github.com/motu-tool/mOTUs/issues/109)). Instead, we can
download and configure the database manually:

```bash
cd <your-motus-database-location>
wget https://zenodo.org/record/5140350/files/db_mOTU_v3.0.1.tar.gz
md5sum db_mOTU_v3.0.1.tar.gz
# expected output: f4fd09fad9b311fb4f21383f6101bfc3
# if you do not see this output, something went wrong and you need to download
# the database again!
tar -zxvf db_mOTU_v3.0.1.tar.gz
```

The version file in the Zenodo repository is not really up-to-date with the
version numbering system of mOTUs, so you will have to adjust it manually:
```bash
cd db_mOTU
sed -i 's/2.6.0/3.0.3/g' db_mOTU_versions
```

### Database for phanta-lite

For running a lite version of [Phanta](https://github.com/bhattlab/phanta), please 
note that Phanta was originally developed as a Snakemake workflow and the port 
to Nextflow is a bit hacky. As of yet, it is pretty unclear if we will want to 
include it in the final analysis. If you want to run it, you will need a 
Kraken2 database. 
[Here](https://github.com/bhattlab/phanta/blob/main/databases.md) is a list of
databases that are made available with the tool. You can download and extract
the tarball following the documentation for the tool.

### CheckM database

You will also need to download the database for checkM. This is available through 
Zenodo as well and can be downloaded and extracted like that:
```bash

cd <your-checkm-database-location>
wget https://zenodo.org/record/7401545/files/checkm_data_2015_01_16.tar.gz
md5sum checkm_data_2015_01_16.tar.gz
# expected output: 631012fa598c43fdeb88c619ad282c4d
# if you do not see this output, something went wrong and you need to download
# the database again!
tar -zxvf checkm_data_2015_01_16.tar.gz
rm checkm_data_2015_01_16.tar.gz
```

## Running the workflows

To run the workflow, you need Java and Nextflow running in your system. On an 
HPC system, those are usually available through the `module load` system. 

```bash
module load java/18.0.2.1
module load nextflow/22.10.5
```

The input for the preprocessing pipeline is the parameter `sample_file` in
the parameter file. The sample file must be organized in the following way:
```bash
sampleID,forward,reverse
SRR9943776,</path/to/raw_data>/SRR9943776_1.fastq.gz,</path/to/raw_data>/SRR9943776_2.fastq.gz
...
```

You can then run each pipeline like this, for example the preprocessing pipeline:
```bash
nextflow run </path/to/awigen/repo>/preprocessing.nf \
	-c </path/to/awigen/repo>/config/run_preprocessing.config \
	-params-file </path/to/your/modified/params/file> \
	-with-trace -with-report
```

The other pipelines take the output from the preprocessing pipeline as input,
for example classification:
```bash
nextflow run </path/to/awigen/repo>/classification.nf \
        -c </path/to/awigen/repo>/config/run_classification.config \
        -params-file </path/to/your/modified/params/file> \
	--input </preprocessing/output/directory>/stats/preprocessed_read.csv \
        -with-trace -with-report
```

> Please note that nextflow takes around 10Gb of memory, so you might have
to run this as a job.

Alternatively, the `sbatch` scripts we used are stored in the `jobs` folder.
