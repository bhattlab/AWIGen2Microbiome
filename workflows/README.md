# Nextflow workflows

This folder contains the workflows for analyzing metagenomic data 
using Nextflow as workflow manager. The workflows are based on those in the 
[bhattlab_workflows Github](https://github.com/bhattlab/bhattlab_workflows)
repository.


## Important preparations

#### Adjusting the `params.yml` file

The `params.yml` file holds all the important settings for running the 
workflows on your system. You will need to adjust the file paths to the 
databases (see below for downloading those to your file system). 
You can also adjust the run parameters for some tools, such as the parameters
to run mOTUs, for example. The `params.yml` file should hold some explanation
for each of the parameters.

#### Indexing the human reference genome

You will need a human reference genome on your file system so that you can 
remove all reads matching to the human genome from your sample, since we are
interested in the bacteria, not the human genome. Since we map with the `bwa`
algorithm, you will need to index the human genome before being able to use
the workflow.
We will download the genome from 
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

#### Downloading the MetaPhlAn4 database

We also need to download the MetaPhlAn4 database, again using the already
prepared singularity container:

```bash
cd <your-metaphlan-database-location>
singularity shell --bind ./:/mnt docker://ghcr.io/jakob-wirbel/micromamba-focal-classification:latest
cd /mnt
metaphlan --install --bowtie2db ./
```

#### Downloading the mOTUs3 database

The mOTU tool needs a custom database of marker genes. Unfortunately, you 
cannot download the database through the tool as of now, but maybe it will be
available as a feautre in the future (see 
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

#### Phanta database

For running [Phanta](https://github.com/bhattlab/phanta), please note that 
Phanta was developed as a Snakemake workflow and the port to Nextflow is a
bit hacky. As of yet, it is pretty unclear if we will want to include it in 
the final analysis. If you want to run it, you will need a Kraken2 database. 
[Here](https://github.com/bhattlab/phanta/blob/main/databases.md) is a list of
databases that are made available with the tool. You can download and extract
the tarball following the documentation for the tool.

## Running the workflows

To run the workflow, you need Java and Nextflow running in your system. On an 
HPC system, those are usually available through the `module load` system:

```bash
module load java/18.0.2.1
module load nextflow/22.10.5
```

On SCG, we have a dedicated script to submit a single job to the cluster:

```bash
ssub -m 6 -t 8 -n nextflow_preprocessing "nextflow run preprocessing.nf -c config/run_preprocessing.config -params-file config/params.yml -with-trace -with-report"
```

The other workflows can then be run after the preprocessing step is done:
```bash
ssub -m 6 -t 8 -n nextflow_classification "nextflow run classification.nf -c config/run_classification.config -params-file config/params.yml -with-trace -with-report --input <path-to-stats-read-csv>"
ssub -m 6 -t 8 -n nextflow_assembly "nextflow run assembly_binning.nf -c config/run_assembly.config -params-file config/params.yml -with-trace -with-report --input <path-to-stats-read-csv>"
```
