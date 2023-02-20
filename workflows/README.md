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

To download and unpack the genome, type:
```bash
cd <your-reference-genome-location>
wget
```

Now, you can index the genome with `bwa`. However, if you want to use the 
singularity container that we use for the workflow (because you might not have
`bwa` installed on your system), you can do so like that:
```bash
cd <your-reference-genome-location>
singularity shell
```



#### Downloading the MetaPhlAn4 database



#### Downloading the mOTUs3 database



## Running the workflows

To run the workflow, you need Java and Nextflow running in your system. On an 
HPC system, those are usually available through the `module load` system:

```
module load java/18.0.2.1
module load nextflow/22.10.5
```

On SCG, we have a dedicated script to submit a single job to the cluster:

```
ssub -m 6 -t 8 -n nextflow_preprocessing "nextflow run preprocessing.nf -c config/run_preprocessing.config -params-file config/params.yml -with-trace -with-report"
```
