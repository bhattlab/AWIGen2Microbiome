# Nextflow workflows

This folder contains the workflows for analyzing metagenomic data using Nextflow as workflow manager.
The workflows are based on those in the [bhattlab_workflows Github](https://github.com/bhattlab/bhattlab_workflows)
repository.


## Important preparations

#### Adjusting the file paths in the `params.yml` file

#### Indecing the human reference genome

#### Downloading the MetaPhlAn4 database

#### Downloading the mOTUs3 database

## Running the workflows

To run the workflow, you need Java and Nextflow running in your system. On a HPC system, those are usually 
available through the `module load` system:

```
module load java/18.0.2.1
module load nextflow/22.10.5
```

On SCG, we have a dedicated script to submit a single job to the cluster:

```
ssub -m 6 -t 8 -n nextflow_preprocessing "nextflow run preprocessing.nf -c config/run.config -params-file config/params.yml -with-trace -with-report"
```
