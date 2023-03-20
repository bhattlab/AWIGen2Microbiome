#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --mem=12G
#SBATCH --account=asbhatt
#SBATCH --partition=batch
#SBATCH --output=/labs/asbhatt/wirbel/SCRATCH/classification_nf.out
#SBATCH --error=/labs/asbhatt/wirbel/SCRATCH/classification_nf.err

module load java/18.0.2.1
module load nextflow/22.10.5

nextflow run /labs/asbhatt/wirbel/AWIgen2/AWIGen2Microbiome/workflows/classification.nf \
	-c /labs/asbhatt/wirbel/AWIgen2/AWIGen2Microbiome/workflows/config/run_classification.config \
	-params-file /labs/asbhatt/wirbel/AWIgen2/AWIGen2Microbiome/workflows/config/params.yml \
	--input /labs/asbhatt/data/bhatt_lab_sequencing/23-03-08_awigen2/stats/preprocessed_reads.csv \
	-with-trace -with-report -resume

