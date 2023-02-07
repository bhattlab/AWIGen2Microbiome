# AWIGen2Microbiome
Repository for the AWI-Gen 2 microbiome project workflows, analysis, and results. 


## Run the workflow

For the preprocessing, run the workflow like that (on SCG):

```
module load java/18.0.2.1
module load nextflow/22.10.5
ssub -m 6 -t 8 -n nextflow "nextflow run main.nf -c config/run.config -params-file config/params.yml -with-trace -with-report"
```


