# Introductory README for the CUT&RUNnextflow GitHub

This repository houses the required materials to run the CUT&RUN pipeline using the https://tower.nf/ webtool.

Alternatively you can acquire the files manually in the directory named Nextflow, which should contain the following: 
- **C&R_pipeline.nf** : This is the effective pipeline, to be invoked by the launch script or manually using\"\$nextflow run C&R_pipeline.nf\". 
- **Launch_C&R_pipeline.sh** : A very simple shell script that launches the pipeline using the SLURM SBATCH-method. 
- **README.md** : A different README file with a quick start guide for basic manual usage of the pipeline. 
- **nextflow.config** : The config file containing the necessary parameters for the pipeline to function, usually you will have to adapt this one to some extent.

These 4 files should be all you need to use this relatively simple pipeline that performs the following steps: 
1. Data organization and sample renaming based on a configurable csv file 
2. Trimming of the data 
3. Mapping of the data 
4. Evaluating the data 
5. Conditionally deduplicating the data 
6. Filtering the data using BEDtools 
7. Creating subdata by splitting at a configurable sequence length
     - Note: This part of the pipeline is outdated but has been retained for historical reasons.
     - Future updates may involve replacing this part of the pipeline with more modern techniques.
9. Creating IGV files based on the data 
10. Performing RPKM normalization on the data 
11. Performing specifiable peakcalling on the data according to a configurable csv file 
12. Annotating the narrow peak data using Homer 
13. Finding motifs in the narrow peak data using Homer
14. Summarising the run in a simple MultiQC

Since Nextflow is capable of asynchronous jobhandling it will submit these jobs whenever the required inputs are supplied.  
This means it can perform various tasks at the same time, as well as flexibly run multiple samples at the same time.
It can also choose to make use of both the hg19 and hg38 reference genome at the same time.

![Diagram of the nextflow pipeline](./Cut_and_run_pipeline_diagram.jpg)
