#!/bin/bash
#SBATCH --mem 4G
#SBATCH -t 24:00:00
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --open-mode=truncate

WORKFLOW="$VSC_DATA_VO/PPOL/resources/sharedData/CUTandRUN_BDC_SLB/Scripts/Scripts_SLB/CR_pipeline.nf"
CONFIG="$VSC_DATA_VO/PPOL/resources/sharedData/CUTandRUN_BDC_SLB/Scripts/Scripts_SLB/nextflow.config"

module purge
module load  Nextflow

# $1 is the runname which should be the name of the input directory. 
# This runname should be found at the end of your inputDir variable so path/to/inputDir/RunName

# 2 is the pipeline mode where you tell it to use hg19, hg38 or both.
# If no argument was supplied it'll use hg19 by default.

# $3 is a reserved slot for any additional nextflow argument you may want to supply. 
# Usually this will be -resume to restart the pipeline from it's last cached position.
nextflow -C ${CONFIG} run ${WORKFLOW} --RunName $1 --pipelineMode ${2:-hg19} $3  -with-apptainer

