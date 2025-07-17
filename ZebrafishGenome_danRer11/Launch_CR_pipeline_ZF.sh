#!/bin/bash
#SBATCH --mem 4G
#SBATCH -t 24:00:00
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --open-mode=truncate

WORKFLOW="$VSC_DATA_VO/PPOL/resources/sharedData/PUBLIC_CUTandRUN_DATASETS/GollMG_HistoneMarks_Zebrafish_GSE178343/Scripts/CR_pipeline_ZF.nf"
CONFIG="$VSC_DATA_VO/PPOL/resources/sharedData/PUBLIC_CUTandRUN_DATASETS/GollMG_HistoneMarks_Zebrafish_GSE178343/Scripts/nextflow_ZF.config"

module purge
module load Nextflow/23.04.2

# $1 is the runname which should be the name of the input directory.
# This runname should be found at the end of your inputDir variable so path/to/inputDir/RunName

# 2 is the pipeline mode where you tell it to use GRCz11
# If no argument was supplied it'll use GRCz11.102 by default.-> remove this!

# $3 is a reserved slot for any additional nextflow argument you may want to supply.
# Usually this will be -resume to restart the pipeline from it's last cached position.
nextflow -C ${CONFIG} run ${WORKFLOW} --RunName $1 -with-apptainer
