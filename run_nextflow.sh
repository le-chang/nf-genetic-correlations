#!/bin/bash
#SBATCH --account=def-gsarah
#SBATCH --job-name=nextflow_run
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --time=47:59:59
#SBATCH --mem=20G

module load StdEnv/2023
module load java/21.0.1
module load nextflow/24.10.2
module load r/4.3.1
module load apptainer/1.3.5

export R_LIBS=~/.local/R/4.3.1/

export NXF_DISABLE_REMOTE_WORKFLOW=true
export NXF_DISABLE_CHECK_LATEST=true

# Use SLURM_SUBMIT_DIR which contains the directory where sbatch was called
WORK_DIR="${SLURM_SUBMIT_DIR}"

# Run nextflow from the submission directory
cd "${WORK_DIR}"
nextflow run main_full.nf -profile beluga -resume