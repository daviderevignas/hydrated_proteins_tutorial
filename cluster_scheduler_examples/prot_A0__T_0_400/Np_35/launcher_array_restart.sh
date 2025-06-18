#!/bin/bash
################
#
# Setting slurm options
#
################
#SBATCH --partition long
#SBATCH --time=13-00:00:00
#SBATCH --array=1-100
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks)
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type ALL
#SBATCH --job-name="MC_job"
#SBATCH --output=/home/revignas/slurm_test/output/%A_%a.out
#SBATCH --error=/home/revignas/slurm_test/errors/%A_%a.err

set -e

# paths to source folder and to output
path_to_master_dir=$(pwd)'/master_dir/*'
path_to_fin_confs=$(pwd)'/fin_confs/'
#path_to_replicas_output=$(pwd)'/replicas_output/'$SLURM_ARRAY_JOB_ID
path_to_replicas_output='/data/biophys/revignas/PKS_calculations/A0_T_0_400_sigmoid/'

# create scratch folder
scratch="/scratch/$USER/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID"
mkdir -p $scratch

# go to scratch folder
cd $scratch

# prepare scratch folder for the execution
cp -r $path_to_master_dir .
#cp $path_to_fin_confs/Final_Configuration_${SLURM_ARRAY_TASK_ID}.dat Final_Configuration.dat
#cp $path_to_fin_confs/Final_Protein_Conformation_${SLURM_ARRAY_TASK_ID}.dat Final_Protein_Conformation.dat
ls -lrt
ifort -O3 PROTEIN_FOLDING_APRIL_2021.f90 -o prot_fold_21_dr

# launch program with auto generated seed
start=`date +%s`
srun ./prot_fold_21_dr $RANDOM
end=`date +%s`
runtime=$((end-start))
echo $runtime >> output_dr.txt
echo $start >> output_dr.txt
echo $end >> output_dr.txt

# copy back results
mkdir -p $path_to_replicas_output'/'$SLURM_ARRAY_JOB_ID
mkdir -p $path_to_replicas_output'/'$SLURM_ARRAY_JOB_ID'/'$SLURM_ARRAY_TASK_ID
rm -f *.mod prot_fold_21_dr protein_lengths.dat protein_moves.dat protein_species.dat
cp "${scratch}/"* $path_to_replicas_output'/'$SLURM_ARRAY_JOB_ID'/'$SLURM_ARRAY_TASK_ID'/'


# Clean up after yourself
cd
rm -rf $scratch

# exit gracefully
exit 0
