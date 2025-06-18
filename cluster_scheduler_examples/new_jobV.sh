#!/bin/bash
# The name of the job, can be anything, simply used when displaying the list of running jobs
#$ -N PROTEIN_FOLDING
# Giving the name of the output log file
#$ -o PROTEIN_FOLDING-$JOB_ID.log
# Combining output/error messages into one file
##$ -j y
# One needs to tell the queue system to use the current directory as the working directory
# Or else the script may fail as it will execute in your top level home directory /home/username
#$ -cwd
# With -V you pass the env variables, it's necessary. And the unset module is needed to remove some errors
#$ -V
# Uncomment the following line if you want to know in which host your job was executed
# echo "Running on " `hostname`
# Now comes the commands to be executed
# Copy exe and required input files to the local disk on the node

# We create a subdirectory in TMP to separate different versions of the same JOB_NAME running in the same node 
replica_execution_dir="$TMPDIR/${JOB_ID}/"

# We create a subdirectory to separate output files generated at differents nodes or different JOB IDs 
replica_out_dir="$SGE_O_WORKDIR/outputs/$HOSTNAME/${JOB_ID}/"

mkdir -p $replica_execution_dir
mkdir -p $replica_output_dir

# Since we are at the compute node, is easier to create the seed and we can already report what number we used
test_seed="$RANDOM"

echo "${test_seed}" > $replica_output_dir/seed_value.txt


prot_fold_exe_file="prot_fold_${HOSTNAME}_${JOB_ID}.exe"
# We store all input filenames in a variable to avoid using '*'
ifort -O3 PROTEIN_FOLDING_APRIL_2021.f90 -o $prot_fold_exe_file
replica_input="PROTEIN_FOLDING_APRIL_2021.f90 aapot_water.dat input_data_folding"
replica_input="${replica_input} input_sequences.dat protein_lengths.dat protein_species.dat protein_moves.dat"
replica_input="${replica_input} target_structures.dat water_parameters $prot_fold_exe_file"


# IMPORTANT VARIABLE, DON'T RUN UNTIL PROPERLY SET
# If we can, we do the same with the files created at runtime 
replica_output_files="*.dat *.out *.txt *.log" 

cp $replica_input $replica_execution_dir

# Change to the execution directory
cd $replica_execution_dir
# Now that we are at the compute node, is easier to create the seed
sed -i "s/_SEED_/${test_seed}/g" input_data_folding 

# And run the exe
./"$prot_fold_exe_file"

output_tarfile="${HOSTNAME}_${JOB_ID}_output.tar"
tar -cf  $output_tarfile $replica_output_files
# Finally, we copy back all important output to the working directory
scp $output_tarfile nodo00:$replica_output_dir
cd $SGE_O_WORKDIR
rm -rf $replica_execution_dir

