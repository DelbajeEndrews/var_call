#!/bin/bash  -v
#SBATCH --partition=SP2
#SBATCH --ntasks=1 		# number of tasks / mpi processes
#SBATCH --cpus-per-task=12	# Number OpenMP Threads per process
#SBATCH -J trim
#SBATCH --time=150:00:00

#OpenMP settings:
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

echo $SLURM_JOB_ID		#ID of job allocation
echo $SLURM_SUBMIT_DIR		#Directory job where was submitted
echo $SLURM_JOB_NODELIST	#File containing allocated hostnames
echo $SLURM_NTASKS		#Total number of cores for job


# Run trimmomatic
for r1 in $(ls isolates/*gz); do
  r2=${r1/_1.fq.gz/_2.fq.gz}
  name=$(echo $r1 | cut -f2 -d "/")
  r1_out=$name"_r1.fq"
  r2_out=$name"_r2.fq"
  r1_out_un=$name"_r1_un.fq"
  r2_out_un=$name"_r2_un.fq"
  java -jar trimmomatic-0.39.jar PE -threads 12 $r1 $r2 $r1_out $r1_out_un $r2_out $r2_out_un ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 MINLEN:60
done
