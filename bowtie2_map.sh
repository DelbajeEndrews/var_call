#!/bin/bash  -v
#SBATCH --partition=SP2
#SBATCH --ntasks=1 		# number of tasks / mpi processes
#SBATCH --cpus-per-task=5	# Number OpenMP Threads per process
#SBATCH -J bowtie2
#SBATCH --time=180:00:00

#OpenMP settings:
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

echo $SLURM_JOB_ID		#ID of job allocation
echo $SLURM_SUBMIT_DIR		#Directory job where was submitted
echo $SLURM_JOB_NODELIST	#File containing allocated hostnames
echo $SLURM_NTASKS		#Total number of cores for job


# source activate maptools

threads=5

bowtie2-build reference/Af293_GCF_000002655.1.fa reference/reference

for r1 in $(ls trimmed_reads/*r1.fq.gz); do
  r2=${r1/r1.fq.gz/r2.fq.gz}
  name=$(basename $r1 _r1.fq.gz)
  bowtie2 -x reference/reference -1 $r1 -2 $r2 -p $threads > $name.sam
  samtools view -bhS $name.sam -@ $threads > $name.bam
  samtools sort $name.bam -o $name"_sorted.bam" -@ $threads
done
