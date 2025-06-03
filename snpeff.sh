#!/bin/bash  -v
#SBATCH --partition=SP2
#SBATCH --ntasks=1 		# number of tasks / mpi processes
#SBATCH --cpus-per-task=1	# Number OpenMP Threads per process
#SBATCH -J snpeff
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



# micromamba activate snpeff

# build snpeff database
mkdir /Storage/user/programs/snpeff-5.2-0/data/Af293
mv Af293_GCF_000002655.1.fna /Storage/user/programs/snpeff-5.2-0/data/Af293/sequences.fa
mv Af293_GCF_000002655.1.gtf /Storage/user/programs/snpeff-5.2-0/data/Af293/genes.gtf
echo "Af293.genome : Af293" >> /Storage/user/programs/snpeff-5.2-0/snpEff.config

snpEff build -gtf22 -v /Storage/user/programs/snpeff-5.2-0/data/Af293 -noCheckCds -noCheckProtein


# Annotate variants
for i in $(ls filtered/*_final.vcf); do
  name=$(basename $i _final.vcf)
  java -jar snpEff.jar Af293 $i > snpeff/$name".vcf"
  mv snpEff_genes.txt snpeff/$name"_genes.txt"
  rm -rf snpEff_summary.html
done
