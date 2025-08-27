#!/bin/bash  -v
#SBATCH --partition=SP2
#SBATCH --ntasks=1 		# number of tasks / mpi processes
#SBATCH --cpus-per-task=1	# Number OpenMP Threads per process
#SBATCH -J gatk
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


# module load Miniconda/3
# source activate micromamba
# micromamba activate gatk4



reference=reference/Af293_GCF_000002655.1.fna

#  create reference
gatk CreateSequenceDictionary -R $reference
samtools faidx $reference


for bam in $(ls MAPs/*_sorted.bam); do
  name=$(basename $i _sorted.bam)
  
#  Format headers of BAM files
  gatk AddOrReplaceReadGroups \
  INPUT=$bam \
  OUTPUT=MAPs/$name"_fix.bam" \
  RGID=H0164.2 \
  RGLB= library_$name \
  RGPL=illumina \
  RGPU=H0164ALXX140820.2 \
  RGSM=$name
  samtools index MAPs/$name"_fix.bam"


# haplotypecaller
  gatk HaplotypeCaller -I 0_MAPs/$name"_fix.bam" -ploidy 1 --reference $reference \
  --annotation QualByDepth \
  --annotation FisherStrand \
  --annotation MappingQuality \
  --annotation MappingQualityRankSumTest \
  --annotation ReadPosRankSumTest \
  --annotation StrandOddsRatio \
  -O haplotypes/$name".vcf" -ERC GVCF
  
  
# genotype
  gatk GenotypeGVCFs -V haplotypes/$name".vcf" -O genotypes/$name"_genotype.vcf" --reference $reference


# filter variants
  gatk VariantFiltration -R $reference \
   -V genotypes/$name"_genotype.vcf" \
   -O filtered/$name"_filt.vcf" \
   --filter-expression "QD < 2.0" --filter-name "QualByDepth" \
   --filter-expression "MQRankSum < -12.5" --filter-name "MappingQualityRankSumTest" \
   --filter-expression "FS > 8.0" --filter-name "FisherStrand" \
   --filter-expression "SOR > 3.0" --filter-name "StrandOddsRatio" \
   --filter-expression "ReadPosRankSum < -2.0" --filter-name "ReadPosRankSumTest" \
   --filter-expression "ReadPosRankSum > 2.0" --filter-name "ReadPosRankSumTest" \
   --filter-expression "MQ < 4.0" --filter-name "MappingQuality"


# exclude low quality variants
  gatk SelectVariants -R $reference -V filtered/$name"_filt.vcf" --exclude-filtered true -O filtered/$name"_final.vcf"
  rm -rf filtered/$name"_filt.vcf"
  rm -rf filtered/$name"_filt.vcf.idx"

done


###### ---------------------------- combined VCF files pipeline for the phylogenetic tree and haplotypes network analyses ---------------------------#########

## Combine VCF files
gatk CombineGVCFs --reference $reference \
--variant haplotypes/LIF3297.vcf \
--variant haplotypes/LIF3365.vcf \
--variant haplotypes/LIF3495.vcf \
--variant haplotypes/LIF3519.vcf \
--variant haplotypes/LIF3545.vcf \
--variant haplotypes/LIF3546.vcf \
--variant haplotypes/LIF3309.vcf \
--variant haplotypes/LIF3492.vcf \
--variant haplotypes/LIF3608.vcf \
-O combined/combined.vcf

gatk GenotypeGVCFs -V combined.vcf -O combined/combined_genotyped.vcf --reference $reference

## filter combined VCF file
gatk VariantFiltration -R $reference \
   -V combined/combined_genotyped.vcf \
   -O combined/combined_filt.vcf \
   --filter-expression "QD < 2.0" --filter-name "QualByDepth" \
   --filter-expression "MQRankSum < -12.5" --filter-name "MappingQualityRankSumTest" \
   --filter-expression "FS > 8.0" --filter-name "FisherStrand" \
   --filter-expression "SOR > 3.0" --filter-name "StrandOddsRatio" \
   --filter-expression "ReadPosRankSum < -2.0" --filter-name "ReadPosRankSumTest" \
   --filter-expression "ReadPosRankSum > 2.0" --filter-name "ReadPosRankSumTest" \
   --filter-expression "MQ < 4.0" --filter-name "MappingQuality"

gatk SelectVariants -R $reference -V combined/$name"_filt.vcf" --exclude-filtered true -O combined/combined_final.vcf
rm -rf combined/combined_filt.vcf
rm -rf combined/combined_filt.vcf.idx

