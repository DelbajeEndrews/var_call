#!/bin/bash  -v
#SBATCH --partition=SP2
#SBATCH --ntasks=1 		# number of tasks / mpi processes
#SBATCH --cpus-per-task=3	# Number OpenMP Threads per process
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

# Combine VFC files (oublic data and serial isolates)

gatk CombineGVCFs --reference $reference \
--variant haplotypes/DRR237582.vcf \
--variant haplotypes/DRR237583.vcf \
--variant haplotypes/ERR769498.vcf \
--variant haplotypes/ERR769499.vcf \
--variant haplotypes/ERR769500.vcf \
--variant haplotypes/ERR769501.vcf \
--variant haplotypes/ERR769502.vcf \
--variant haplotypes/ERR769503.vcf \
--variant haplotypes/ERR769504.vcf \
--variant haplotypes/ERR769505.vcf \
--variant haplotypes/ERR769506.vcf \
--variant haplotypes/ERR769507.vcf \
--variant haplotypes/ERR769508.vcf \
--variant haplotypes/ERR769509.vcf \
--variant haplotypes/ERR769510.vcf \
--variant haplotypes/ERR769511.vcf \
--variant haplotypes/ERR769512.vcf \
--variant haplotypes/ERR769513.vcf \
--variant haplotypes/ERR769514.vcf \
--variant haplotypes/ERR769515.vcf \
--variant haplotypes/ERR769516.vcf \
--variant haplotypes/ERR769517.vcf \
--variant haplotypes/ERR769518.vcf \
--variant haplotypes/ERR769519.vcf \
--variant haplotypes/ERR769520.vcf \
--variant haplotypes/ERR769521.vcf \
--variant haplotypes/ERR9791571.vcf \
--variant haplotypes/ERR9791572.vcf \
--variant haplotypes/ERR9791573.vcf \
--variant haplotypes/ERR9791574.vcf \
--variant haplotypes/ERR9791575.vcf \
--variant haplotypes/ERR9791576.vcf \
--variant haplotypes/ERR9791577.vcf \
--variant haplotypes/ERR9791578.vcf \
--variant haplotypes/ERR9791579.vcf \
--variant haplotypes/ERR9791580.vcf \
--variant haplotypes/ERR9791581.vcf \
--variant haplotypes/ERR9791582.vcf \
--variant haplotypes/ERR9791583.vcf \
--variant haplotypes/ERR9791584.vcf \
--variant haplotypes/ERR9791585.vcf \
--variant haplotypes/ERR9791586.vcf \
--variant haplotypes/ERR9791587.vcf \
--variant haplotypes/ERR9791588.vcf \
--variant haplotypes/ERR9791589.vcf \
--variant haplotypes/ERR9791590.vcf \
--variant haplotypes/ERR9791591.vcf \
--variant haplotypes/ERR9791592.vcf \
--variant haplotypes/ERR9791593.vcf \
--variant haplotypes/ERR9791594.vcf \
--variant haplotypes/ERR9791595.vcf \
--variant haplotypes/ERR9791596.vcf \
--variant haplotypes/ERR9791597.vcf \
--variant haplotypes/ERR9791598.vcf \
--variant haplotypes/ERR9791599.vcf \
--variant haplotypes/ERR9791600.vcf \
--variant haplotypes/ERR9791601.vcf \
--variant haplotypes/ERR9791602.vcf \
--variant haplotypes/ERR9791603.vcf \
--variant haplotypes/ERR9791604.vcf \
--variant haplotypes/ERR9791605.vcf \
--variant haplotypes/ERR9791606.vcf \
--variant haplotypes/ERR9791607.vcf \
--variant haplotypes/ERR9791608.vcf \
--variant haplotypes/ERR9791609.vcf \
--variant haplotypes/ERR9791610.vcf \
--variant haplotypes/ERR9791611.vcf \
--variant haplotypes/ERR9791612.vcf \
--variant haplotypes/ERR9791613.vcf \
--variant haplotypes/ERR9791614.vcf \
--variant haplotypes/ERR9791615.vcf \
--variant haplotypes/ERR9791616.vcf \
--variant haplotypes/ERR9791617.vcf \
--variant haplotypes/ERR9791618.vcf \
--variant haplotypes/ERR9791619.vcf \
--variant haplotypes/ERR9791620.vcf \
--variant haplotypes/ERR9791621.vcf \
--variant haplotypes/ERR9791622.vcf \
--variant haplotypes/ERR9791623.vcf \
--variant haplotypes/ERR9791624.vcf \
--variant haplotypes/ERR9791625.vcf \
--variant haplotypes/ERR9791626.vcf \
--variant haplotypes/ERR9791627.vcf \
--variant haplotypes/ERR9791628.vcf \
--variant haplotypes/ERR9791629.vcf \
--variant haplotypes/ERR9791630.vcf \
--variant haplotypes/ERR9791631.vcf \
--variant haplotypes/ERR9791632.vcf \
--variant haplotypes/ERR9791633.vcf \
--variant haplotypes/ERR9791634.vcf \
--variant haplotypes/ERR9791635.vcf \
--variant haplotypes/ERR9791636.vcf \
--variant haplotypes/ERR9791637.vcf \
--variant haplotypes/ERR9791638.vcf \
--variant haplotypes/ERR9791639.vcf \
--variant haplotypes/ERR9791640.vcf \
--variant haplotypes/ERR9791641.vcf \
--variant haplotypes/ERR9791642.vcf \
--variant haplotypes/ERR9791643.vcf \
--variant haplotypes/ERR9791644.vcf \
--variant haplotypes/ERR9791645.vcf \
--variant haplotypes/ERR9791646.vcf \
--variant haplotypes/ERR9791650.vcf \
--variant haplotypes/ERR9791651.vcf \
--variant haplotypes/ERR9791652.vcf \
--variant haplotypes/ERR9791653.vcf \
--variant haplotypes/ERR9791664.vcf \
--variant haplotypes/ERR9791665.vcf \
--variant haplotypes/ERR9791666.vcf \
--variant haplotypes/ERR9791667.vcf \
--variant haplotypes/ERR9791668.vcf \
--variant haplotypes/ERR9791669.vcf \
--variant haplotypes/ERR9791670.vcf \
--variant haplotypes/ERR9791671.vcf \
--variant haplotypes/ERR9791672.vcf \
--variant haplotypes/ERR9791673.vcf \
--variant haplotypes/ERR9791674.vcf \
--variant haplotypes/ERR9791675.vcf \
--variant haplotypes/ERR9791676.vcf \
--variant haplotypes/ERR9791677.vcf \
--variant haplotypes/ERR9791678.vcf \
--variant haplotypes/ERR9791679.vcf \
--variant haplotypes/ERR9791680.vcf \
--variant haplotypes/ERR9791681.vcf \
--variant haplotypes/ERR9791682.vcf \
--variant haplotypes/ERR9791683.vcf \
--variant haplotypes/ERR9791684.vcf \
--variant haplotypes/ERR9791685.vcf \
--variant haplotypes/ERR9791686.vcf \
--variant haplotypes/ERR9791687.vcf \
--variant haplotypes/ERR9791688.vcf \
--variant haplotypes/ERR9791689.vcf \
--variant haplotypes/ERR9791690.vcf \
--variant haplotypes/ERR9791691.vcf \
--variant haplotypes/ERR9791692.vcf \
--variant haplotypes/ERR9791693.vcf \
--variant haplotypes/ERR9791694.vcf \
--variant haplotypes/ERR9791695.vcf \
--variant haplotypes/ERR9791696.vcf \
--variant haplotypes/ERR9791697.vcf \
--variant haplotypes/ERR9791698.vcf \
--variant haplotypes/ERR9791699.vcf \
--variant haplotypes/ERR9791700.vcf \
--variant haplotypes/ERR9791701.vcf \
--variant haplotypes/ERR9791702.vcf \
--variant haplotypes/ERR9791703.vcf \
--variant haplotypes/ERR9791704.vcf \
--variant haplotypes/ERR9791705.vcf \
--variant haplotypes/ERR9791706.vcf \
--variant haplotypes/ERR9791707.vcf \
--variant haplotypes/ERR9791708.vcf \
--variant haplotypes/ERR9791709.vcf \
--variant haplotypes/ERR9791710.vcf \
--variant haplotypes/ERR9791711.vcf \
--variant haplotypes/ERR9791712.vcf \
--variant haplotypes/ERR9791713.vcf \
--variant haplotypes/ERR9791714.vcf \
--variant haplotypes/ERR9791715.vcf \
--variant haplotypes/ERR9791716.vcf \
--variant haplotypes/ERR9791717.vcf \
--variant haplotypes/ERR9791718.vcf \
--variant haplotypes/ERR9791719.vcf \
--variant haplotypes/ERR9791720.vcf \
--variant haplotypes/ERR9791721.vcf \
--variant haplotypes/ERR9791722.vcf \
--variant haplotypes/ERR9791723.vcf \
--variant haplotypes/ERR9791724.vcf \
--variant haplotypes/ERR9791725.vcf \
--variant haplotypes/ERR9791726.vcf \
--variant haplotypes/ERR9791727.vcf \
--variant haplotypes/ERR9791728.vcf \
--variant haplotypes/ERR9791729.vcf \
--variant haplotypes/ERR9791730.vcf \
--variant haplotypes/ERR9791731.vcf \
--variant haplotypes/ERR9791732.vcf \
--variant haplotypes/ERR9791733.vcf \
--variant haplotypes/ERR9791734.vcf \
--variant haplotypes/ERR9791735.vcf \
--variant haplotypes/ERR9791736.vcf \
--variant haplotypes/ERR9791737.vcf \
--variant haplotypes/ERR9791738.vcf \
--variant haplotypes/ERR9791739.vcf \
--variant haplotypes/ERR9791740.vcf \
--variant haplotypes/ERR9791741.vcf \
--variant haplotypes/ERR9791743.vcf \
--variant haplotypes/ERR9791744.vcf \
--variant haplotypes/ERR9791745.vcf \
--variant haplotypes/ERR9791746.vcf \
--variant haplotypes/ERR9791747.vcf \
--variant haplotypes/ERR9791748.vcf \
--variant haplotypes/ERR9791749.vcf \
--variant haplotypes/ERR9791750.vcf \
--variant haplotypes/ERR9791751.vcf \
--variant haplotypes/ERR9791752.vcf \
--variant haplotypes/ERR9791753.vcf \
--variant haplotypes/ERR9791754.vcf \
--variant haplotypes/ERR9791755.vcf \
--variant haplotypes/ERR9791756.vcf \
--variant haplotypes/ERR9791757.vcf \
--variant haplotypes/ERR9791758.vcf \
--variant haplotypes/ERR9791759.vcf \
--variant haplotypes/ERR9791760.vcf \
--variant haplotypes/ERR9791761.vcf \
--variant haplotypes/ERR9791762.vcf \
--variant haplotypes/ERR9791763.vcf \
--variant haplotypes/ERR9791764.vcf \
--variant haplotypes/ERR9791765.vcf \
--variant haplotypes/ERR9791766.vcf \
--variant haplotypes/ERR9791767.vcf \
--variant haplotypes/ERR9791768.vcf \
--variant haplotypes/ERR9791769.vcf \
--variant haplotypes/ERR9791772.vcf \
--variant haplotypes/ERR9791775.vcf \
--variant haplotypes/ERR9791776.vcf \
--variant haplotypes/ERR9791777.vcf \
--variant haplotypes/ERR9791778.vcf \
--variant haplotypes/ERR9791779.vcf \
--variant haplotypes/ERR9791780.vcf \
--variant haplotypes/ERR9791781.vcf \
--variant haplotypes/ERR9791782.vcf \
--variant haplotypes/SRR10592629.vcf \
--variant haplotypes/SRR10592631.vcf \
--variant haplotypes/SRR10592632.vcf \
--variant haplotypes/SRR10592633.vcf \
--variant haplotypes/SRR11977809.vcf \
--variant haplotypes/SRR11977810.vcf \
--variant haplotypes/SRR11977811.vcf \
--variant haplotypes/SRR11977812.vcf \
--variant haplotypes/SRR11977813.vcf \
--variant haplotypes/SRR11977814.vcf \
--variant haplotypes/SRR11977815.vcf \
--variant haplotypes/SRR11977816.vcf \
--variant haplotypes/SRR11977817.vcf \
--variant haplotypes/SRR11977818.vcf \
--variant haplotypes/SRR11977819.vcf \
--variant haplotypes/SRR11977820.vcf \
--variant haplotypes/SRR11977821.vcf \
--variant haplotypes/SRR11977822.vcf \
--variant haplotypes/SRR11977823.vcf \
--variant haplotypes/SRR11977824.vcf \
--variant haplotypes/SRR11977825.vcf \
--variant haplotypes/SRR11977826.vcf \
--variant haplotypes/SRR11977827.vcf \
--variant haplotypes/SRR11977828.vcf \
--variant haplotypes/SRR11977829.vcf \
--variant haplotypes/SRR11977830.vcf \
--variant haplotypes/SRR11977831.vcf \
--variant haplotypes/SRR11977832.vcf \
--variant haplotypes/SRR11977833.vcf \
--variant haplotypes/SRR11977834.vcf \
--variant haplotypes/SRR11977835.vcf \
--variant haplotypes/SRR11977836.vcf \
--variant haplotypes/SRR11977837.vcf \
--variant haplotypes/SRR11977838.vcf \
--variant haplotypes/SRR11977839.vcf \
--variant haplotypes/SRR11977840.vcf \
--variant haplotypes/SRR11977841.vcf \
--variant haplotypes/SRR11977842.vcf \
--variant haplotypes/SRR11977843.vcf \
--variant haplotypes/SRR11977844.vcf \
--variant haplotypes/SRR11977845.vcf \
--variant haplotypes/SRR11977846.vcf \
--variant haplotypes/SRR11977847.vcf \
--variant haplotypes/SRR11977848.vcf \
--variant haplotypes/SRR11977849.vcf \
--variant haplotypes/SRR11977850.vcf \
--variant haplotypes/SRR11977851.vcf \
--variant haplotypes/SRR11977852.vcf \
--variant haplotypes/SRR11977853.vcf \
--variant haplotypes/SRR11977854.vcf \
--variant haplotypes/SRR11977855.vcf \
--variant haplotypes/SRR11977856.vcf \
--variant haplotypes/SRR11977857.vcf \
--variant haplotypes/SRR11977858.vcf \
--variant haplotypes/SRR11977859.vcf \
--variant haplotypes/SRR11977860.vcf \
--variant haplotypes/SRR11977861.vcf \
--variant haplotypes/SRR11977862.vcf \
--variant haplotypes/SRR11977863.vcf \
--variant haplotypes/SRR11977864.vcf \
--variant haplotypes/SRR11977865.vcf \
--variant haplotypes/SRR11977866.vcf \
--variant haplotypes/SRR12190174.vcf \
--variant haplotypes/SRR12190175.vcf \
--variant haplotypes/SRR12190176.vcf \
--variant haplotypes/SRR12190177.vcf \
--variant haplotypes/SRR12190178.vcf \
--variant haplotypes/SRR12190179.vcf \
--variant haplotypes/SRR12190180.vcf \
--variant haplotypes/SRR12190181.vcf \
--variant haplotypes/SRR12190182.vcf \
--variant haplotypes/SRR12190183.vcf \
--variant haplotypes/SRR17660940.vcf \
--variant haplotypes/SRR17660941.vcf \
--variant haplotypes/SRR17660942.vcf \
--variant haplotypes/SRR17660943.vcf \
--variant haplotypes/SRR17660944.vcf \
--variant haplotypes/SRR17660945.vcf \
--variant haplotypes/SRR17660946.vcf \
--variant haplotypes/SRR17660948.vcf \
--variant haplotypes/SRR17660950.vcf \
--variant haplotypes/SRR17660951.vcf \
--variant haplotypes/SRR17660952.vcf \
--variant haplotypes/SRR17660953.vcf \
--variant haplotypes/SRR17660954.vcf \
--variant haplotypes/SRR17660955.vcf \
--variant haplotypes/SRR17660956.vcf \
--variant haplotypes/SRR24308377.vcf \
--variant haplotypes/SRR24318406.vcf \
--variant haplotypes/SRR24318407.vcf \
--variant haplotypes/SRR24318408.vcf \
--variant haplotypes/SRR24318409.vcf \
--variant haplotypes/SRR24318410.vcf \
--variant haplotypes/SRR24318411.vcf \
--variant haplotypes/SRR24318412.vcf \
--variant haplotypes/SRR24318413.vcf \
--variant haplotypes/SRR9265316.vcf \
--variant haplotypes/SRR9265318.vcf \
--variant haplotypes/IFM_59355-1.vcf \
--variant haplotypes/IFM_59355-2.vcf \
--variant haplotypes/IFM_59356-1.vcf \
--variant haplotypes/IFM_59356-2.vcf \
--variant haplotypes/IFM_59356-3.vcf \
--variant haplotypes/IFM_59361-1.vcf \
--variant haplotypes/IFM_59361-2.vcf \
--variant haplotypes/IFM_60237.vcf \
--variant haplotypes/V130-15.vcf \
--variant haplotypes/V130-14.vcf \
--variant haplotypes/V130-18.vcf \
--variant haplotypes/V130-54.vcf \
--variant haplotypes/V157-39.vcf \
--variant haplotypes/V157-40.vcf \
--variant haplotypes/V157-47.vcf \
--variant haplotypes/V157-48.vcf \
--variant haplotypes/V157-62.vcf \
--variant haplotypes/V157-59.vcf \
--variant haplotypes/V157-60.vcf \
--variant haplotypes/V157-61.vcf \
--variant haplotypes/V157-80.vcf \
--variant haplotypes/LIF3297.vcf \
--variant haplotypes/LIF3365.vcf \
--variant haplotypes/LIF3495.vcf \
--variant haplotypes/LIF3519.vcf \
--variant haplotypes/LIF3545.vcf \
--variant haplotypes/LIF3546.vcf \
--variant haplotypes/LIF3309.vcf \
--variant haplotypes/LIF3492.vcf \
--variant haplotypes/LIF3608.vcf \
-O combined.vcf

gatk GenotypeGVCFs -V combined.vcf -O combined_genotyped.vcf --reference $reference

## filter combined VCF file
gatk VariantFiltration -R $reference \
   -V combined_genotyped.vcf \
   -O combined_filt.vcf \
   --filter-expression "QD < 2.0" --filter-name "QualByDepth" \
   --filter-expression "MQRankSum < -12.5" --filter-name "MappingQualityRankSumTest" \
   --filter-expression "FS > 8.0" --filter-name "FisherStrand" \
   --filter-expression "SOR > 3.0" --filter-name "StrandOddsRatio" \
   --filter-expression "ReadPosRankSum < -2.0" --filter-name "ReadPosRankSumTest" \
   --filter-expression "ReadPosRankSum > 2.0" --filter-name "ReadPosRankSumTest" \
   --filter-expression "MQ < 4.0" --filter-name "MappingQuality"

gatk SelectVariants -R $reference -V combined/$name"_filt.vcf" --exclude-filtered true -O combined/combined_final.vcf
rm -rf combined_filt.vcf
rm -rf combined_filt.vcf.idx

# convert SNP data (VCF) to alignment (FASTA)
python3 /Storage/user/programs/vcf2phylip.py -i combined_filt.vcf --output-folder out --fasta

# Run phylogenetic analysis with IQtree (the best model tested was GTR+F+I+G4+ASC)
iqtree -s combined.min4.varsites.fasta -B 1000 --seqtype DNA -m TEST -T 3
