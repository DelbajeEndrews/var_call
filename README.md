### Variant calling for _Aspergillus fumigatus_

The scripts in this repository were used in the manuscript titled **Increased azole persistence precedes azole resistance acquisition of the human fungal pathogen _Aspergillus fumigatus_** for variant calling.

This pipeline was applied to sequencing data from **9 newly sequenced serial isolates** and **339 publicly available whole-genome sequences (WGS)**.

Summary of the workflow (with associated scripts in parentheses):
1. **Quality control** of Illumina reads using Trimmomatic (_trimmomatic.sh_);
2. **Read mapping** to the _Af293_ reference genome with Bowtie2 (_bowtie2_map.sh_);
3. **Variant calling** with GATK4 (_GATK4.sh_):
   - Used the modules _HaplotypeCaller_, _GenotypeGVCFs_, and _VariantFiltration_ to generate individual VCF files per isolate.
   - Combined _HaplotypeCaller_ outputs with _CombineGVCFs_, followed by joint genotyping (_GenotypeGVCFs_) and filtering (_VariantFiltration_) for phylogenomic analysis.
6. **Variant annotation** of each separated VCF file with SnpEff (_snpeff.sh_);
7. **Phylogenomic analysis** based on SNPs using the combined VCF file, with tree inference performed in IQ-TREE (_snp_phylogenomic_tree.sh_).
