#!/bin/sh
#SBATCH --mail-user=audrey.m.arner@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --time=1:30:00
#SBATCH --mem=25GB
#SBATCH -o prepare_EHH_22Oct25_%A_%a.out 
#SBATCH -e prepare_EHH_22Oct25_%A_%a.err  
#SBATCH --array=1-22

#############################################
## script created for R9 by Audrey Arner   ## 
#############################################

chr=`sed -n ${SLURM_ARRAY_TASK_ID}p /nobackup/lea_lab/arneram/OA_WGS_chrs.txt`

#gsutil cp gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_chr${chr}.filtered.SNV_INDEL.phased.shapeit5.bcf /data/lea_lab/arneram/gnomAD

module purge
module load StdEnv/2023 gcc/12.3 bcftools/1.19

#get only Han chinese samples and turn into a vcf
#bcftools view -S /data/lea_lab/arneram/han_1KG_hgdp_samplenames.txt -o /data/lea_lab/arneram/gnomAD/hgdp1kgp_chr${chr}.han.phased.shapeit5.vcf /data/lea_lab/arneram/gnomAD/hgdp1kgp_chr${chr}.filtered.SNV_INDEL.phased.shapeit5.bcf

#need to liftover to hg19
module load picard/3.1.0

in_chain=/data/lea_lab/shared/hg38ToHg19.over.chain
in_genome=/data/lea_lab/shared/genomes/hg19.fa

# define output
out_reject=/nobackup/lea_lab/arneram/rehh_files/high_cov_1KGP_HGDP.reject_chr${chr}.vcf
out_lifted_vcf=/nobackup/lea_lab/arneram/rehh_files/high_cov_1KGP_HGDP.hg19_chr${chr}.vcf

in_vcf=/data/lea_lab/arneram/gnomAD/hgdp1kgp_chr${chr}.han.phased.shapeit5.vcf

#java -Xmx8g -jar $EBROOTPICARD/picard.jar LiftoverVcf \
   #I=$in_vcf \
   #OUTPUT=$out_lifted_vcf \
   #CHAIN=$in_chain \
   #REJECT=$out_reject \
   #R=$in_genome WARN_ON_MISSING_CONTIG=true

#java -jar $EBROOTPICARD/picard.jar SortVcf \
      #I=$out_lifted_vcf O=/nobackup/lea_lab/arneram/rehh_files/high_cov_1KGP_HGDP.hg19_chr${chr}.sort.vcf

#rm $in_vcf
#rm $out_lifted_vcf

#now I want to get rid of the "chr"
#sed 's/^chr//' /nobackup/lea_lab/arneram/rehh_files/high_cov_1KGP_HGDP.hg19_chr${chr}.sort.vcf > /nobackup/lea_lab/arneram/rehh_files/high_cov_1KGP_HGDP.hg19_chr${chr}.sort2.vcf

#compress file
#bgzip /nobackup/lea_lab/arneram/rehh_files/high_cov_1KGP_HGDP.hg19_chr${chr}.sort2.vcf
#tabix -p vcf /nobackup/lea_lab/arneram/rehh_files/high_cov_1KGP_HGDP.hg19_chr${chr}.sort2.vcf.gz

#make an index for the high coverage files
tabix -p vcf /nobackup/lea_lab/arneram/OA_high_cov_final/rehh_files/high_cov.hg19_chr${chr}.phased_withAA.vcf.gz

#get only variants that are in both files
bcftools isec -n=2 -w1 \
  /nobackup/lea_lab/arneram/rehh_files/high_cov_1KGP_HGDP.hg19_chr${chr}.sort2.vcf.gz \
  /nobackup/lea_lab/arneram/OA_high_cov_final/rehh_files/high_cov.hg19_chr${chr}.phased_withAA.vcf.gz \
  -O z -o /nobackup/lea_lab/arneram/OA_high_cov_final/rehh_files/1KGP_HGDP.hg19_chr${chr}_overlap.vcf.gz


module load StdEnv/2023  gcc/12.3
module load vcftools/0.1.16
module load samtools/1.20

#add ancestral variant
gunzip -c /nobackup/lea_lab/arneram/OA_high_cov_final/rehh_files/1KGP_HGDP.hg19_chr${chr}_overlap.vcf.gz > /nobackup/lea_lab/arneram/OA_high_cov_final/rehh_files/1KGP_HGDP.hg19_chr${chr}_overlap.vcf

cat /nobackup/lea_lab/arneram/OA_high_cov_final/rehh_files/1KGP_HGDP.hg19_chr${chr}_overlap.vcf | fill-aa -a /data/lea_lab/shared/human_ancestor_GRCh37_e59/human_ancestor_${chr}.fa.gz | bgzip -c > /nobackup/lea_lab/arneram/OA_high_cov_final/rehh_files/1KGP_HGDP.hg19_chr${chr}_overlap_withAA.vcf.gz

bcftools view -H /nobackup/lea_lab/arneram/OA_high_cov_final/rehh_files/1KGP_HGDP.hg19_chr${chr}_overlap_withAA.vcf.gz | wc -l

bcftools view -i 'MAF>=0.01' \
  /nobackup/lea_lab/arneram/OA_high_cov_final/rehh_files/1KGP_HGDP.hg19_chr${chr}_overlap_withAA.vcf.gz \
  -Oz -o /nobackup/lea_lab/arneram/OA_high_cov_final/rehh_files/1KGP_HGDP.hg19_chr${chr}_overlap_withAA_MAF01.vcf.gz

bcftools view -H /nobackup/lea_lab/arneram/OA_high_cov_final/1KGP_HGDP.hg19_chr${chr}_overlap_withAA_MAF01.vcf.gz | wc -l

#any last things to remove
rm /nobackup/lea_lab/arneram/rehh_files/high_cov_1KGP_HGDP.hg19_chr${chr}.sort.vcf
rm /nobackup/lea_lab/arneram/OA_high_cov_final/rehh_files/1KGP_HGDP.hg19_chr${chr}_overlap.vcf

