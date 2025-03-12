#!/bin/sh -l
#SBATCH --mail-user=audrey.m.arner@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=100GB
#SBATCH -o filter_hc_WGS_24Feb25_try2.out
#SBATCH -e filter_hc_WGS_24Feb25_try2.err

#############################################
## script created for R9 by Audrey Arner   ## 
#############################################

module purge
module load StdEnv/2023 gcc/12.3 bcftools/1.19

bcftools concat -a -Oz -o /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_27Jan25_OA_allCHR.vcf.gz /nobackup/lea_lab/arneram/OA_high_cov_Jan25/*.all_high_cov.vcf.gz

module purge
module load StdEnv/2020
module load plink/1.9b_6.21-x86_64 

# check for problematic samples
plink --vcf /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_27Jan25_OA_allCHR.vcf.gz --snps-only --missing --out /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_missing
plink --vcf /nobackup/lea_lab/arneram/OA_WGS_2Dec24/low_cov_27Jan25_OA_allCHR.vcf.gz --snps-only --missing --out /nobackup/lea_lab/arneram/OA_WGS_2Dec24/low_missing

##########
# USE HARD FILTERING FIRST TO REMOVE STUFF WE DEFINITELY WON'T ANALYZE
##########

# region filtering
cpgs=/data/lea_lab/shared/cpgs/CpG_hg38_tab.bed
mask=/data/lea_lab/shared/human_masked/20160622.allChr.pilot_mask.bed
super=/data/lea_lab/shared/human_superdups/hg38_superdups.bed

#already done
#cat $cpgs $super > /data/lea_lab/shared/temp.bed #already done in previous WGS round

#awk '$4=(FNR FS)' /data/lea_lab/shared/temp.bed > /data/lea_lab/shared/temp_v2.bed #done in previous WGS round

#awk '$1 !~ /_/' /data/lea_lab/shared/temp_v2.bed > /data/lea_lab/shared/temp_v3.bed #done in previous WGS round

# GATK reccomends filtering for excess hets before VQSR, HWE filtering should do something similar: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
plink --vcf /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_27Jan25_OA_allCHR.vcf.gz --snps-only --hwe 0.000001 --make-just-bim --extract range $mask --maf 0.01 --exclude range /data/lea_lab/shared/temp_v3.bed --allow-extra-chr --out /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs

#this will give any with missingness over .1
awk '$6 > 0.10 {print $1, $2}' /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_missing.imiss > /nobackup/lea_lab/arneram/OA_WGS_2Dec24/many_missing.txt 

# turn into bed
awk '{OFS=""; print "chr",$1,"\t",$4-1,"\t",$4}' /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs.bim >  /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs_loc1.bed

awk '{OFS=""; print "chr",$1,"\t",$4,"\t",$4+1}' /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs.bim >  /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs_loc2.bed

module purge
module load StdEnv/2023 tabix/0.2.6

tabix -p vcf /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_27Jan25_OA_allCHR.vcf.gz

# remove low MAF and problem regions, remove 2 problem samples
module purge
module load StdEnv/2023 gcc/12.3 bcftools/1.19

#run above first. Look at samples to see if any are terrible
in_gatk=/nobackup/lea_lab/petersrm/mSTARR/gatk-4.1.4.0

regions=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs_loc1.bed
bcftools view -Oz -o /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs1.vcf.gz --samples ^tid_1305,tid_743 --regions-file $regions /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_27Jan25_OA_allCHR.vcf.gz

regions=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs_loc2.bed
bcftools view -Oz -o /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs2.vcf.gz --samples ^tid_1305,tid_743 --regions-file $regions /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_27Jan25_OA_allCHR.vcf.gz

# filter and sort
bcftools sort -Oz -o /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs1._sort.vcf.gz /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs1.vcf.gz

module purge
module load StdEnv/2023 picard/3.1.0 gcc/12.3 samtools/1.20

$in_gatk/gatk IndexFeatureFile -F /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs1._sort.vcf.gz

##########
# USE GATK MIXTURE MODEL FILTERING
##########

module load r/4.4.0 
echo 'calculating VQSLOD tranches for SNPs...'

in_vcf_gz=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_24Feb25_WGS_allCHR.SNPs1._sort.vcf.gz

in_genome=/data/lea_lab/arneram/hg38.fa
in_resourceDIR=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets

out_r_plots=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_WGS_allCHR.SNP1.plots.R
out_snps_recal=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_WGS_allCHR.SNP1.recal
out_snps_tranc=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_WGS_allCHR.SNP1.tranches
out_vcf_gz=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_WGS_allCHR.SNP1.vqsr.vcf.gz

$in_gatk/gatk --java-options "-Xmx60g -Xms60g" VariantRecalibrator \
            -V $in_vcf_gz \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 95.0 -tranche 90.0 \
            -an QD -an MQ -an ReadPosRankSum -an FS -an SOR -an DP \
            -mode SNP \
            --max-gaussians 6 \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /nobackup/lea_lab/arneram/resources/hapmap_3.3.hg38.vcf.gz \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 /nobackup/lea_lab/arneram/resources/1000G_omni2.5.hg38.vcf.gz \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 /nobackup/lea_lab/arneram/resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /data/lea_lab/arneram/dbsnp/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
            -O $out_snps_recal \
            --tranches-file $out_snps_tranc \
            -R $in_genome \
            --rscript-file $out_r_plots

echo 'filtering SNPs with VQSLOD...'

$in_gatk/gatk --java-options "-Xmx60g -Xms60g" ApplyVQSR \
            -V $in_vcf_gz \
            --recal-file $out_snps_recal \
            --tranches-file $out_snps_tranc \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            --exclude-filtered true \
            -mode SNP \
            -O $out_vcf_gz \
            -R $in_genome

module purge
module load StdEnv/2020
module load plink/2.00a3.6

# create one VCF per chromosome
for chr in {1..22}; do plink2 --vcf /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_WGS_allCHR.SNP1.vqsr.vcf.gz --var-filter --recode vcf --out /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_WGS_allCHR.SNP1.chr${chr}.vcf.gz --chr ${chr}; echo ${chr}; done

echo 'donezo!'

