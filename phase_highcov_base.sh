#!/bin/sh
#SBATCH --mail-user=audrey.m.arner@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --time=16:00:00
#SBATCH --mem=50GB
#SBATCH -o phase_redo_%A_%a.out 
#SBATCH -e phase_redo_%A_%a.err  
#SBATCH --array=1-22

#############################################
## script created for R9 by Audrey Arner   ## 
#############################################

chr=`sed -n ${SLURM_ARRAY_TASK_ID}p /nobackup/lea_lab/arneram/OA_WGS_chrs_phase.txt`

# liftover to hg19
module load StdEnv/2023
module load picard/3.1.0
in_chain=/data/lea_lab/shared/hg38ToHg19.over.chain
in_genome=/data/lea_lab/shared/genomes/hg19.fa

# define output
out_reject=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_OA_SNP1.reject_chr${chr}.vcf
out_lifted_vcf=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_OA_SNP1.hg19_chr${chr}.vcf

grep -v '*' /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_WGS_allCHR.SNP1.chr${chr}.vcf.gz.vcf > /nobackup/lea_lab/arneram/OA_WGS_2Dec24/temp_${chr}.vcf

in_vcf=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/temp_${chr}.vcf

cat $in_vcf | awk '{if ($1 !~ /^#/) print "chr"$0; else print $0}' > /nobackup/lea_lab/arneram/OA_WGS_2Dec24/temp2_${chr}.vcf

java -Xmx8g -jar $EBROOTPICARD/picard.jar LiftoverVcf \
   I=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/temp2_${chr}.vcf \
   OUTPUT=$out_lifted_vcf \
   CHAIN=$in_chain \
   REJECT=$out_reject \
   R=$in_genome WARN_ON_MISSING_CONTIG=true

java -jar $EBROOTPICARD/picard.jar SortVcf \
      I=$out_lifted_vcf O=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_OA_SNP1.hg19_chr${chr}.sort.vcf

module purge
module load StdEnv/2023 gcc/12.3 bcftools/1.19

#sort vcf file
bcftools view -H /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_OA_SNP1.hg19_chr${chr}.sort.vcf | wc -l

module load StdEnv/2020
module load plink/2.00a3.6 

#grab only SNPs, some other filtering
plink2 --recode vcf bgz --vcf /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_OA_SNP1.hg19_chr${chr}.sort.vcf --allow-extra-chr --snps-only --max-alleles 2 --min-alleles 2 --maf 0.01 --keep-allele-order --out /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_OA_SNP1.hg19_chr${chr}.sort --geno 0.05 --chr ${chr}

# phase 
map=/data/lea_lab/shared/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt
phased_haps=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/OA_SNP1.hg19_chr${chr}.phased.haps
phased_sample=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/OA_SNP1.hg19_chr${chr}.phased.sample

module purge
module load StdEnv/2023 gcc/12.3 bcftools/1.19

#remove bad tid
bcftools view -s ^tid_743 -o /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_OA_SNP1.hg19_chr${chr}.sort2.vcf.gz /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_OA_SNP1.hg19_chr${chr}.sort.vcf.gz

module purge
module load StdEnv/2020
module load shapeit/2.r904

shapeit --input-vcf /nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov_OA_SNP1.hg19_chr${chr}.sort2.vcf.gz \
        --input-map $map \
        --output-max $phased_haps $phased_sample 
