#!/bin/bash

#SBATCH --mail-user=audrey.m.arner@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --time=0:30:00
#SBATCH --mem=20GB
#SBATCH -o iHS_make_%A_%a.out
#SBATCH -e iHS_make_%A_%a.err
#SBATCH --array=1-22

#############################################
## script created for R9 by Audrey Arner   ## 
#############################################

# NOTE: chr and CHROMNUM were replaced with the chromosome number
CHROMNUM=`sed -n ${SLURM_ARRAY_TASK_ID}p /nobackup/lea_lab/arneram/OA_WGS_chrs.txt`

######
# WGS data, high coverage, phased haplotypes - convert from haplotype to VCF format
######

phased_haps=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/OA_SNP1.hg19_chr${CHROMNUM}.phased
phased_haps_vcf=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/high_cov.SNP1.hg19_chr${CHROMNUM}.phased.vcf

module load StdEnv/2020
module load shapeit/2.r904

shapeit -convert --input-haps $phased_haps \
        --output-vcf $phased_haps_vcf --thread=4

#######
# add ancestral allele info
#######
#already done
# wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2

cd /nobackup/lea_lab/arneram/rehh_files

module load StdEnv/2023  gcc/12.3
module load vcftools/0.1.16
module load samtools/1.20

sed 's,^>.*,>'$CHROMNUM',' /data/lea_lab/shared/human_ancestor_GRCh37_e59/human_ancestor_${CHROMNUM}.fa | bgzip > /data/lea_lab/shared/human_ancestor_GRCh37_e59/human_ancestor_${CHROMNUM}.fa.gz

samtools faidx /data/lea_lab/shared/human_ancestor_GRCh37_e59/human_ancestor_${CHROMNUM}.fa.gz

cat $phased_haps_vcf | fill-aa -a /data/lea_lab/shared/human_ancestor_GRCh37_e59/human_ancestor_$CHROMNUM.fa.gz | bgzip -c > /nobackup/lea_lab/arneram/rehh_files/high_cov.hg19_chr$CHROMNUM.phased_withAA.vcf.gz

