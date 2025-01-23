#!/bin/sh
#SBATCH --mail-user=audrey.m.arner@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=125GB
#SBATCH -o Bam_processing_round2_%A_%a.out
#SBATCH -e Bam_processing_round2_%A_%a.err
#SBATCH --array=1-214%50

#############################################
## script created by Audrey Arner 		   ## 
#############################################

# define input
barcode=`sed -n ${SLURM_ARRAY_TASK_ID}p /data/lea_lab/arneram/OA_files/OA_WGS_IDs_round2.txt`
echo ${barcode}

# define output
out_bam1=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/$barcode.hg38.bam

out_bam_sorted1=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/$barcode.hg38.sorted1.bam

out_bam_duplicates=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/$barcode.hg38.dup.bam
out_txt_dupmetrics=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/$barcode.hg38.dup.txt
out_bam_readgroups=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/$barcode.hg38.reg.bam
out_table_recalibr=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/$barcode.rbqs.table
out_bam_recalibrat=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/$barcode.rbqs.bam
out_g_vcf=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/$barcode.hg38.g.vcf

in_genome=/data/lea_lab/arneram/hg38.fa
in_variants=/data/lea_lab/arneram/dbsnp/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf

in_gatk=/nobackup/lea_lab/petersrm/mSTARR/gatk-4.1.4.0

# Step   I: Identify duplicated reads with sorted bam or sam files
# Step  II: Recalibrate base quality scores
# Step III: Calling GATK variants and export VCF

# Step I
module purge
module load GCC/11.3.0
module load SAMtools/1.18

echo 'samtools sorting...' # sort bam file

samtools sort -m 3G -o $out_bam_sorted1 $out_bam1

module load picard/2.18.27

echo 'picard duplications...' # run picard to mark duplicates

java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$out_bam_sorted1 O=$out_bam_duplicates M=$out_txt_dupmetrics

# Step II
echo 'adding read groups...'

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$out_bam_duplicates O=$out_bam_readgroups SO=coordinate RGLB=$barcode RGPL=illumina RGPU=orangasli RGSM=$barcode

echo 'pre-indexing...'

samtools index $out_bam_readgroups --threads 4

echo 'recalibrating base quality scores...' # sort bam file

module purge
module load GCC/8.2.0 SAMtools/1.9 picard/2.18.27

$in_gatk/gatk --java-options "-Xmx80g -Xms80g" BaseRecalibrator  -I $out_bam_readgroups -R $in_genome --known-sites $in_variants -O $out_table_recalibr

$in_gatk/gatk --java-options "-Xmx80g -Xms80g" ApplyBQSR -R $in_genome -I $out_bam_readgroups -bqsr $out_table_recalibr -O $out_bam_recalibrat

echo 'post-indexing...'

module purge
module load GCC/11.3.0
module load SAMtools/1.18

samtools index $out_bam_recalibrat --threads 4

# Step III
echo 'calling GATK variants as VCF and zip...' # sort bam file

module purge
module load GCC/8.2.0 SAMtools/1.9 picard/2.18.27

$in_gatk/gatk --java-options "-Xmx80g -Xms80g" HaplotypeCaller -R $in_genome -I $out_bam_recalibrat -O $out_g_vcf -ERC GVCF --max-alternate-alleles 2

module purge
module load GCC/5.4.0-2.26 tabix/0.2.6

bgzip /nobackup/lea_lab/arneram/OA_WGS_2Dec24/VANP_TID0071_1_PB_WBC_C1_IDPFT_A14778_22GYWWLT4_ATGTTGTTGG_L007.hg38.g.vcf
tabix -f -p vcf /nobackup/lea_lab/arneram/OA_WGS_2Dec24/VANP_TID0287_1_PB_WBC_C1_IDPFT_A14766_22GYWWLT4_GGTAGAATTA_L007.hg38.g.vcf.gz

zcat /nobackup/lea_lab/arneram/OA_WGS_2Dec24/${barcode}.hg38.g.vcf.gz | echo "$barcode $(wc -l)" >> /nobackup/lea_lab/arneram/OA_WGS_2Dec24/variants_gcf_round2.txt

