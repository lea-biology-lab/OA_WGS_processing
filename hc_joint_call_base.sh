#!/bin/sh -l
#SBATCH --mail-user=audrey.m.arner@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=100GB
#SBATCH -o high_joint_call_22Feb25_%A_%a.out
#SBATCH -e high_joint_call_22Feb25_%A_%a.err

#############################################
## script created for R9 by Audrey Arner   ## 
#############################################

INTERVAL=`sed -n ${SLURM_ARRAY_TASK_ID}p /data/lea_lab/arneram/50MB_chr_redo_22Feb25`

rm -r /nobackup/lea_lab/arneram/OA_high_cov_Jan25/$INTERVAL &>/dev/null 
rm /nobackup/lea_lab/arneram/OA_high_cov_Jan25/${INTERVAL}* &>/dev/null 

chrom=$INTERVAL
in_genome=/data/lea_lab/arneram/hg38.fa
in_database=/nobackup/lea_lab/arneram/OA_high_cov_Jan25/$chrom
in_samples=/data/lea_lab/arneram/OA_files/high_cov.sample_map.txt
out_vcf_gz=/nobackup/lea_lab/arneram/OA_high_cov_Jan25/$chrom.all_high_cov.vcf.gz
tmp_dir=/nobackup/lea_lab/arneram/OA_high_cov_Jan25/tmp_files

in_gatk=/nobackup/lea_lab/petersrm/mSTARR/gatk-4.1.4.0

module purge
module load StdEnv/2023 picard/3.1.0 gcc/12.3 samtools/1.20

#merge GVCFs
$in_gatk/gatk --java-options "-Xmx80g -Xms80g" GenomicsDBImport --sample-name-map $in_samples --genomicsdb-workspace-path $in_database --tmp-dir=$tmp_dir -L $chrom --batch-size 15

# call genotypes

echo 'calling joint genotypes...'

echo $chrom > temp.${chrom}.list

$in_gatk/gatk --java-options "-Xmx80g" GenotypeGVCFs \
            -R $in_genome \
            -L temp.${chrom}.list \
            -V gendb://$in_database \
            -O $out_vcf_gz \
            --tmp-dir=$tmp_dir --max-alternate-alleles 2

echo 'donezo!'
