#!/bin/sh
#SBATCH --mail-user=audrey.m.arner@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=100GB
#SBATCH -o Trim_map_round2_%A_%a.out
#SBATCH -e Trim_map_round2_%A_%a.err
#SBATCH --array=1-214%50

#############################################
## script created by Audrey Arner 		   ## 
#############################################

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p /data/lea_lab/arneram/OA_files/OA_WGS_IDs_round2.txt`

R1=/data/lea_lab/archive_raw_fastq/OrangAsli_WGS_NovaSeqX_13Nov24/${sampleID}_R1_001.fastq.gz
R2=/data/lea_lab/archive_raw_fastq/OrangAsli_WGS_NovaSeqX_13Nov24/${sampleID}_R2_001.fastq.gz
R1_trim=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/${sampleID}.R1.trim.fastq.gz
R2_trim=/nobackup/lea_lab/arneram/OA_WGS_2Dec24/${sampleID}.R2.trim.fastq.gz

module purge
module load GCC/11.3.0
module load SAMtools/1.18

gunzip -cd $R1 | echo "$sampleID $(wc -l)" >> /nobackup/lea_lab/arneram/OA_WGS_scripts/total_reads_round2.txt

module purge
module load GCC/6.4.0-2.28
module load cutadapt/1.16-Python-3.6.3
module load BWA/0.7.17

cutadapt --nextseq-trim 20 -e 0.05 --overlap 2 -a CTGTCTCTTATACACATCT -a ATGTGTATAAGAGACA -A CTGTCTCTTATACACATCT -A ATGTGTATAAGAGACA -a "G{101}" --minimum-length=20 --trim-n -o $R1_trim -p $R2_trim $R1 $R2

echo 'trimming done'

#total mapped reads
bwa mem -t 10 /data/lea_lab/arneram/hg38.fa $R1_trim $R2_trim > /nobackup/lea_lab/arneram/OA_WGS_2Dec24/${sampleID}.hg38.sam

module purge
module load GCC/11.3.0
module load SAMtools/1.18

samtools view /nobackup/lea_lab/arneram/OA_WGS_2Dec24/${sampleID}.hg38.sam | echo "$sampleID $(wc -l)" >> /nobackup/lea_lab/arneram/OA_WGS_scripts/mapped_reads_round2.txt

echo 'mapping done'

#uniquely mapped reads
samtools view -Sbq 1 /nobackup/lea_lab/arneram/OA_WGS_2Dec24/${sampleID}.hg38.sam > /nobackup/lea_lab/arneram/OA_WGS_2Dec24/${sampleID}.hg38.bam

samtools view /nobackup/lea_lab/arneram/OA_WGS_2Dec24/${sampleID}.hg38.bam | echo "$sampleID $(wc -l)" >> /nobackup/lea_lab/arneram/OA_WGS_scripts/unique_reads_round2.txt

samtools view /nobackup/lea_lab/arneram/OA_WGS_2Dec24/${sampleID}.hg38.bam | grep 'chrY' | echo "$sampleID $(wc -l)" >> /nobackup/lea_lab/arneram/OA_WGS_scripts/chrY_reads_round2.txt

echo ${sampleID}

