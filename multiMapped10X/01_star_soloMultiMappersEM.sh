#!/bin/bash

#SBATCH --job-name=Star_count%j  # Job name
#SBATCH --nodes=1                                    # Run all processes on a single node	
#SBATCH --ntasks=2                                   # Run a single task		
#SBATCH --cpus-per-task=4                            # Number of CPU cores per task
#SBATCH --mem=60gb                                   # Job memory request
#SBATCH --time=05:00:00                              # Time limit hrs:min:sec
#SBATCH --partition=short 
#SBATCH --output=Star_count%j.log                
#SBATCH --array=1-9

pwd
# load starsolo and other software
module load gcc/6.2.0 star/2.7.9a samtools/1.10

# set genome and data directory path(removed temporally)
# 
genome_dir="~/genome1"
data_dir_lane12="raw_data/s12"
data_dir_lane3="raw_data/s3"
data_dir_lane45="raw_data/s45"
whitelist_file="~/genome1/3M-february-2018.txt"
date 
mkdir -p Data
cd Data

# list sample name for job array
samples={}
sampleID="${samples[${SLURM_ARRAY_TASK_ID}]}"

mkdir -p "$sampleID"
cd "$sampleID"

STAR --runThreadN 12 \
    --genomeDir $genome_dir \
    --soloType CB_UMI_Simple --soloCBwhitelist $whitelist_file --soloUMIlen 12 \
    --readFilesIn $data_dir_lane12/"$sampleID"_L001_R2_001.fastq.gz,$data_dir_lane12/"$sampleID"_L002_R2_001.fastq.gz,$data_dir_lane3/"$sampleID"_L003_R2_001.fastq.gz,$data_dir_lane45/"$sampleID"_L004_R2_001.fastq.gz,$data_dir_lane45/"$sampleID"_L005_R2_001.fastq.gz \
    $data_dir_lane12/"$sampleID"_L001_R1_001.fastq.gz,$data_dir_lane12/"$sampleID"_L002_R1_001.fastq.gz,$data_dir_lane3/"$sampleID"_L003_R1_001.fastq.gz,$data_dir_lane45/"$sampleID"_L004_R1_001.fastq.gz,$data_dir_lane45/"$sampleID"_L005_R1_001.fastq.gz \
    --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
    --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
    --soloCellFilter EmptyDrops_CR --soloFeatures Gene GeneFull SJ Velocyto \
    --twopassMode Basic --outFilterType BySJout --alignSJoverhangMin 5 --alignSJDBoverhangMin 1 --limitOutSJcollapsed 2000000\
    --outFilterMismatchNmax 999 --alignIntronMin 5 --alignIntronMax 10000 --alignMatesGapMax 100000 \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --outFilterMultimapNmax 1000 --winAnchorMultimapNmax 100 --outSAMmultNmax 1 --outMultimapperOrder Random \
    --seedSearchStartLmaxOverLread 0.5 --genomeSAindexNbases 10 --sjdbOverhang 89 --soloMultiMappers Uniform EM

sleep 5s

# comment about parameter in command above
# line 2 : specify single-cell library chemistry
# line 4 : specify format of input/output
# line 5 : Matching CellRanger 4.x.x and 5.x.x results | https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
# line 6 : specify cell filtering and feature/gene option
# line 7 : specify splice junction paramters
# line 8 : specify mapping paramterss
# line 9 : specify tag in bam
# line 10-11 : specify multi-mapping paramters

# get reads per barcode 
samtools view -@ 8 Aligned.sortedByCoord.out.bam | cut -f 21 | cut -d ':' -f 3 | sort | uniq -c | sort -nr| grep '[ATCG]\{16\}' > readCount.txt

date








