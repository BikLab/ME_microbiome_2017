#!/bin/bash

#SBATCH --job-name="demultiplex-w-QIIME"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=taruna.aggarwal@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e demultiplex-w-QIIME.err-%N
#SBATCH -o demultiplex-w-QIIME.out-%N


module unload python
module load qiime

# Step 0 - set file paths

RAW=/rhome/taruna/shared/taruna/memb/data-raw/uc-davis/16S/fastqs
COMBINED_BARCODES=/rhome/taruna/shared/taruna/memb/data-raw/uc-davis/16S/combined-barcodes
COMBINED_FASTQS=/rhome/taruna/shared/taruna/memb/data-raw/uc-davis/16S/combined-fastqs
MAP=/rhome/taruna/shared/taruna/memb/qiime-files/mapping
SPLIT=/rhome/taruna/shared/taruna/memb/data-raw/uc-davis/16S/split-libraries-output

# Step 1 - extract barcodes from each index file (there should be two) and combine them into a single file using QIIME's extract_barcodes.py 
extract_barcodes.py \
--fastq1 $RAW/16S-bac_memb1_S0_L001_I2_001_ForBC.fastq.gz \
--fastq2 $RAW/16S-bac_memb1_S0_L001_I1_001_RevBC.fastq.gz \
--bc1_len 15 \
--bc2_len 15 \
--input_type barcode_paired_end \
-o $COMBINED_BARCODES


# Step 2 - Modify the mapping file so that it matches the fastQ generated in Step 1. This can be done in excel or using a custom script that is yet to be written

# Step 3 - Remove the trailing " #:N:0" from R1, R2 and the combined barcodes fastQ. For example:

zcat $RAW/16S-bac_memb1_S0_L001_R1_001.fastq.gz | sed 's/ 1:N:0:0//g' > $RAW/16S-bac_memb1_S0_L001_R1_001_filtered.fastq
zcat $RAW/16S-bac_memb1_S0_L001_R2_001.fastq.gz | sed 's/ 2:N:0:0//g' > $RAW/16S-bac_memb1_S0_L001_R2_001_filtered.fastq
sed 's/ 2:N:0:0//g' $COMBINED_BARCODES/16S-combined-barcodes.fastq > $COMBINED_BARCODES/16S-combined-barcodes-filtered.fastq

# Step 4 - Join the R1 and R2 fastQ files using join_paired_ends.py with -b option

join_paired_ends.py \
	-f $RAW/16S-bac_memb1_S0_L001_R1_001_filtered.fastq \
	-r $RAW/16S-bac_memb1_S0_L001_R2_001_filtered.fastq \
	-b $COMBINED_BARCODES/16S-combined-barcodes-filtered.fastq \
	-o $COMBINED_FASTQS/ \
	-j 10 -p 15


# Step 4 - Demultiplex the files using split_libraries_fastq.py

split_libraries_fastq.py -i $COMBINED_FASTQS/fastqjoin.join.fastq -b $COMBINED_FASTQS/fastqjoin.join_barcodes.fastq -m $MAP/16S-bac-QIIME-mapping-MEmicrobiome-combined-RC.txt -o $SPLIT --rev_comp_mapping_barcodes -q 19 -r 5 -p 0.70 --barcode_type 30
