#!/bin/bash

#SBATCH --job-name="demultiplex-16S"
#SeukTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=taruna.aggarwal@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e demultiplex-16S.err-%N
#SBATCH -o demultiplex-16S.out-%N

# before running this script, make sure you reverse complement all your reverse primers

perl \
/rhome/taruna/shared/taruna/memb/16S/scripts/Demul_trim_prep_flipedmerge2.1.pl --reverse \
--trim-file /rhome/taruna/shared/taruna/pkgs/trimmomatic-0.36/adapters/memb1-adapters-fixed.fa \
--trim-tool /rhome/taruna/shared/taruna/pkgs/trimmomatic-0.36 \
/rhome/taruna/shared/taruna/memb/16S/data-raw/uc-davis/16S/fastqs \
/rhome/taruna/shared/taruna/memb/16S/qiime-files/mapping/16S-bac-QIIME-mapping-MEmicrobiome_RC.txt \
/rhome/taruna/shared/taruna/memb/16S/data-clean/ \
16S-bac-memb1

# this script generates a qiime ready directory. Cat all the *.M.* files into a single file that should be used for OTU picking 