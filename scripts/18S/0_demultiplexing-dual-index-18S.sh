#!/bin/bash

#SBATCH --job-name="demultiplex-18S"
#SeukTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=taruna.aggarwal@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e demultiplex-18S.err-%N
#SBATCH -o demultiplex-18S.out-%N

# before running this script, make sure you reverse complement all your reverse primers

perl \
/rhome/taruna/shared/taruna/memb/18S/scripts/Demul_trim_prep_flipedmerge2.1.pl --reverse \
--trim-file /rhome/taruna/shared/taruna/pkgs/trimmomatic-0.36/adapters/memb1-adapters-fixed.fa \
--trim-tool /rhome/taruna/shared/taruna/pkgs/trimmomatic-0.36 \
/rhome/taruna/shared/taruna/memb/18S/data-raw/uc-davis/18S/fastqs \
/rhome/taruna/shared/taruna/memb/18S/qiime-files/mapping/18S-euk-QIIME-mapping-MEmicrobiome_RC.txt \
/rhome/taruna/shared/taruna/memb/18S/data-clean/ \
18S-euk-memb1

# this script generates a qiime ready directory. Cat all the *.M.* files into a single file that should be used for OTU picking
 