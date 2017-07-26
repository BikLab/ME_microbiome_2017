#!/bin/bash

#SBATCH --job-name="pickOTUopenRef-18s"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=taruna.aggarwal@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e pickOTUopenRef-18s.err-%N
#SBATCH -o pickOTUopenRef-18s.out-%N


module unload python
module load qiime

# These paths are the same in the entire script
#MAP=
REF=/rhome/taruna/shared/taruna/GOMRI/ref-dbs/SILVA_128_QIIME_release/rep_set/rep_set_18S_only/99/99_otus_18S.fasta
TAXA=/rhome/taruna/shared/taruna/GOMRI/ref-dbs/SILVA_128_QIIME_release/taxonomy/18S_only/99/consensus_taxonomy_7_levels.txt
REFALIGN=/rhome/taruna/shared/taruna/GOMRI/ref-dbs/SILVA_128_QIIME_release/rep_set_aligned/99/99_otus_aligned.fasta

# Change these paths according to files
INPUT=/rhome/taruna/shared/taruna/memb/18S/data-clean/18S-euk-memb1.fa
OUT=/rhome/taruna/shared/taruna/memb/18S/analysis-results/uclust_openref99_MEMB1_18s_25July2017
PARA99=/rhome/taruna/shared/taruna/memb/18S/qiime-files/parameters/18S_openref99_rdp_silva128.txt

# Step1: pick otu using 99 similarity score
pick_open_reference_otus.py \
    -r $REF \
    -i $INPUT \
    -o $OUT \
    -p $PARA99 \
    -s 0.10 \
    --prefilter_percent_id 0.0 \
    --suppress_align_and_tree   
# -s subsamples at 10% of discarded sequence reads
