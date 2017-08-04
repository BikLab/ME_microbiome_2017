#!/bin/bash

#SBATCH --job-name="pickOTUopenRef-16s"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=taruna.aggarwal@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e pickOTUopenRef-16s.err-%N
#SBATCH -o pickOTUopenRef-16s.out-%N


module unload python
module load qiime

# These paths are the same in the entire script

REF=/rhome/taruna/shared/taruna/dbs/greengenes_13_8/rep_set/97_otus.fasta
TAXA=/rhome/taruna/shared/taruna/dbs/greengenes_13_8/taxonomy/97_otu_taxonomy.txt
REFALIGN=/rhome/taruna/shared/taruna/dbs/greengenes_13_8/rep_set_aligned/97_otus.fasta

# Change these paths according to files
INPUT=/rhome/taruna/shared/taruna/memb/16S/data-clean/16S-bac-memb1.fa
OUT=/rhome/taruna/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017
PARA97=/rhome/taruna/shared/taruna/memb/16S/qiime-files/parameters/16S_openref97_rdp_gg_13_8.txt

# Step1: pick otu using 97 similarity score
pick_open_reference_otus.py \
    -r $REF \
    -i $INPUT \
    -o $OUT \
    -p $PARA97 \
    -s 0.10 \
    --prefilter_percent_id 0.0 \
    --suppress_align_and_tree   
# -s subsamples at 10% of discarded sequence reads

