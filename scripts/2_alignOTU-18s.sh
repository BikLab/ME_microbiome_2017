#!/bin/bash

#SBATCH --job-name="alignOTU-18s"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=taruna.aggarwal@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e alignOTU-18s.err-%N
#SBATCH -o alignOTU-18s.out-%N

module unload python
module load qiime


#### Step2: align sequences and remove gaps from the resulting alignment 

#MAP=
REF=/rhome/taruna/shared/taruna/GOMRI/ref-dbs/SILVA_128_QIIME_release/rep_set/rep_set_18S_only/99/99_otus_18S.fasta
TAXA=/rhome/taruna/shared/taruna/GOMRI/ref-dbs/SILVA_128_QIIME_release/taxonomy/18S_only/99/consensus_taxonomy_7_levels.txt
REFALIGN=/rhome/taruna/shared/taruna/GOMRI/ref-dbs/SILVA_128_QIIME_release/rep_set_aligned/99/99_otus_aligned.fasta
OTUS=/rhome/taruna/shared/taruna/memb/18S/analysis-results/uclust_openref99_MEMB1_18s_25July2017

# align seqs with Pynast
align_seqs.py -i $OTUS/rep_set.fna \
	-o $OTUS/pynast_aligned_seqs \
	-t $REFALIGN \
	--alignment_method pynast \
	--pairwise_alignment_method uclust \
	--min_percent_id 70.0