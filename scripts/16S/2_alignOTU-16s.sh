#!/bin/bash

#SBATCH --job-name="alignOTU-16s-arc"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=taruna.aggarwal@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e alignOTU-16s.err-%N
#SBATCH -o alignOTU-16s.out-%N

module unload python
module load qiime


#### Step2: align sequences and remove gaps from the resulting alignment 


MAP=/rhome/taruna/shared/taruna/memb/16S/qiime-files/mapping/16S-bac-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt 
REF=/rhome/taruna/shared/taruna/dbs/greengenes_13_8/rep_set/97_otus.fasta
TAXA=/rhome/taruna/shared/taruna/dbs/greengenes_13_8/taxonomy/97_otu_taxonomy.txt
REFALIGN=/rhome/taruna/shared/taruna/dbs/greengenes_13_8/rep_set_aligned/97_otus.fasta
OTUS=/rhome/taruna/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017


# align seqs with Pynast
align_seqs.py -i $OTUS/rep_set.fna \
	-o $OTUS/pynast_aligned_seqs \
	-t $REFALIGN \
	--alignment_method pynast \
	--pairwise_alignment_method uclust \
	--min_percent_id 70.0
	
# removes gaps from the Pynast aligned seqs 
filter_alignment.py -i $OTUS/pynast_aligned_seqs/rep_set_aligned.fasta \
	-o $OTUS/pynast_aligned_seqs \
	--suppress_lane_mask_filter


# make a tree before filtering the BIOM table and the aligned rep sets fasta
make_phylogeny.py -i $OTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta \
	-o $OTUS/trees/rep_set.tre \
	--tree_method fasttree
