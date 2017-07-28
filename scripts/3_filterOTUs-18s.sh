#!/bin/bash

#SBATCH --job-name="filter-18s-OTUs"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=taruna.aggarwal@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e filter-18s-OTU.err-%N
#SBATCH -o filter-18s-OTUs.out-%N

module unload python
module load qiime

MAP=/rhome/taruna/shared/taruna/memb/18S/qiime-files/mapping/18S-euk-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt
REF=/rhome/taruna/shared/taruna/GOMRI/ref-dbs/SILVA_128_QIIME_release/rep_set/rep_set_18S_only/99/99_otus_18S.fasta
TAXA=/rhome/taruna/shared/taruna/GOMRI/ref-dbs/SILVA_128_QIIME_release/taxonomy/18S_only/99/consensus_taxonomy_7_levels.txt
REFALIGN=/rhome/taruna/shared/taruna/GOMRI/ref-dbs/SILVA_128_QIIME_release/rep_set_aligned/99/99_otus_aligned.fasta
OTUS=/rhome/taruna/shared/taruna/memb/18S/analysis-results/uclust_openref99_MEMB1_18s_mergedONLY_26July2017
BLANKSAMPLES=/rhome/taruna/shared/taruna/memb/16S/qiime-files/mapping/MEMB-blank-SampleIDs.txt

# removes gaps from the Pynast aligned seqs 
filter_alignment.py -i $OTUS/pynast_aligned_seqs/rep_set_aligned.fasta \
	-o $OTUS/pynast_aligned_seqs \
	--suppress_lane_mask_filter
	
# filter BIOM table for D_3__Metazoa (Animalia) only
filter_taxa_from_otu_table.py \
	-i $OTUS/otu_table_mc2_w_tax.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_metazoaONLY.biom \
	-p "D_3__Metazoa (Animalia)"


# filter the otu_table_mc2_w_tax_metazoaONLY.biom to remove samples that have no OTU counts
filter_samples_from_otu_table.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_metazoaONLY.biom -o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_metazoaONLY_allNemaSamples.biom -n 1


filter_otus_from_otu_table.py -i $OTUS/otu_table_mc2_w_tax.biom -o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_NOmetazoa.biom -e $OTUS/filtered-biom-tables/list-of-otus-MetazoaONLY.txt

