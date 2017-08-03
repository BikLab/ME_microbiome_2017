#!/bin/bash

#SBATCH --job-name="filter-16s-OTUs"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=taruna.aggarwal@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e filter-16s-OTU.err-%N
#SBATCH -o filter-16s-OTUs.out-%N

module unload python
module load qiime

# set file paths 
MAP=/rhome/taruna/shared/taruna/memb/16S/qiime-files/mapping/16S-bac-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt 
REF=/rhome/taruna/shared/taruna/dbs/greengenes_13_8/rep_set/97_otus.fasta
TAXA=/rhome/taruna/shared/taruna/dbs/greengenes_13_8/taxonomy/97_otu_taxonomy.txt
REFALIGN=/rhome/taruna/shared/taruna/dbs/greengenes_13_8/rep_set_aligned/97_otus.fasta
OTUS=/rhome/taruna/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017


# removes gaps from the Pynast aligned seqs 
filter_alignment.py -i $OTUS/pynast_aligned_seqs/rep_set_aligned.fasta \
	-o $OTUS/pynast_aligned_seqs \
	--suppress_lane_mask_filter

# make a tree before filtering the BIOM table and the aligned rep sets fasta
make_phylogeny.py -i $OTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta \
	-o $OTUS/trees/rep_set.tre \
	--tree_method fasttree

# create a BIOM table with Blanks only
filter_samples_from_otu_table.py \
	-i $OTUS/otu_table_mc2_w_tax.biom \
	-o $OTUS/blank-biom-tables/otu_table_mc2_w_tax_BlanksONLY.biom \
	-m $MAP \
	-s "Habitat:Blank"

# filter out OTUs with zero results bc we only want the OTUs with positive counts from the Blanks
filter_otus_from_otu_table.py \
	-i $OTUS/blank-biom-tables/otu_table_mc2_w_tax_BlanksONLY.biom \
	-o $OTUS/blank-biom-tables/otu_table_mc2_w_tax_BlanksONLY_noZero.biom \
	-n 1

# create a tab separated version of the filtered BIOM table
biom convert \
	-i $OTUS/blank-biom-tables/otu_table_mc2_w_tax_BlanksONLY_noZero.biom \
	-o $OTUS/blank-biom-tables/otu_table_mc2_w_tax_BlanksONLY_noZero.txt \
	--to-tsv

# filter out OTUs from the original BIOM table that were determined to be present in the Blank samples
filter_otus_from_otu_table.py \
	-i $OTUS/otu_table_mc2_w_tax.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_noBlanks.biom \
	-e $OTUS/blank-biom-tables/otu_table_mc2_w_tax_BlanksONLY_noZero.txt

# summarize BIOM table with no blanks
biom summarize-table \
	-i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_noBlanks.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_noBlanks.biom.summary


# remove pynast failures from BIOM table with tax
filter_otus_from_otu_table.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_noBlanks.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_noBlanks_NOpynastfail.biom \
	-e $OTUS/pynast_aligned_seqs/rep_set_failures.fasta

# remove low abundance OTUs from the pynast failures filtered BIOM table
filter_otus_from_otu_table.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_noBlanks_NOpynastfail.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_noBlanks_NOpynastfail_NOlowabundance0.005.biom \
	--min_count_fraction 0.005
	
# summarize the filtered BIOM table
biom summarize-table -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_noBlanks_NOpynastfail_NOlowabundance0.005.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_noBlanks_NOpynastfail_NOlowabundance0.005.biom.summary
	
# filter fasta file of aligned rep set sequences to only keep OTUs in final abundance filtered BIOM file
filter_fasta.py -f $OTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta \
	-o $OTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered_noBlanks_NOpynastfail_NOlowabundance0.005.fasta \
	-b $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_noBlanks_NOpynastfail_NOlowabundance0.005.biom
	
# make a tree before filtering the BIOM table and the aligned rep sets fasta
make_phylogeny.py -i $OTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta \
	-o $OTUS/trees/rep_set.tre \
	--tree_method fasttree
