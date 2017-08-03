#!/bin/bash

#SBATCH --job-name="noBadNems-PCoAs"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=holly.bik@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e noBadNems-PCoAs.err-%N
#SBATCH -o noBadNems-PCoAs.out-%N

module unload python
module load qiime

# set file paths 

BACMAP=/rhome/hbik/shared/taruna/memb/16S/qiime-files/mapping/16S-bac-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt 
BACOTUS=/rhome/hbik/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017
BACBIOM=/rhome/hbik/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.biom

EUKMAP=/rhome/hbik/shared/taruna/memb/18S/qiime-files/mapping/18S-euk-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt
EUKBIOM=/rhome/hbik/shared/taruna/memb/18S/analysis-results/uclust_openref99_MEMB1_18s_mergedONLY_26July2017/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.biom
EUKOTUS=/rhome/hbik/shared/taruna/memb/18S/analysis-results/uclust_openref99_MEMB1_18s_mergedONLY_26July2017

############# Unifrac PCoAs for 16S-bac #######################

# filter fasta file of aligned rep set sequences to only keep OTUs in final abundance filtered BIOM file
filter_fasta.py -f $BACOTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta \
	-o $BACOTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.fasta \
	-b $BACBIOM
	
# make a tree before filtering the BIOM table and the aligned rep sets fasta
make_phylogeny.py -i $BACOTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.fasta \
	-o $BACOTUS/trees/rep_set_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.tre \
	--tree_method fasttree

#Unifrac Analyses at 100, 500, and 1000 seqs/sample rarefaction
beta_diversity_through_plots.py -i $BACBIOM \
	-m $BACMAP \
	-o $BACOTUS/div-analyses/NObadnems-unifrac-100 \
	-t $BACOTUS/trees/rep_set_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.tre \
	-e 100

beta_diversity_through_plots.py -i $BACBIOM \
	-m $BACMAP \
	-o $BACOTUS/div-analyses/NObadnems-unifrac-500 \
	-t $BACOTUS/trees/rep_set_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.tre \
	-e 500
	
beta_diversity_through_plots.py -i $BACBIOM \
	-m $BACMAP \
	-o $BACOTUS/div-analyses/NObadnems-unifrac-1000 \
	-t $BACOTUS/trees/rep_set_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.tre \
	-e 1000

############# Unifrac PCoAs for 18S-euk #######################

# filter fasta file of aligned rep set sequences to only keep OTUs in final abundance filtered BIOM file
filter_fasta.py -f $EUKOTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta \
	-o $EUKOTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.fasta \
	-b $EUKBIOM
	
# make a tree before filtering the BIOM table and the aligned rep sets fasta
make_phylogeny.py -i $EUKOTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.fasta \
	-o $EUKOTUS/trees/rep_set_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.tre \
	--tree_method fasttree

#Unifrac Analyses at 100, 500, and 1000 seqs/sample rarefaction
beta_diversity_through_plots.py -i $EUKBIOM \
	-m $EUKMAP \
	-o $EUKOTUS/div-analyses/NObadnems-unifrac-100 \
	-t $EUKOTUS/trees/rep_set_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.tre \
	-e 100

beta_diversity_through_plots.py -i $EUKBIOM \
	-m $EUKMAP \
	-o $EUKOTUS/div-analyses/NObadnems-unifrac-500 \
	-t $EUKOTUS/trees/rep_set_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.tre \
	-e 500
	
beta_diversity_through_plots.py -i $EUKBIOM \
	-m $EUKMAP \
	-o $EUKOTUS/div-analyses/NObadnems-unifrac-1000 \
	-t $EUKOTUS/trees/rep_set_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.tre \
	-e 1000