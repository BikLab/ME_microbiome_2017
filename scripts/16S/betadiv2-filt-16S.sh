#!/bin/bash

#SBATCH --job-name="betadiv2-filt-16s-OTUs"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=holly.bik@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e betadiv2-filt-16s-OTUs.err-%N
#SBATCH -o betadiv2-filt-16s-OTUs.out-%N

module unload python
module load qiime

# set file paths 
MAP=/rhome/hbik/shared/taruna/memb/16S/qiime-files/mapping/16S-bac-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt 
REF=/rhome/hbik/shared/taruna/dbs/greengenes_13_8/rep_set/97_otus.fasta
TAXA=/rhome/hbik/shared/taruna/dbs/greengenes_13_8/taxonomy/97_otu_taxonomy.txt
REFALIGN=/rhome/hbik/shared/taruna/dbs/greengenes_13_8/rep_set_aligned/97_otus.fasta

OTUS=/rhome/hbik/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017
PARAMS=/rhome/hbik/shared/taruna/memb/scripts-holly
DIV=/rhome/hbik/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017/div-analyses
	
# filter fasta file of aligned rep set sequences to only keep OTUs in final abundance filtered BIOM file
filter_fasta.py -f $OTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta \
	-o $OTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered_BlankOTUsRemoved_BlankSamplesRemoved.fasta \
	-b $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom
	
# make a tree before filtering the BIOM table and the aligned rep sets fasta
make_phylogeny.py -i $OTUS/pynast_aligned_seqs/rep_set_aligned_pfiltered_BlankOTUsRemoved_BlankSamplesRemoved.fasta \
	-o $OTUS/trees/rep_set_BlankOTUsRemoved_BlankSamplesRemoved.tre \
	--tree_method fasttree

#Unifrac Analyses at 500, 1000, and 10,000 seqs/sample rarefaction
beta_diversity_through_plots.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-m $MAP \
	-o $DIV/BlankOTUsRemoved-unifrac-500 \
	-t $OTUS/trees/rep_set_BlankOTUsRemoved_BlankSamplesRemoved.tre \
	-e 500

beta_diversity_through_plots.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-m $MAP \
	-o $DIV/BlankOTUsRemoved-unifrac-1000 \
	-t $OTUS/trees/rep_set_BlankOTUsRemoved_BlankSamplesRemoved.tre \
	-e 1000
	
beta_diversity_through_plots.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-m $MAP \
	-o $DIV/BlankOTUsRemoved-unifrac-10000 \
	-t $OTUS/trees/rep_set_BlankOTUsRemoved_BlankSamplesRemoved.tre \
	-e 10000

#Bray curtis, Canberra, and Jaccard analyses at 500, 1000, and 10,000 seqs/sample rarefaction
beta_diversity_through_plots.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-m $MAP \
	-o $DIV/BlankOTUsRemoved-bc-1000 \
	-p $PARAMS/params_braycurtis.txt \
	-e 1000
	
beta_diversity_through_plots.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-m $MAP \
	-o $DIV/BlankOTUsRemoved-canberra-1000 \
	-p $PARAMS/params_canberra.txt \
	-e 1000

beta_diversity_through_plots.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-m $MAP \
	-o $DIV/BlankOTUsRemoved-jaccard-1000 \
	-p $PARAMS/params_jaccard.txt \
	-e 1000

beta_diversity_through_plots.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-m $MAP \
	-o $DIV/BlankOTUsRemoved-bc-10000 \
	-p $PARAMS/params_braycurtis.txt \
	-e 10000
	
beta_diversity_through_plots.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-m $MAP \
	-o $DIV/BlankOTUsRemoved-canberra-10000 \
	-p $PARAMS/params_canberra.txt \
	-e 10000

beta_diversity_through_plots.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-m $MAP \
	-o $DIV/BlankOTUsRemoved-jaccard-10000 \
	-p $PARAMS/params_jaccard.txt \
	-e 10000
