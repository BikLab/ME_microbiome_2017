#!/bin/bash

#SBATCH --job-name="filter-18s-OTUs"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=holly.bik@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e filter-18s-OTU.err-%N
#SBATCH -o filter-18s-OTUs.out-%N

module unload python
module load qiime

# set file paths 

BACMAP=/rhome/taruna/shared/taruna/memb/16S/qiime-files/mapping/16S-bac-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt 
BACOTUS=/rhome/taruna/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017

EUKMAP=/rhome/hbik/shared/taruna/memb/18S/qiime-files/mapping/18S-euk-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt
EUKOTUS=/rhome/hbik/shared/taruna/memb/18S/analysis-results/uclust_openref99_MEMB1_18s_mergedONLY_26July2017

BADNEMS=

############# Removing Blank SampleIDs, then filtering by Abundance (More Strict Filtering) #######################

# create a BIOM table with Blanks only
filter_samples_from_otu_table.py \
	-i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems.biom \
	-m $MAP \
	--sample_id_fp $BADNEMS \
	--negate_sample_id_fp
	
# Filter out OTU ids that have zero counts, as we only want the OTUs with positive counts remaining:
filter_otus_from_otu_table.py \
	-i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.biom \
	-n 1


	