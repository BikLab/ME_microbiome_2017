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
BLANKSAMPLES=/rhome/taruna/shared/taruna/memb/16S/qiime-files/mapping/MEMB-blank-SampleIDs.txt

############## Summarize and add metadata to unfiltered BIOM TABLE ###########################

# summarize BIOM table with no blanks
biom summarize-table \
	-i $OTUS/otu_table_mc2_w_tax.biom \
	-o $OTUS/otu_table_mc2_w_tax_SUMMARY.txt

#Convert BIOM table to "classic" tab-delimited format
biom convert -i $OTUS/otu_table_mc2_w_tax.biom \
	-o $OTUS/otu_table_mc2_w_tax.txt --to-tsv \
	--header-key taxonomy

#Convert classic table to JSON BIOM 1.0 format
biom convert -i $OTUS/otu_table_mc2_w_tax.txt \
	-o $OTUS/otu_table_mc2_w_tax-JSON.biom \
	--to-json \
	--table-type="OTU table" \
	--process-obs-metadata taxonomy

#Add metadata to JSON BIOM 1.0 table
biom add-metadata -i $OTUS/otu_table_mc2_w_tax-JSON.biom \
	-o $OTUS/otu_table_mc2_w_tax-JSON-metadata.biom \
	-m $MAP


#################### Abundance Filtering ONLY (Loosest Filtering) ############################

# remove low abundance OTUs from the unfiltered BIOM table
filter_otus_from_otu_table.py -i $OTUS/otu_table_mc2_w_tax.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_af.biom \
	--min_count_fraction 0.00005

# summarize BIOM table with no blanks
biom summarize-table \
	-i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_af.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_af_SUMMARY.txt

#Convert BIOM table to "classic" tab-delimited format
biom convert -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_af.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_af.txt --to-tsv \
	--header-key taxonomy

#Convert classic table to JSON BIOM 1.0 format
biom convert -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_af.txt \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_af-JSON.biom \
	--to-json \
	--table-type="OTU table" \
	--process-obs-metadata taxonomy

#Add metadata to JSON BIOM 1.0 table
biom add-metadata -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_af-JSON.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_af-JSON-metadata.biom \
	-m $MAP
	

############# Removing Blank SampleIDs, then filtering by Abundance (More Strict Filtering) #######################

# create a BIOM table with Blanks only
filter_samples_from_otu_table.py \
	-i $OTUS/otu_table_mc2_w_tax.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved.biom \
	-m $MAP \
	--sample_id_fp $BLANKSAMPLES \
	--negate_sample_id_fp
	
# Filter out OTU ids that have zero counts, as we only want the OTUs with positive counts remaining:
filter_otus_from_otu_table.py \
	-i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved_noZeros.biom \
	-n 1

# remove low abundance OTUs from the unfiltered BIOM table
filter_otus_from_otu_table.py -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved_noZeros.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved_noZeros_af.biom \
	--min_count_fraction 0.00005

# summarize BIOM table 
biom summarize-table \
	-i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved_noZeros_af.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved_noZeros_af_SUMMARY.txt

#Convert BIOM table to "classic" tab-delimited format
biom convert -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved_noZeros_af.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved_noZeros_af.txt --to-tsv \
	--header-key taxonomy

#Convert classic table to JSON BIOM 1.0 format
biom convert -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved_noZeros_af.txt \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved_noZeros_af-JSON.biom \
	--to-json \
	--table-type="OTU table" \
	--process-obs-metadata taxonomy

#Add metadata to JSON BIOM 1.0 table
biom add-metadata -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved_noZeros_af-JSON.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankSamplesRemoved_noZeros_af-JSON-metadata.biom \
	-m $MAP


################## Removing ALL OTUs found in Blank/Control Samples (Most Stringent Filtering) ##############################

# create a BIOM table with Blanks only
filter_samples_from_otu_table.py \
	-i $OTUS/otu_table_mc2_w_tax.biom \
	-o $OTUS/blank-biom-tables/otu_table_mc2_w_tax_BlanksONLY.biom \
	-m $MAP \
	--sample_id_fp $BLANKSAMPLES

# Filter out OTU ids that have zero counts, as we only want the OTUs with positive counts from the Control_Blank samples:
filter_otus_from_otu_table.py \
	-i $OTUS/blank-biom-tables/otu_table_mc2_w_tax_BlanksONLY.biom \
	-o $OTUS/blank-biom-tables/otu_table_mc2_w_tax_BlanksONLY_noZeros.biom \
	-n 1

# Then create a tab separated version of this OTU table (use flag to retain taxonomy):
biom convert \
	-i $OTUS/blank-biom-tables/otu_table_mc2_w_tax_BlanksONLY_noZeros.biom  \
	-o $OTUS/blank-biom-tables/otu_table_mc2_w_tax_BlanksONLY_noZeros.txt  \
	--to-tsv \
	--header-key taxonomy
	
# filter out OTUs from the original BIOM table that were determined to be present in the Blank samples
filter_otus_from_otu_table.py \
	-i $OTUS/otu_table_mc2_w_tax.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved.biom \
	-e $OTUS/blank-biom-tables/otu_table_mc2_w_tax_BlanksONLY_noZeros.txt

# Blank samples now have zero sequences associated with them. Remove these to get a final OTU table:
filter_samples_from_otu_table.py \
	-i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-n 1

# summarize BIOM table with no blanks
biom summarize-table \
	-i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_SUMMARY.txt

#Convert BIOM table to "classic" tab-delimited format
biom convert -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.txt --to-tsv \
	--header-key taxonomy

#Convert classic table to JSON BIOM 1.0 format
biom convert -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.txt \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved-JSON.biom \
	--to-json \
	--table-type="OTU table" \
	--process-obs-metadata taxonomy

#Add metadata to JSON BIOM 1.0 table
biom add-metadata -i $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved-JSON.biom \
	-o $OTUS/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved-JSON-metadata.biom \
	-m $MAP





	
