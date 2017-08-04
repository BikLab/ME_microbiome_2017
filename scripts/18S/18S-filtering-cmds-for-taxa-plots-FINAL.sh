#!/bin/bash

filter_samples_from_otu_table.py -i otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.biom \
-o otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_NemsRemoved_above500.biom \
--sample_id_fp /rhome/taruna/shared/taruna/memb/18S/analysis-results/uclust_openref99_MEMB1_18s_mergedONLY_26July2017/lists-for-filtering/NemsRemoved-SampleIDs-above-500-reads.txt --negate_sample_id_fp



filter_samples_from_otu_table.py -i otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.biom \
-o otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_MetazoaRemoved-aove-500.biom \
--sample_id_fp ../lists-for-filtering/MetazoaRemoved-SampleIDs-above-500-reads.txt --negate_sample_id_fp
