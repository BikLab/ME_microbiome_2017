#!/bin/bash

#SBATCH --job-name="tax-summ-18S-16S"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=holly.bik@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e tax-summ-18S-16S.err-%N
#SBATCH -o tax-summ-18S-16S.out-%N

module unload python
module load qiime

# set file paths 
16SOTUS=/rhome/hbik/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017
16STAXSUM=/rhome/hbik/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017/tax-summary

18SOTUS=/rhome/hbik/shared/taruna/memb/18S/analysis-results/uclust_openref99_MEMB1_18s_mergedONLY_26July2017
18STAXSUM=/rhome/hbik/shared/taruna/memb/18S/analysis-results/uclust_openref99_MEMB1_18s_mergedONLY_26July2017/tax-summary

#Summarize taxonomy levels and make bar chart plots for 16S-bac data
summarize_taxa.py -i $16SOTUS/otu_table_mc2_w_tax.biom \
	-o $16STAXSUM \
	-L 3,4,5,6,7

plot_taxa_summary.py -i $16STAXSUM/otu_table_mc2_w_tax_L3.txt \
	-c bar \
	-o $16STAXSUM/tax_plots/

plot_taxa_summary.py -i $16STAXSUM/otu_table_mc2_w_tax_L4.txt \
	-c bar \
	-o $16STAXSUM/tax_plots/

plot_taxa_summary.py -i $16STAXSUM/otu_table_mc2_w_tax_L5.txt \
	-c bar \
	-o $16STAXSUM/tax_plots/

plot_taxa_summary.py -i $16STAXSUM/otu_table_mc2_w_tax_L6.txt \
	-c bar \
	-o $16STAXSUM/tax_plots/
	
plot_taxa_summary.py -i $16STAXSUM/otu_table_mc2_w_tax_L7.txt \
	-c bar \
	-o $16STAXSUM/tax_plots/

#Summarize taxonomy levels and make bar chart plots for 18S-euk data
summarize_taxa.py -i $18SOTUS/otu_table_mc2_w_tax.biom \
	-o $18STAXSUM \
	-L 3,4,5,6,7

plot_taxa_summary.py -i $18STAXSUM/otu_table_mc2_w_tax_L3.txt \
	-c bar \
	-o $18STAXSUM/tax_plots/

plot_taxa_summary.py -i $18STAXSUM/otu_table_mc2_w_tax_L4.txt \
	-c bar \
	-o $18STAXSUM/tax_plots/

plot_taxa_summary.py -i $18STAXSUM/otu_table_mc2_w_tax_L5.txt \
	-c bar \
	-o $18STAXSUM/tax_plots/

plot_taxa_summary.py -i $18STAXSUM/otu_table_mc2_w_tax_L6.txt \
	-c bar \
	-o $18STAXSUM/tax_plots/
	
plot_taxa_summary.py -i $18STAXSUM/otu_table_mc2_w_tax_L7.txt \
	-c bar \
	-o $18STAXSUM/tax_plots/