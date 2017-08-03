#!/bin/bash

#SBATCH --job-name="betadiv-unfilt-16S"
#SBACTH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=holly.bik@ucr.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e betadiv-16s-OTU.err-%N
#SBATCH -o betadiv-16s-OTUs.out-%N

module unload python
module load qiime

# set file paths 
MAP=/rhome/hbik/shared/taruna/memb/16S/qiime-files/mapping/16S-bac-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt 
REF=/rhome/hbik/shared/taruna/dbs/greengenes_13_8/rep_set/97_otus.fasta
TAXA=/rhome/hbik/shared/taruna/dbs/greengenes_13_8/taxonomy/97_otu_taxonomy.txt
REFALIGN=/rhome/hbik/shared/taruna/dbs/greengenes_13_8/rep_set_aligned/97_otus.fasta

OTUS=/rhome/hbik/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017
DIV=/rhome/hbik/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017/div-analyses

TREE=/rhome/hbik/shared/taruna/memb/16S/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017/trees/rep_set_unfiltered.tre

# Change rarefaction values below (-e flag) based on sample with lowest number of reads in the BIOM table summary
# Values listed below were chosen based on the BIOM table summaries generated in the previous script

beta_diversity_through_plots.py -i $OTUS/otu_table_mc2_w_tax.biom \
	-m $MAP \
	-o $DIV/unfilt-beta-div-1000 \
	-t $TREE \
	-e 1000

alpha_rarefaction.py -i $OTUS/otu_table_mc2_w_tax.biom \
	-m $MAP \
	-o $DIV/unfilt-alpha-rare-1000 \
	-t $TREE \
	-e 1000	
	
beta_diversity_through_plots.py -i $OTUS/otu_table_mc2_w_tax.biom \
	-m $MAP \
	-o $DIV/unfilt-beta-div-10000 \
	-t $TREE \
	-e 10000

alpha_rarefaction.py -i $OTUS/otu_table_mc2_w_tax.biom \
	-m $MAP \
	-o $DIV/unfilt-alpha-rare-10000 \
	-t $TREE \
	-e 10000	

exit 0;
