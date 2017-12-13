source('http://bioconductor.org/biocLite.R')
biocLite("BiocUpgrade")
biocLite('phyloseq')
library(devtools)
library(phyloseq)
library(ggplot2)
library(biomformat)
library(wesanderson)
library(dplyr)

#import mapping file into phyloseq
map <- import_qiime_sample_data("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/18S/18S-euk-QIIME-mapping-MEmicrobiome-collapsedBYOceanRegion-15Nov17.txt")

#import biom tables into phyloseq
otu_metazoa <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/18S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_NemsRemoved_above-500-kept-collapsedBYOceanRegion.biom")


#merge mapping file and otu tables into a phyloseq object 
physeq_metazoa = merge_phyloseq(otu_metazoa,map)


#print merge phyloseq object to make sure the merging was successful
print(physeq_metazoa)


#transform the abundance values into relative abundance
physeq_metazoa_trans = transform_sample_counts(physeq_metazoa, function(x) 100 * x/sum(x))



#run some commands on the phyloseq object 
rank_names(physeq_metazoa_trans)
sample_variables(physeq_metazoa_trans)

######################### Remove unobserved OTUs - Rank2 ######################### 

any(taxa_sums(physeq_metazoa_trans) == 0)
sum(taxa_sums(physeq_metazoa_trans) == 0)
physeq_metazoa_trans_trimmed = prune_taxa(taxa_sums(physeq_metazoa_trans) > 0, physeq_metazoa_trans)
physeq_metazoa_trans_trimmed_glom <- tax_glom(physeq_metazoa_trans_trimmed, taxrank = "Rank5")
pb_meta <- plot_bar(physeq_metazoa_trans_trimmed, x = "Sample", y ="Abundance", fill = "Rank5")
pb_meta
gb_meta <- pb_meta+geom_bar(aes(color=Rank5, fill=Rank5), stat = "identity", position="stack")
gb_meta
