library(devtools)
library(phyloseq)
library(ggplot2)
library(biomformat)
library(wesanderson)

#import mapping file into phyloseq
map <- import_qiime_sample_data("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/18S/18S-euk-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt")

#import biom tables into phyloseq
otu_metazoaRemoved <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/18S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_MetazoaRemoved-above-500.biom")
otu_nemsRemoved <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/18S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_NemsRemoved_above-500.biom")

#merge mapping file and otu tables into a phyloseq object 
physeq_metazoaRemoved = merge_phyloseq(otu_metazoaRemoved,map)
physeq_nemsRemoved = merge_phyloseq(otu_nemsRemoved, map)

#print merge phyloseq object to make sure the merging was successful
print(physeq_metazoaRemoved)
print(physeq_nemsRemoved)

#transform the abundance values into relative abundance
physeq_metazoaRemoved_trans = transform_sample_counts(physeq_metazoaRemoved, function(x) 100 * x/sum(x))
physeq_nemsRemoved_trans = transform_sample_counts(physeq_nemsRemoved, function(x) 100 * x/sum(x))


#run some commands on the phyloseq object 
rank_names(physeq_metazoaRemoved_trans)

############################# Rank2 ############################# 
#plot based on raw Rank2 (Phylum) data

pdf('~/Desktop/18S-metazoaRemoved_Rank2.pdf')
plot_bar(physeq_metazoaRemoved_trans, fill="Rank2")
dev.off()

pdf('~/Desktop/18S-nemsRemoved_Rank2.pdf')
plot_bar(physeq_nemsRemoved_trans, fill="Rank2")
dev.off()



######################### Remove unobserved OTUs - Rank2 ######################### 
any(taxa_sums(physeq_metazoaRemoved_trans) == 0)
sum(taxa_sums(physeq_metazoaRemoved_trans) == 0)

physeq_metazoaRemoved_trans_trimmed = prune_taxa(taxa_sums(physeq_metazoaRemoved_trans) > 0, physeq_metazoaRemoved_trans)
plot_bar(physeq_metazoaRemoved_trans_trimmed, x = "Sample", y ="Abundance", fill = "Rank2")


############################# Remove unobserved OTUs - top taxa - Rank2 #############################

#extract top 10 OTUs - I don't know if this is really doing top ten taxa. Whatever the legends come out like, I'd say use that number being the top X.
TopTenOTUs = names(sort(taxa_sums(physeq_metazoaRemoved_trans_trimmed), TRUE)[1:20])
head(TopTenOTUs)
physeqtop10_metazoaRemoved = prune_species(TopTenOTUs, physeq_metazoaRemoved_trans_trimmed)
plot_bar(physeqtop10_metazoaRemoved, x = "Sample", y ="Abundance", fill = "Rank3")

