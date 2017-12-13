source('http://bioconductor.org/biocLite.R')
biocLite("BiocUpgrade")
biocLite('phyloseq')
library(devtools)
library(phyloseq)
library(ggplot2)
library(biomformat)
library(wesanderson)
library(ape)

#import mapping file into phyloseq
map <- import_qiime_sample_data("~/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/18S/18S-euk-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt")

#import biom tables into phyloseq
otu_metazoaRemoved <- import_biom("~/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/18S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_NemsRemoved_above-500-kept-B1200A.biom")
otu_nemsRemoved <- import_biom("~/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/18S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_NemsRemoved_above-500-kept.biom")

#import tree files
tree_metazoa <- read.tree("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/18S/rep_set_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros.tre")
tree_metazoa

#merge mapping file and otu tables into a phyloseq object 
physeq_metazoaRemoved = merge_phyloseq(otu_metazoaRemoved,map,tree_metazoa)
physeq_nemsRemoved = merge_phyloseq(otu_nemsRemoved, map)

#print merge phyloseq object to make sure the merging was successful
print(physeq_metazoaRemoved)
print(physeq_nemsRemoved)

#transform the abundance values into relative abundance
physeq_metazoaRemoved_trans = transform_sample_counts(physeq_metazoaRemoved, function(x) 100 * x/sum(x))
physeq_nemsRemoved_trans = transform_sample_counts(physeq_nemsRemoved, function(x) 100 * x/sum(x))


#run some commands on the phyloseq object 
rank_names(physeq_metazoaRemoved_trans)
sample_variables(physeq_metazoaRemoved_trans)

############################# Rank2 ############################# 
#plot based on raw Rank2 (Phylum) data

pdf('~/Desktop/18S-metazoaRemoved_Rank2.pdf')
plot_bar(physeq_metazoaRemoved_trans, fill="Rank3")
dev.off()

pdf('~/Desktop/18S-nemsRemoved_Rank2.pdf')
plot_bar(physeq_nemsRemoved_trans, fill="Rank3")
dev.off()



######################### Remove unobserved OTUs - Rank2 ######################### 
any(taxa_sums(physeq_metazoaRemoved_trans) == 0)
sum(taxa_sums(physeq_metazoaRemoved_trans) == 0)
physeq_metazoaRemoved_trans_trimmed = prune_taxa(taxa_sums(physeq_metazoaRemoved_trans) > 0, physeq_metazoaRemoved_trans)
physeq_metazoaRemoved_trans_trimmed_glom <- tax_glom(physeq_metazoaRemoved_trans_trimmed, taxrank = "Rank3")
pb_meta <- plot_bar(physeq_metazoaRemoved_trans_trimmed_glom, x = "Sample", y ="Abundance", fill = "Rank3") #facet_grid=~OceanRegion)
gb_meta <- pb_meta+geom_bar(aes(color=Rank3, fill=Rank3), stat = "identity", position="stack")
gb_meta


meta_ord <- ordinate(physeq_metazoaRemoved_trans_trimmed, "NMDS", "bray")
meta_ord_plot <- plot_ordination(physeq_metazoaRemoved_trans_trimmed, meta_ord, color="OceanRegion", shape="FeedingGroup" )

any(taxa_sums(physeq_nemsRemoved_trans) == 0)
sum(taxa_sums(physeq_nemsRemoved_trans) == 0)
physeq_nemsRemoved_trans_trimmed = prune_taxa(taxa_sums(physeq_nemsRemoved_trans) > 0, physeq_nemsRemoved_trans)
physeq_nemsRemoved_trans_trimmed_glom <- tax_glom(physeq_nemsRemoved_trans_trimmed, taxrank = "Rank5")
pb_nems <- plot_bar(physeq_nemsRemoved_trans_trimmed_glom, x = "Sample", y ="Abundance", fill = "Rank5")
gb_nems <- pb_nems+geom_bar(aes(color=Rank5, fill=Rank5), stat = "identity", position="stack")
gb_nems


############################# Remove unobserved OTUs - top taxa - Rank2 #############################

#extract top 10 OTUs - I don't know if this is really doing top ten taxa. Whatever the legends come out like, I'd say use that number being the top X.
TopTenOTUs = names(sort(taxa_sums(physeq_metazoaRemoved_trans_trimmed_glom), TRUE)[1:10])
head(TopTenOTUs)
physeqtop10_metazoaRemoved = prune_taxa(TopTenOTUs, physeq_metazoaRemoved_trans_trimmed_glom)
rowSums(otu_table(physeqtop10_metazoaRemoved))
pb_meta10 <- plot_bar(physeqtop10_metazoaRemoved, x = "OceanRegion", y ="Abundance", fill = "Rank2", facet_grid=~Habitat)
gb_meta10 <- pb_meta10+geom_bar(aes(color=Rank2, fill=Rank2), stat = "identity", position="stack")
gb_meta10
