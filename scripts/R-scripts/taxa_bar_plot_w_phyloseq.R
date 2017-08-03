library(devtools)
library(phyloseq)
library(ggplot2)
library(biomformat)
library(wesanderson)

#import mapping file and biom tables into phyloseq
map <- import_qiime_sample_data("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/16S-bac-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt")
otu_1 <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved.biom")
otu_2 <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_by_SampleSite.biom")
biomformat::read_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_by_SampleSite.biom")
# merge mapping file and biom table into a single phyloseq object
physeq1 = merge_phyloseq(otu_1,map)
physeq2 = merge_phyloseq(otu_2,map)

# print merge phyloseq object to make sure the merging was successful
print(physeq1)
print(physeq2)

# run some commands on the phyloseq object 
sample_variables(physeq1)
levels(sample_data(physeq1)$OceanRegion)
rank_names(physeq1)
names(physeq1) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")


plot_bar(physeq1, fill="Rank2", facet_grid=~OceanRegion)
plot_bar(physeq1, x = "SampleSite", y = "Abundance", fill = "Rank2")


# make basic plots
title = "Figure 1 16S Taxonomic Distribution based by Sample Sites"
plot_bar(physeq1, "SampleSite", fill = "Rank2", title = title)
plot_bar(physeq2, fill = "Rank2", title = title)



# extract top 10 OTUs
TopTenOTUs = names(sort(taxa_sums(physeq2), TRUE)[1:20])
head(TopTenOTUs)
physeqtop10 = prune_species(TopTenOTUs, physeq2)
plot_bar(physeqtop10, x = "Sample", y ="Abundance", fill = "Rank2")
plot + geom_bar(aes(color=Rank2, fill=Rank2), stat = "identity", position = "stack") 
+   scale_fill_manual(values = wes_palette(21, name = "Zissou", type = "continuous")) 
+   scale_color_manual(values = wes_palette(21, name = "Zissou", type = "continuous"))

#g__Sphingomonas

