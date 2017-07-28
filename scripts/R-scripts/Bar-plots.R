library("phyloseq")
packageVersion("phyloseq")
library("ggplot2")
packageVersion("ggplot2")
library(biomformat)
library(devtools)
install.packages("wesanderson")
library(wesanderson)
theme_set(theme_bw())

biomfile <- "/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2_w_tax_by_SampleSite.biom"
biom <- import_biom(biomfile)
sdfile <- "/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/16S-bac-QIIME-mapping-MEmicrobiome-FINAL-26Jul17_by_SampleSite.txt"
sample_metadata <- import_qiime_sample_data(sdfile)
memb <- merge_phyloseq(biom, sample_metadata)
memb
#memb <- transform_sample_counts(memb, function(x) 100 * x/ sum(x) )
#memb <- transform_sample_counts(memb, fun)
memb
# Run the following commands to look at a few things about the phyloseq object
rank_names(memb)
sample_variables(memb)
levels(sample_data(memb)$Habitat)
any(taxa_sums(memb) == 0)
tax_table(memb)[1:3, 1:7]

plot = plot_bar(memb, x = "Sample", y = "Abundance", fill = "Rank2")
plot

plot + geom_bar(aes(color=Rank2, fill=Rank2), stat = "identity", position = "stack") +
  scale_fill_manual(values = wes_palette(name="FantasticFox")) + 
  scale_color_manual(values = wes_palette(name="FantasticFox"))




# readsumsdf = data.frame(nreads = sort(taxa_sums(memb), TRUE), sorted = 1:ntaxa(memb), 
#                         type = "OTUs")
# readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(memb), 
#                                                         TRUE), sorted = 1:nsamples(memb), type = "Samples"))
# title = "Total number of reads"
# p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
# p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")