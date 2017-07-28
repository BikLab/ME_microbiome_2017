install_github(karthik/wesanderson)
library(devtools)
library(phyloseq)
library(ggplot2)
library(biomformat)
library(wesanderson)

#Sys.setlocale(locale="C")

### import qiime files 
# use final OTU table containing metadata
# use sperate taxanomy file
# no tree file
# merge files into phyloseq object

map <- import_qiime_sample_data("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/16S-bac-QIIME-mapping-MEmicrobiome-combined-RC.txt")
taxonomy <- import_qiime_sample_data("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/rep_set_tax_assignments_fixed.txt")
dat <- read_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2.biom")
#transform_sample_counts(function(x) {x/sum(x)} )
otu_table <- as.data.frame(as.matrix(biom_data(dat)))

OTU = otu_table(otu_table, taxa_are_rows =TRUE)

taxa <- as.matrix(taxonomy)
TAX = tax_table(taxa)

physeq = merge_phyloseq(OTU,TAX,map)
#physeq <- transform_sample_counts(physeq_merge, function(x) {x / sum(x)} )

print(physeq)

physeq_pruned = prune_taxa("Phylum", 20)
physeq.merged = merge_taxa(physeq, taxa_names(physeq)[1:5])
print(physeq.merged)
rank_names(physeq)

par(mar=c(4,4.5,2,1))
par(oma=c(0,0,0,0) )
#plot_bar(physeq, fill = "Class")

### relative abundance graph
plot = plot_bar(physeq, x = "Sample", y = "Abundance", fill = "Phylum")
plot

plot + geom_bar(aes(color=Phylum, fill=Phylum), stat = "identity", position = "stack") +
  scale_fill_manual(values = wes_palette(n=63, name="FantasticFox", type = "continuous")) + 
  scale_color_manual(values = wes_palette(n=63, name="FantasticFox", type = "continuous"))
