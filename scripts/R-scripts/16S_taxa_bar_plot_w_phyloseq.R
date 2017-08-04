library(devtools)
library(phyloseq)
library(ggplot2)
library(biomformat)
library(wesanderson)

#import mapping file into phyloseq
map <- import_qiime_sample_data("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/16S-bac-QIIME-mapping-MEmicrobiome-FINAL-26Jul17.txt")

#import biom tables into phyloseq
otu_Cervonema <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_Cervonema-nems.biom")
otu_Chromadorida_Desmodorida <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_Chromadorida-Desmodorida-nems.biom")
otu_Chromadoridae <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_Chromadoridae-nems.biom")
otu_Desmoscolecidae <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_Desmoscolecidae-nems.biom")
otu_Enoplida <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_Enoplida-nems.biom")
otu_Monhysterida <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_Monhysterida-nems.biom")
otu_Oxystominidae <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_Oxystominidae-nems.biom")
otu_Sabatieria_Setosabatieria <- import_biom("/Users/Taruna/Dropbox/bik_lab/MEmicrobiome/files-4-phyloseq/16S/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_Sabatieria-Setosabatieria.biom")

#merge mapping file and otu tables into a phyloseq object 
physeq_Cervonema = merge_phyloseq(otu_Cervonema,map)
physeq_Chromadorida_Desmodorida = merge_phyloseq(otu_Chromadorida_Desmodorida, map)
physeq_Chromadoridae = merge_phyloseq(otu_Chromadoridae, map)
physeq_Desmoscolecidae = merge_phyloseq(otu_Desmoscolecidae, map)
physeq_Enoplida = merge_phyloseq(otu_Enoplida, map)
physeq_Monhysterida = merge_phyloseq(otu_Monhysterida, map)
physeq_Oxystominidae = merge_phyloseq(otu_Oxystominidae, map)
physeq_Sabatieria_Setosabatieria = merge_phyloseq(otu_Sabatieria_Setosabatieria, map)

physeq_Cervonema
#print merge phyloseq object to make sure the merging was successful
#print(physeq_Cervonema)

#transform the abundance values into relative abundance
physeq_Cervonema_trans = transform_sample_counts(physeq_Cervonema, function(x) 100 * x/sum(x))
physeq_Chromadorida_Desmodorida_trans = transform_sample_counts(physeq_Chromadorida_Desmodorida, function(x) 100 * x/sum(x))
physeq_Chromadoridae_trans = transform_sample_counts(physeq_Chromadoridae, function(x) 100 * x/sum(x))
physeq_Desmoscolecidae_trans = transform_sample_counts(physeq_Desmoscolecidae, function(x) 100 * x/sum(x))
physeq_Enoplida_trans = transform_sample_counts(physeq_Enoplida, function(x) 100 * x/sum(x))
physeq_Monhysterida_trans = transform_sample_counts(physeq_Monhysterida, function(x) 100 * x/sum(x))
physeq_Oxystominidae_trans = transform_sample_counts(physeq_Oxystominidae, function(x) 100 * x/sum(x))
physeq_Sabatieria_Setosabatieria_trans = transform_sample_counts(physeq_Sabatieria_Setosabatieria, function(x) 100 * x/sum(x))


#run some commands on the phyloseq object 
sample_variables(physeq_Cervonema)
levels(sample_data(physeq_Cervonema)$OceanRegion)
rank_names(physeq_Cervonema)

############################# Rank2 ############################# 
#plot based on raw Rank2 (Phylum) data

pdf('~/Desktop/Cervonema_Rank2.pdf')
plot_bar(physeq_Cervonema_trans, fill="Rank2")
dev.off()

pdf('~/Desktop/Chromadorida_Desmodorida_Rank2.pdf')
plot_bar(physeq_Chromadorida_Desmodorida_trans, fill="Rank2")
dev.off()

pdf('~/Desktop/Chromadoridae_Rank2.pdf')
plot_bar(physeq_Chromadoridae_trans, fill="Rank2")
dev.off()

pdf('~/Desktop/Desmoscolecidae_Rank2.pdf')
plot_bar(physeq_Desmoscolecidae_trans, fill="Rank2")
dev.off()

pdf('~/Desktop/Enoplida_Rank2.pdf')
plot_bar(physeq_Enoplida_trans, fill="Rank2")
dev.off()

pdf('~/Desktop/Monhysterida_Rank2.pdf')
plot_bar(physeq_Monhysterida_trans, fill="Rank2")
dev.off()

pdf('~/Desktop/Oxystominidae_Rank2.pdf')
plot_bar(physeq_Oxystominidae_trans, fill="Rank2")
dev.off()

pdf('~/Desktop/Sabatieria_Setosabatieria_Rank2.pdf')
plot_bar(physeq_Sabatieria_Setosabatieria_trans, fill="Rank2")
dev.off()



######################### Remove unobserved OTUs - Rank2 ######################### 
any(taxa_sums(physeq_Cervonema_trans) == 0)
sum(taxa_sums(physeq_Cervonema_trans) == 0)

physeq_Cervonema_trans_trimmed = prune_taxa(taxa_sums(physeq_Cervonema_trans) > 0, physeq_Cervonema_trans)
plot_bar(physeq_Cervonema_trans_trimmed, x = "Sample", y ="Abundance", fill = "Rank2")


############################# Remove unobserved OTUs - top taxa - Rank2 #############################

#extract top 10 OTUs - I don't know if this is really doing top ten taxa. Whatever the legends come out like, I'd say use that number being the top X.
TopTenOTUs = names(sort(taxa_sums(physeq_Cervonema_trans_trimmed), TRUE)[1:20])
head(TopTenOTUs)
physeqtop10_cervonema = prune_species(TopTenOTUs, physeq_Cervonema_trans_trimmed)
plot_bar(physeqtop10_cervonema, x = "Sample", y ="Abundance", fill = "Rank2")

