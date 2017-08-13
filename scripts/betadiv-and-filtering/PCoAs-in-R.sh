###### 16S PCoA - 500seqs/sample - 97% OTUs #####

#Load R packages
library("phyloseq")
library("ggplot2")

#Load OTU Table and QIIME mapping file 

otufile = "/Users/hollybik/Dropbox/Projects/MEMB/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017/filtered-biom-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros-JSON.biom"

mapfile = "/Users/hollybik/Dropbox/Projects/MEMB/qiime-files/mapping-files/16S-bac-QIIME-mapping-MEmicrobiome-FINAL-26Jul17-ForR.txt"

# import the data into phyloseq
qiimedata = import_biom(BIOMfilename=otufile)
mapdata = import_qiime(mapfilename=mapfile)

# make a Phyloseq object, this is required for most phyloseq functions
memb = merge_phyloseq(qiimedata, mapdata)

# load in the distance matrix and make it the correct format
memb.unifrac = read.table("/Users/hollybik/Dropbox/Projects/MEMB/analysis-results/uclust_openref97_MEMB1_16s_mergedONLY_26July2017/div-analyses/NObadnems-unifrac-500/unweighted_unifrac_dm.txt", header=TRUE)
memb.unifrac = as.dist(as(memb.unifrac, "matrix"))

# ordinate the distance matrix
# this is Phyloseqs version of 'principal_coordinates.py' but you have a lot more options than just principal coordinates
ord.memb.unifrac = ordinate(memb, "PCoA", memb.unifrac)

# plot the ordination, this is the basic command
memb_unweighted_unifrac = plot_ordination(memb, ord.memb.unifrac)


# Format data for ggplot
memb_unweighted_unifrac = plot_ordination(memb, ord.memb.unifrac, color="FeedingGroup", shape="OceanRegion")

#Make Final PCoA Plot
memb_unweighted_unifrac + geom_point(size=5) + scale_colour_brewer(palette="Paired") + scale_shape_manual(values=c(16,17,8)) + labs(title="16S rRNA OTUs Amplified from Single Nematodes") + theme(plot.title = element_text(size=22))

#Resize Points and Title Size for better export image - FINAL IMAGE USED
memb_unweighted_unifrac + geom_point(size=3) + scale_colour_brewer(palette="Paired") + scale_shape_manual(values=c(16,17,8)) + labs(title="16S rRNA OTUs Amplified from Single Nematodes") + theme(plot.title = element_text(size=16))

# Format data for ggplot - by taxon
memb_unweighted_unifrac_NemFam = plot_ordination(memb, ord.memb.unifrac, color="NematodeFamily", shape="OceanRegion")

#Resize Points and Title Size for better export image - FINAL IMAGE USED
memb_unweighted_unifrac_NemFam + geom_point(size=3) + scale_colour_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5"))+ scale_shape_manual(values=c(16,17,8)) + labs(title="16S rRNA OTUs Amplified from Single Nematodes") + theme(plot.title = element_text(size=16))

#Nematode Family Shapes
memb_unweighted_unifrac_NemFamswapshapes = plot_ordination(memb, ord.memb.unifrac, color="OceanRegion", shape="NematodeFamily")

memb_unweighted_unifrac_NemFamswapshapes + geom_point(size=3) + scale_colour_brewer(palette="Set2") + scale_shape_manual(values=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)) + labs(title="16S rRNA OTUs Amplified from Single Nematodes") + theme(plot.title = element_text(size=12))

#Nematode Family + Feeding Group
memb_unweighted_unifrac_NemFamFG = plot_ordination(memb, ord.memb.unifrac, color="FeedingGroup", shape="NematodeFamily")

memb_unweighted_unifrac_NemFamFG + geom_point(size=2) + scale_colour_brewer(palette="Paired") + scale_shape_manual(values=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)) + labs(title="16S rRNA OTUs Amplified from Single Nematodes") + theme(plot.title = element_text(size=12)) + theme(legend.text=element_text(size=10)) + theme(legend.key.size = unit(0.5, "cm"))


###### 18S PCoA - 500seqs/sample - 99% OTUs #####

#Load R packages
library("phyloseq")
library("ggplot2")

#Load OTU Table and QIIME mapping file 

otufile = "/Users/hollybik/Dropbox/Projects/MEMB/analysis-results/manually-filtered-18S-tables/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros-ALLmetazoaRemoved-JSON.biom"

mapfile = "/Users/hollybik/Dropbox/Projects/MEMB/qiime-files/mapping-files/18S-euk-QIIME-mapping-MEmicrobiome-FINAL-26Jul17-ForR.txt"

# import the data into phyloseq
qiimedata = import_biom(BIOMfilename=otufile)
mapdata = import_qiime(mapfilename=mapfile)

# make a Phyloseq object, this is required for most phyloseq functions
memb = merge_phyloseq(qiimedata, mapdata)

# load in the distance matrix and make it the correct format
memb.unifrac = read.table("/Users/hollybik/Dropbox/Projects/MEMB/analysis-results/uclust_openref99_MEMB1_18s_mergedONLY_26July2017/div-analyses/NObadnems-unifrac-500/unweighted_unifrac_dm.txt", header=TRUE)
memb.unifrac = as.dist(as(memb.unifrac, "matrix"))

# ordinate the distance matrix
# this is Phyloseqs version of 'principal_coordinates.py' but you have a lot more options than just principal coordinates
ord.memb.unifrac = ordinate(memb, "PCoA", memb.unifrac)

# plot the ordination, this is the basic command
memb_unweighted_unifrac = plot_ordination(memb, ord.memb.unifrac)

# But the above command doesn't seem to plot anything, so you have to add other commands
# Look through Phyloseq and ggplot2 documentation to find options for customizing

# Phyloseq ones Julia uses the most are 
# color and shape, where you can make the points different based on categories in the mapping file ie:
memb_unweighted_unifrac = plot_ordination(memb, ord.memb.unifrac, color="FeedingGroup", shape="OceanRegion")

#Make Final PCoA Plot
memb_unweighted_unifrac + geom_point(size=5) + scale_colour_brewer(palette="Paired") + scale_shape_manual(values=c(16,17,8)) + labs(title="18S rRNA OTUs Amplified from Single Nematodes - Metazoan OTUs Removed") + theme(plot.title = element_text(size=22))

#Resize Points and Title Size for better export image - FINAL IMAGE USED
memb_unweighted_unifrac + geom_point(size=3) + scale_colour_brewer(palette="Accent") + scale_shape_manual(values=c(16,17,8)) + labs(title="Non-Metazoan 18S rRNA OTUs from Single Nematodes") + theme(plot.title = element_text(size=12))

#Replot using Nematode ID as color scheme
memb_unweighted_unifrac_NemID = plot_ordination(memb, ord.memb.unifrac, color="NematodeID", shape="OceanRegion")

#Replot using Nematode family - FINAL Figure
memb_unweighted_unifrac_NemIDfam + geom_point(size=2) + scale_color_hue(l = 70, c = 200) + scale_shape_manual(values=c(16,17,8)) + labs(title="Non-Metazoan 18S rRNA OTUs from Single Nematodes") + theme(plot.title = element_text(size=10))




