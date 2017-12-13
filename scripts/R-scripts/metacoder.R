install.packages("metacoder")
library(metacoder)
input <- readLines("~/Dropbox/bik_lab/MEmicrobiome/metacoder/otu_table_mc2_w_tax_BlankOTUsRemoved_BlankSamplesRemoved_NObadnems_noZeros_Desmoscolecidae-196-nems-n2.txt")
input
data <- extract_taxonomy(input, regex = "^(.*)\\t(.*)", key = c(id = "obs_info", "class"),
                         class_sep = ";")
data
taxon_data(data)
heat_tree(data, node_size = n_obs, node_label = name, node_color = n_obs)
set.seed(8)
plot <- heat_tree(data, title = "MEMB.nem.196", node_size = n_obs, edge_color = n_supertaxa,
          node_label = name, node_color = n_obs,
          node_color_range = c("cyan", "magenta", "green"),
          edge_color_range   = c("#555555", "#EEEEEE"),
          initial_layout = "reingold-tilford", layout = "davidson-harel",
          overlap_avoidance = 1)

plot
tiff(file = "~/Desktop/Taxonomy-4-MEMB.nem.196.tiff", 
     width=8, height=8, units = "in", res = 300)
plot
dev.off()
