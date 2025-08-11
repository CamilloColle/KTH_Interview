library("Rcpp")
library("dada2") 
library("phyloseq") 
library("ggplot2") 
library("readr") 
library("reshape2") 
library("RColorBrewer") 
library("ggpubr") 
library("permute")
library("lattice")
library("vegan")
library("qiime2R")

###INPUT FILES
setwd("C:/Users/ccoll/Desktop/Thesis/thesis_master/comparison")
metadata <- "metadata2.txt"

gg_table_qza <- "gg_table.qza"
gg_tree_qza <- "gg_rooted-tree.qza"
gg_taxa_qza <- "gg_taxonomy.qza"

ss_table_qza <- "ss_table.qza"
ss_tree_qza <- "ss_rooted-tree.qza"
ss_taxa_qza <- "ss_taxonomy.qza"

gt_table_qza <- "gt_table.qza"
gt_tree_qza <- "gt_rooted-tree.qza"
gt_taxa_qza <- "gt_taxonomy.qza"

design<-read.delim(metadata)

######################################
#TAXONOMY

#features
ASVs_gg <- read_qza(gg_table_qza)
seqtab.nochim_gg <- as.matrix(t(ASVs_gg$data))

ASVs_ss <- read_qza(ss_table_qza)
seqtab.nochim_ss <- as.matrix(t(ASVs_ss$data))

ASVs_gt <- read_qza(gt_table_qza)
seqtab.nochim_gt <- as.matrix(t(ASVs_gt$data))

#discard ctrl samples
seqtab.nochim_gg <- seqtab.nochim_gg[!grepl('ctrl', rownames(seqtab.nochim_gg)),]
seqtab.nochim_ss <- seqtab.nochim_ss[!grepl('ctrl', rownames(seqtab.nochim_ss)),]
seqtab.nochim_gt <- seqtab.nochim_gt[!grepl('ctrl', rownames(seqtab.nochim_gt)),]

#taxa
taxa_gg <- read_qza(gg_taxa_qza)
taxonomy_gg <- as.matrix(do.call(rbind, strsplit(as.character(taxa_gg$data$Taxon), "; ")))
colnames(taxonomy_gg) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(taxonomy_gg) <- taxa_gg$data$Feature.ID
#write.csv(taxonomy_gg, file = "gg_taxonomy.csv")

taxa_ss <- read_qza(ss_taxa_qza)
taxonomy_ss <- as.matrix(do.call(rbind, strsplit(as.character(taxa_ss$data$Taxon), "; ")))
colnames(taxonomy_ss) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(taxonomy_ss) <- taxa_ss$data$Feature.ID
#write.csv(taxonomy_ss, file = "ss_taxonomy.csv")

taxa_gt <- read_qza(gt_taxa_qza)
taxonomy_gt <- as.matrix(do.call(rbind, strsplit(as.character(taxa_gt$data$Taxon), ";")))
colnames(taxonomy_gt) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(taxonomy_gt) <- taxa_gt$data$Feature.ID
#write.csv(taxonomy_gt, file = "gt_taxonomy.csv")

#Transpose seqtab.nochim and make it show abundance as proportions
t.seqtab.nochim_gg <- as.data.frame(t(seqtab.nochim_gg))
t.seqtab.nochim_gg <- as.data.frame(prop.table(as.matrix(t.seqtab.nochim_gg), margin=2))

t.seqtab.nochim_ss <- as.data.frame(t(seqtab.nochim_ss))
t.seqtab.nochim_ss <- as.data.frame(prop.table(as.matrix(t.seqtab.nochim_ss), margin=2))

t.seqtab.nochim_gt <- as.data.frame(t(seqtab.nochim_gt))
t.seqtab.nochim_gt <- as.data.frame(prop.table(as.matrix(t.seqtab.nochim_gt), margin=2))

#'All columns should therefore add up to 1. The colSums function can be used to check.
colSums(t.seqtab.nochim_gg)
colSums(t.seqtab.nochim_ss)
colSums(t.seqtab.nochim_gt)

#Check how many columns and rows there are in the data
dim(design)

dim(seqtab.nochim_gg)
dim(taxonomy_gg)

dim(seqtab.nochim_ss)
dim(taxonomy_ss)

dim(seqtab.nochim_gt)
dim(taxonomy_gt)

#create one single table from taxonomy and sequence variant table
taxtab_gg <- merge(taxonomy_gg, t.seqtab.nochim_gg, by="row.names")

taxtab_ss <- merge(taxonomy_ss, t.seqtab.nochim_ss, by="row.names")

taxtab_gt <- merge(taxonomy_gt, t.seqtab.nochim_gt, by="row.names")

#'For tables, when you want to get rid of a column, you need to use '<-NULL'...
taxtab_gg$Row.names<-NULL
taxtab_ss$Row.names<-NULL
taxtab_gt$Row.names<-NULL

#enforce each species name has the same prefix
taxtab_gg$Species <- ifelse(grepl("^s__", taxtab_gg$Species), taxtab_gg$Species, "s__")
taxtab_ss$Species <- ifelse(grepl("^s__", taxtab_ss$Species), taxtab_ss$Species, "s__")
taxtab_gt$Species <- ifelse(grepl("^s__", taxtab_gt$Species), taxtab_gt$Species, "s__")

#Aggregate the table by Species
myfun<- function(x) {if (class(x)%in%c("numeric", "integer")) {sum(x)} 
  else {if (length(unique(x))==1) {x[1]}
    else {NA}}}

taxtab_gg <- aggregate(taxtab_gg, by=list(taxtab_gg$Species), FUN=myfun)
taxtab_ss <- aggregate(taxtab_ss, by=list(taxtab_ss$Species), FUN=myfun)
taxtab_gt <- aggregate(taxtab_gt, by=list(taxtab_gt$Species), FUN=myfun)


#Replace the NAs and emtpy rows in the table with unclassified_something
taxtab_gg$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_gg$Phylum<-as.character(taxtab_gg$Phylum)
taxtab_gg$Class<-as.character(taxtab_gg$Class)
taxtab_gg$Order<-as.character(taxtab_gg$Order)
taxtab_gg$Family<-as.character(taxtab_gg$Family)
taxtab_gg$Genus<-as.character(taxtab_gg$Genus)
taxtab_gg$Species<-as.character(taxtab_gg$Species)

taxtab_gg$Species[which(is.na(taxtab_gg$Species))]<-"unclassified species"
taxtab_gg$Species[which(taxtab_gg$Species == "s__")]<-"unclassified species"


taxtab_ss$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_ss$Phylum<-as.character(taxtab_ss$Phylum)
taxtab_ss$Class<-as.character(taxtab_ss$Class)
taxtab_ss$Order<-as.character(taxtab_ss$Order)
taxtab_ss$Family<-as.character(taxtab_ss$Family)
taxtab_ss$Genus<-as.character(taxtab_ss$Genus)
taxtab_ss$Species<-as.character(taxtab_ss$Species)

taxtab_ss$Species[which(is.na(taxtab_ss$Species))]<-"unclassified species"
taxtab_ss$Species[which(taxtab_ss$Species == "s__")]<-"unclassified species"

taxtab_gt$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_gt$Phylum<-as.character(taxtab_gt$Phylum)
taxtab_gt$Class<-as.character(taxtab_gt$Class)
taxtab_gt$Order<-as.character(taxtab_gt$Order)
taxtab_gt$Family<-as.character(taxtab_gt$Family)
taxtab_gt$Genus<-as.character(taxtab_gt$Genus)
taxtab_gt$Species<-as.character(taxtab_gt$Species)

taxtab_gt$Species[which(is.na(taxtab_gt$Species))]<-"unclassified species"
taxtab_gt$Species[which(taxtab_gt$Species == "s__")]<-"unclassified species"


#Simplify the table keeping only the Species information
species.tab_gg <- taxtab_gg[, 7:ncol(taxtab_gg)]
row.names(species.tab_gg)<-taxtab_gg$Species

species.tab_ss <- taxtab_ss[, 7:ncol(taxtab_ss)]
row.names(species.tab_ss)<-taxtab_ss$Species

species.tab_gt <- taxtab_gt[, 7:ncol(taxtab_gt)]
row.names(species.tab_gt)<-taxtab_gt$Species

#'At this stage, you can also remove any unnecessary columns. 
#'We no longer need species or species for instance, as the species is now the row name of the table.
species.tab_gg$Species<-NULL
species.tab_ss$Species<-NULL
species.tab_gt$Species<-NULL

#full_genera <- full_join(species.tab_gg, species.tab_gt, species.tab_ss, by = "Species")

# species.tab_gg$avg <- rowMeans(species.tab_gg)
# species.tab_ss$avg <- rowMeans(species.tab_ss)
# species.tab_gt$avg <- rowMeans(species.tab_gt)
# 
# species.tab_gg$avg[which(rownames(species.tab_gg) == "unclassified species")]
# species.tab_ss$avg[which(rownames(species.tab_ss) == "unclassified species")]
# species.tab_gt$avg[which(rownames(species.tab_gt) == "unclassified species")]

#uncomment to save these tables
#write.csv(genus.tab_gg, "genus_tab_gg.csv")
#write.csv(genus.tab_ss, "genus_tab_ss.csv")
#write.csv(genus.tab_gt, "genus_tab_gt.csv")


#Order taxa by mean abundance
species.tab_gg <-species.tab_gg[order(rowMeans(species.tab_gg), decreasing = T),]
species.tab_ss <-species.tab_ss[order(rowMeans(species.tab_ss), decreasing = T),]
species.tab_gt <-species.tab_gt[order(rowMeans(species.tab_gt), decreasing = T),]

#Group together genera under 1% mean representation and check
species.short_gg <- species.tab_gg[rowMeans(species.tab_gg)>=0.01,]
others_gg <- colSums(species.tab_gg[rowMeans(species.tab_gg)<0.01, ])
species.short_gg[nrow(species.short_gg)+1,] <- others_gg
row.names(species.short_gg)[nrow(species.short_gg)] <- "under_1%"
colSums(species.short_gg)

species.short_ss <- species.tab_ss[rowMeans(species.tab_ss)>=0.01,]
others_ss <- colSums(species.tab_ss[rowMeans(species.tab_ss)<0.01, ])
species.short_ss[nrow(species.short_ss)+1,] <- others_ss
row.names(species.short_ss)[nrow(species.short_ss)] <- "under_1%"
colSums(species.short_ss)

species.short_gt <- species.tab_gt[rowMeans(species.tab_gt)>=0.01,]
others_gt <- colSums(species.tab_gt[rowMeans(species.tab_gt)<0.01, ])
species.short_gt[nrow(species.short_gt)+1,] <- others_gt
row.names(species.short_gt)[nrow(species.short_gt)] <- "under_1%"
colSums(species.short_gt)


#Put species names in a column
species.short_gg$species <- row.names(species.short_gg)

species.short_ss$species <- row.names(species.short_ss)

species.short_gt$species <- row.names(species.short_gt)


#Transform the table to a long format for stacked barplot
library(tidyr)
library(dplyr)
species.long_gg <- gather(species.short_gg, key = sample, value = Proportion, -species)

species.long_ss <- gather(species.short_ss, key = sample, value = Proportion, -species)

species.long_gt <- gather(species.short_gt, key = sample, value = Proportion, -species)

#Order the species in the same way they were ordered in the previous table (by abundance)
species.long_gg$species <- factor(species.long_gg$species, levels = species.short_gg$species)

species.long_ss$species <- factor(species.long_ss$species, levels = species.short_ss$species)

species.long_gt$species <- factor(species.long_gt$species, levels = species.short_gt$species)

#Now we can merge the design table with species.long. 
final.tab_gg <- inner_join(species.long_gg, design, by=c("sample" = "Subject"))

#we then manually add a "Database" column corresponding to the respective database
final.tab_gg$Database <- "GreenGenes2"

#uncomment to save table
#write.csv(final.tab_gg, file = "gg_taxa_table.csv")

#repeat for SILVA and GTDB
final.tab_ss <- inner_join(species.long_ss, design, by=c("sample" = "Subject"))
final.tab_ss$Database <- "SILVA"
#write.csv(final.tab_ss, file = "ss_taxa_table.csv")

final.tab_gt <- inner_join(species.long_gt, design, by=c("sample" = "Subject"))
final.tab_gt$Database <- "GTDB"
#write.csv(final.tab_gt, file = "gt_taxa_table.csv")

# Combining the DataFrames
final.tab <- rbind(final.tab_gg, final.tab_ss, final.tab_gt)

#get rid of prefixes
final.tab$species <- gsub('^s__', '', final.tab$species)

#uncomment save the table that will be plotted
#write.csv(final.tab, file = "combined_species_table.csv")

###PLOTS
library(forcats)

##Define order of the columns

# First, calculate the sum of 'Proportion' for 'unclassified species' within each 'sample' and 'Database'.
ranked_samples <- final.tab %>%
  filter(species == "unclassified species") %>%
  group_by(Database, sample) %>%
  summarise(Total_Proportion = sum(Proportion)) %>%
  ungroup() %>%
  arrange(Database, desc(Total_Proportion)) %>%
  mutate(Rank = row_number())

# Next, join this back to your original data to get the ranking for each sample.
final.tab <- final.tab %>%
  left_join(ranked_samples, by = c("Database", "sample"))

#relevel categories to determine where to position them in the plot
final.tab$species <- fct_relevel(final.tab$species, "unclassified species")
final.tab$species <- fct_rev(final.tab$species)


#Black Plot
colors <- setNames(rep("black", 48), levels(final.tab$species))
colors["unclassified species"] <- "red"  # Specify the "unclassified species" group to be red

faceted_plot_black <- ggplot(final.tab, aes(x = reorder(sample, -Rank), y = Proportion, fill = species)) +
  geom_bar(stat = "identity", , position = "stack", color = "white") +
  xlab("Samples") +
  scale_fill_manual(values = colors) +
  theme_classic2() +
  theme(plot.margin=unit(c(0,0,0,0), "cm"),
        text = element_text(size = 13), 
        legend.position = "bottom", 
        axis.text.x = element_blank()) +
  facet_grid(~ Database, space = "free", scales = "free", switch = "x") +
  guides(fill = guide_legend(ncol = 3))

faceted_plot_black

#save plot
ggsave(faceted_plot_black, file = "species_barplot_faceted_black.png", dpi=300, height = 200, width=370, units = "mm")


#Colored Plot - UNUSED
# 
# library(scales) # For hue_pal()
# levels <- sort(unique(final.tab$species))
# 
# # Generate a set of colors that ggplot2 might use by default
# default_colors <- hue_pal()(length(levels))
# 
# # Now override the color for "under_1%"
# colors <- setNames(default_colors, levels)
# colors["unclassified species"] <- "red"
# 
# # Now, plot with the samples ordered by this rank. The samples will be ordered within each facet.
# faceted_plot_color <- ggplot(final.tab, aes(x = reorder(sample, -Rank), y = Proportion, fill = species)) +
#   geom_bar(stat = "identity", , position = "stack", color = "white") +
#   xlab("Samples") +
#   scale_fill_manual(values = colors)+
#   theme_classic2() +
#   theme(plot.margin = common_margins,
#         text = element_text(size = 13), 
#         legend.position = "bottom", 
#         axis.text.x = element_blank()) +
#   facet_grid(~ Database, space = "free", scales = "free", switch = "x") +
#   guides(fill = guide_legend(ncol = 3))
# 
# ggsave(faceted_plot_color, file = "phyla_barplot_faceted_color.png", dpi=300, height = 200, width=370, units = "mm")
# 
