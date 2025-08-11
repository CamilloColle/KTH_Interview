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
metadata <- "metadata.txt"

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

#enforce each phylum name has the same prefix
taxtab_gg$Phylum <- ifelse(grepl("^p__", taxtab_gg$Phylum), taxtab_gg$Phylum, "p__")
taxtab_ss$Phylum <- ifelse(grepl("^p__", taxtab_ss$Phylum), taxtab_ss$Phylum, "p__")
taxtab_gt$Phylum <- ifelse(grepl("^p__", taxtab_gt$Phylum), taxtab_gt$Phylum, "p__")

#Aggregate the table by phylum
myfun<- function(x) {if (class(x)%in%c("numeric", "integer")) {sum(x)} 
  else {if (length(unique(x))==1) {x[1]}
    else {NA}}}

taxtab_gg <- aggregate(taxtab_gg, by=list(taxtab_gg$Phylum), FUN=myfun)
taxtab_ss <- aggregate(taxtab_ss, by=list(taxtab_ss$Phylum), FUN=myfun)
taxtab_gt <- aggregate(taxtab_gt, by=list(taxtab_gt$Phylum), FUN=myfun)


#Replace the NAs and empty rows in the table with unclassified_something
taxtab_gg$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_gg$Phylum<-as.character(taxtab_gg$Phylum)
taxtab_gg$Class<-as.character(taxtab_gg$Class)
taxtab_gg$Order<-as.character(taxtab_gg$Order)
taxtab_gg$Family<-as.character(taxtab_gg$Family)
taxtab_gg$Genus<-as.character(taxtab_gg$Genus)
taxtab_gg$Species<-as.character(taxtab_gg$Species)

taxtab_gg$Phylum[which(is.na(taxtab_gg$Phylum))]<-"unclassified phylum"
taxtab_gg$Phylum[which(taxtab_gg$Phylum == "p__")]<-"unclassified phylum"


taxtab_ss$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_ss$Phylum<-as.character(taxtab_ss$Phylum)
taxtab_ss$Class<-as.character(taxtab_ss$Class)
taxtab_ss$Order<-as.character(taxtab_ss$Order)
taxtab_ss$Family<-as.character(taxtab_ss$Family)
taxtab_ss$Genus<-as.character(taxtab_ss$Genus)
taxtab_ss$Species<-as.character(taxtab_ss$Species)

taxtab_ss$Phylum[which(is.na(taxtab_ss$Phylum))]<-"unclassified phylum"
taxtab_ss$Phylum[which(taxtab_ss$Phylum == "p__")]<-"unclassified phylum"

taxtab_gt$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_gt$Phylum<-as.character(taxtab_gt$Phylum)
taxtab_gt$Class<-as.character(taxtab_gt$Class)
taxtab_gt$Order<-as.character(taxtab_gt$Order)
taxtab_gt$Family<-as.character(taxtab_gt$Family)
taxtab_gt$Genus<-as.character(taxtab_gt$Genus)
taxtab_gt$Species<-as.character(taxtab_gt$Species)

taxtab_gt$Phylum[which(is.na(taxtab_gt$Phylum))]<-"unclassified phylum"
taxtab_gt$Phylum[which(taxtab_gt$Phylum == "p__")]<-"unclassified phylum"


#Simplify the table keeping only the phylum information
phylum.tab_gg <- subset(taxtab_gg, select = -c(Group.1, Class, Order, Family, Genus, Species))
row.names(phylum.tab_gg) <- phylum.tab_gg$Phylum

phylum.tab_ss <- subset(taxtab_ss, select = -c(Group.1, Class, Order, Family, Genus, Species))
row.names(phylum.tab_ss)<-phylum.tab_ss$Phylum

phylum.tab_gt <- subset(taxtab_gt, select = -c(Group.1, Class, Order, Family, Genus, Species))
row.names(phylum.tab_gt)<-phylum.tab_gt$Phylum

#'At this stage, you can also remove any unnecessary columns. 
#'We no longer need genus or species for instance, as the genus is now the row name of the table.
phylum.tab_gg$Phylum<-NULL
phylum.tab_ss$Phylum<-NULL
phylum.tab_gt$Phylum<-NULL


#Order taxa by mean abundance
phylum.tab_gg <-phylum.tab_gg[order(rowMeans(phylum.tab_gg), decreasing = T),]
phylum.tab_ss <-phylum.tab_ss[order(rowMeans(phylum.tab_ss), decreasing = T),]
phylum.tab_gt <-phylum.tab_gt[order(rowMeans(phylum.tab_gt), decreasing = T),]

colSums(phylum.tab_gg)
colSums(phylum.tab_ss)
colSums(phylum.tab_gt)

#compute average unclassified rate
# phylum.tab_gg$avg <- rowMeans(phylum.tab_gg)
# phylum.tab_ss$avg <- rowMeans(phylum.tab_ss)
# phylum.tab_gt$avg <- rowMeans(phylum.tab_gt)
# 
# phylum.tab_gg$avg[which(rownames(phylum.tab_gg) == "unclassified phylum")]
# phylum.tab_ss$avg[which(rownames(phylum.tab_ss) == "unclassified phylum")]
# phylum.tab_gt$avg[which(rownames(phylum.tab_gt) == "unclassified phylum")]

#uncomment to save these tables
#write.csv(phylum.tab_gg, "phylum_tab_gg.csv")
#write.csv(phylum.tab_ss, "phylum_tab_ss.csv")
#write.csv(phylum.tab_gt, "phylum_tab_gt.csv")

#Group together genera under 1% mean representation and check
phylum.short_gg <- phylum.tab_gg[rowMeans(phylum.tab_gg)>=0.01,]
others_gg <- colSums(phylum.tab_gg[rowMeans(phylum.tab_gg)<0.01, ])
phylum.short_gg[nrow(phylum.short_gg)+1,] <- others_gg
row.names(phylum.short_gg)[nrow(phylum.short_gg)] <- "under_1%"
colSums(phylum.short_gg)

phylum.short_ss <- phylum.tab_ss[rowMeans(phylum.tab_ss)>=0.01,]
others_ss <- colSums(phylum.tab_ss[rowMeans(phylum.tab_ss)<0.01, ])
phylum.short_ss[nrow(phylum.short_ss)+1,] <- others_ss
row.names(phylum.short_ss)[nrow(phylum.short_ss)] <- "under_1%"
colSums(phylum.short_ss)

phylum.short_gt <- phylum.tab_gt[rowMeans(phylum.tab_gt)>=0.01,]
others_gt <- colSums(phylum.tab_gt[rowMeans(phylum.tab_gt)<0.01, ])
phylum.short_gt[nrow(phylum.short_gt)+1,] <- others_gt
row.names(phylum.short_gt)[nrow(phylum.short_gt)] <- "under_1%"
colSums(phylum.short_gt)


#Put phylum names in a column
phylum.short_gg$phylum <- row.names(phylum.short_gg)

phylum.short_ss$phylum <- row.names(phylum.short_ss)

phylum.short_gt$phylum <- row.names(phylum.short_gt)


#Transform the table to a long format for stacked barplot
library(tidyr)
library(dplyr)
phylum.long_gg <- gather(phylum.short_gg, key = sample, value = Proportion, -phylum)

phylum.long_ss <- gather(phylum.short_ss, key = sample, value = Proportion, -phylum)

phylum.long_gt <- gather(phylum.short_gt, key = sample, value = Proportion, -phylum)

#Order the phyla in the same way they were ordered in the previous table (by abundance)
phylum.long_gg$phylum <- factor(phylum.long_gg$phylum, levels = phylum.short_gg$phylum)

phylum.long_ss$phylum <- factor(phylum.long_ss$phylum, levels = phylum.short_ss$phylum)

phylum.long_gt$phylum <- factor(phylum.long_gt$phylum, levels = phylum.short_gt$phylum)

#'The design file will now be merged with the phylum.long file. 
#'For the 'by = ' parameter, type the column name you want to use to merge the two tables. In this case, sample. 

#Now we can merge the design table with phylum.long. 
final.tab_gg <- inner_join(phylum.long_gg, design, by=c("sample" = "Subject"))

#we then manually add a "Database" column corresponding to the respective database
final.tab_gg$Database <- "GreenGenes2"

#uncomment to save table
#write.csv(final.tab_gg, file = "gg_taxa_table.csv")

#repeat for SILVA and GTDB
final.tab_ss <- inner_join(phylum.long_ss, design, by=c("sample" = "Subject"))
final.tab_ss$Database <- "SILVA"
#write.csv(final.tab_ss, file = "ss_taxa_table.csv")

final.tab_gt <- inner_join(phylum.long_gt, design, by=c("sample" = "Subject"))
final.tab_gt$Database <- "GTDB"
#write.csv(final.tab_gt, file = "gt_taxa_table.csv")

# Combining the DataFrames
final.tab.raw <- rbind(final.tab_gg, final.tab_ss, final.tab_gt)

# Combining the different Frimicutes entries
firmicutes_summed <- final.tab.raw %>%
  filter(grepl("^p__Firmicutes", phylum)) %>%
  group_by(sample, Database) %>%
  summarise(phylum = "p__Firmicutes", Proportion = sum(Proportion), .groups = 'drop')

# Remove original 'p_Firmicutes' entries from the data
data_without_firmicutes <- final.tab.raw %>%
  filter(!grepl("^p__Firmicutes", phylum))

# Combine the 'p_Firmicutes' summed data with the rest of the data
final.tab <- bind_rows(data_without_firmicutes, firmicutes_summed)

#get rid of prefixes
final.tab$phylum <- gsub('^p__', '', final.tab$phylum)

#uncomment save the table that will be plotted
#write.csv(final.tab, file = "combined_phyla_table.csv")

###PLOTS
library(forcats)
common_margins <- margin(t = 1, r = 10, b = 2, l = 10, unit = "mm")

##Define order of the columns

# First, calculate the sum of 'Proportion' for 'unclassified phylum' within each 'sample' and 'Database'.
ranked_samples <- final.tab %>%
  filter(phylum == "unclassified phylum") %>%
  group_by(Database, sample) %>%
  summarise(Total_Proportion = sum(Proportion)) %>%
  ungroup() %>%
  arrange(Database, desc(Total_Proportion)) %>%
  mutate(Rank = row_number())

# Next, join this back to your original data to get the ranking for each sample.
final.tab <- final.tab %>%
  left_join(ranked_samples, by = c("Database", "sample"))

#relevel categories to determine where to position them in the plot
final.tab$phylum <- fct_relevel(final.tab$phylum, "unclassified phylum")
final.tab$phylum <- fct_rev(final.tab$phylum)


#Black Plot
colors <- setNames(rep("black", 48), levels(final.tab$phylum))
colors["unclassified phylum"] <- "red"  # Specify the "unclassified phylum" group to be red

faceted_plot_black <- ggplot(final.tab, aes(x = sample, y = Proportion, fill = phylum)) +
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
ggsave(faceted_plot_black, file = "phylum_barplot_faceted_black.png", dpi=300, height = 200, width=370, units = "mm")


#Colored Plot - UNUSED
# 
# library(scales) # For hue_pal()
# levels <- sort(unique(final.tab$phylum))
# 
# # Generate a set of colors that ggplot2 might use by default
# default_colors <- hue_pal()(length(levels))
# 
# # Now override the color for "under_1%"
# colors <- setNames(default_colors, levels)
# colors["unclassified phylum"] <- "red"
# 
# # Now, plot with the samples ordered by this rank. The samples will be ordered within each facet.
# faceted_plot_color <- ggplot(final.tab, aes(x = reorder(sample, -Rank), y = Proportion, fill = phylum)) +
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
