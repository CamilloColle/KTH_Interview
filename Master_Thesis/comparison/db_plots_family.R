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

#enforce each family name has the same prefix
taxtab_gg$Family <- ifelse(grepl("^f__", taxtab_gg$Family), taxtab_gg$Family, "f__")
taxtab_ss$Family <- ifelse(grepl("^f__", taxtab_ss$Family), taxtab_ss$Family, "f__")
taxtab_gt$Family <- ifelse(grepl("^f__", taxtab_gt$Family), taxtab_gt$Family, "f__")

#Aggregate the table by family
myfun<- function(x) {if (class(x)%in%c("numeric", "integer")) {sum(x)} 
  else {if (length(unique(x))==1) {x[1]}
    else {NA}}}

taxtab_gg <- aggregate(taxtab_gg, by=list(taxtab_gg$Family), FUN=myfun)
taxtab_ss <- aggregate(taxtab_ss, by=list(taxtab_ss$Family), FUN=myfun)
taxtab_gt <- aggregate(taxtab_gt, by=list(taxtab_gt$Family), FUN=myfun)


#Replace the NAs and empty rows in the table with unclassified_something
taxtab_gg$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_gg$Phylum<-as.character(taxtab_gg$Phylum)
taxtab_gg$Class<-as.character(taxtab_gg$Class)
taxtab_gg$Order<-as.character(taxtab_gg$Order)
taxtab_gg$Family<-as.character(taxtab_gg$Family)
taxtab_gg$Genus<-as.character(taxtab_gg$Genus)
taxtab_gg$Species<-as.character(taxtab_gg$Species)

taxtab_gg$Family[which(is.na(taxtab_gg$Family))]<-"unclassified family"
taxtab_gg$Family[which(taxtab_gg$Family == "f__")]<-"unclassified family"


taxtab_ss$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_ss$Phylum<-as.character(taxtab_ss$Phylum)
taxtab_ss$Class<-as.character(taxtab_ss$Class)
taxtab_ss$Order<-as.character(taxtab_ss$Order)
taxtab_ss$Family<-as.character(taxtab_ss$Family)
taxtab_ss$Genus<-as.character(taxtab_ss$Genus)
taxtab_ss$Species<-as.character(taxtab_ss$Species)

taxtab_ss$Family[which(is.na(taxtab_ss$Family))]<-"unclassified family"
taxtab_ss$Family[which(taxtab_ss$Family == "f__")]<-"unclassified family"


taxtab_gt$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_gt$Phylum<-as.character(taxtab_gt$Phylum)
taxtab_gt$Class<-as.character(taxtab_gt$Class)
taxtab_gt$Order<-as.character(taxtab_gt$Order)
taxtab_gt$Family<-as.character(taxtab_gt$Family)
taxtab_gt$Genus<-as.character(taxtab_gt$Genus)
taxtab_gt$Species<-as.character(taxtab_gt$Species)


taxtab_gt$Family[which(is.na(taxtab_gt$Family))]<-"unclassified family"
taxtab_gt$Family[which(taxtab_gt$Family == "f__")]<-"unclassified family"


#Simplify the table keeping only the family information
family.tab_gg <- subset(taxtab_gg, select = -c(1,2,3,4,6,7)) #Phylum, Class, Order, Genus, Species
row.names(family.tab_gg) <- family.tab_gg$Family

family.tab_ss <- subset(taxtab_ss, select = -c(1,2,3,4,6,7)) #Phylum, Class, Order, Genus, Species
row.names(family.tab_ss)<-family.tab_ss$Family

family.tab_gt <- subset(taxtab_gt, select = -c(1,2,3,4,6,7)) #Phylum, Class, Order, Genus, Species
row.names(family.tab_gt)<-family.tab_gt$Family

#'At this stage, you can also remove any unnecessary columns. 
#'We no longer need genus or species for instance, as the genus is now the row name of the table.
family.tab_gg$Family<-NULL
family.tab_ss$Family<-NULL
family.tab_gt$Family<-NULL

#compute average unclassified rate
# family.tab_gg$avg <- rowMeans(family.tab_gg)
# family.tab_ss$avg <- rowMeans(family.tab_ss)
# family.tab_gt$avg <- rowMeans(family.tab_gt)
# 
# family.tab_gg$avg[which(rownames(family.tab_gg) == "unclassified family")]
# family.tab_ss$avg[which(rownames(family.tab_ss) == "unclassified family")]
# family.tab_ss$avg[which(rownames(family.tab_gt) == "unclassified family")]

#uncomment to save these tables
#write.csv(family.tab_gg, "family_tab_gg.csv")
#write.csv(family.tab_ss, "family_tab_ss.csv")
#write.csv(family.tab_gt, "family_tab_gt.csv")


#Order taxa by mean abundance
family.tab_gg <-family.tab_gg[order(rowMeans(family.tab_gg), decreasing = T),]
family.tab_ss <-family.tab_ss[order(rowMeans(family.tab_ss), decreasing = T),]
family.tab_gt <-family.tab_gt[order(rowMeans(family.tab_gt), decreasing = T),]

#verify that columns sum up to 1
colSums(family.tab_gg)
colSums(family.tab_ss)
colSums(family.tab_gt)

#Group together genera under 1% mean representation and check
family.short_gg <- family.tab_gg[rowMeans(family.tab_gg)>=0.01,]
others_gg <- colSums(family.tab_gg[rowMeans(family.tab_gg)<0.01, ])
family.short_gg[nrow(family.short_gg)+1,] <- others_gg
row.names(family.short_gg)[nrow(family.short_gg)] <- "under_1%"
colSums(family.short_gg)

family.short_ss <- family.tab_ss[rowMeans(family.tab_ss)>=0.01,]
others_ss <- colSums(family.tab_ss[rowMeans(family.tab_ss)<0.01, ])
family.short_ss[nrow(family.short_ss)+1,] <- others_ss
row.names(family.short_ss)[nrow(family.short_ss)] <- "under_1%"
colSums(family.short_ss)

family.short_gt <- family.tab_gt[rowMeans(family.tab_gt)>=0.01,]
others_gt <- colSums(family.tab_gt[rowMeans(family.tab_gt)<0.01, ])
family.short_gt[nrow(family.short_gt)+1,] <- others_gt
row.names(family.short_gt)[nrow(family.short_gt)] <- "under_1%"
colSums(family.short_gt)


#Put family names in a column
family.short_gg$family <- row.names(family.short_gg)

family.short_ss$family <- row.names(family.short_ss)

family.short_gt$family <- row.names(family.short_gt)


#Transform the table to a long format for stacked barplot
library(tidyr)
library(dplyr)
family.long_gg <- gather(family.short_gg, key = sample, value = Proportion, -family)

family.long_ss <- gather(family.short_ss, key = sample, value = Proportion, -family)

family.long_gt <- gather(family.short_gt, key = sample, value = Proportion, -family)

#Order the families in the same way they were ordered in the previous table (by abundance)
family.long_gg$family <- factor(family.long_gg$family, levels = family.short_gg$family)

family.long_ss$family <- factor(family.long_ss$family, levels = family.short_ss$family)

family.long_gt$family <- factor(family.long_gt$family, levels = family.short_gt$family)

#'The design file will now be merged with the family.long file. 
#'For the 'by = ' parameter, type the column name you want to use to merge the two tables. In this case, sample. 

#Now we can merge the design table with family.long. 
final.tab_gg <- inner_join(family.long_gg, design, by=c("sample" = "Subject"))

#we then manually add a "Database" column corresponding to the respective database
final.tab_gg$Database <- "GreenGenes2"

#uncomment to save table
#write.csv(final.tab_gg, file = "gg_taxa_table.csv")

#repeat for SILVA and GTDB
final.tab_ss <- inner_join(family.long_ss, design, by=c("sample" = "Subject"))
final.tab_ss$Database <- "SILVA"
#write.csv(final.tab_ss, file = "ss_taxa_table.csv")

final.tab_gt <- inner_join(family.long_gt, design, by=c("sample" = "Subject"))
final.tab_gt$Database <- "GTDB"
#write.csv(final.tab_gt, file = "gt_taxa_table.csv")

# Combining the DataFrames by row binding
final.tab <- rbind(final.tab_gg, final.tab_ss, final.tab_gt)

#get rid of the prefix
final.tab$family <- gsub('^f__', '', final.tab$family)

#uncomment save the table that will be plotted
#write.csv(final.tab, file = "combined_families_table.csv")

###PLOTS
library(forcats)
##Define order of the columns

# First, calculate the sum of 'Proportion' for 'unclassified family' within each 'sample' and 'Database'.
ranked_samples <- final.tab %>%
  filter(family == "unclassified family") %>%
  group_by(Database, sample) %>%
  summarise(Total_Proportion = sum(Proportion)) %>%
  ungroup() %>%
  arrange(Database, desc(Total_Proportion)) %>%
  mutate(Rank = row_number())

# Next, join this back to your original data to get the ranking for each sample.
final.tab <- final.tab %>%
  left_join(ranked_samples, by = c("Database", "sample"))

#relevel categories to determine where to position them in the plot
final.tab$family <- fct_relevel(final.tab$family, "unclassified family")
final.tab$family <- fct_rev(final.tab$family)


#Black Plot
colors <- setNames(rep("black", 48), levels(final.tab$family))
colors["unclassified family"] <- "red"  # Specify the "unclassified family" group to be red

faceted_plot_black <- ggplot(final.tab, aes(x = reorder(sample, -Rank), y = Proportion, fill = family)) +
  geom_bar(stat = "identity", , position = "stack", color = "white") +
  xlab("Samples") +
  scale_fill_manual(values = colors) +
  theme_classic2() +
  theme(text = element_text(size = 13), 
        legend.position = "bottom", 
        axis.text.x = element_blank()) +
  facet_grid(~ Database, space = "free", scales = "free", switch = "x") +
  guides(fill = guide_legend(ncol = 3))

faceted_plot_black

#save plot
ggsave(faceted_plot_black, file = "families_barplot_faceted_black.png", dpi=300, height = 200, width=370, units = "mm")


#Colored Plot - UNUSED

# library(scales) # For hue_pal()
# levels <- sort(unique(final.tab$family))
# 
# # Generate a set of colors that ggplot2 might use by default
# default_colors <- hue_pal()(length(levels))
# 
# # Now override the color for "under_0.1%"
# colors <- setNames(default_colors, levels)
# colors["unclassified family"] <- "red"
# 
# # Now, plot with the samples ordered by this rank. The samples will be ordered within each facet.
# faceted_plot_color <- ggplot(final.tab, aes(x = reorder(sample, -Rank), y = Proportion, fill = family)) +
#   geom_bar(stat = "identity", , position = "stack", color = "white") +
#   xlab("Samples") +
#   scale_fill_manual(values = colors)+
#   theme_classic2() +
#   theme(text = element_text(size = 13), 
#         legend.position = "bottom", 
#         axis.text.x = element_blank()) +
#   facet_grid(~ Database, space = "free", scales = "free", switch = "x") +
#   guides(fill = guide_legend(ncol = 3))
# 
# ggsave(faceted_plot_color, file = "families_barplot_faceted_color.png", dpi=300, height = 200, width=370, units = "mm")
# 

