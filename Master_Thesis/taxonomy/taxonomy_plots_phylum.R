library("Rcpp")
library("dada2") 
library("phyloseq") 
library("ggplot2") 
library("readr") 
library("dplyr")
library("tidyr")
library("reshape2") 
library("RColorBrewer") 
library("ggpubr") 
library("permute")
library("lattice")
library("vegan")
library("qiime2R")
library("forcats")
library("scales") # For hue_pal()
library("xtable")


###INPUT FILES
setwd("C:/Users/ccoll/Desktop/Thesis/thesis_master/comparison")
metadata <- "metadata2.txt"

ss_table_qza <- "ss_table.qza"
ss_tree_qza <- "ss_rooted-tree.qza"
ss_taxa_qza <- "ss_taxonomy.qza"

design<-read.delim(metadata)

######################################
#TAXONOMY

#features
ASVs_ss <- read_qza(ss_table_qza)
seqtab.nochim_ss <- as.matrix(t(ASVs_ss$data))

#discard ctrl samples
seqtab.nochim_ss <- seqtab.nochim_ss[!grepl('ctrl', rownames(seqtab.nochim_ss)),]

#taxa
taxa_ss <- read_qza(ss_taxa_qza)
taxonomy_ss <- as.matrix(do.call(rbind, strsplit(as.character(taxa_ss$data$Taxon), "; ")))
colnames(taxonomy_ss) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(taxonomy_ss) <- taxa_ss$data$Feature.ID
#write.csv(taxonomy_ss, file = "ss_taxonomy.csv")

#Transpose seqtab.nochim and make it show abundance as proportions
t.seqtab.nochim_ss <- as.data.frame(t(seqtab.nochim_ss))
t.seqtab.nochim_ss <- as.data.frame(prop.table(as.matrix(t.seqtab.nochim_ss), margin=2))

#this line was used to isolate outliers for more in detail inspection
#t.seqtab.nochim_ss <- t.seqtab.nochim_ss[, c('05BE05', '16IE02', '02BE32', '16IE11', '02BE21', '16IE15')]

#'All columns should therefore add up to 1. The colSums function can be used to check.
colSums(t.seqtab.nochim_ss)

#create one single table from taxonomy and sequence variant table
taxtab_ss <- merge(taxonomy_ss, t.seqtab.nochim_ss, by="row.names")

#'For tables, when you want to get rid of a column, you need to use '<-NULL'...
taxtab_ss$Row.names<-NULL

# replace unclassified entries with the same value so that they're aggregated together later
taxtab_ss$Phylum <- ifelse(grepl("^p__", taxtab_ss$Phylum), taxtab_ss$Phylum, "p__")

#Aggregate the table by phylum
myfun<- function(x) {if (class(x)%in%c("numeric", "integer")) {sum(x)} 
  else {if (length(unique(x))==1) {x[1]}
    else {NA}}}

taxtab_ss <- aggregate(taxtab_ss, by=list(taxtab_ss$Phylum), FUN=myfun)


#Replace the NAs in the table with unclassified_something
taxtab_ss$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_ss$Phylum<-as.character(taxtab_ss$Phylum)
taxtab_ss$Class<-as.character(taxtab_ss$Class)
taxtab_ss$Order<-as.character(taxtab_ss$Order)
taxtab_ss$Family<-as.character(taxtab_ss$Family)
taxtab_ss$Genus<-as.character(taxtab_ss$Genus)
taxtab_ss$Species<-as.character(taxtab_ss$Species)

taxtab_ss$Phylum[which(is.na(taxtab_ss$Phylum))]<-"unclassified phylum"
taxtab_ss$Phylum[which(taxtab_ss$Phylum == "p__")]<-"unclassified phylum"


#Simplify the table keeping only the phylum information
phylum.tab_ss <- subset(taxtab_ss, select = -c(Group.1, Class, Order, Family, Genus, Species))
row.names(phylum.tab_ss)<-phylum.tab_ss$Phylum

#'At this stage, you can also remove any unnecessary columns. 
#'We no longer need genus or species for instance, as the genus is now the row name of the table.
phylum.tab_ss$Phylum<-NULL

#Order taxa by mean abundance
phylum.tab_ss <-phylum.tab_ss[order(rowMeans(phylum.tab_ss), decreasing = T),]

colSums(phylum.tab_ss)

#Group together genera under 1% mean representation and check
phylum.short_ss <- phylum.tab_ss[rowMeans(phylum.tab_ss)>=0.01,]
others_ss <- colSums(phylum.tab_ss[rowMeans(phylum.tab_ss)<0.01, ])
phylum.short_ss[nrow(phylum.short_ss)+1,] <- others_ss
row.names(phylum.short_ss)[nrow(phylum.short_ss)] <- "under_1%"
colSums(phylum.short_ss)

#Put phylum names in a column
phylum.short_ss$phylum <- row.names(phylum.short_ss)
phylum.short_ss$phylum <- gsub('^p__', '', phylum.short_ss$phylum)

#write.table(phylum.short_ss, 'phylum_int.tsv', sep = '\t')

#Transform the table to a long format
phylum.long_ss <- gather(phylum.short_ss, key = sample, value = Proportion, -phylum)

#Order the phyla in the same way they were ordered in the previous table (by abundance)
phylum.long_ss$phylum <- factor(phylum.long_ss$phylum, levels = phylum.short_ss$phylum)

#'The design file will now be merged with the phylum.long file. 
#'For the 'by = ' parameter, type the column name you want to use to merge the two tables. In this case, sample. 

#Now we can merge the design table with phylum.long. 
final.tab <- inner_join(phylum.long_ss, design, by=c("sample" = "Subject"))

#write.csv(final.tab, file = "ss_taxa_table.csv")


#write.csv(final.tab, file = "combined_phyla_table.csv")

#get rid of prefix
final.tab$phylum <- gsub('^p__', '', final.tab$phylum)

#rename variables/groups
final.tab$Inside.outside <- gsub('Inside', 'Indoors', final.tab$Inside.outside)
final.tab$Inside.outside <- gsub('Outside', 'Outdoors', final.tab$Inside.outside)
final.tab$Bedding.material <- gsub('Wood chips and straw', 'Wood chips\n and straw', final.tab$Bedding.material)

###PLOTS

#manually adjust margins if need be
common_margins <- margin(t = 1, r = 10, b = 2, l = 10, unit = "mm")

##Define order of the columns

# First, calculate the sum of 'Proportion' for 'under_1%' within each 'sample' and 'Database'.
ranked_samples <- final.tab %>%
  filter(phylum == "unclassified phylum") %>%
  group_by(sample) %>%
  summarise(Total_Proportion = sum(Proportion)) %>%
  ungroup() %>%
  arrange(sample, desc(Total_Proportion)) %>%
  mutate(Rank = row_number())

# Next, join this back to your original data to get the ranking for each sample.

final.tab$phylum <- fct_relevel(final.tab$phylum, "unclassified phylum")
final.tab$phylum <- fct_rev(final.tab$phylum)


#Colored Plot
levels <- sort(unique(final.tab$phylum))

# Generate a set of colors that ggplot2 might use by default
default_colors <- hue_pal()(length(levels))

# Now override the color for "under_1%"
colors <- setNames(default_colors, levels)
colors["unclassified phylum"] <- "red"

# Now, plot with the samples ordered by this rank. The samples will be ordered within each facet.
faceted_plot_color <- ggplot(final.tab, aes(x =sample, y = Proportion, fill = phylum)) +
  geom_bar(stat = "identity", , position = "stack", color = "white") +
  xlab("Samples") +
  scale_fill_manual(values = colors)+
  theme_classic2() +
  theme(plot.margin = common_margins,
        text = element_text(size = 13), 
        legend.position = "bottom", 
        axis.text.x = element_text(angle = 90),
        #axis.text.x = element_blank()) #remove x axis labels
  ) +
  facet_grid(~ Inside.outside, space = "free", scales = "free", switch = "x")
  #guides(fill = guide_legend(ncol = 3)) #enforce number of columns
  
faceted_plot_color

#save plot
ggsave(faceted_plot_color, file = "phyla_barplot_outliers.png", dpi=300, height = 200, width=370, units = "mm")


###PLOT AVERAGE ABUNDANCES PER GROUP
avg.table <- final.tab %>%
  group_by(phylum, Bedding.material) %>%
  summarise(Average_Abundance = mean(Proportion)) %>%
  ungroup()

#define groups order
#uncomment respective variable

# avg.table$Age_category <- factor(avg.table$Age_category, 
#                                  levels = c("Lactation", "Nursery", "Growing", "Finishing", "Mature"))
# 
# avg.table$Inside.outside <- factor(avg.table$Inside.outside, 
#                                    levels = c("Indoors", "Outdoors", "Both"))

#produce colors
levels <- sort(unique(avg.table$phylum))

# Generate a set of colors that ggplot2 might use by default
default_colors <- hue_pal()(length(levels))
colors <- setNames(default_colors, levels)#in this case no need to shuffle since we have few enough colors

# Now, plot with the samples ordered by this rank. The samples will be ordered within each facet.
faceted_plot_avg <- ggplot(avg.table, aes(x = Bedding.material, y = Average_Abundance, fill = phylum)) +
  geom_bar(stat = "identity", , position = "stack", color = "white") +
  xlab("Bedding Material") +
  ylab("Average Proportional Abundance") +
  scale_fill_manual(values = colors)+
  theme_classic2() +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(angle = 315, hjust = 0, face = "bold"),
        legend.position = "right") 
  #guides(fill = guide_legend(ncol = 3)) #enforce number of columns

faceted_plot_avg

#save plot
ggsave(faceted_plot_avg, file = "phylum_barplot_avg_bedding.png", dpi=300, height = 200, width=370, units = "mm")

#write tex table to include in Latex report 
# colnames(avg.table) <- c('Phylum', 'Indoors or Outdoors', 'Average Proportion')
# print(xtable(avg.table, type = "latex"), file = "avg_phylum_in_out.tex")
