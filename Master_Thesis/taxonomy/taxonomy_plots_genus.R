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
setwd("C:/Users/ccoll/Desktop/Thesis/thesis_master/comparison_and_taxonomy")
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

#create one single table from taxonomy and sequence variant table
taxtab_ss <- merge(taxonomy_ss, t.seqtab.nochim_ss, by="row.names")

#'For tables, when you want to get rid of a column, you need to use '<-NULL'...
taxtab_ss$Row.names<-NULL


# replace unclassified entries with the same value so that they're aggregated together later
taxtab_ss$Genus <- ifelse(grepl("^g__", taxtab_ss$Genus), taxtab_ss$Genus, "g__")


#Aggregate the table by genus
myfun<- function(x) {if (class(x)%in%c("numeric", "integer")) {sum(x)} 
  else {if (length(unique(x))==1) {x[1]}
    else {NA}}}

taxtab_ss <- aggregate(taxtab_ss, by=list(taxtab_ss$Genus), FUN=myfun)

#Replace the NAs in the table with unclassified_something
taxtab_ss$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_ss$Phylum<-as.character(taxtab_ss$Phylum)
taxtab_ss$Class<-as.character(taxtab_ss$Class)
taxtab_ss$Order<-as.character(taxtab_ss$Order)
taxtab_ss$Family<-as.character(taxtab_ss$Family)
taxtab_ss$Genus<-as.character(taxtab_ss$Genus)
taxtab_ss$Species<-as.character(taxtab_ss$Species)

taxtab_ss$Genus[which(is.na(taxtab_ss$Genus))] <- "unclassified genus"
taxtab_ss$Genus[which(taxtab_ss$Genus == "g__")] <- "unclassified genus"


#Simplify the table keeping only the genus information
genus.tab_ss <- taxtab_ss[, 6:ncol(taxtab_ss)]
row.names(genus.tab_ss)<-taxtab_ss$Genus

#'At this stage, you can also remove any unnecessary columns. 
#'We no longer need genus or species for instance, as the genus is now the row name of the table.
genus.tab_ss$Genus<-NULL
genus.tab_ss$Species<-NULL


#Order taxa by mean abundance
genus.tab_ss <-genus.tab_ss[order(rowMeans(genus.tab_ss), decreasing = T),]
colSums(genus.tab_ss)

#Group together genera under 1% mean representation and check
genus.short_ss <- genus.tab_ss[rowMeans(genus.tab_ss)>=0.01,]
others_ss <- colSums(genus.tab_ss[rowMeans(genus.tab_ss)<0.01, ])
genus.short_ss[nrow(genus.short_ss)+1,] <- others_ss
row.names(genus.short_ss)[nrow(genus.short_ss)] <- "under_1%"
colSums(genus.short_ss)

#Put genus names in a column
genus.short_ss$genus <- row.names(genus.short_ss)

#write.table(genus.short_ss, 'genus_int.tsv', sep = '\t')

#Transform the table to a long format
genus.long_ss <- gather(genus.short_ss, key = sample, value = Proportion, -genus)

#Order the genera in the same way they were ordered in the previous table (by abundance)
genus.long_ss$genus <- factor(genus.long_ss$genus, levels = genus.short_ss$genus)

#Now we can merge the design table with genus.long. 
final.tab <- inner_join(genus.long_ss, design, by=c("sample" = "Subject"))

#write.csv(final.tab, file = "combined_genera_table.csv")

#get rid of prefix
final.tab$genus <- gsub('^g__', '', final.tab$genus)

#rename variables/groups
final.tab$Inside.outside <- gsub('Inside', 'Indoors', final.tab$Inside.outside)
final.tab$Inside.outside <- gsub('Outside', 'Outdoors', final.tab$Inside.outside)
final.tab$Bedding.material <- gsub('Wood chips and straw', 'Wood chips\n and straw', final.tab$Bedding.material)

###PLOTS
##Define order of the columns

# First, calculate the sum of 'Proportion' for 'under_1%' within each 'sample' and 'Database'.
ranked_samples <- final.tab %>%
  filter(genus == "unclassified genus") %>%
  group_by(sample) %>%
  summarise(Total_Proportion = sum(Proportion)) %>%
  ungroup() %>%
  arrange(sample, desc(Total_Proportion)) %>%
  mutate(Rank = row_number())

# Next, join this back to your original data to get the ranking for each sample.
final.tab <- final.tab %>%
  left_join(ranked_samples, by = "sample")

#relevel categories to determine where to position them in the plot
final.tab$genus <- fct_relevel(final.tab$genus, "unclassified genus")
final.tab$genus <- fct_rev(final.tab$genus)

#order groups
final.tab$Age_category <- factor(final.tab$Age_category, 
                                 levels = c("Lactation", "Nursery", "Growing", "Finishing", "Mature"))

final.tab$Inside.outside <- factor(final.tab$Inside.outside, 
                                 levels = c("Indoors", "Outdoors", "Both"))

# #Black Plot - UNUSED
# colors <- setNames(rep("black", 48), levels(final.tab$genus))
# colors["unclassified genus"] <- "red"  # Specify the "under_1%" group to be red
# 
# faceted_plot_black <- ggplot(final.tab, aes(x = reorder(sample, -Rank), y = Proportion, fill = genus)) +
#   geom_bar(stat = "identity", , position = "stack", color = "white") +
#   xlab("Samples") +
#   scale_fill_manual(values = colors) +
#   theme_classic2() +
#   theme(text = element_text(size = 13), 
#         legend.position = "bottom", 
#         axis.text.x = element_blank()) +
#   facet_grid(~ Age_category, space = "free", scales = "free", switch = "x") +
#   guides(fill = guide_legend(ncol = 3))
# 
# ggsave(faceted_plot_black, file = "silva_genera_barplot_black.png", dpi=300, height = 200, width=370, units = "mm")
# 

#Colored Plot
levels <- sort(unique(final.tab$genus))

# Generate a set of colors that ggplot2 might use by default
default_colors <- hue_pal()(length(levels))

# Now override the color for "under_1%"
colors <- setNames(default_colors, levels)
colors["unclassified genus"] <- "red"

#set seed and shuffle colors to avoid similar colors being adjacent
set.seed(50)
colors <- setNames(default_colors, sample(levels))

# Now, plot with the samples ordered by this rank. The samples will be ordered within each facet.
faceted_plot_color <- ggplot(final.tab, aes(x = reorder(sample, -Rank), y = Proportion, fill = genus)) +
  geom_bar(stat = "identity", , position = "stack", color = "white") +
  xlab("Samples") +
  scale_fill_manual(values = colors)+
  theme_classic2() +
  theme(text = element_text(size = 13), 
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "bottom") +
  facet_grid(~ Inside.outside, space = "free", scales = "free", switch = "x")
  #guides(fill = guide_legend(ncol = 3)) #enforce number of columns
faceted_plot_color

#save plot
ggsave(faceted_plot_color, file = "silva_genus_barplot_outliers.png", dpi=300, height = 200, width=370, units = "mm")


###PLOT AVERAGE ABUNDANCES PER GROUP
avg.table <- final.tab %>%
  group_by(genus, Age_category) %>% #change to the variable you want to compute avg for
  summarise(Average_Abundance = mean(Proportion)) %>%
  ungroup()


#define groups order
#uncomment respective variable
avg.table$Age_category <- factor(avg.table$Age_category, 
                                 levels = c("Lactation", "Nursery", "Growing", "Finishing", "Mature"))

# avg.table$Inside.outside <- factor(avg.table$Inside.outside, 
#                                    levels = c("Indoors", "Outdoors", "Both"))

#produce colors
levels <- unique(avg.table$genus)
levels <- sample(levels)

# Generate a set of colors that ggplot2 might use by default
default_colors <- hue_pal()(length(levels))

#set seed and shuffle colors to avoid similar colors being adjacent
set.seed(116)
colors <- setNames(default_colors, sample(levels))

#plot average barplot
faceted_plot_avg <- ggplot(avg.table, aes(x = Age_category, y = Average_Abundance, fill = genus)) +
  geom_bar(stat = "identity", , position = "stack", color = "white") +
  xlab("Age Category") +
  ylab("Average Proportional Abundance") +
  scale_fill_manual(values = colors)+
  theme_classic2() +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(angle = 315, hjust = 0, face = "bold"),
        legend.position = "right")
  #guides(fill = guide_legend(ncol = 3)) #enforce number of columns

faceted_plot_avg

#saveplot
ggsave(faceted_plot_avg, file = "genus_barplot_avg_age_cat_ugly.png", dpi=300, height = 200, width=370, units = "mm")

#write tex table to include in Latex report
# colnames(avg.table) <- c('Genus', 'Age Category', 'Average Proportion')
# in_out_order <- c("Indoors", "Outdoors", "Both")
# age_order <- c("Lactation", "Nursery", "Growing", "Finishing", "Mature")
# 
# avg.table$Inside.outside <- factor(avg.table$Inside.outside, levels=(in_out_order))
# 
# df_wide <- spread(avg.table, key = Inside.outside, value = Average_Abundance)
# 
# print(xtable(df_wide, type = "latex"), file = "avg_genus_in_out.tex")
