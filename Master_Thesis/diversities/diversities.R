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
library("sjstats")


###INPUT FILES
#in this case diversity metrics are computed for results coming from SILVA
setwd("C:/Users/ccoll/Desktop/Thesis/thesis_master/comparison")
metadata <- "metadata2.txt"

#gg_table_qza <- "gg_table.qza"
#gg_tree_qza <- "gg_rooted-tree.qza"
#gg_taxa_qza <- "gg_taxonomy.qza"

ss_table_qza <- "ss_table.qza"
ss_tree_qza <- "ss_rooted-tree.qza"
ss_taxa_qza <- "ss_taxonomy.qza"

#gt_table_qza <- "gt_table.qza"
#gt_tree_qza <- "gt_rooted-tree.qza"
#gt_taxa_qza <- "gt_taxonomy.qza"

design<-read.delim(metadata)

###ps objects

# ps_gg <- qza_to_phyloseq(
#   features = gg_table_qza,
#   taxonomy = gg_taxa_qza,
#   #  tree = gg_tree_qza,
#   metadata = metadata)

ps_ss <- qza_to_phyloseq(
  features = ss_table_qza,
  taxonomy = ss_taxa_qza,
  tree = ss_tree_qza,
  metadata = metadata)
# 
# ps_gt <- qza_to_phyloseq(
#   features = gt_table_qza,
#   taxonomy = gt_taxa_qza,
#   #  tree = gt_tree_qza,
#   metadata = metadata)


#select the ps object of the correspnding database
ps <- ps_ss


###TESTING FOR NORMALITY
otu <- as.data.frame(otu_table(ps))

#subset randomly because the shapiro.test function has a limit of 5000 entries
otu_sub <- otu[sample(rownames(otu), 5000),]


for(i in colnames(otu_sub)){
  print(shapiro.test(as.numeric(otu_sub[,i])))
}


################################################################################
###DIVERSITY INDICES

##ALPHA DIVERSITY

#'The plot_richness package offers several alpha diversity indices to choose from. 
#'These are: Observed, Chao1, ACE, Shannon, Simpson, InvSimpson and Fisher.

#extract sample data information for exploration
sam <- sample_data(ps)

#renaming variables for visualisation
sample_data(ps)$Inside.outside <- gsub("Inside", "Indoors", sample_data(ps)$Inside.outside)
sample_data(ps)$Inside.outside <- gsub("Outside", "Outdoors", sample_data(ps)$Inside.outside)
sample_data(ps)$Bedding.material <- gsub('Wood chips and straw', 'Wood chips\n and straw', sample_data(ps)$Bedding.material)


#alpha-diversity boxplots boxplots
#customised palette
#palette <- c("#BBF5DF", "#F7D1C2")

alpha1 <- plot_richness(ps, x="Age_category", measures = c("Chao1", "Shannon", "Simpson"))+
  geom_boxplot() + 
  stat_compare_means(method = "kruskal") +
  theme(text = element_text(size = 18), 
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold")
        #legend.position = "right"
        ) 
alpha1

#define order of x axis
in_out_order <- c("Indoors", "Outdoors", "Both")
age_order <- c("Lactation", "Nursery", "Growing", "Finishing", "Mature", "NA")

#sort x axis based on previously defined order
#uncomment to reorder respective variable


# alpha1$data$Age_category <- as.character(alpha1$data$Age_category)
# alpha1$data$Age_category <- factor(alpha1$data$Age_category, levels=age_order)

# alpha1$data$Inside.outside <- as.character(alpha1$data$Inside.outside)
# alpha1$data$Inside.outside <- factor(alpha1$data$Inside.outside, levels=in_out_order)



alpha1

#save the plot
ggsave(alpha1, file = "alpha_age_cat.png", dpi=300, height = 200, width=370, units = "mm")


################################################################################
##BETA DIVERSITY


sample_data(ps)$samples <- rownames(sample_data(ps))
sam <- sample_data(ps)


#Log transformation
ps_log <- transform_sample_counts(ps, function(x) log(1+x))

#'To plot a PCoA of the data the argument "MDS" (multi-dimensional scaling) is used. Other measures can be used instead to plot different graphs if required, e.g. NMDS.
#'The argument "bray" is the default distance method (Bray-Curtis). If necessary, this can also be changed. There are many alternatives e.g. unifrac, jsd (Jensen-Shannon divergence), jaccard etc.

#PCoA
#Bray-Curtis

pcoa_bray <- ordinate(ps_log, "MDS", "unifrac")

pcoa_plot <- plot_ordination(ps_log, pcoa_bray, type = "samples", 
                             #shape = "Location",
                             color = "Age_category", label = "samples") +
                             geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 1, level = 0.95) +
  theme(text = element_text(size = 18),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold")
        #legend.position = "right"
  )

#uncomment to reorder respective variable

# pcoa_plot$data$Inside.outside <- as.character(pcoa_plot$data$Inside.outside)
# pcoa_plot$data$Inside.outside <- factor(pcoa_plot$data$Inside.outside, levels=in_out_order)

# pcoa_plot$data$Age_category <- as.character(pcoa_plot$data$Age_category)
# pcoa_plot$data$Age_category <- factor(pcoa_plot$data$Age_category, levels=age_order)

pcoa_plot
ggsave(pcoa_plot, file = "beta_location_unif.png", dpi=600, height = 200, width=370, units = "mm")


#Unifrac

#create a different ps_log object to avoid further recomputing
ps_log_unif <- ps_log

sample_data(ps_log_unif)[sample_data(ps_log_unif) == 0] <- 0.0001
sample_data(ps_log_unif)[is.na(sample_data(ps_log_unif))] <- 0.0001

ps_log_unif <- transform_sample_counts(ps, function(x) log(1+x))

#Unweighted
pcoa_unif <- ordinate(ps_log_unif, "MDS", "unifrac")


pcoa_plot_unif <- plot_ordination(ps_log_unif, pcoa_unif, type="samples", 
                              shape = "Bedding.material",
                              color= "Bedding.material", label = "samples")+
  geom_point(size = 3)+
  stat_ellipse(type = "t", linetype = 1, level = 0.95)

#pcoa_plot2$data$Inside.outside <- as.character(pcoa_plot2$data$Inside.outside)
#pcoa_plot2$data$Inside.outside <- factor(pcoa_plot2$data$Inside.outside, levels=in_out_order)

pcoa_plot_unif
ggsave(pcoa_plot_unif, file = "beta_bedding_unif.png", dpi=300, height = 200, width=370, units = "mm")


#Weighted
pcoa_wunif <- ordinate(ps_log_unif, "MDS", "wunifrac")

pcoa_plot_wunif <- plot_ordination(ps_log_wunif, pcoa_unif, type="samples", 
                              shape = "Bedding.material",
                              color= "Bedding.material", label = "samples")+
  geom_point(size = 3)+
  stat_ellipse(type = "t", linetype = 1, level = 0.95)

#pcoa_plot2$data$Inside.outside <- as.character(pcoa_plot2$data$Inside.outside)
#pcoa_plot2$data$Inside.outside <- factor(pcoa_plot2$data$Inside.outside, levels=in_out_order)

pcoa_plot_wunif
ggsave(pcoa_plot_wunif, file = "beta_bedding_wunif.png", dpi=300, height = 200, width=370, units = "mm")



#'The geom_point function can be used to change the size of the dots on the graph. Or if desired, the shapes/colours of points can also be changed.

#'Finally, to determine statistical significance, a PERMANOVA test can be performed. You will need to check that you use the same distance measures as for the PCoA, e.g. bray in this case.

#PERMANOVA
#extract from phyloseq the tables

seqtab_log<-as(otu_table(ps_log), "matrix")
taxonomy_log<-as(tax_table(ps_log), "matrix")
design_log<-as(sample_data(ps_log), "data.frame")

#dissimilarity matrix
dist.bray<-vegdist(seqtab_log, method = "bray")

#define variable
groups <- as.factor(design_log$Age_category)
levels(groups)

# PERMANOVA test
adonis <- adonis2(dist.bray ~ groups, design_log, permutations = 1000)

