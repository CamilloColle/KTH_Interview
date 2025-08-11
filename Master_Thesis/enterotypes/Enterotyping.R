library(cluster)
library(MASS)
library(clusterSim)
library(ade4)
library(dplyr)
library(ggplot2)
library(factoextra)
library(phyloseq)

#SETUP
main_dir <- "C:/Users/ccoll/Desktop/Thesis/thesis_master/"
setwd(main_dir)

#INPUT DATA
#in this case the input is the complete genus table (from SILVA)
clusterdata <- read.table(paste0(main_dir, "enterotypes/genus_tab_full.tsv"), sep = '\t')

#clean up input data as necessary

#remove X from column names
colnames(clusterdata) <- gsub("^X", "", colnames(clusterdata))

#remove ctrl samples
clusterdata <- clusterdata[,!grepl('ctrl', colnames(clusterdata))]

#backup clean data to avoid redoing the cleaning
clusterdata_backup <- clusterdata

#define needed functions
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}


#compute distances
data.denoised=noise.removal(clusterdata, percent=0.01)

data.dist=dist.JSD(clusterdata)

#compute ideal number of clusters
fviz_nbclust(t(clusterdata), FUNcluster = pam, method="silhouette") +
  theme_classic()

#compute clusters
data.cluster=pam.clustering(data.dist, k=2)

#verify clusters via silhouette
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
obs.silhouette #0.174


#PCoA - alternative approach
obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), label = names(data.dist), grid=F,sub="Principal coordiante analysis")

#extract clustering information
#make cluster and saample vectors
cluster_vec <- data.cluster
sample_vec <- rownames(t(clusterdata))

#build df from defined vectors
cluster_df <- data.frame(sample_vec, cluster_vec)
colnames(cluster_df) <- c("sample", "cluster")

#VISUALISATIONS

#import ps object
ps <- readRDS(paste0(main_dir, "enterotypes/ps_latest.RDS"))

#include rownames as a column
sample_data(ps)$sample <- rownames(sample_data(ps))

#extract sample data to inspect
sam <- sample_data(ps)

#include clustering information in ps object
sample_data(ps)$enterotype <- as.factor(data.cluster)

#sam_cluster <- inner_join(sam, cluster_df, by=("samples" = "sample"))

#Log transformation
ps_log <- transform_sample_counts(ps, function(x) log(1+x))


#now use PCoA (BC index) to visualise the clustering
pcoa_bray <- ordinate(ps_log, "MDS", "bray")

pcoa_plot <- plot_ordination(ps_log, pcoa_bray, type = "samples", shape = "enterotype" ,
                             color = "enterotype", label = "Age_category") +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 1, level = 0.95)

pcoa_plot

#save plot
ggsave(pcoa_plot, file = "beta_entero_full_bc.png", dpi=300, height = 200, width=370, units = "mm")


# #PERMANOVA
# library(vegan)
# seqtab_log<-as(otu_table(ps_log), "matrix")
# taxonomy_log<-as(tax_table(ps_log), "matrix")
# design_log<-as(sample_data(ps_log), "data.frame")
# 
# #dissimilarity matrix
# dist.bray<-vegdist(seqtab_log, method = "bray")
# 
# #define variable
# groups <- as.factor(design_log$Age_category)
# levels(groups)
# 
# # PERMANOVA test
# adonis <- adonis2(dist.bray ~ groups, design_log, permutations = 1000)


###Merge data frames for subsetting based on variables
taxa_tab <- read.table('enterotypes/merged_genus_table.tsv', sep = '\t')

final.tab <- inner_join(taxa_tab, cluster_df, by = "sample")

#rename values for better visualisation
final.tab$cluster[which(final.tab$cluster == '1')] <- 'Enterotype 1'
final.tab$cluster[which(final.tab$cluster == '2')] <- 'Enterotype 2'

#remove ctrl samples in case they are still present (should not be)
final.tab <- final.tab[!grepl('ctrl', final.tab$sample),]

#get rid of prefix
final.tab$genus <- gsub('^g__', '', final.tab$genus)

#save table
#write.csv(final.tab, file = "combined_genera_table.csv")


###PLOTS
library(forcats)
##Define order of the columns

# First, calculate the sum of 'Proportion' for 'unclassified genus' within each 'sample'.
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

#reorder groups to determine position in the visualisation
final.tab$genus <- fct_relevel(final.tab$genus, "unclassified genus")
final.tab$genus <- fct_rev(final.tab$genus)

#Colored Plot
library(scales) # For hue_pal()
library(ggpubr) 
levels <- sort(unique(final.tab$genus))

# Generate a set of colors that ggplot2 might use by default
default_colors <- hue_pal()(length(levels))

# Now override the color for "unclassified genus"
colors <- setNames(default_colors, levels)
colors["unclassified genus"] <- "red"

# Now, plot with the samples ordered by this rank. The samples will be ordered within each facet.
faceted_plot_color <- ggplot(final.tab, aes(x = reorder(sample, -Rank), y = Proportion, fill = genus)) +
  geom_bar(stat = "identity", , position = "stack", color = "white") +
  xlab("Samples") +
  scale_fill_manual(values = colors)+
  theme_classic2() +
  theme(text = element_text(size = 13), 
        axis.text.x = element_blank(),
        legend.position = "bottom") +
  facet_grid(~ cluster, space = "free", scales = "free", switch = "x")
  #guides(fill = guide_legend(ncol = 3)) #enforce number of columns

faceted_plot_color

#save plot
#ggsave(faceted_plot_color, file = "enterotypes/genera_faceted_entero_full.png", dpi=300, height = 200, width=370, units = "mm")


#DEFINING SUBSETS - this step is needed as input for the Re-enterotyping script
#subsetting samples by enterotype
metadata <- read.table(paste0(main_dir, 'input_files/metadata2.txt'), sep = '\t', header = T)

#select  columns used for subsetting
metadata <- metadata[,c("Subject", "Age_category")]

#join with prvious table
meta_joined <- inner_join(metadata, cluster_df, by =  c('Subject' = 'sample'))

#define enterotype-based subsets
cluster_1 <- meta_joined[which(meta_joined$cluster == 1),]
cluster_2 <- meta_joined[which(meta_joined$cluster == 2),]

#define age category-based subsets
cluster_young <- meta_joined[which(meta_joined$Age_category %in% c("Lactation", "Nursery")),]
cluster_mature <- meta_joined[which(!meta_joined$Age_category %in% c("Lactation", "Nursery")),]

#lactation/non-lactation subsets
#cluster_lactation <- meta_joined[which(meta_joined$Age_category %in% c("Lactation")),]
#cluster_non_lactation <- meta_joined[which(!meta_joined$Age_category %in% c("Lactation")),]

###VISUALISE YOUNG/MATURE

#First PCoA for both
#YOUNG
ps_young <- ps
sample_data(ps_young) <- sample_data(ps_young)[which(sample_data(ps_young)$sample %in% cluster_young$Subject),]
sam_young <- sample_data(ps_young)

#Log transformation
ps_log_y <- transform_sample_counts(ps_young, function(x) log(1+x))

pcoa_bray_y <- ordinate(ps_log_y, "MDS", "bray")

pcoa_plot_young <- plot_ordination(ps_log_y, pcoa_bray_y, type = "samples", shape = "enterotype" ,
                             color = "enterotype", label = "sample") +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 1, level = 0.95)

pcoa_plot_young
ggsave(pcoa_plot_young, file = "enterotypes/beta_entero_young_bc.png", dpi=300, height = 200, width=370, units = "mm")


#MATURE
ps_mature <- ps
sample_data(ps_mature) <- sample_data(ps_mature)[which(sample_data(ps_mature)$sample %in% cluster_mature$Subject),]
sam_mature <- sample_data(ps_mature)

#Log transformation
ps_log_m <- transform_sample_counts(ps_mature, function(x) log(1+x))

pcoa_bray_m <- ordinate(ps_log_m, "MDS", "bray")

pcoa_plot_mature <- plot_ordination(ps_log_m, pcoa_bray_m, type = "samples", shape = "enterotype" ,
                                   color = "enterotype", label = "sample") +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 1, level = 0.95)

pcoa_plot_mature
ggsave(pcoa_plot_mature, file = "enterotypes/beta_entero_mature_bc.png", dpi=300, height = 200, width=370, units = "mm")



#Then taxonomy barplots for both
library(forcats)
library(scales) # For hue_pal()
library(ggpubr) 


##YOUNG
final.tab.y <- inner_join(taxa_tab, cluster_young, by = c("sample" = "Subject"))

#rename values for better visualisation
final.tab.y$cluster[which(final.tab.y$cluster == '1')] <- 'Enterotype 1'
final.tab.y$cluster[which(final.tab.y$cluster == '2')] <- 'Enterotype 2'

# First, calculate the sum of 'Proportion' for 'unclassified genus' within each 'sample'.
ranked_samples_y <- final.tab.y %>%
  filter(genus == "unclassified genus") %>%
  group_by(sample) %>%
  summarise(Total_Proportion = sum(Proportion)) %>%
  ungroup() %>%
  arrange(sample, desc(Total_Proportion)) %>%
  mutate(Rank = row_number())

# Next, join this back to your original data to get the ranking for each sample.
final.tab.y <- final.tab.y %>%
  left_join(ranked_samples_y, by = "sample")

final.tab.y$genus <- fct_relevel(final.tab.y$genus, "unclassified genus")
final.tab.y$genus <- fct_rev(final.tab.y$genus)

#get rid of prefix
final.tab.y$genus <- gsub('^g__', '', final.tab.y$genus)

#Colored Plot
levels_y <- sort(unique(final.tab.y$genus))

# Generate a set of colors that ggplot2 might use by default
default_colors_y <- hue_pal()(length(levels_y))

# Now override the color for "unclassified genus"
colors_y <- setNames(default_colors_y, levels_y)
colors_y["unclassified genus"] <- "red"

# Now, plot with the samples ordered by this rank. The samples will be ordered within each facet.
faceted_plot_young <- ggplot(final.tab.y, aes(x = reorder(sample, -Rank), y = Proportion, fill = genus)) +
  geom_bar(stat = "identity", , position = "stack", color = "white") +
  xlab("Samples") +
  scale_fill_manual(values = colors_y)+
  theme_classic2() +
  theme(text = element_text(size = 13), 
        axis.text.x = element_blank(),
        legend.position = "bottom") +
  facet_grid(~ cluster, space = "free", scales = "free", switch = "x")
  #guides(fill = guide_legend(ncol = 3)) #enforce number of columns

faceted_plot_young

ggsave(faceted_plot_young, file = "enterotypes/genera_faceted_entero_young.png", dpi=300, height = 200, width=370, units = "mm")


##MATURE
final.tab.m <- inner_join(taxa_tab, cluster_mature, by = c("sample" = "Subject"))

#rename values for better visualisation
final.tab.m$cluster[which(final.tab.m$cluster == '1')] <- 'Enterotype 1'
final.tab.m$cluster[which(final.tab.m$cluster == '2')] <- 'Enterotype 2'

# First, calculate the sum of 'Proportion' for 'unclassified genus' within each 'sample'.
ranked_samples_m <- final.tab.m %>%
  filter(genus == "unclassified genus") %>%
  group_by(sample) %>%
  summarise(Total_Proportion = sum(Proportion)) %>%
  ungroup() %>%
  arrange(sample, desc(Total_Proportion)) %>%
  mutate(Rank = row_number())

# Next, join this back to your original data to get the ranking for each sample.
final.tab.m <- final.tab.m %>%
  left_join(ranked_samples_m, by = "sample")

final.tab.m$genus <- fct_relevel(final.tab.m$genus, "unclassified genus")
final.tab.m$genus <- fct_rev(final.tab.m$genus)

#get rid of prefix
final.tab.m$genus <- gsub('^g__', '', final.tab.m$genus)

#Colored Plot
levels_m <- sort(unique(final.tab.m$genus))

# Generate a set of colors that ggplot2 might use by default
default_colors_m <- hue_pal()(length(levels_m))

# Now override the color for "unclassified genus"
colors_m <- setNames(default_colors_m, levels_m)
colors_m["unclassified genus"] <- "red"

# Now, plot with the samples ordered by this rank. The samples will be ordered within each facet.
faceted_plot_mature <- ggplot(final.tab.m, aes(x = reorder(sample, -Rank), y = Proportion, fill = genus)) +
  geom_bar(stat = "identity", , position = "stack", color = "white") +
  xlab("Samples") +
  scale_fill_manual(values = colors_m)+
  theme_classic2() +
  theme(text = element_text(size = 13), 
        axis.text.x = element_blank(),
        legend.position = "bottom") +
  facet_grid(~ cluster, space = "free", scales = "free", switch = "x")
  #guides(fill = guide_legend(ncol = 3)) #enforce number of columns

faceted_plot_mature

ggsave(faceted_plot_mature, file = "enterotypes/genera_faceted_entero_mature.png", dpi=300, height = 200, width=370, units = "mm")






