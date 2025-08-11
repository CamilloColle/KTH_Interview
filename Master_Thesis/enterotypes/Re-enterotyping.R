#Dependent from variables in Enterotyping.R.
#First run that code up to the definition of cluster subsets (line 218-234) and
#use those as input for this step

##INPUT STEP
#define subset vector
subset_vector <- cluster_young$Subject #in this case 'cluster_XXX' is defined in the previous scipt, 
                                        #replace with the subset you want to compute enterotypes for

#reset initial df
clusterdata <- clusterdata_backup

#subset cluster data
clusterdata <- clusterdata[, which(colnames(clusterdata) %in% subset_vector)]

#compute distance matrix and assign clusters
data.dist=dist.JSD(clusterdata)

fviz_nbclust(t(clusterdata), FUNcluster = pam, method="silhouette") +
  theme_classic()

data.cluster=pam.clustering(data.dist, k=2)


#verify clusters via silhouette
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
obs.silhouette

#denoised 1 = 0.083
#1 = 0.099
#2 = 0.113

#young = 0.164
#mature = 0.204

#denoised_young = 0.164
#denoised_mature = 0.205


#build merging df
cluster_vec <- data.cluster
sample_vec <- colnames(clusterdata)
cluster_df <- data.frame(sample_vec, cluster_vec)
colnames(cluster_df) <- c("sample", "cluster")

#merge dfs

ps <- readRDS(paste0(main_dir, "enterotypes/ps_latest.RDS"))
sample_data(ps)$sample <- rownames(sample_data(ps))
sam <- sample_data(ps)

#subset ps object
ps_sub <- ps
sample_data(ps_sub) <- sample_data(ps_sub)[which(sample_data(ps_sub)$sample %in% cluster_df$sample),]

#insert enterotype column
sample_data(ps_sub)$enterotype <- as.factor(data.cluster)
sam_sub <- sample_data(ps_sub)

#at this stage the ps object related to the subset can be saved as an RDS object
#saveRDS(ps_sub, file = "enterotypes/re-enterotyping/ps_enterotype_1.RDS")


###PCoA VISUALISATION
#Log transformation
ps_log <- transform_sample_counts(ps_sub, function(x) log(1+x))

pcoa_bray <- ordinate(ps_log, "MDS", "bray")

pcoa_plot <- plot_ordination(ps_log, pcoa_bray, type = "samples", shape = "enterotype" ,
                               color = "enterotype", label = "sample") +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 1, level = 0.95) +
  scale_color_manual(values = c("#1AA70A", "#5D3A9B")) #alterantive color options

pcoa_plot

#save plot
ggsave(pcoa_plot, file = "enterotypes/re-enterotyping/beta_entero_mature2_bc.png", dpi=300, height = 200, width=370, units = "mm")


###BARPLOT
#taxonomy barplots for young/mature subgroups
library(forcats)
library(scales) # For hue_pal()
library(ggpubr) 

final.tab <- inner_join(taxa_tab, cluster_df, by = "sample")

#rename values for better visualisation
final.tab$cluster[which(final.tab$cluster == '1')] <- 'Enterotype 1'
final.tab$cluster[which(final.tab$cluster == '2')] <- 'Enterotype 2'

#again run to make sure ctrl samples are gone
final.tab <- final.tab[!grepl('ctrl', final.tab$sample),]

#save table
#write.csv(final.tab, file = "combined_genera_table.csv")

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

final.tab$genus <- fct_relevel(final.tab$genus, "unclassified genus")
final.tab$genus <- fct_rev(final.tab$genus)

#get rid of prefix
final.tab$genus <- gsub('^g__', '', final.tab$genus)

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
ggsave(faceted_plot_color, file = "enterotypes/re-enterotyping/genera_faceted_entero_non_lactation.png", dpi=300, height = 200, width=370, units = "mm")

