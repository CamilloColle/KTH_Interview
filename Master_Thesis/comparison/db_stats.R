library("dplyr")
library("tidyr")
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


########################
#Stats Wilcoxon/Krsukal-Wallis

#setup
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


#TAXONOMY

#features
ASVs_gg <- read_qza(gg_table_qza)
seqtab.nochim_gg <- as.matrix(t(ASVs_gg$data))

ASVs_ss <- read_qza(ss_table_qza)
seqtab.nochim_ss <- as.matrix(t(ASVs_ss$data))

ASVs_gt <- read_qza(gt_table_qza)
seqtab.nochim_gt <- as.matrix(t(ASVs_gt$data))

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

# replace unclassified entries with the same value so that they're aggregated together later
taxtab_gg$Genus <- ifelse(grepl("^g__", taxtab_gg$Genus), taxtab_gg$Genus, "g__")
taxtab_ss$Genus <- ifelse(grepl("^g__", taxtab_ss$Genus), taxtab_ss$Genus, "g__")
taxtab_gt$Genus <- ifelse(grepl("^g__", taxtab_gt$Genus), taxtab_gt$Genus, "g__")

#Aggregate the table by genus
myfun<- function(x) {if (class(x)%in%c("numeric", "integer")) {sum(x)} 
  else {if (length(unique(x))==1) {x[1]}
    else {NA}}}

taxtab_gg <- aggregate(taxtab_gg, by=list(taxtab_gg$Genus), FUN=myfun)
taxtab_ss <- aggregate(taxtab_ss, by=list(taxtab_ss$Genus), FUN=myfun)
taxtab_gt <- aggregate(taxtab_gt, by=list(taxtab_gt$Genus), FUN=myfun)


#Replace the NAs in the table with unclassified_something
taxtab_gg$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_gg$Phylum<-as.character(taxtab_gg$Phylum)
taxtab_gg$Class<-as.character(taxtab_gg$Class)
taxtab_gg$Order<-as.character(taxtab_gg$Order)
taxtab_gg$Family<-as.character(taxtab_gg$Family)
taxtab_gg$Genus<-as.character(taxtab_gg$Genus)
taxtab_gg$Species<-as.character(taxtab_gg$Species)

taxtab_gg$Genus[which(is.na(taxtab_gg$Genus))] <- "unclassified genus"
taxtab_gg$Genus[which(taxtab_gg$Genus == "g__")] <- "unclassified genus"


taxtab_ss$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_ss$Phylum<-as.character(taxtab_ss$Phylum)
taxtab_ss$Class<-as.character(taxtab_ss$Class)
taxtab_ss$Order<-as.character(taxtab_ss$Order)
taxtab_ss$Family<-as.character(taxtab_ss$Family)
taxtab_ss$Genus<-as.character(taxtab_ss$Genus)
taxtab_ss$Species<-as.character(taxtab_ss$Species)

taxtab_ss$Genus[which(is.na(taxtab_ss$Genus))] <- "unclassified genus"
taxtab_ss$Genus[which(taxtab_ss$Genus == "g__")] <- "unclassified genus"


taxtab_gt$Kingdom<-NULL #(Only use if the table also has kingdom)
taxtab_gt$Phylum<-as.character(taxtab_gt$Phylum)
taxtab_gt$Class<-as.character(taxtab_gt$Class)
taxtab_gt$Order<-as.character(taxtab_gt$Order)
taxtab_gt$Family<-as.character(taxtab_gt$Family)
taxtab_gt$Genus<-as.character(taxtab_gt$Genus)
taxtab_gt$Species<-as.character(taxtab_gt$Species)

taxtab_gt$Genus[which(is.na(taxtab_gt$Genus))] <- "unclassified genus"
taxtab_gt$Genus[which(taxtab_gt$Genus == "g__")] <- "unclassified genus"


#Simplify the table keeping only the genus information
genus.tab_gg <- taxtab_gg[, 6:ncol(taxtab_gg)]
row.names(genus.tab_gg)<-taxtab_gg$Genus

genus.tab_ss <- taxtab_ss[, 6:ncol(taxtab_ss)]
row.names(genus.tab_ss)<-taxtab_ss$Genus

genus.tab_gt <- taxtab_gt[, 6:ncol(taxtab_gt)]
row.names(genus.tab_gt)<-taxtab_gt$Genus

#'At this stage, you can also remove any unnecessary columns. 
#'We no longer need genus or species for instance, as the genus is now the row name of the table.
genus.tab_gg$Genus<-NULL
genus.tab_gg$Species<-NULL

genus.tab_ss$Genus<-NULL
genus.tab_ss$Species<-NULL

genus.tab_gt$Genus<-NULL
genus.tab_gt$Species<-NULL

#Order taxa by mean abundance
genus.tab_gg <-genus.tab_gg[order(rowMeans(genus.tab_gg), decreasing = T),]
genus.tab_ss <-genus.tab_ss[order(rowMeans(genus.tab_ss), decreasing = T),]
genus.tab_gt <-genus.tab_gt[order(rowMeans(genus.tab_gt), decreasing = T),]

#first transpose
genus_gg_t <- as.data.frame(t(genus.tab_gg))
genus_ss_t <- as.data.frame(t(genus.tab_ss))
genus_gt_t <- as.data.frame(t(genus.tab_gt))

#add Database column to each
genus_gg_t$Database <- 'GreenGenes'
genus_ss_t$Database <- 'SILVA'
genus_gt_t$Database <- 'GTDB'

#add suffix to samples
genus_gg_t$sample_db <- paste(row.names(genus_gg_t), "GreenGenes", sep="_")
genus_ss_t$sample_db <- paste(row.names(genus_ss_t), "SILVA", sep="_")
genus_gt_t$sample_db <- paste(row.names(genus_gt_t), "GTDB", sep="_")

#full join dataframes
combined_data <- full_join(genus_gg_t, genus_ss_t, by = "sample_db", 
                           suffix = c("_GreenGenes",  "_SILVA"))
combined_data <- full_join(combined_data, genus_gt_t, by = "sample_db", 
                           suffix = c("", "_GTDB"))

#get rid of partial Database columns
combined_data$Database_GreenGenes <- NULL
combined_data$Database_SILVA <- NULL
combined_data$Database <- NULL

#store Database and sample columns in vectors
sample_db <- combined_data$sample
db_vec <- rep(c("GreenGenes", "SILVA", "GTDB"), times=c(230, 230, 230)) 

#standardize columns
# Adjust the pattern based on the actual suffixes in your dataset
standardized_col_names <- gsub("_GreenGenes|_SILVA|_GTDB$", "", 
                               names(combined_data))

# Apply the standardized column names to the dataframe
names(combined_data) <- standardized_col_names

#reinsert sample and Database columns
combined_data <- cbind(db_vec, combined_data)

#rename columns
names(combined_data)[names(combined_data) == 'sample'] <- 'sample_db'
names(combined_data)[names(combined_data) == 'db_vec'] <- 'Database'

#pivot data to merge overlapping taxa columns
long_data <- pivot_longer(combined_data, cols = -c(sample_db, Database), 
                          names_to = "taxon", values_to = "abundance")

# Aggregate the data by summing the abundances of each taxon for each sample
aggregated_data <- long_data %>%
  group_by(sample_db, taxon) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = 'drop')

# Reshape the data back to a wide format
final_data <- pivot_wider(aggregated_data, names_from = taxon, 
                          values_from = abundance)

#store dataframe in different variable
genus.met <- final_data

#save table
#write.csv(genus.met, "merged_genus_tab.csv")


#reinsert and rename the Database column
genus.met <- cbind(db_vec, genus.met)            
names(genus.met)[names(genus.met) == 'db_vec'] <- 'Database'

#write.csv(genus.met, "merged_genus_tab.csv")

################################################################################

#join with metadata to test other factors
#genus.met <- inner_join(genus.met, design, by = c("sample" = "Subject"))

#Test for all samples 
library(tibble)
results<-data.frame(Genus=character(length = ncol(genus.met)-1),
                    kruskal_p=numeric(length = ncol(genus.met)-1))

for (i in 3:ncol(genus.met)){
  results[i-1,1]<-colnames(genus.met)[i]
  kruskal.test(genus.met[,i]~genus.met$Database)->res
  results[i-1,2]<-res$p.value
}

results<-results[(results$Genus=="")==F,]

# Apply Bonferroni correction
sig.results<-results[results$kruskal_p< 0.05/(nrow(results)-1),] #Bonferroni

#create new object for plotting with no prefix
sig.results.clean <- sig.results
sig.results.clean$Genus <- gsub("^g__", "", sig.results.clean$Genus)

#Plot results
kruskal_plot <- ggplot(sig.results.clean, aes(x = reorder(Genus, -kruskal_p), y = kruskal_p)) +
  geom_point() + xlab("Genus") + ylab("Adjusted p-value") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12), #adjust text size if needed
        legend.position="none") # Rotate x labels for better visibility

kruskal_plot

#save plot
#ggsave(kruskal_plot, file = "kruskal_p_adj.png", dpi=300, height = 200, width=370, units = "mm")

#save results in tables
#write.csv(results, "kruskal_results.csv")
#write.csv(sig.results, "kruskal_results_sig.csv")

###Post-Hoc pairwise test
library(dunn.test)

#copy significant results from KW test
dunn_results <- sig.results

# Add needed column in the results table 
dunn_results$GG_GTDB <- NA  
dunn_results$GG_SILVA <- NA 
dunn_results$GTDB_SILVA <- NA  


# Perform Dunn's test for significant taxa
# Dunn's test is performed on every genus and if the adjusted p-value for the 
# specific test is significant then the Z score is added to the table
for(taxon in dunn_results$Genus){
  taxon_data <- genus.met[, c("Database", taxon)]
  dunn_res <- dunn.test(taxon_data[,2], g=taxon_data$Database, method="bonferroni")
  if(dunn_res$P.adjusted[1]<0.05){
    dunn_results[dunn_results$Genus == taxon, "GG_GTDB"] <- dunn_res$Z[1]
  }
  if(dunn_res$P.adjusted[2]<0.05){
    dunn_results[dunn_results$Genus == taxon, "GG_SILVA"] <- dunn_res$Z[2]
  }
  if(dunn_res$P.adjusted[3]<0.05){
    dunn_results[dunn_results$Genus == taxon, "GTDB_SILVA"] <- dunn_res$Z[3]
  }
}

#write.csv(dunn_results, file = "dunn_results.csv")

#plots
data_long <- dunn_results %>%
  pivot_longer(cols = c("GG_GTDB", "GG_SILVA", "GTDB_SILVA"), names_to = "variable", values_to = "value") %>%
  mutate(symbol = case_when(
    value > 0 ~ "+",
    value < 0 ~ "-",
    TRUE ~ ""
  ))

#get rid of prefix
data_long$Genus <- gsub("^g__", "", data_long$Genus)


# Plotting
dunn_plot <- ggplot(data_long, aes(x = Genus, y = variable, label = symbol)) +
  geom_text(aes(color = symbol), size = 6) +
  scale_color_manual(values = c("+" = "green4", "-" = "red2", "NA" = "transparent")) +
  theme_bw() +
  labs(x = "Genus", y = "Variable", title = "Positive and Negative values for each Genus") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12), #adjust text size if needed
        legend.position="none") # Rotate x labels for better visibility

dunn_plot

#save plot
#ggsave(dunn_plot, file = "dunn_plot.png", dpi=300, height = 200, width=370, units = "mm")


#select the most relevant genera
taxa_sel <- c("g__Anaerostipes", "g__Anaerovibrio", "g__Bacteroides_F", "g__Cellulosilyticum",
              "g__Desulfovibrio", "g__Dialister", "g__Fusobacterium_A", "g__Fusobacterium_B",
              "g__Lachnospira", "g__Lactobacillus", "g__Succinivibrio", "g__Treponema_D",
              "g__Turicibacter", "g__Veillonella")

dunn_subset <- dunn_results[dunn_results$Genus%in%taxa_sel,]

#get rid of prefix
dunn_subset$Genus <- gsub("^g__", "", dunn_subset$Genus)

#save table
#write.csv(dunn_subset, file = "dunn_subset.csv")

data_long_subset <- dunn_subset %>%
  pivot_longer(cols = c("GG_GTDB", "GG_SILVA", "GTDB_SILVA"), names_to = "variable", values_to = "value") %>%
  mutate(symbol = case_when(
    value > 0 ~ "+",
    value < 0 ~ "–",
    TRUE ~ ""
  ))

#plot Dunn subset
dunn_plot_sel <- ggplot(data_long_subset, aes(x = Genus, y = variable, label = symbol)) +
  geom_text(aes(color = symbol), size = 12) +
    scale_color_manual(values = c("+" = "green4", "–" = "red2", "NA" = "transparent")) +
  theme_bw() +
  labs(x = "Genus", y = "Variable", title = "Positive and Negative values for each Genus") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 15, face = "bold"),
        legend.position="none") # Rotate x labels for better visibility

dunn_plot_sel

#save plot
ggsave(dunn_plot_sel, file = "dunn_plot_sel.png", dpi=300, height = 200, width=370, units = "mm")

