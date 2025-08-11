library(circlize)
library(dplyr)
library(qiime2R)
library(RColorBrewer)
library(gridBase)
library(ComplexHeatmap)
library(reshape2)
library(ggplot2)

###SETUP
main_dir <- "C:/Users/ccoll/Desktop/Thesis/thesis_master/"
setwd(paste0(main_dir, "correlations"))

###INPUT FILES - the input for this file consists of the metadata and the 'genus.short' object
# produced in the 'taxonomy_plots_genus' file
metadata <- paste0(main_dir, "input_files/metadata2.txt")
metadata <- read.csv(metadata, sep = "\t")

genus.short <- read.table("genus_short.tsv",  sep = '\t')
colnames(genus.short) <- gsub('^X', '', colnames(genus.short))

# Clean metadata to keep only columns of interest
metadata_clean <- metadata[,c("Subject",
                              "Crude.protein", "Crude.fibre", "Crude.fat", "Carbohydrates",
                              "Calcium", "Phosphorous", "Sodum", "Lysine", "Methionine",
                              "SCFA.", "BA.", "VA.","MCFA.", "CA.", "EA.", "BCFA.", "Water.content....")]


# Transpose genus tab
genus.int_t <- as.data.frame(t(genus.short))

# Select interesting genera from the genus table if needed
# genus.int_t <- genus.int_t[,which(colnames(genus.int_t) %in% c("Lactobacillus", "Megasphaera", 
#                                                    "Prevotella", "Rikenellaceae_RC9_gut_group",
#                                                    "Treponema", "Clostridium_sensu_stricto_1",
#                                                    "UCG-005","Ruminococcus"))]

# Include genus names as a column for merging
genus.int_t$Subject <- rownames(genus.int_t)


# Merge taxonomy table and metdata
corr_tab <- inner_join(metadata_clean, genus.int_t, by=c("Subject"))

# Remove subjects for computation of correlations
corr_tab$Subject <- NULL

# Convert empty entries to NAs
corr_tab <- data.frame(lapply(corr_tab, function(x) ifelse(x == "", NA, x)))

# Get rid of incomplete lines
corr_tab <- corr_tab[complete.cases(corr_tab),]

# Check type of columns (should all be numeric)
sapply(corr_tab, class)

# Convert columns to numeric (in case they are not) - first define function then apply it to df
convert_to_numeric <- function(x) {
  if (is.character(x)) {
    as.numeric(gsub(",", ".", x))
  } else {
    x
  }
}
corr_tab[] <- lapply(corr_tab, convert_to_numeric)

# Check that entries are all numeric
sapply(corr_tab, class)

# This line can be used to remove any unwanted pathway entries from the table
# In this case we only want the actual genera to be shown
corr_tab <- corr_tab[,-which(colnames(corr_tab) %in% c("unclassified.genus", "uncultured", "under_1."))]

# Adjust variable names
colnames(corr_tab)[colnames(corr_tab) == "SCFA."] <- "SCFA %"
colnames(corr_tab)[colnames(corr_tab) == "Water.content...."] <- "Water Content %"
colnames(corr_tab)[colnames(corr_tab) == "BA."] <- "Butyric Acid %"
colnames(corr_tab)[colnames(corr_tab) == "VA."] <- "Valeric Acid %"
colnames(corr_tab)[colnames(corr_tab) == "MCFA."] <- "MCFA %"
colnames(corr_tab)[colnames(corr_tab) == "CA."] <- "Caprohic Acid %"
colnames(corr_tab)[colnames(corr_tab) == "EA."] <- "Enanthic Acid %"
colnames(corr_tab)[colnames(corr_tab) == "BCFA."] <- "BCFA %"

colnames(corr_tab)[colnames(corr_tab) == "Crude.protein"] <- "Crude Protein"
colnames(corr_tab)[colnames(corr_tab) == "Crude.fibre"] <- "Crude Fibre"
colnames(corr_tab)[colnames(corr_tab) == "Crude.fat"] <- "Crude Fat"
colnames(corr_tab)[colnames(corr_tab) == "Sodum"] <- "Sodium"

colnames(corr_tab)[colnames(corr_tab) == "Christensenellaceae_R.7_group"] <- "Christensenellaceae R.7"
colnames(corr_tab)[colnames(corr_tab) == "Rikenellaceae_RC9_gut_group"] <- "Rikenellaceae RC9"
colnames(corr_tab)[colnames(corr_tab) == "Lachnospiraceae_NK4A136_group"] <- "Lachnospiraceae NK4A136"


# Compute correlations
cor_matrix_full <- rcorr(as.matrix(corr_tab), type="pearson")
cor_matrix <- cor_matrix_full$r

# Extract p-values of computed correlations for verification
cor_p <- cor_matrix_full$P


# Prepare the correlation matrix for the chord diagram
# diag(cor_matrix) <- NA #get rid of autocorrelations
# cor_matrix[lower.tri(cor_matrix)] <- NA  # Remove redundant lower half to avoid duplication
# cor_matrix <- as.data.frame(as.table(cor_matrix)) #transform matrix
# cor_matrix <- na.omit(cor_matrix) #get rid of NAs

# Define and apply function to reorder variables based on strength of correlation
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
cormat <- reorder_cormat(cor_matrix)

# Melt matrix for plotting
melted_cormat <- melt(cor_matrix)
head(melted_cormat)

# Make Heatmap
ggheatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#DC3220", high = "#005AB5", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 6, hjust = 1)) +
  coord_fixed() +
  xlab("") +
  ylab("")

ggheatmap

# Save plot
ggsave(ggheatmap, file = "corr_heatmap.png", dpi=300, height = 200, width=370, units = "mm")





