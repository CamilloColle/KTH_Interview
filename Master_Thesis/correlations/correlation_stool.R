library(circlize)
library(dplyr)
library(qiime2R)
library(RColorBrewer)
library(gridBase)
library(ComplexHeatmap)
library(tidyr)
library(Hmisc)

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
metadata_clean <- metadata[,c("Subject", "SCFA.", "BA.", "VA.",
                              "MCFA.", "CA.", "EA.", "BCFA.", "Water.content....")]

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

colnames(corr_tab)[colnames(corr_tab) == "Christensenellaceae_R.7_group"] <- "Christensenellaceae R.7"
colnames(corr_tab)[colnames(corr_tab) == "Rikenellaceae_RC9_gut_group"] <- "Rikenellaceae RC9"
colnames(corr_tab)[colnames(corr_tab) == "Lachnospiraceae_NK4A136_group"] <- "Lachnospiraceae NK4A136"


# Compute correlations
cor_matrix_full <- rcorr(as.matrix(corr_tab), type="pearson")
cor_matrix <- cor_matrix_full$r

# Extract p-values of computed correlations for verification
cor_p <- cor_matrix_full$P

# Extract most relevant correlations
# For this we first create a subsets of the correlation matrix with only the 
# genus information on the rows and the metadata on the columns
corr_sub <- cor_matrix[c(1:8),-c(1:8)]


# Convert the dataframe to a long format with row and column identifiers
long_df <- as.data.frame(corr_sub) %>%
  mutate(row_number = row_number()) %>%
  pivot_longer(
    cols = -row_number,
    names_to = "column",
    values_to = "value"
  )

# Calculate absolute values and find the maximum per column
top_values_per_column <- long_df %>%
  mutate(abs_value = abs(value)) %>%
  group_by(column) %>%
  dplyr::summarize(
    max_abs_value = max(abs_value),
    row_at_max = row_number[which.max(abs_value)],
    value_at_max = value[which.max(abs_value)],
    .groups = "drop"
  )

# Sort and select the top n unique columns
top_columns <- top_values_per_column %>%
  arrange(desc(max_abs_value)) %>%
  slice_head(n = 8)

# Check the results
print(top_columns)

# Define vector for selected variables
sub_vec <- c(rownames(corr_sub), top_columns$column)

# Use the vector to subset the correlation and p-value matrices
cor_matrix <- cor_matrix[sub_vec, sub_vec]
cor_p <- cor_p[sub_vec, sub_vec]

# Prepare the correlation matrix for the chord diagram
diag(cor_matrix) <- NA #get rid of autocorrelations
cor_matrix[lower.tri(cor_matrix)] <- NA  # Remove redundant lower half to avoid duplication
cor_matrix <- as.data.frame(as.table(cor_matrix)) #transform matrix
cor_matrix <- na.omit(cor_matrix) #get rid of NAs

# Function to determine the color based on correlation
get_color <- function(value) {
  if (value > 0) {
    # Blue for positive
    color <- grDevices::adjustcolor("#005AB5", alpha.f = abs(value))
  } else {
    # Red for negative
    color <- grDevices::adjustcolor("#DC3220", alpha.f = abs(value))
  }
  return(color)
}

# Apply color function to correlation matrix
cor_matrix$color <- apply(cor_matrix, 1, function(x) get_color(as.numeric(x[3])))

# Define sectors for order
sectors <- unique(c(cor_matrix$Var1, cor_matrix$Var2))

# Define the groups and the relative number of entries
df.groups <- c(rep("Stool", 8), rep("Taxa", 8))
names(df.groups) <- sectors

# Check defined groups/names association
df.groups

# Define colors for the groups. Note that this must be done manually, meaning if we have a "Stool" group
# with 8 entries and a "Taxa" group with 8 entries, we will need to create a colors vector which has the color
# for each of the groups repeated for the same number of times, i.e. 8 and 8 in this case
sec_cols <- c(rep("#C55A11", 8), rep("#660066", 8))
names(sec_cols) <- sectors

# Associate group names with colors
group_cols <- c(Stool = "#C55A11", Taxa = "#660066")


# Create the chord diagram

# This line is used to initialise the file that will contain the plot
# Comment out if you want the plot to appear in the plot window on the right instead
png("chord_stool.png", units="in", width=8, height=5, res=300)

# Reset the plot space in case something was there already
circos.clear()

# Define font size (cex) and margins
par(cex = 0.78, mar = c(0, 0, 0, 8))

# Start building the chord diagram
chordDiagramFromDataFrame(cor_matrix %>% select(Var1, Var2), 
                          order = sectors,
                          grid.col = sec_cols, 
                          col = cor_matrix$color,
                          group = df.groups,
                          annotationTrack = c("grid"), # here you can add 'name' next to 'grid' to also display the labels of the entries
                          #annotationTrackHeight = convert_height(c(3, 2), "mm"),
                          preAllocateTracks = 1
)

# Add labels
for(si in 1:length(get.all.sector.index())) {
  xlim = get.cell.meta.data("xlim", sector.index = get.all.sector.index()[si], track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = get.all.sector.index()[si], track.index = 2)
  circos.text(mean(xlim), mean(ylim), si, sector.index = get.all.sector.index()[si], 
              col = "white", track.index = 2, cex = .7)
}


# Add sector labels
for (i in unique(df.groups)) {
  highlight.sector(names(df.groups[df.groups==i]), track.index = 1, 
                   col = sec_cols[df.groups==i], text = i, 
                   text.vjust = -1, padding = c(-0.95, 0, -.4, 0))
}


# LEGENDS
# Sectors
lgd_stool = Legend(title = "Stool", 
                   labels =  names(df.groups[df.groups == "Stool"]),
                   type = "points", 
                   pch = as.character(1:8), #change this range according to the number of entries in each sector (in this case 8)
                   legend_gp = gpar(col = "white", cex = 0.7),
                   background = group_cols[names(group_cols) == "Stool"])

lgd_taxa = Legend(title = "Taxa", 
                  labels =  names(df.groups[df.groups == "Taxa"]),
                  type = "points", 
                  pch = as.character(9:16), #change this range according to the number of entries in each sector (in this case 8)
                  legend_gp = gpar(col = "white", cex = 0.7),
                  background = group_cols[names(group_cols) == "Taxa"])


# Links
col_fun = colorRamp2(c(-1, 0, 1), c("#DC3220", "white", "#005AB5"))
lgd_links = Legend(at = c(-1, -0.5, 0, 0.5, 1), col_fun = col_fun, 
                   title_position = "topleft", title = "Links")

# Join and draw legends, the draw function is repeated twice because a bug 
# sometimes causes the legens to be drawn incompletely
lgd_list_vertical = packLegend(lgd_stool, lgd_taxa, lgd_links)
lgd_list_vertical

draw(lgd_list_vertical, x = unit(0.88, 'npc'), y = unit(0.5, "npc"), just = "centre")
draw(lgd_list_vertical, x = unit(0.88, 'npc'), y = unit(0.5, "npc"), just = "centre")

# Save and close the file
dev.off()










