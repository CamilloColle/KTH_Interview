library(Tax4Fun2)
library(qiime2R)
library(dplyr)
library(tidyr)
library(circlize)
library(RColorBrewer)
library(gridBase)
library(ComplexHeatmap)
library(Hmisc)

# SETUP
setwd("C:/Users/ccoll/Desktop/Thesis/thesis_master/f_prediction")

# INPUT FILES
func <- as.data.frame((read.table("picrust_table.tsv"))) #picrust pathway abundance output
genus_t <- read.table('genus_t.tsv')

# Clean the input
colnames(func) <- func[1,]
func <- func[-1,]
rownames(func) <- func[,1]
func <- func[,-1]
func <- t(func)

# Select interesting genera from the genus table
genus_t <- genus_t[,which(colnames(genus_t) %in% c("Lactobacillus", "Megasphaera", 
                                                       "Prevotella", "Rikenellaceae_RC9_gut_group",
                                                   "Treponema", "Streptococcus",
                                                   "UCG.005","Ruminococcus"))]

# Merge the two input tables
corr_tab <- merge(genus_t, func, by = 'row.names')

###COMPUTE CORRELATIONS

# Remove subjects for computation of correlations
rownames(corr_tab) <- corr_tab$Row.names
corr_tab$Row.names <- NULL
corr_tab$Subject <- NULL

# Convert empty entries to NAs
corr_tab <- data.frame(lapply(corr_tab, function(x) ifelse(x == "", NA, x)))

# Get rid of incomplete lines
corr_tab <- corr_tab[complete.cases(corr_tab),]

# Check type of columns
sapply(corr_tab, class)

# Convert columns to numeric (first define the function then apply it)
convert_to_numeric <- function(x) {
  if (is.character(x)) {
    as.numeric(gsub(",", ".", x))
  } else {
    x
  }
}
corr_tab[] <- lapply(corr_tab, convert_to_numeric)


# This line can be used to remove any unwanted pathway entries from the table
# In this case two entries were removed because the similarity in pathways would have resulted in redundant information
corr_tab <- corr_tab[,-which(colnames(corr_tab) %in% c("EC.4.2.3.4", "EC.2.5.1.19"))]

# Rename entries for better visualisations
colnames(corr_tab)[colnames(corr_tab) == "Christensenellaceae_R.7_group"] <- "Christensenellaceae R.7"
colnames(corr_tab)[colnames(corr_tab) == "Rikenellaceae_RC9_gut_group"] <- "Rikenellaceae RC9"

colnames(corr_tab)[colnames(corr_tab) == "EC.1.1.3.21"] <- "G3P oxidation"
colnames(corr_tab)[colnames(corr_tab) == "EC.1.1.1.219"] <- "leucodelphinidin biosynthesis"
colnames(corr_tab)[colnames(corr_tab) == "EC.2.4.2.6"] <- "nucleoside deoxyribosyltransf."
colnames(corr_tab)[colnames(corr_tab) == "EC.3.4.11.5"] <- "proline hydrolysis"
colnames(corr_tab)[colnames(corr_tab) == "EC.2.6.1.55"] <- "taurine degradation"
colnames(corr_tab)[colnames(corr_tab) == "EC.4.2.3.5"] <- "chorismate biosynthesis"
colnames(corr_tab)[colnames(corr_tab) == "EC.2.5.1.2"] <- "thiamine degradation"
colnames(corr_tab)[colnames(corr_tab) == "EC.2.7.7.8"] <- "polynucleotide phosphorylation"


# Compute correlations
cor_matrix_full <- rcorr(as.matrix(corr_tab), type="pearson")
cor_matrix <- cor_matrix_full$r

# Extract p-values of computed correlations for verification
cor_p <- cor_matrix_full$P

# Extract most relevant correlations
# For this we first create two subsets of the correlation matrix, 
# one with only the genus information on the rows and the pathways on the columns and one the other way around
corr_sub1 <- cor_matrix[c(1:8),-c(1:8)]
corr_sub2 <- cor_matrix[-c(1:8),c(1:8)]

# Then define a function that sorts the entries based on correlation and takes the top n (in this case 8)
top_corr <- function(corr_sub, n = 8){
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
    slice_head(n = n)
  
  # Print the results
  return(top_columns$column)
}


# Make vector for selected variables by applying the above function
sub_vec <- c(top_corr(corr_sub1), top_corr(corr_sub2))

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
df.groups <- c(rep("Function", 8), rep("Taxa", 8))
names(df.groups) <- sectors

# Check defined groups/names association
df.groups

# Define colors for the groups. Note that this must be done manually, meaning if we have a "Function" group
# with 8 entries and a "Taxa" group with 8 entries, we will need to create a colors vector which has the color
# for each of the groups repeated for the same number of times, i.e. 8 and 8 in this case
sec_cols <- c(rep("#006666", 8), rep("#660066", 8))
names(sec_cols) <- sectors

# Associate group names with colors
group_cols <- c(Function = "#006666", Taxa = "#660066")


# Create the chord diagram

# This line is used to initialise the file that will contain the plot
# Comment out if you want the plot to appear in the plot window on the right instead
png("func_picrust.png", units="in", width=8, height=5, res=300)

# Reset the plot space in case something was there already
circos.clear()

# Define font size (cex) and margins
par(cex = 0.78, mar = c(0, 0, 0, 10))

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

# Legend sectors
lgd_function = Legend(title = "Function", 
                      labels =  names(df.groups[df.groups == "Function"]),
                      type = "points", 
                      pch = as.character(1:8), #change this range according to the number of entries in each sector (in this case 8)
                      legend_gp = gpar(col = "white", cex = 0.7),
                      background = group_cols[names(group_cols) == "Function"])

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
lgd_list_vertical = packLegend(lgd_function, lgd_taxa, lgd_links)

draw(lgd_list_vertical, x = unit(0.86, 'npc'), y = unit(0.5, "npc"), just = "centre")
draw(lgd_list_vertical, x = unit(0.86, 'npc'), y = unit(0.5, "npc"), just = "centre")

# Save and close the file
dev.off()








