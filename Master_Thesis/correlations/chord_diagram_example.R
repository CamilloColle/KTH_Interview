
library(circlize)
library(RColorBrewer)


# Prepare the correlation matrix for the chord diagram
diag(cor_matrix) <- NA
cor_matrix[lower.tri(cor_matrix)] <- NA  # Remove redundant lower half to avoid duplication
cor_matrix <- as.data.frame(as.table(cor_matrix))
cor_matrix <- na.omit(cor_matrix)

# Function to determine the color based on correlation
get_color <- function(value) {
  if (value > 0) {
    color <- rgb(0, 0, 1, abs(value))  # Blue for positive
  } else {
    color <- rgb(1, 0, 0, abs(value))  # Red for negative
  }
  return(color)
}

cor_matrix$color <- apply(cor_matrix, 1, function(x) get_color(as.numeric(x[3])))

cor_matrix$Var1 <- gsub("Water.content....", "WC",cor_matrix$Var1)
cor_matrix$Var2 <- gsub("Water.content....", "WC",cor_matrix$Var2)

cor_matrix$Var1 <- gsub("SCFA.", "SCFA%",cor_matrix$Var1)
cor_matrix$Var2 <- gsub("SCFA.", "SCFA%",cor_matrix$Var2)

#define sectors
sectors <- unique(c(cor_matrix$Var1, cor_matrix$Var2))

# Define the groups and their colors
df.groups <- c(rep("Stool", 2), rep("Diet", 3), rep("Taxa", 8))
names(df.groups) <- sectors
df.groups

sec_cols <- c(rep("#C55A11", 2), rep("#548235", 3), rep("#660066", 8))
names(sec_cols) <- sectors

group_cols <- c(Stool = "#C55A11", Diet = "#548235", Taxa = "#660066")


# Create the chord diagram
circos.clear()
chordDiagramFromDataFrame(cor_matrix %>% select(Var1, Var2), 
                          order = sectors,
                          grid.col = sec_cols, 
                          col = cor_matrix$color,
                          group = df.groups,
                          annotationTrack = c("grid", "axis"),
                          preAllocateTracks = 1
                          )

#add labels
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), mean(ylim), substr(si, 1, 10), sector.index = si, 
              col = "white", track.index = 2, cex = .7)
}


#add sector labels
for (i in unique(df.groups)) {
  highlight.sector(names(df.groups[df.groups==i]), track.index = 1, 
                   col = sec_cols[df.groups==i], text = i, 
                   text.vjust = -1, padding = c(-1, 0, -.5, 0))
}

