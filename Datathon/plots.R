library(tidyverse) # For data manipulation and plotting
library(sf)        # For handling spatial data
library(ggplot2)   # For plotting

# Load the dataset
data <- read.csv("final_rawcounts_merged.csv")

# Calculate the ratio
# Note: Ensure to handle division by zero or missing values appropriately
data$ratio <- with(data, aed_count / DOA_cardiac_count)

# Load Belgium postal code map
belgium_map <- st_read("4326/postaldistricts.shp")

# Ensure the postal code column in the map data matches the format in your dataset

# Merge the ratio data with the map data based on postal codes
# Adjust the column names accordingly if they differ
merged_data <- belgium_map %>% 
  left_join(data, by = c("postal_code_column_in_map" = "postal_code_column_in_dataset"))

# Plot the map
ggplot(data = merged_data) +
  geom_sf(aes(fill = ratio), color = NA) + 
  scale_fill_gradient(low = "lightblue", high = "purple", na.value = "white", 
                      name = "AED to DOA Cardiac Ratio") +
  labs(title = "Belgium Postal Codes by AED to DOA Cardiac Ratio",
       subtitle = "Darker colors indicate higher ratios") +
  theme_minimal()
