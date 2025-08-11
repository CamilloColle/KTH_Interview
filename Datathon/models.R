data <- read.csv("final_rawcounts_merged.csv", header = T)
data <- data[,-1]
rownames(data) <- data$postal_code

data_sub <- data[,-1]

# Load necessary libraries
library(keras)
library(mlbench)
library(dplyr)
library(tidyr)
library(magrittr)
library(neuralnet)
library(tidyverse)
library(tidymodels)

# Generate all possible combinations of interactions among variables
interactions <- data_sub %>%
  colnames() %>%
  combn(2, simplify = FALSE) %>%
  lapply(function(cols) paste(cols, collapse = "*")) %>%
  unlist()

# Add main effects
variables <- colnames(data_sub)[-1]
all_vars <- c(variables, interactions)

# Generate all possible combinations of models
model_formulas <- lapply(seq_along(all_vars), function(i) {
  combn(all_vars, i, simplify = FALSE) %>%
    lapply(function(vars) paste(vars, collapse = "+"))
  }) %>%
    unlist(recursive = FALSE)

# Fit models and evaluate statistical significance
results <- lapply(model_formulas, function(formula) {
  model <- lm(as.formula(paste("mean_DOA_cardiac ~", formula)), data = data_sub)
  anova(model)
})

# Print results
names(results) <- model_formulas
lapply(results, print)


###ANN

set.seed(1)
data_rows <- floor(0.80 * nrow(data_sub))
train_indices <- sample(c(1:nrow(data_sub)), data_rows)
train_data <- data_sub[train_indices,]
test_data <- data_sub[-train_indices,]

model = neuralnet(
  DOA_cardiac_count~avg_intervention_time+cardiac_intervention_count+
    DOA_all_count+interventions_count+aed_count+population,
  data=train_data,
  hidden=c(4,2),
  linear.output = FALSE
)













