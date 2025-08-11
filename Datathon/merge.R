library(dplyr)

full_df = read.csv("final_raw_count.csv")
pop = read.csv("Data/population_by_municipality.csv")

#drop unused columns
full_df = full_df[,-c(1)]
pop = pop[,-c(1,5)]

#merge dataframes
final_df_with_pop = inner_join(full_df, pop, by = c("postal_code" = "Post.code"))
names(final_df_with_pop)[names(final_df_with_pop) == 'CD_REFNIS'] <- 'municipality'
names(final_df_with_pop)[names(final_df_with_pop) == 'TOTAL'] <- 'population'


#group postal codes
final_df <- final_df_with_pop[,-8] # Remove the municipality column


final_df <- final_df %>%
  group_by(postal_code) %>%
  summarise(
    avg_intervention_time = mean(V1),
    cardiac_intervention_count = mean(cardiac_intervention_count),
    DOA_cardiac_count = mean(DOA_cardiac_count),
    DOA_all_count = mean(DOA_all_count), 
    interventions_count = mean(interventions_count),
    aed_count = mean(aed_count),
    population = mean(population)
  )  



#write to csv
write.csv(final_df, file = "final_rawcounts_merged.csv")

