rm(list = ls())
library(here)

threshold_level = "13_threshold_analysis"
for (time_point in c("t0", "t1", "t2", "t3", "t4")){
  
  load(here(threshold_level, time_point, "5_merged_data", "Population_Data.RData"))
  load(here(threshold_level, time_point, "5_merged_data", "Morph_mvt.RData"))
  
  write.csv(pop_output,
            file = here("results", "populations", paste0(time_point, ".csv")),
            row.names = F)
  
  write.csv(morph_mvt,
            file = here("results", "individuals", paste0(time_point, ".csv")),
            row.names = F)
  
  rm(pop_output)
  rm(morph_mvt)
}
