# Example of the use of species identification script (Author: Lynn Goavert, adapted by Emanuele Giacomuzzo)

#!!!! Open this script from the project, otherwise the package "here" doesn't work !!!!!

#This is an example of how to use the species identification script. 
#The data comes from the .... 
#The training dataset is the dataset with inside the monocultures. 
#The mixed cultures data is the videos which I've taken during t0.


#Remove all the parameters present in your environemnt.
rm(list = ls())

#Install packages if they are not already install on your computer.
#library(devtools)
#install.packages("e1071",dependencies = T)
#install.packages("devtools",dependencies = T)
#install_github("pennekampster/bemovi", ref="master")

library(bemovi)
library(e1071)
library("here")
library("tidyverse")

time_points_in_experiment = c("t0")

for (time_point in time_points_in_experiment) {
  
  #### Step 1: define the parameters. ####
  
  #a) folder names and paths
  video.description.folder = "0_video_description/"
  video.description.file = "video_description.txt"
  merged.data.folder = "5_merged_data/"
  monocultures_folder_path = here("40_threshold_analysis", "training", "")
  mixed_cultures_folder_path = here("40_threshold_analysis", time_point, "")
  
  #b) parameters used in the video analysis script
  fps = 25
  nsv = 5
  measured_volume = 34.4
  pixel_to_scale = 4.05
  filter_min_net_disp = 25
  filter_min_duration = 1
  filter_detection_freq = 0.1
  filter_median_step_length = 3
  
  #### Step 2: get the training data. ####
  #To do so, turn the trajectories of the monocultures data into the training data by:
  #a) loading them
  #b) filtering them
  #c) summarising them into individuals data (morph_mvt)
  #d) get rid of the individuals without some columns that could not be calculated
  
  #a) load trajectories of monocultures
  load(paste0(
    monocultures_folder_path,
    merged.data.folder,
    "Master.RData"
  ))
  
  trajectory.data_monocultures = trajectory.data
  rm(trajectory.data)
  
  #b) filter trajectories of monocultures. Filter according to the same parameters you used in the videos analysis script.
  trajectory.data_monocultures.filtered = filter_data(
    trajectory.data_monocultures,
    filter_min_net_disp,
    filter_min_duration,
    filter_detection_freq,
    filter_median_step_length
  )
  
  #c) summarise trajectories of monocultures into individuals data (morph_mvt)  
  morph_mvt = summarize_trajectories(
    data = trajectory.data_monocultures.filtered,
    calculate.median = FALSE,
    write = TRUE,
    to.data = monocultures_folder_path,
    merged.data.folder = merged.data.folder
  ) %>%
    mutate(comment = NULL)
  
  #d) get rid of the individuals without some columns that could not be calculated
  training_data = morph_mvt[complete.cases(morph_mvt), ]
  
  #### Step 3: use the training model to construct the model. ####
  #a) construct the model
  #b) construct confusion matrix
  #c) get species names (???)
  
  #a) construct model
  svm1 = svm(
    factor(species) ~
      mean_grey +
      sd_grey +
      mean_area +
      sd_area +
      mean_perimeter +
      mean_turning +
      sd_turning +
      sd_perimeter +
      mean_major +
      sd_major +
      mean_minor +
      sd_minor +
      mean_ar +
      sd_ar +
      duration +
      max_net  +
      net_disp +
      net_speed +
      gross_disp +
      max_step +
      min_step +
      sd_step +
      sd_gross_speed +
      max_gross_speed +
      min_gross_speed ,
    data = training_data,
    probability = T,
    na.action = na.pass
  )
  
  #b) construct confusion matrix
  confusion.matrix = table(svm1$fitted, training_data$species)
  confusion.matrix.nd = confusion.matrix
  diag(confusion.matrix.nd) = 0
  svm1$confusion = cbind(confusion.matrix,
                         class.error = rowSums(confusion.matrix.nd) / rowSums(confusion.matrix))
  
  print(paste("Confusion matrix of time point", time_point))
  print(svm1$confusion)
  
  #c) get species names
  species.names = unique(trajectory.data_monocultures$species)
  
  #### Step 4: get the predictor data. ####
  #To do so, turn the trajectories of the mixed cultures data (time point) into the predictor data by:
  #a) loading them
  #b) filtering them
  #c) summarising them into individuals data (morph_mvt)
  #d) get rid of the individuals without some columns that could not be calculated

  #a) load mixed cultures trajectories
  load(paste0(
    mixed_cultures_folder_path,
    merged.data.folder,
    "Master.RData"
  ))
  
  trajectory.data_mixed = trajectory.data
  rm(trajectory.data)
  
  #b) filter mixed cultures trajectories (same parameters as in the videos analysis script)
  trajectory.data_mixed.filtered = filter_data(
    trajectory.data_mixed,
    filter_min_net_disp,
    filter_min_duration,
    filter_detection_freq,
    filter_median_step_length
  )
  
  #c) summarise mixed cultures trajectories into individuals data (morph_mvt)
  morph_mvt = summarize_trajectories(
    data = trajectory.data_mixed.filtered,
    calculate.median = FALSE,
    write = TRUE,
    to.data = mixed_cultures_folder_path,
    merged.data.folder = merged.data.folder) %>%
    mutate(comment = NULL)
  
  #d) get rid of the individuals without some columns that could not be calculated
  data.to.predict = morph_mvt[complete.cases(morph_mvt),]
  
  #### Step 5: identify species in the predictor data. ####
  #a) predict species ID of individuals
  #b) add predicted species ID to the dataset
  #c) summarise individuals into populations
  
  #a) predict species ID of individuals
  p.id = predict(object = svm1,
                 data.to.predict,
                 type = "response")
  
  #b) add predicted species ID to the dataset
  data.to.predict$predicted_species = as.character(p.id)
  
  #c) summarise individuals into populations
  pop.data = summarize_populations(
    traj.data = trajectory.data_mixed.filtered,
    sum.data = morph_mvt,
    write = TRUE,
    to.data = mixed_cultures_folder_path,
    merged.data.folder = merged.data.folder,
    video.description.folder = video.description.folder,
    video.description.file = video.description.file,
    total_frames = fps * nsv
  )
  
  #### Step 6: get species densities. ####
  #a) construct a function that creates a matrix with the density of each species in each sample.
  #b) use the function
  
  #a) construct a function that creates a matrix with the density of each species in each sample.
  species.density = function(sample_output,
                             indiv_predicted,
                             species_names,
                             total_frames,
                             mv = measured_volume) {
    samples = unique(indiv_predicted$file)
    
    sp.dens = matrix(data = 0,
                     nrow = nrow(sample_output),
                     ncol = length(species_names))
    
    colnames(sp.dens) = species_names
    
    #For each sample
    for (i in 1:length(samples)) {
      
      #Get all the individuals from that sample
      indiv = subset(indiv_predicted, file == samples[i])
      
      #Get all the species detected in that sample
      spec = unique(indiv$predicted_species)
      
      #For each species identified in the sample
      for (j in 1:length(spec)) {
        
        #Get all the individuals of that species detected in that sample
        all.indiv.sp = subset(indiv, predicted_species == spec[j])
        
        #Calculate the density of that species as the average number of individuals / measured volume
        dens = sum(all.indiv.sp$N_frames) / total_frames / mv
        
        #Add this density to the density of that species for that sample
        sp.dens[which(sample_output$file == as.character(samples[i])), 
                which(species_names == spec[j])] = dens
      }
    }
    
    return(cbind(sample_output, sp.dens))
    
  }
  
  #b) use the function
  output = species.density(
    pop.data,
    data.to.predict,
    species.names,
    total_frames = fps * nsv,
    mv = measured_volume
  )
  
  #### Step 7: save your output. ####
  
  non_species_columns <- output %>%
    select(-all_of(species.names))
  
  species_columns <- output %>%
    select(all_of(species.names)) %>%
    rename_with(~paste0(., "_indiv_per_volume"))
  
  output <- cbind(non_species_columns,
                  species_columns)
  
    
  file_name = paste0("species_ID_", time_point, ".csv")
  
  write.csv(output,
            here("40_threshold_analysis", "species_ID_results", file_name))
  
  rm(output)
  
}