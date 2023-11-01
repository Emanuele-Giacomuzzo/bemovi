#Originally written by Lynn Govaert, then modified by Emanuele Giacomuzzo. 
#See https://github.com/Emanuele-Giacomuzzo/bemovi/blob/master/species_id_example.html for an explanation with example of this script.
#The script is not commented to be kept clean.

rm(list = ls())
#library(devtools)
#install.packages("e1071",dependencies = T)
#install.packages("devtools",dependencies = T)
#install_github("pennekampster/bemovi", ref="master")

library(bemovi)
library(e1071)
library("here")
library("tidyverse")

time_points_in_experiment = c("t0", "t1", "t2", "t3", "t4")
threshold_levels = c(13, 40)

for (threshold in threshold_levels) {
  for (time_point in time_points_in_experiment) {
  
    #### Step 1: define the parameters. ####
    
      video.description.folder = "0_video_description/"
    video.description.file = "video_description.txt"
    merged_data_folder = "5_merged_data/"
    training_folder_path = here(paste0(threshold, "_threshold_analysis"), "training", "")
    time_point_folder_path = here(paste0(threshold, "_threshold_analysis"), time_point, "") 
    
    fps = 25
    seconds_per_video = 5
    total_frames = fps * seconds_per_video
    measured_volume = 34.4 #volume measured for each video
    pixel_to_scale = 4.05
    filter_min_net_disp = 25 
    filter_min_duration = 1
    filter_detection_freq = 0.1
    filter_median_step_length = 3
    
    #### Step 2: get the training data. ####
    
    load(paste0(training_folder_path,
                merged_data_folder,
                "Master.RData")) 
    
    trajectory_data_training = trajectory.data
    rm(trajectory.data)
    trajectory_data_training[1]
    
    trajectory_data_training_filtered = filter_data(
      raw_data = trajectory_data_training,
      net_filter = filter_min_net_disp,
      duration_filter = filter_min_duration,
      detect_filter = filter_detection_freq,
      median_step_filter = filter_median_step_length
    )
    
    morph_mvt_training = summarize_trajectories(
      data = trajectory_data_training_filtered,
      calculate.median = FALSE,
      write = TRUE,
      to.data = training_folder_path,
      merged.data.folder = merged_data_folder
    ) %>%
      mutate(comment = NULL)
    
    training_data = morph_mvt_training[complete.cases(morph_mvt_training),]
    
    #### Step 3: get the predictor data. ####
    
    load(paste0(
      time_point_folder_path,
      merged_data_folder,
      "Master.RData"
    ))
    trajectory_data_time_point = trajectory.data
    rm(trajectory.data)
    
    trajectory_data_time_point_filtered = filter_data(
      raw_data = trajectory_data_time_point,
      net_filter = filter_min_net_disp,
      duration_filter = filter_min_duration,
      detect_filter = filter_detection_freq,
      median_step_filter = filter_median_step_length
    )
    
    morph_mvt_time_point = summarize_trajectories(
      data = trajectory_data_time_point_filtered,
      calculate.median = FALSE,
      write = TRUE,
      to.data = time_point_folder_path,
      merged.data.folder = merged_data_folder
    ) %>%
      mutate(comment = NULL)
    
    data_to_predict = morph_mvt_time_point[complete.cases(morph_mvt_time_point), ]
    
    #### Step 4: use the training model to construct the model. ####
    
    model = svm(
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
    
    confusion_matrix = table(model$fitted, 
                             training_data$species)
    
    confusion_matrix_with_diagonal_zeros = confusion_matrix
    diag(confusion_matrix_with_diagonal_zeros) = 0
    
    model$confusion = cbind(confusion_matrix,
                            class_error = rowSums(confusion_matrix_with_diagonal_zeros) / rowSums(confusion_matrix))
    
    print(paste("Confusion matrix of time point", time_point, "with threshold", threshold))
    
    model$confusion %>%
      as.data.frame() %>%
      mutate(class_error_percentage = class_error * 100,
             class_error = NULL)
    
    #### Step 5: identify species in the predictor data. ####
  
    predicted_IDs = predict(object = model,
                            data_to_predict,
                            type = "response")
    predicted_IDs[1:10]
    
    if(nrow(data_to_predict) != length(predicted_IDs)) stop()
    
    data_to_predict$predicted_species = as.character(predicted_IDs)
    
    #### Step 6: get species densities. ####
    
    species_names = unique(trajectory_data_training$species)
    
    species.density = function(pop.data,
                               data_to_predict,
                               species_names,
                               total_frames,
                               measured_volume) {
      
      samples = unique(data_to_predict$file)
      
      species_density_matrix = matrix(
        data = 0,
        nrow = nrow(pop.data),
        ncol = length(species_names)
      )
      
      colnames(species_density_matrix) = species_names
      
      for (sample_i in 1:length(samples)) {
        
        individuals_in_sample = subset(data_to_predict, 
                                       file == samples[sample_i])
        
        species_in_sample = unique(individuals_in_sample$predicted_species)
        
        n_species_in_sample = length(species_in_sample)
        
        for (species_i in 1:n_species_in_sample) {
          
          indiv_of_species_i = subset(individuals_in_sample, 
                                      predicted_species == species_in_sample[species_i])
          
          density_species_i = sum(indiv_of_species_i$N_frames) / total_frames / measured_volume
          
          species_density_matrix[
            which(pop.data$file == as.character(samples[sample_i])),
            which(species_names == species_in_sample[species_i])
          ] = 
            density_species_i
        }
      }
      
      return(cbind(pop.data, 
                   species_density_matrix))
      
    }
    
    pop.data = summarize_populations(
      traj.data = trajectory_data_time_point_filtered,
      sum.data = morph_mvt_time_point,
      write = TRUE,
      to.data = time_point_folder_path,
      merged.data.folder = merged_data_folder,
      video.description.folder = video.description.folder,
      video.description.file = video.description.file,
      total_frames = total_frames
    )
    
    output = species.density(
      pop.data = pop.data,
      data_to_predict = data_to_predict,
      species_names = species_names,
      total_frames = total_frames,
      measured_volume = measured_volume
    )
    
    #### Step 7: save your output. ####
    
    non_species_columns <- output %>%
      select(-all_of(species_names))
    
    species_columns <- output %>%
      select(all_of(species_names)) %>%
      rename_with( ~ paste0(., "_indiv_per_volume"))
    
    output <- cbind(non_species_columns,
                    species_columns)
    
    file_name = paste0(time_point, ".csv")
    
    write.csv(output,
              here("results", paste0("species_ID_", threshold, "_threshold"), file_name))
    
    rm(output)
 
  }
}
