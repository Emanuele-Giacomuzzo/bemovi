rm(list = ls())

#install.packages("e1071",dependencies = T)
#install.packages("devtools",dependencies = T)
#install_github("pennekampster/bemovi", ref="master")

library(devtools)
library(bemovi)
library(e1071)

video.description.folder = "0_video_description/"
video.description.file = "video_description.txt"
merged.data.folder = "5_merged_data/"
to.data.mono = "/media/mendel-himself/ID_058_Ema1/training_002_Lynn_script/Mendel/threshold_10/"
to.data.mixt = "/media/mendel-himself/ID_058_Ema1/training_002_Lynn_script/Mendel/threshold_10/"
to.data = to.data.mixt

fps = 25 
nsv = 20 
measured_volume = 34.4 
pixel_to_scale = 4.05 
filter_min_net_disp = 20 
filter_min_duration = 0.5 
filter_detection_freq = 0.1 
filter_median_step_length = 3

load(paste0(to.data.mono, 
            merged.data.folder, 
            "Master.RData"))

trajectory.data.filtered = filter_data(trajectory.data, 
                                       filter_min_net_disp, 
                                       filter_min_duration, 
                                       filter_detection_freq, 
                                       filter_median_step_length)

morph_mvt = summarize_trajectories(trajectory.data.filtered, 
                                   write = T,
                                   to.data, 
                                   calculate.median = F, 
                                   merged.data.folder)

training_data = morph_mvt[complete.cases(morph_mvt), ]

svm1 = svm(factor(Community) ~   
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
           na.action = na.pass)

confusion.matrix = table(svm1$fitted,training_data$Community)
confusion.matrix.nd = confusion.matrix
diag(confusion.matrix.nd) = 0
svm1$confusion = cbind(confusion.matrix,class.error=rowSums(confusion.matrix.nd)/rowSums(confusion.matrix))
svm1$confusion

species.names = c("P","C")
load(paste0(to.data.mixt, merged.data.folder, "Master.RData"))
trajectory.data.filtered = filter_data(trajectory.data, filter_min_net_disp, filter_min_duration, filter_detection_freq, filter_median_step_length)
morph_mvt = summarize_trajectories(trajectory.data.filtered, write = T, to.data, calculate.median=F, merged.data.folder)[,which(colnames(morph_mvt)!="Col_manual")]

data.to.predict = morph_mvt[complete.cases(morph_mvt),]

p.id = predict(svm1, data.to.predict, type="response")     
data.to.predict$predicted_species = as.character(p.id)

pop.data = summarize_populations(trajectory.data.filtered, morph_mvt, write=T, to.data, merged.data.folder, video.description.folder, video.description.file,total_frame = fps * nsv)

species.density = function(sample_output,indiv_predicted,species_names,total_frames,mv = measured_volume){
  
  samples = unique(indiv_predicted$file)
  
  sp.dens = matrix(0,nrow(sample_output),length(species_names))
  colnames(sp.dens) = species_names
  
  for(i in 1:length(samples)){
    
    indiv = subset(indiv_predicted,file == samples[i])
    
    spec = unique(indiv$predicted_species)
    
    for(j in 1:length(spec)){
      
      all.indiv.sp = subset(indiv,predicted_species == spec[j])
      
      
      dens = sum(all.indiv.sp$N_frames)/total_frames/mv
      sp.dens[which(sample_output$file == as.character(samples[i])) ,which(species_names == spec[j])] = dens}
  }
  
  return(cbind(sample_output,sp.dens)) 

  }

output = species.density(pop.data,data.to.predict,species.names,total_frames = fps*nsv, mv = measured_volume)
