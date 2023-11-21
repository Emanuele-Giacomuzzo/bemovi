rm(list=ls())
setwd('/media/mendel-himself/ID_061_Ema2/training_004_with_Felix/t1')

# library(devtools)
# install_github("femoerman/bemovi", ref="master")
library(bemovi)
library(parallel)
library(doParallel)
library(foreach)

memory.alloc <- 240000 #-needs_to_be_specified
memory.per.identifier <- 40000 #-needs_to_be_specified
memory.per.linker <- 5000 #-needs_to_be_specified
memory.per.overlay <- 60000 #-needs_to_be_specified

tools.path <- "/home/mendel-himself/bemovi_tools/" #-needs_to_be_specified
to.particlelinker <- tools.path

to.data <- paste(getwd(),"/",sep="")
video.description.folder <- "0_video_description/"
video.description.file <- "video_description.txt"
raw.video.folder <- "1_raw/"
raw.avi.folder <- "1a_raw_avi/"
metadata.folder <- "1b_raw_meta/"
particle.data.folder <- "2_particle_data/"
trajectory.data.folder <- "3_trajectory_data/"
temp.overlay.folder <- "4a_temp_overlays/"
overlay.folder <- "4_overlays/"
merged.data.folder <- "5_merged_data/"
ijmacs.folder <- "ijmacs/"


fps <- 25 #-needs_to_be_specified
total_frames <- 125 #-needs_to_be_specified
width=2048 #-needs_to_be_specified
height=2048 #-needs_to_be_specified
measured_volume <- 34.4 # for Leica M205 C with 1.6 fold magnification, sample height 0.5 mm and Hamamatsu Orca Flash 4
#measured_volume <- 14.9 # for Nikon SMZ1500 with 2 fold magnification, sample height 0.5 mm and Canon 5D Mark III

pixel_to_scale <- 4.05 # for Leica M205 C with 1.6 fold magnification, sample height 0.5 mm and Hamamatsu Orca Flash 4
#pixel_to_scale <- 3.79 # for Nikon SMZ1500 with 2 fold magnification, sample height 0.5 mm and Canon 5D Mark III

video.format <- "cxd" #-needs_to_be_specified

difference.lag <- 10
thresholds <- c(13,255) # don't change the second value

particle_min_size <- 10
particle_max_size <- 1000
trajectory_link_range <- 3
trajectory_displacement <- 16
filter_min_net_disp <- 25
filter_min_duration <- 1
filter_detection_freq <- 0.1
filter_median_step_length <- 3

check_tools_folder(tools.path)

system(paste0("chmod a+x ", tools.path, "bftools/bf.sh"))
system(paste0("chmod a+x ", tools.path, "bftools/bfconvert"))
system(paste0("chmod a+x ", tools.path, "bftools/showinf"))

convert_to_avi(to.data, raw.video.folder, raw.avi.folder, metadata.folder, tools.path, fps, video.format)

# check_video_file_names(to.data,raw.avi.folder,video.description.folder,video.description.file)
# check_threshold_values(to.data, raw.avi.folder, ijmacs.folder, 2, difference.lag, thresholds, tools.path,  memory.alloc)

locate_and_measure_particles(to.data, raw.avi.folder, particle.data.folder, difference.lag, min_size = particle_min_size, 
                             max_size = particle_max_size, thresholds=thresholds, tools.path, 
                             memory=memory.alloc, memory.per.identifier=memory.per.identifier, max.cores=detectCores()-1)

link_particles(to.data, particle.data.folder, trajectory.data.folder, linkrange = trajectory_link_range, disp = trajectory_displacement, 
               start_vid = 1, memory = memory.alloc, memory_per_linkerProcess = memory.per.linker, raw.avi.folder, max.cores=detectCores()-1, max_time = 1)

merge_data(to.data, particle.data.folder, trajectory.data.folder, video.description.folder, video.description.file, merged.data.folder)

load(paste0(to.data, merged.data.folder, "Master.RData"))

trajectory.data.filtered <- filter_data(trajectory.data, filter_min_net_disp, filter_min_duration, filter_detection_freq, filter_median_step_length)

morph_mvt <- summarize_trajectories(trajectory.data.filtered, calculate.median=F, write = T, to.data, merged.data.folder)

summarize_populations(trajectory.data.filtered, morph_mvt, write=T, to.data, merged.data.folder, video.description.folder, video.description.file, total_frames)

create.subtitle.overlays(to.data, traj.data=trajectory.data.filtered, raw.video.folder, raw.avi.folder, temp.overlay.folder, overlay.folder, fps,
                         vid.length=total_frames/fps, width, height, tools.path = tools.path, overlay.type="number", video.format)

create_overlays(traj.data = trajectory.data.filtered, to.data = to.data, merged.data.folder = merged.data.folder, raw.video.folder = raw.avi.folder, temp.overlay.folder = "4a_temp_overlays_old/",
                overlay.folder ="4_overlays_old/", width = width, height = height, difference.lag = difference.lag, type = "traj", predict_spec = F, contrast.enhancement = 1, 
                IJ.path = "/home/mendel-himself/bemovi_tools", memory = memory.alloc, max.cores = detectCores()-1, memory.per.overlay = memory.per.overlay)
