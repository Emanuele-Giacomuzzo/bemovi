####################################################################
# R script for analysing video files with BEMOVI (www.bemovi.info)
#
# Felix Moerman
# Computer = Mendel
# 12.01.2022 adapted by
# Samuel Huerlemann, Emanuele Giacomuzzo
######################################################################

# library(devtools)
# install_github("femoerman/bemovi", ref="master")
library(bemovi)
library(parallel)
library(doParallel)
library(foreach)

rm(list = ls())

memory.alloc <- 240000 # Define memory to be allocated
memory.per.identifier <- 40000 # Define memory to be allocated
memory.per.linker <- 5000 # Define memory to be allocated
memory.per.overlay <- 60000 # Define memory to be allocated
tools.path <-
  "/home/mendel-himself/bemovi_tools/" # set paths to tools folder and particle linker - needs to be specified
to.particlelinker <- tools.path
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
maximum_time_allotted_time_per_video <- 24

######################################################################
# VIDEO PARAMETERS
fps <-
  25   # video frame rate (in frames per second) - needs to be specified
total_frames <-
  125 # length of video (in frames) - needs to be specified
width = 2048 #Width of the videos in pixels - needs to be specified
height = 2048 #-needs_to_be_specified
measured_volume <-
  34.4 # measued µL during sampling for Leica M205 C with 1.6 fold magnification, sample height 0.5 mm and Hamamatsu Orca Flash 4 - needs to be specified
#measured_volume <- 14.9 # measued µL during sampling for Nikon SMZ1500 with 2 fold magnification, sample height 0.5 mm and Canon 5D Mark III
pixel_to_scale <-
  4.05 # µm of a pixel for Leica M205 C with 1.6 fold magnification, sample height 0.5 mm and Hamamatsu Orca Flash 4
#pixel_to_scale <- 3.79 # µm of a pixel for Nikon SMZ1500 with 2 fold magnification, sample height 0.5 mm and Canon 5D Mark III
video.format <-
  "cxd" # format of the videos taken  (one of "avi","cxd","mov","tiff") - needs to be specified
difference.lag <- 50 #not sure
min_threshold = 13 #40 for Ble 13 for others
thresholds <-
  c(min_threshold, 255) # threshold for what is considered background and what is not in ImageJ - the first number should be adjusted, the second should not


# MORE PARAMETERS (USUALLY NOT CHANGED)
######################################################################
# FILTERING PARAMETERS
# optimized for Perfex Pro 10 stereomicrocope with Perfex SC38800 (IDS UI-3880LE-M-GL) camera
# tested stereomicroscopes: Perfex Pro 10, Nikon SMZ1500, Leica M205 C
# tested cameras: Perfex SC38800, Canon 5D Mark III, Hamamatsu Orca Flash 4
# tested species: Tet, Col, Pau, Pca, Eug, Chi, Ble, Ceph, Lox, Spi
particle_min_size <- 10 #Smallest particle filtered (area in pixel)
particle_max_size <- 6000 #Largest particle filtered (area in pixel)
trajectory_link_range <-
  3 # number of adjacent frames to be considered for linking particles
trajectory_displacement <-
  16 # maximum distance a particle can move between two frames
filter_min_net_disp <-
  25 # these values are in the units defined by the parameters above: fps (seconds), measured_volume (microliters) and pixel_to_scale (micometers)
filter_min_duration <-
  1 # these values are in the units defined by the parameters above: fps (seconds), measured_volume (microliters) and pixel_to_scale (micometers)
filter_detection_freq <-
  0.1 # these values are in the units defined by the parameters above: fps (seconds), measured_volume (microliters) and pixel_to_scale (micometers)
filter_median_step_length <-
  3 # these values are in the units defined by the parameters above: fps (seconds), measured_volume (microliters) and pixel_to_scale (micometers)

# directories and file names
to.data <- "/media/mendel-himself/ID_061_Ema2/AHSize/main_analysis/t0/"

######################################################################
# VIDEO ANALYSIS

#Check if all tools are installed, and if not install them
check_tools_folder(tools.path)

#Ensure computer has permission to run bftools
system(paste0("chmod a+x ", tools.path, "bftools/bf.sh"))
system(paste0("chmod a+x ", tools.path, "bftools/bfconvert"))
system(paste0("chmod a+x ", tools.path, "bftools/showinf"))

# Convert files to compressed avi (takes approx. 2.25 minutes per video)
convert_to_avi(
  to.data,
  raw.video.folder,
  raw.avi.folder,
  metadata.folder,
  tools.path,
  fps,
  video.format
)


# TESTING

# check file format and naming
# check_video_file_names(to.data,raw.avi.folder,video.description.folder,video.description.file)

# check whether the thresholds make sense (set "dark backgroud" and "red")
#check_threshold_values(to.data, raw.avi.folder, ijmacs.folder, 2, difference.lag, thresholds, tools.path,  memory.alloc)

# identify particles
locate_and_measure_particles(
  to.data,
  raw.avi.folder,
  particle.data.folder,
  difference.lag,
  min_size = particle_min_size,
  max_size = particle_max_size,
  thresholds = thresholds,
  tools.path,
  memory = memory.alloc,
  memory.per.identifier = memory.per.identifier,
  max.cores = detectCores() - 1
)

# link the particles
link_particles(
  to.data,
  particle.data.folder,
  trajectory.data.folder,
  linkrange = trajectory_link_range,
  disp = trajectory_displacement,
  start_vid = 1,
  memory = memory.alloc,
  memory_per_linkerProcess = memory.per.linker,
  raw.avi.folder,
  max.cores = detectCores() - 1,
  max_time = maximum_time_allotted_time_per_video
)

# merge info from description file and data
merge_data(
  to.data,
  particle.data.folder,
  trajectory.data.folder,
  video.description.folder,
  video.description.file,
  merged.data.folder
)

# load the merged data
load(paste0(to.data, merged.data.folder, "Master.RData"))

# filter data: minimum net displacement, their duration, the detection frequency and the median step length
trajectory.data.filtered <-
  filter_data(
    trajectory.data,
    filter_min_net_disp,
    filter_min_duration,
    filter_detection_freq,
    filter_median_step_length
  )

# summarize trajectory data to individual-based data
morph_mvt <-
  summarize_trajectories(
    trajectory.data.filtered,
    calculate.median = F,
    write = T,
    to.data,
    merged.data.folder
  )

# get sample level info
summarize_populations(
  trajectory.data.filtered,
  morph_mvt,
  write = T,
  to.data,
  merged.data.folder,
  video.description.folder,
  video.description.file,
  total_frames
)

# create overlays for validation
create.subtitle.overlays(
  to.data,
  traj.data = trajectory.data.filtered,
  raw.video.folder,
  raw.avi.folder,
  temp.overlay.folder,
  overlay.folder,
  fps,
  vid.length = total_frames / fps,
  width,
  height,
  tools.path = tools.path,
  overlay.type = "number",
  video.format
)

# Create overlays (old method)
create_overlays(
  traj.data = trajectory.data.filtered,
  to.data = to.data,
  merged.data.folder = merged.data.folder,
  raw.video.folder = raw.avi.folder,
  temp.overlay.folder = "4a_temp_overlays_old/",
  overlay.folder = "4_overlays_old/",
  width = width,
  height = height,
  difference.lag = difference.lag,
  type = "traj",
  predict_spec = F,
  contrast.enhancement = 1,
  IJ.path = "/home/mendel-himself/bemovi_tools",
  memory = memory.alloc,
  max.cores = detectCores() - 1,
  memory.per.overlay = memory.per.overlay
)

system("rm -r ijmacs")