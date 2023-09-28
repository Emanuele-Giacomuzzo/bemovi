rm(list=ls())
files <- list.files("/media/mendel-himself/ID_061_Ema2/AHSize/Scripts/", full.names = T)
for (i in files){
  system(paste0("Rscript ", i))
}