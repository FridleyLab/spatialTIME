data.dir<-"/Volumes/Lab_Fridley/IHC/Peres/data/"
setwd(data.dir)

library(spatialIHC)
library(tidyverse)

load("TMA.clin_08172020.Rdata")
summary = summary_tma %>% mutate(image.tag = `Image Tag`,
                                 suid = as.character(suid))
clin.data = clin.data %>% mutate(suid = as.character(suid))
#################################################################################
tma.data_small = list()
count = 0
for(i in 1:length(tma.data)){
  if(tma.data[[i]]$image.tag[1] %in% c('Peres_P1_AACES_TMA 2017_[1,F].tif', 
                                       'Peres_P1_AACES_TMA 2017_[10,D].tif',
                                       'Peres_P1_AACES_TMA 2017_[11,R].tif',
                                       'Peres_P1_AACES_TMA 2017_[14,I].tif',
                                       'Peres_P1_AACES_TMA 2017_[14,O].tif')){
    tma.data_small =  rlist::list.append(tma.data_small,tma.data[[i]])
    count = 1 + count
    print(count)
  }
}


tma_mif <- create_mif(clinical_data = clin.data,
                      sample_data = summary,
                      spatial_list = tma.data_small,
                      patient_id = 'suid', 
                      sample_id = 'image.tag',
                      clean_columns = FALSE
)

tma_mif[["derived"]][["spatial_plots"]] <- plot_immunoflo(tma_mif, plot_title = 'image.tag', mnames = mnames, mlabels = mlabels, cell_type = "Classifier Label")

num_perms = 5
r_range = seq(0,200,5)
tma_mif[["derived"]][["ripleys"]] <- ripleys_k(tma_mif, id = "image.tag", 
                                                mnames = mnames,
                                                num_permutations = num_perms,
                                                r_range = r_range, 
                                                edge_correction = 'translation',
                                                keep_perm_dis = FALSE, 
                                               mlabels = mlabels)


# This final dataset should have number of spatial objects * number of markers * number of r values

final = inner_join(tma_mif$clinical, tma_mif$sample) %>% 
  inner_join(tma_mif[["derived"]][["ripleys"]])

length(mnames)
length(r_range)
length(tma_mif$spatial)
length(mnames) * length(r_range) *length(tma_mif$spatial)

dim(final)


