#compare distributions with biased and unbiased background points

library(tidyverse)
library(rgdal)
library(raster)
rasterOptions(tmpdir = "E:/backups/routinne_backup/transfer for backup/gis/raster_temp/")
library(tmap)
library(dismo)

#for win

setwd("D:/Aniruddha/IISc_postdoc/")
load("sdm_v1/frogs/all_frogs_enmeval_biased_unbiased_compare_13-6-2023.RData")

#load predictions with biased background
temp_env <- new.env()
load("sdm/all_frogs_enmeval//all_frogs_sdm_enmevla_18-4-2023.RData",envir = temp_env)
pred_bias <- get("pred_current",envir = temp_env)
eval_bias <- get("eval_current_avg",envir = temp_env)
thresh_bias <- get("thresh_current_avg",envir = temp_env)
frogs_final_sp <- get("frogs_final_sp",envir = temp_env)
wg_bnd_atree <- get("wg_bnd_atree",envir = temp_env)
area_dat_bias <- get("area_dat",envir = temp_env)
thresholds_out_bias <- get("thresholds_out",envir = temp_env)
rm(temp_env)

#load predictions with random background
temp_env <- new.env()
load("sdm/all_frogs_enmeval_without_bias/all_frogs_enmeval_without_bias_14-6-2023.RData",envir = temp_env)
pred_nobias <- get("pred_current",envir = temp_env)
eval_nobias <- get("eval_current_avg",envir = temp_env)
thresh_nobias <- get("thresh_current_avg",envir = temp_env)
area_dat_nobias <- get("area_dat",envir = temp_env)
thresholds_out_nobias <- get("thresholds_out",envir = temp_env)
rm(temp_env)

#get data for mapping from sdm outputs
#predictions with biased background

#name:character;name of speices 
#pts:spatialpoints. Occurrence points
#pred:raster. predicted distribution
#threshold: named numeric. named vector of estimated thresholds
#threshold_metric: character. threshold metric to be used. must be in names of
#                 threshold
#eval: object of class 'eval' in 'package:dismo', or named vector of evaluation metrics
#metric: evaluation metric to be used.


dat_bias <- list()
dat_nobias <- list()

for(i in 1:nrow(thresh_bias)){
  family <- thresh_bias$Family[i]
  spp <- thresh_bias$Binomial_species[i]
  l1 <- list()
  l1$name <- spp
  l1$pts <- frogs_final_sp[frogs_final_sp$Binomial_species == spp,]
  l1$pred <- pred_bias[[spp]]
  l1$threshold <- thresh_bias[i,c(-1,-2)]
  l1$threshold_metric <- "Maximum.training.sensitivity.plus.specificity"
  l1$eval <- eval_bias[i,]
  l1$metric <- "Training.AUC"
  dat_bias[[spp]] <- l1
  
  l1 <- list()
  l1$name <- spp
  l1$pts <- frogs_final_sp[frogs_final_sp$Binomial_species == spp,]
  l1$pred <- pred_nobias[[spp]]
  l1$threshold <- thresh_nobias[i,c(-1,-2)]
  l1$threshold_metric <- "Maximum.training.sensitivity.plus.specificity"
  l1$eval <- eval_nobias[i,]
  l1$metric <- "Training.AUC"
  dat_nobias[[spp]] <- l1
}


maps_bias <- lapply(dat_bias,
                    function(x){
                      with(x,sdm_pred_plots_fun1(name = name,pts = pts,
                                                 pred = pred,
                                                 threshold = threshold,
                                                 threshold_metric = threshold_metric,
                                                 eval = eval,metric = metric,
                                                 binary = F))
                      
                    })

maps_bias_out <- lapply(maps_bias,function(x) tmap_arrange(x$pred,x$inset))

maps_nobias <- lapply(dat_nobias,
                      function(x){
                        with(x,sdm_pred_plots_fun1(name = name,pts = pts,
                                                   pred = pred,
                                                   threshold = threshold,
                                                   threshold_metric = threshold_metric,
                                                   eval = eval,metric = metric,
                                                   binary = F))
                        
                      })

maps_nobias_out <- lapply(maps_nobias,function(x) tmap_arrange(x$pred,x$inset))


maps_pred_compare <- list()

for(i in names(maps_bias)){
  maps_pred_compare[[i]] <- tmap_arrange(
    tm_graticules(lines = F)+
      maps_bias[[i]]$pred,
    tm_graticules(lines = F)+
      maps_nobias[[i]]$pred
  )
}

##make maps to compare predictions with biased and random background points 
## side by side. 
ext2 <- extent(wg_bnd_atree)
ext2_asp <- (ext2@xmax - ext2@xmin)/(ext2@ymax - ext2@ymin)
w <- 0.1
h <- w/ext2_asp
vp <- grid::viewport(x=0.95,y= 0.7,width = w,height = h,
                     just = c("right", "top"))

vp1 <- grid::viewport(x = 0.46,y = 0.7, width = w,height = h,
                      just=c("right", "top"))


pdf("sdm/all_frogs_enmeval/all_frogs_enmeval_compare_biased_unbiased_outputs_14-6-2023.pdf",
    width = 12.3,height = 11.4,onefile = T)

for(i in 1:length(maps_bias)){
  print(tmap_arrange(maps_bias[[i]]$main+
                       tm_layout(main.title.size = 0.8,
                                 legend.outside = T,
                                 legend.position = c(0.2,0.75)),
                     maps_nobias[[i]]$main+
                       tm_layout(main.title.size = 0.8,
                                 legend.outside = T,
                                 legend.position = c(0.2,0.75)),
                     ncol = 2
  )
  )
  print(maps_bias[[i]]$inset,vp = vp1)
  print(maps_nobias[[i]]$inset,vp = vp)
}

dev.off()

par_dflt <- par()
par(mfcol = c(1,2))


##These maps are visually examined by me, Kartik and other experts - Vijay, Varun,
#and Dinesh (but mostly dinesh!!) and each distribtion is scored as better in 
#biased background or ranodm. 

## Furhter cocde combines the 'thresholds_out' objects from biased and random 
#background outputs, that contain binary files, and are the strarting point for 
#most analysis, will be combined to contain only the suitable outputs among
#biased and random backgrounds

##make a observed distribution summary file to score the outputs

temp1 <- lapply(thresholds_out,
                function(x){
                  tryCatch({
                    pts <- x$pts
                    zones <-  wg_bnd_atree_massif
                    pts_zones <- unique(raster::intersect(pts,zones)@data$MASSIF)
                    
                    if(length(pts_zones) == 0){
                      data.frame(Nzone = NA,Szone = NA,N_bnd = NA,S_bnd = NA,
                                 min_elev = NA,num_zones = NA) 
                    }else{
                      
                      zones_dat <- zones[zones$MASSIF%in%pts_zones,]
                      
                      Nzone <- zones_dat@data[which.max(zones_dat$ymax),"MASSIF"]
                      Szone <- zones_dat@data[which.min(zones_dat$ymin),"MASSIF"]
                      N_bnd <- zones_dat@data[which.max(zones_dat$ymax),"nzone"]
                      S_bnd <- zones_dat@data[which.min(zones_dat$ymin),"szone"]
                      
                      elev <- raster::extract(dem_crop,pts)
                      
                      min_elev <- min(elev)
                      
                      data.frame(Nzone = Nzone,Szone = Szone,N_bnd = N_bnd,S_bnd = S_bnd,
                                 min_elev = min_elev,num_zones = nrow(zones_dat))}
                  },
                  error = function(e){
                    data.frame(Nzone = NA,Szone = NA,N_bnd = NA,S_bnd = NA,
                               min_elev = NA,num_zones = NA)
                    
                  })
                  
                })

temp1 <- do.call(rbind,temp1)

temp1 <- data.frame(frogs_final$Family[match(rownames(temp1),frogs_final$Binomial_species)],
                    rownames(temp1),temp1)
write.csv(temp1,"sdm/all_frogs_enmeval/spp_zones_25-4-2023.csv")

#read the current version of this file. 

spp_zones <- read.csv("sdm/all_frogs_enmeval/spp_zones_21-4-2024.csv")

spp_nobias <- spp_zones$spp[spp_zones$bias == "no"]

thresholds_out <- thresholds_out_bias

thresholds_out[spp_nobias] <- thresholds_out_nobias[spp_nobias]

area_dat <- area_dat_bias

area_dat[area_dat$Binomial_species%in%spp_nobias,] <- area_dat_nobias[area_dat$Binomial_species%in%spp_nobias,]

thresh_current_avg <- thresh_bias
thresh_current_avg <- thresh_nobias[spp_nobias]

#thresholds_out in this workspace will now contain the binary predictions
#correctly from biased and random backgrounds

####save.image####
#save.image("sdm_v1/frogs/all_frogs_enmeval_biased_unbiased_compare_13-6-2023.RData")