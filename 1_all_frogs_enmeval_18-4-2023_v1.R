###Distribution models frogs in Western Ghats###

##using enmeval and biased background

#This script starts from species occurrence data along with data extracted 
#for each occurrence point and a set of background locations with bioclim data
#it runs ENMEvaluate within a function enmeval_fit. 
#The outputs from enmeval_fit as as follows:
#1. The specified path is selected for staring outputs
#2. Folders for each family
#3. One folder for each species
#4. Whitin the speices folder
##-1.output of ENMevaluate: _enmeval
##-2.selected model: model parameters_sltd
##-3.lamdas file: lambdas file for the selected model
##-4.ENMeval results file: _res.csv
##-5.predicted distribution: _preds_avg.asc ##This is actuallt predicted with all
##occurrences, not average acrooss cv replicates
##-6.max sum bin output: mss_bin.asc
##-7.mintrain bin output: mintrain_bin.asc
##-8.predictions clipped within the lowest contour: _binzones_clip.asc
##-9.predictions clipped at barrier contours: _binzones_mask.asc
##-10.folder /masked: contains outputs masked by wg
##-11.plots: contains response curves and var imp
##-11.projections: contains projections for future scenarios

#IN this script outputs up-to 4.5 and 4.9 to 4.11 are created
## final output is 'thresholds_out' object which is a list with species name, 
## binary outputs and occurrence records for each species

#for all species models with biased background as well as random background are 
# prepared and outputs are compared for over-prediction and most suitable 
# output will be used for further analysis. 

##All families

options(java.parameters = "-Xmx20240m")
library(tidyverse)
library(dplyr)
library(rgdal)
library(raster)
rasterOptions(tmpdir = "E:/backups/routinne_backup/transfer for backup/gis/raster_temp/")
library(tmap)
library(ENMeval)

#for win
#setwd("D:/My_Documents/OneDrive - Indian Institute of Science/My_Documents/IISc_postdoc/")
setwd("D:/Aniruddha/IISc_postdoc/")
load("sdm/all_frogs_enmeval/all_frogs_sdm_enmevla_18-4-2023.RData")


source("sdm/sdm_functions.R")
source("D:/My_Documents/OneDrive - Indian Institute of Science/My_Documents/R/tutorials/sdm/dismo/maxent_in_R/workshop_maxent_R-master/code/Appendix2_prepPara.R")

source("D:/My_Documents/OneDrive - Indian Institute of Science/My_Documents/R/ani_misc_functions.R")

crs_wgs84 <- CRS("+proj=longlat +ellps=WGS84")

#####@@data@@####

temp_env <- new.env()
load("raw_data_dumps/herp_data_final_8-5-2020.RData",envir = temp_env)
frogs_final <- get("frogs_final",envir = temp_env)
frogs_final_sp <- get("frogs_final_sp",envir = temp_env)
wg_bnd <- get("wg_bnd",envir = temp_env)
wg_bnd_atree <- get("wg_bnd_atree",envir = temp_env)
wg_bnd_atree_massif <- get("wg_bnd_atree_massif",envir = temp_env)
rm(temp_env)

wg_bnd_merge <- readOGR("gis/wg_merged",layer = "wg_cepf_atree_MCF_merged")


frogs_final <- data.frame(id = 1:nrow(frogs_final),
                          frogs_final)

frogs_final_sp@data <- data.frame(id = 1:nrow(frogs_final_sp@data),
                                  frogs_final_sp@data)

#writeRaster(bias_grid_allpoints,"sdm/input/bias_files/bias_grid_allpoints.asc")

temp_env <- new.env()
load("gis/rasters_wg_3-10-2020.RData",envir = temp_env)
myexpl_uncor1 <- get("myexpl_uncor1",envir = temp_env)

rm(temp_env)

temp_env <- new.env()
load("sdm/all_frogs/all_frogs_sdm_17-12-2021.RData",envir = temp_env)
back_rnd_swd <- get("back_rnd_swd",envir = temp_env)
bias_grids_dat <- get("bias_grids_dat",envir = temp_env)
maxent_swd <- get("maxent_swd",envir = temp_env)
rm(temp_env)

occ <- table(maxent_swd$Binomial_species)

occ_5 <- occ[occ>=5]

maxent_swd_occ5 <- maxent_swd[maxent_swd$Binomial_species%in%names(occ_5),]

maxent_swd_occ5$Binomial_species <- gsub(" ","_",
                                         maxent_swd_occ5$Binomial_species)

frogs_final$Binomial_species <- gsub(" ","_",
                                     frogs_final$Binomial_species)

frogs_final_sp$Binomial_species <- gsub(" ","_",
                                        frogs_final_sp$Binomial_species)

#contours from dem
cntr <- readOGR("gis/wg_massif_divisions",layer = "dem_srtm_contour")

cntr$level <- as.numeric(cntr$level)

cntr_single <- disaggregate(cntr)
cntr_single$id <- 1:nrow(cntr_single)

cntr_poly <- SpatialPolygons(lapply(slot(cntr_single,"lines"), function(x) {
  Polygons(list(Polygon(slot(slot(x, "Lines")[[1]], "coords"))),
           ID=slot(x, "ID"))
}),
proj4string=CRS("+proj=longlat +ellps=WGS84"))

cntr_poly <- SpatialPolygonsDataFrame(cntr_poly,data = cntr_single@data)


####ENMevaluate####
path <- "E:/backups/all_frogs_enmeval_18-4-2023"

dir.create(path)
spp <- unique(maxent_swd_occ5$Binomial_species)

#enmeval_fit(spp[1],path = path)

mods_current <- lapply(spp,enmeval_fit,path = path)

mods_current1 <- do.call(rbind,mods_current[!sapply(mods_current,is.null)])

mods_current2 <- lapply(sapply(paste0(mods_current1$Binomial_species,"_bioclim_enmeval"),
                               function(x){
                                 list.files("E:/backups/all_frogs_enmeval_18-4-2023/",pattern = x,full.names = T,recursive = T)
                               }),load)


#single model could not be selected with the above criteria for 7 species 
# so selecting the model with minimum feature classes

spp <- mods_current1$Binomial_species[1]
lapply(mods_current1$Binomial_species,
       function(spp){
         
         name <- paste0(spp,"_bioclim_enmeval")
         
         p1 <-unlist(list.files("E:/backups/all_frogs_enmeval_18-4-2023",
                                pattern = name,full.names = T,recursive = T))
         
         t1 <- unlist(strsplit(p1,"/"))
         
         path_spp <- paste(t1[1:(length(t1)-1)],collapse = "/")
         
         load(unlist(list.files("E:/backups/all_frogs_enmeval_18-4-2023/",
                                pattern = name,full.names = T,recursive = T)))
         
         sel <- mod_out@results %>% 
           filter(auc.train > 0.6) %>%
           filter(or.mtp.avg == min(or.mtp.avg,na.rm = T)) %>%
           filter(auc.diff.avg == min(auc.diff.avg,na.rm = T)) %>%
           filter(ncoef == min(ncoef,na.rm = T)) %>% 
           filter(rm == max(as.numeric(as.character(rm)))) %>% 
           mutate(n_f = sapply(strsplit(as.character(fc),""),length)) %>% 
           filter(n_f == min(n_f,na.rm = T)) %>% 
           data.frame(Binomial_species = x,.)
         
         mod_sel_id <- which(names(mod_out@models)%in%sel$tune.args)
         mod_sel <- mod_out@models[[mod_sel_id]]
         
         mod_name <- paste(spp,"bioclim_preds",sel$tune.args,'sltd',sep = "_")
         save(mod_sel,file = paste0(path_spp,"/",mod_name))
         
         l1 <- mod_outputs(mod_sel)
         
         writeRaster(l1$mod_pred,
                     paste0(path_spp,"/",paste(spp,sel$tune.args,"bioclim_preds",'avg.asc',
                                               sep = "_")),overwrite = T)
         writeRaster(l1$mod_bin,
                     paste0(path_spp,"/",paste(spp,sel$tune.args,"bioclim_preds","ses_bin.asc",
                                               sep = "_")),overwrite = T)
         
         dir.create(paste0(path_spp,"/plots/"))
         
         jpeg(file=paste0(path_spp,"/plots/",paste(spp,"bioclim_preds","percon.jpeg",
                                                   sep = "_")), height=8, width = 11, units = 'in', res=300)
         print(l1$p_con)
         dev.off()
         
         jpeg(file=paste0(path_spp,"/plots/",paste(spp,"bioclim_preds","perimp.jpeg",
                                                   sep = "_")), height=8, width = 11, units = 'in', res=300)
         print(l1$p_imp)
         dev.off()
         
         jpeg(file=paste0(path_spp,"/plots/",paste(spp,"bioclim_preds","response.jpeg",
                                                   sep = "_")), height=8, width = 11, units = 'in', res=300)
         dismo::response(mod_sel)
         dev.off()
         
       })



dsn = "E:/backups/all_frogs_enmeval_18-4-2023"
#load selected models

mods_name <- list.files(dsn,pattern = "_sltd",recursive = T)

f1 <- function(x) {
  a <- load(x)
  get(a)
}

mods_current_sltd <- fetch_files(dsn,target = "_sltd",fun = f1)




####fetch pred eval####
##get evaluation and predictions from maxent output

pred_current <- fetch_pred(
  dsn = "E:/backups/all_frogs_enmeval_18-4-2023/",
  target ="_bioclim_preds_avg.asc",recursive = T)

# 
# 
# pred_had_26_50 <- fetch_pred(
#   dsn = "E:/backups/routinne_backup/transfer for backup/gis/sdm_8-10-2021/",
#   target = "/species_proj_layers_avg.asc")

eval_current <- t(do.call(cbind,lapply(mods_current_sltd,
                                       function(x){
                                         x@results
                                       })))

rownames(eval_current) <- names(mods_current_sltd)

eval_current_avg <- eval_current

eval_current_avg <- data.frame(
  maxent_swd_occ5[
    match(rownames(eval_current_avg),maxent_swd_occ5$Binomial_species),
    c("Family","Binomial_species")],
  eval_current_avg)

####make data for mapping####


thresh_current_avg <- eval_current_avg[,c(1,2,grep("Cloglog.threshold",
                                                   names(eval_current_avg)))]

names(thresh_current_avg) <- gsub("\\.Cloglog.threshold","",
                                  names(thresh_current_avg))



#make a list with species name, binary outputs and occurrence records for 
#each species. 

area_dat <- data.frame()
#map_list <- list()
#i <- 1
thresholds_out <- list()
#i <- 6
pb <- txtProgressBar(min = 1,max = nrow(thresh_current_avg),style = 3)
for( i in 1:nrow(thresh_current_avg)){
  fam <- thresh_current_avg[i,"Family"]
  spp <- thresh_current_avg[i,"Binomial_species"]
  
  threshold_min <- thresh_current_avg[i,"Minimum.training.presence"]
  threshold_max <- max(thresh_current_avg[
    i,"Maximum.training.sensitivity.plus.specificity"])
  
  pred <- pred_current[[thresh_current_avg[i,"Binomial_species"]]]
  
  bin_min <- pred
  
  bin_min[bin_min<threshold_min] <- NA
  
  bin_min[bin_min>=threshold_min] <- 1
  
  #bin_min_vec <- rasterToPolygons(bin_min,dissolve = T)
  
  bin_min_area <- sum(getValues(raster::area(bin_min,na.rm = T)),na.rm = T)
  
  #min_mased <- raster::mask(pred,bin_min)
  
  bin_max <- pred
  
  bin_max[bin_max<threshold_max] <- NA
  bin_max[bin_max>=threshold_max] <- 1
  
  #  bin_max_vec <- rasterToPolygons(bin_max,dissolve = T)
  
  bin_max_area <- sum(getValues(raster::area(bin_max,na.rm = T)),na.rm = T)
  
  
  pts <- frogs_final_sp[
    frogs_final_sp$Binomial_species ==
      thresh_current_avg[i,"Binomial_species"],]
  
  mcp <- adehabitatHR::mcp(as(pts,"SpatialPoints"))
  
  mcp_area <- raster::area(mcp)/1000000
  
  AUC <- eval_current_avg[i,grep("Test.AUC",names(eval_current_avg))]
  
  area_dat <- rbind(area_dat,
                    
                    data.frame(Family = fam,
                               Binomial_species = spp,
                               AUC = AUC,
                               threshold_min = threshold_min,
                               threshold_max = threshold_max,
                               threshold_min_area = bin_min_area,
                               theshold_max_area = bin_max_area,
                               mcp_area = mcp_area)
                    
  )
  
  outpath <- paste0(dsn,"/",fam,"/",spp)
  
  writeRaster(bin_max,paste0(dsn,"/",fam,"/",spp,"/",spp,"_mss_bin.asc"),overwrite = TRUE)
  bin_max <- raster(paste0(dsn,"/",fam,"/",spp,"/",spp,"_mss_bin.asc"))
  
  writeRaster(bin_min,paste0(dsn,"/",fam,"/",spp,"/",spp,"_mintrain_bin.asc"),overwrite = TRUE)
  bin_min <- raster(paste0(dsn,"/",fam,"/",spp,"/",spp,"_mintrain_bin.asc"))
  
  list_out <- list(	spp = spp,
                    bin_min = bin_min,
                    bin_max = bin_max,
                    pts = pts)
  thresholds_out[[spp]] <- list_out
  
  setTxtProgressBar(pb,i)
  
  
}

####save.image####
save.image("sdm/all_frogs_enmeval/all_frogs_sdm_enmevla_18-4-2023.RData")