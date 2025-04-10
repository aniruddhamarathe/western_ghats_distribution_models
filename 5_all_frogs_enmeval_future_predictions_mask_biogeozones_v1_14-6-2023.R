##Future predictions for distribution models.
#Predict distribution models for the three of the most pessimistic model in 
#SSP585 and most optimistic models in SSP126

#see initial description in the '/sdm_v1/all_frogs_enmeval_18-4-2024_v1.R' for 
#the outputs description

#this script will genetate numbers 8:10 of the outputs for predictions with 
#future climate conditions. 

# The main outputs from this script are:
#   1. bin_2070_zones_area: it contains area predicted in each future scenario
                            # after masking by western ghats, and area after 
                            # limiting distributions with biogeo boudaries by the
                            # three methods

##This script uses a function to predict an existing maxent output stored in a 
#folder and not through a dismo object. So the '.lambda' files have to be 
#exported to the corresponding folders first. This can be avoided by using 
#predict on the enmeval model outputs object. 

temp_env <- new.env()
load("sdm_v1/frogs/all_frogs_enmeval_biased_unbiased_compare_13-6-2023.RData",envir = temo_env)
thresh_current_avg <- get("thresh_current_avg",envir = temp_env)
spp_zones <- get("spp_zones",envir = temp_env)
rm(temp_env)



temp_env <- new.env()
load("all_frogs_enmeval_mask_biogeozones_v1_14-6-2024.RData",envir = temp_env)
thresholds_out <- get("thresholds_out",envir = temp_env)
rm(temp_env)

####@@lambda files@@####
#write the .lambda file to folders as the predict function used here takes the 
# .lmbdas from disk and not from model object
#           full.names = T,recursive = F)

lapply(thresholds_out,function(x){
  
  a <- unlist(strsplit(x$bin_max@file@name,"[\\]"))
  b <- paste(a[1:5],collapse = "\\")
  

  outname <- paste0(list.files(b,pattern = "_sltd$"),".lambdas")
  if(!file.exists(paste0(path_spp,"/",outname))){
    writeLines(mods_current_sltd[[x]]@lambdas,paste0(b,"/",outname))
  }
  
})


####@@predictions@@####
#predict distribution models for three most optimistic models and combine
#results.

lapply(thresholds_out, 
       function(a){
         O <- unlist(strsplit(a$bin_max@file@name,"[\\]"))
         a <- paste(a[1:5],collapse = "\\")
         r1 <- predict_maxent_spp(a,preds = "E:/backups/future_preds/INM-CM5-0_ssp126_2061-2080",overwrite = F)
         r2 <- predict_maxent_spp(a,preds = "E:/backups/future_preds/GFDL-ESM4_ssp126_2061-2080",overwrite = F)
         r3 <- predict_maxent_spp(a,preds = "E:/backups/future_preds/MIROC6_ssp126_2061-2080",overwrite = F)
         s <- stack(r1,r2,r3)
         
         r4 <- stackApply(s,c(1,1,1),median,na.rm = T)
         
         spp_name <- tail(strsplit(a,"/")[[1]],n = 1)
         
         
         outname <- paste0(spp_name,"_ensemble_ssp126min_2061-2080_med.asc")
         
         writeRaster(r4,paste0(a,"/projections/",outname),overwrite = T)
         setTxtProgressBar(pb,which(spp == a))
         
       })


#predict distribution models for three most pessimistic models and combine
#results. 

pb <- txtProgressBar(min = 0,max = length(spp),style = 3)
lapply(thresholds_out, 
       function(a){
         O <- unlist(strsplit(a$bin_max@file@name,"[\\]"))
         a <- paste(a[1:5],collapse = "\\")
         
         r1 <- predict_maxent_spp(a,preds = "E:/backups/future_preds/HadGEM3-GC31-LL_ssp585_2061-2080",overwrite = F)
         r2 <- predict_maxent_spp(a,preds = "E:/backups/future_preds/ACCESS-CM2_ssp585_2061-2080",overwrite = F)
         r3 <- predict_maxent_spp(a,preds = "E:/backups/future_preds/UKESM1-0-LL_ssp585_2061-2080",overwrite = F)
         s <- stack(r1,r2,r3)
         
         r4 <- stackApply(s,c(1,1,1),median,na.rm = T)
         
         spp_name <- tail(strsplit(a,"/")[[1]],n = 1)
         
         
         outname <- paste0(spp_name,"_ensemble_ssp585max_2061-2080_med.asc")
         
         writeRaster(r4,paste0(a,"/projections/",outname),overwrite = T)
         setTxtProgressBar(pb,which(spp == a))
         
       })

####@@mask and biogeozones@@####

##make a dataframe to hold species names, preditor data, and time period

spp_scnr_2070 <- base::cbind.data.frame(
  Family = rep(thresh_current_avg$Family,each = 2),
  Binomial_species = rep(thresh_current_avg$Binomial_species, each = 2),
  mod = "ensemble",
  ssp = c("ssp126min","ssp585max"),
  time = "2061-2080",
  cat = c("min","max"))

spp_scnr_2070$name <- with(spp_scnr_2070,paste(mod,ssp,time,sep = "_"))

spp_scnr_2070$area <- NA



##mask future predictions within WG
pb <- txtProgressBar(min = 1,max = nrow(spp_scnr_2070),style = 3)


for( i in 1:nrow(spp_scnr_2070)){
  fam <- spp_scnr_2070[i,"Family"]
  spp <- spp_scnr_2070[i,"Binomial_species"]
  
  a <- unlist(strsplit(thresholds_out[[spp]]$bin_max@file@name,"[\\]"))
  b <- paste(a[1:5],collapse = "\\")
  
  r <- get_pred(b,target = paste0(spp_scnr_2070$name[i],
                                  "_bin.asc"),
                recursive = T)
  
  outname <- gsub("_bin.asc","_mask.asc",r@file@name)
  
  a <- unlist(strsplit(outname,"[\\]"))
  
  b <- paste(a[1:5],collapse = "\\")
  
  outpath <- paste(b,"masked",a[length(a)],sep = "\\")
  
  #if(!file.exists(outpath)){
  
  
  pred_mask <- raster::mask(r,wg_bnd_merge)
  
  
  writeRaster(pred_mask,outpath,overwrite = T)
  
  
  pred_mask <- raster(outpath)
  
  crs(pred_mask) <- crs_wgs84
  
  x <- projectRaster(pred_mask,crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  
  
  spp_scnr_2070$area[i] <- raster::cellStats(x,stat = 'sum',na.rm = T)*0.801520
  
  
  #bin_2070_mask[[spp]] <- pred_mask
  # }else{
  #   pred_mask <- raster(outpath)
  #   x <- projectRaster(pred_mask,crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  #   
  #   
  #   spp_scnr_2070$area[i] <- raster::cellStats(x,stat = 'sum',na.rm = T)*0.801520
  #   
  # }
  setTxtProgressBar(pb,i)
}

## limit future predictions within biogeographic limits
#get gis objects 
temp_env <- new.env()
load("gis/wg_massif_divisions/wg_massif_divisions_25-2-2022.RData",
     envir = temp_env)
dem_crop <- get("dem_crop",envir = temp_env)
cntr_poly <- get("cntr_poly",envir = temp_env)
cntr_single <- get("cntr_single",envir = temp_env)
clip_sdm_wgbiogeo_mask <- get("clip_sdm_wgbiogeo_mask",envir = temp_env)
clip_sdm_wgbiogeo_occelev <- get("clip_sdm_wgbiogeo_occelev",envir = temp_env)
rm(temp_env)

bin_2070_zones<- 
  apply(spp_scnr_2070,1,
        function(x){
          tryCatch(
            {
              fam <- unlist(x[1])
              spp <- unlist(x[2])
              name <- unlist(x[7])
              
              filename <- paste(spp,
                                name,
                                "mask.asc",sep = "_"
              )
              
              a <- unlist(strsplit(thresholds_out[[spp]]$bin_max@file@name,"[\\]"))
              b <- paste(a[1:5],collapse = "\\")
              
              path <- paste(b,"masked",filename,sep = "/")
              
              outname_zones <- paste(spp,
                                     name,
                                     "zones.asc",sep = "_"
              )
              outpath_zones <- paste(b,"masked",outname_zones,sep = "/")
              
              outname_clip <- paste(spp,
                                    name,
                                    "clip.asc",sep = "_"
              )
              
              outpath_clip <- paste(b,"masked",outname_clip,sep = "/")
              
              outname_bins <- paste(spp,
                                    name,
                                    "bins.asc",sep = "_"
              )
              
              outpath_bins <- paste(b,"masked",outname_bins,sep = "/")
              
              bin <- raster(path)
              out_zones <- clip_sdm_wgbiogeo_mask(bin = bin,
                                                  pts = thresholds_out[[spp]]$pts,
                                                  zones = wg_bnd_merge_massif)
              
              c <- projectRaster(out_zones$out,crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
              
              
              x$area_masked <- raster::cellStats(c,stat = 'sum',na.rm = T)*0.801520
              
              writeRaster(out_zones$out,outpath_zones,overwrite = T)
              
              
              out_clip <- clip_sdm_wgbiogeo_occelev(bin = bin,
                                                    pts = thresholds_out[[spp]]$pts)
              
              c <- projectRaster(out_clip$out,crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
              
              
              x$area_clipped <- raster::cellStats(c,stat = 'sum',na.rm = T)*0.801520
              
              
              writeRaster(out_clip$out,outpath_clip,overwrite = T)
              
              x
              
              if(spp%in%spp_zones$spp[spp_zones$bins == "good"]){
                out_bins <- clip_sdm_wgbiogeo_mask(bin = bin,
                                                    pts = thresholds_out[[spp]]$pts,
                                                    zones = wg_bnd_merge_massif)
                
                c <- projectRaster(out_bins$out,crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
                x$area_bins <- raster::cellStats(c,stat = 'sum',na.rm = T)*0.801520
                
                writeRaster(out_bins$out,outpath_bins,overwrite = T)
              }
              
            },
            error = function(x){
              message(x)
              return(NA)}
          )
          
        })


bin_2070_zones_area <- do.call(rbind,
                               lapply(bin_2070_zones[!is.na(bin_2070_zones)],
                                      data.frame)
)

bin_2070_zones_area$area <- as.numeric(bin_2070_zones_area$area)  

####save.image####
#save.image("sdm_v1/frogs/all_frogs_enmeval_future_predictions_mask_biogeozones_v1_14-6-2023.RData")