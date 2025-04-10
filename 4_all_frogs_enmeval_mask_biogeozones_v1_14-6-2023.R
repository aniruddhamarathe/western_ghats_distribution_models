## Mask all outputs to wg boudary
##limit the predicted distribution extent within the potential biogeographic
# boundaries 

#see initial description in the '/sdm_v1/all_frogs_enmeval_18-4-2024_v1.R' for 
#the outputs description

#this script will genetate numbers 8:10 of the outputs. 

# The main outputs from this script are:
#   1. bin_current_mask: Binary outputs of models under current conditions. 
#   2. all_frogs_zones: Binary outputs within biogeographic limits, after selecting 
#                       between the three methods. 

##this script will start from the combined 'thresholds_out' object so that 
# correct predictions with respect to biased or random background are used. 

temp_env <- new.env()
load("sdm_v1/frogs/all_frogs_enmeval_biased_unbiased_compare_13-6-2023.RData",envir = temp_env)
thresholds_out <- get("thresholds_out",envir = temp_env)
rm(temp_env)

temp_env <- new.env()
load("gis/wg_massif_divisions/wg_massif_divisions_25-2-2022.RData",
     envir = temp_env)
dem_crop <- get("dem_crop",envir = temp_env)
#wg_bnd_merge_massif <- get("wg_bnd_merge_massif",envir = temp_env)
rm(temp_env)
wg_bnd_merge_massif <- readOGR("gis/wg_massif_divisions",layer = "wg_bnd_merge_massif")

####functions####
zones_compare_f <- function(x){
  raw <- thresholds_out[[x]]
  clipped <- all_frogs_zones_clip[[x]]
  clipped_out <- trim(clipped$out)
  masked <- all_frogs_zones_mask[[x]]
  masked_out <- trim(masked$out)
  
  
  #ext <- as(extent(out),"SpatialPolygons")
  
  (map_masked <- tm_graticules(lines = F)+
      tm_shape(trim(raw$bin_max))+
      tm_raster(palette = "Set1",labels = "out",title = "")+
      tm_shape(masked_out)+
      tm_raster(palette = "Set2",labels = "in",title = "")+
      tm_shape(masked$barrier)+
      tm_lines()+
      tm_shape(raw$pts)+
      tm_dots(col = "yellow",size = 0.3)+
      tm_shape(raw$pts)+
      tm_dots(shape = 1, size = 0.3)+
      tm_shape(wg_bnd_atree_massif)+
      tm_borders()+
      #tm_text("MASSIF",just = "left")+
      tm_layout(main.title = x,legend.outside = T,
                legend.text.size = 1.1))
  
  map_clipped <- tm_graticules(lines = F)+
    tm_shape(trim(raw$bin_max))+
    tm_raster(palette = "Set1",labels = "out",title = "")+
    tm_shape(clipped_out)+
    tm_raster(palette = "Set2",labels = "in",title = "")+
    tm_shape(clipped$barrier)+
    tm_lines()+
    tm_shape(raw$pts)+
    tm_dots(col = "yellow",size = 0.2)+
    tm_shape(raw$pts)+
    tm_dots(shape = 1, size = 0.2)+
    tm_shape(wg_bnd_atree_massif)+
    tm_borders()+
    #tm_text("MASSIF",just = "left")+
    tm_layout(main.title = x,legend.outside = T,
              legend.text.size = 1.1)
  
  tmap_arrange(map_masked,map_clipped,ncol = 2)
}


####mask the sdm outputs to western ghats boudary

bin_current_mask <- lapply(thresholds_out,
                           function(x){
                             m <- mask(x$bin_max,wg_bnd_merge)
                             file <- x$bin_max@file@name
                             t1 <- unlist(strsplit(file,"[\\]"))
                             path <- paste(c(t1[1:(length(t1))-1],"masked"),
                                           collapse = "\\")
                             filename <- gsub("_bin","_mask",t1[length(t1)])
                             
                             if(!dir.exists(path)){
                               dir.create(path)
                             }
                             
                             writeRaster(m,paste(path,filename,sep = "\\"),overwrite = T)
                             
                             raster(paste(path,filename,sep = "\\"))
                             
                           })


##The predictions are clipped using three methods 
#1. mask by the actual contour at the infered biogeographic divide.
#2. mask by the lowest contour in occurrence points
#3. when 1 or 2 dont produce sensible results, mask by the latitudinal bin
#   in which the occurrences are present. 


all_frogs_zones_mask <- 
  lapply(thresholds_out,
         function(x){
           tryCatch(
             {
               out <- clip_sdm_wgbiogeo_mask(bin = x$bin_max,
                                             pts = x$pts,zones = wg_bnd_merge_massif)
               writeRaster(out$out,paste0(grep(x$spp,spp_paths,value = T),"/",x$spp,
                                          "_mss_binzones_mask.asc"),overwrite = T)
               out$out <- raster(paste0(grep(x$spp,spp_paths,value = T),"/",x$spp,
                                        "_mss_binzones_mask.asc"))
               out
               
             },
             error = function(x){
               message(x)
               return(NA)}
           )
         })

all_frogs_zones_clip <- 
  lapply(thresholds_out,
         function(x){
           tryCatch(
             {
               out <- clip_sdm_wgbiogeo_occelev(bin = x$bin_max,
                                                pts = x$pts)
               writeRaster(out$out,paste0(grep(x$spp,spp_paths,value = T),"/",x$spp,
                                          "_mss_binzones_clip.asc"),overwrite = T)
               out$out <- raster(paste0(grep(x$spp,spp_paths,value = T),"/",x$spp,
                                        "_mss_binzones_clip.asc"))
               out
               
             },
             error = function(x){
               message(x)
               return(NA)}
           )
         })

all_frogs_zones_clip_mask <- lapply(all_frogs_zones_clip,
                                    function(x){
                                      m <- mask(x$out,wg_bnd_merge)
                                      file <- x$out@file@name
                                      t1 <- unlist(strsplit(file,"[\\]"))
                                      path <- paste(c(t1[1:(length(t1))-1],"masked"),
                                                    collapse = "\\")
                                      filename <- gsub("binzones_clip","clip",t1[length(t1)])
                                      
                                      if(!dir.exists(path)){
                                        dir.create(path)
                                      }
                                      
                                      writeRaster(m,paste(path,filename,sep = "\\"),
                                                  overwrite = T)
                                      
                                      x$out <- raster(paste(path,filename,sep = "\\"))
                                      
                                      x
                                    })

##prepare output maps to compare the two methods and score these as good and 
#bad in the spp_zones file

all_frogs_zones_compare <- 
  lapply(names(all_frogs_zones_mask),
         function(x){
           tryCatch(
             zones_compare_f(x),
             error = function(e){ x }
           )
         }
  ) 

pdf("sdm/all_frogs_enmeval/all_frogs_sdm_withbias_masked_clipped_compare.pdf",
    onefile = T)
all_frogs_zones_compare
dev.off()

##when either the mask or clip outputs are not considered good, mask the outputs 
#using latitudinal bins

#score the predictions as good and bad based on 
#1. if there are abrupt cuts 
#2. if the occurrences are lower than the barrier elevation at boudary
#3. overall if it looks ok. 
#This is later crosschecked by experts

spp_bins <- spp_zones$spp[spp_zones$bins == "good"|spp_zones$bins == "Good"]

all_frogs_zones_bins <- 
  lapply(thresholds_out[spp_bins],
         function(x){
           tryCatch(
             {
               out <- clip_sdm_wgbiogeo(bin = x$bin_max,
                                        pts = x$pts,zones = wg_bnd_merge_massif)
               writeRaster(out$out,paste0(grep(x$spp,spp_paths,value = T),"/masked/",x$spp,
                                          "_mss_bins.asc"),overwrite = T)
               out$out <- raster(paste0(grep(x$spp,spp_paths,value = T),"/masked/",x$spp,
                                        "_mss_bins.asc"))
               out
               
             },
             error = function(x){
               message(x)
               return(NA)}
           )
         })


## combine outputs into a single list 


all_frogs_zones <- all_frogs_zones_mask[names(ssp126_zones_mask)%in%
                                          spp_zones$spp[spp_zones$masked == "good"]]

all_frogs_zones <- append(all_frogs_zones,
                          all_frogs_zones_clip_mask[
                            names(all_frogs_zones_clip_mask)%in%
                              spp_zones$spp[spp_zones$masked == "bad"&
                                              spp_zones$clipped == "good"]
                          ])

all_frogs_zones <- append(all_frogs_zones,
                          all_frogs_zones_bins[
                            names(all_frogs_zones_clip_mask)%in%
                              spp_zones$spp[spp_zones$masked == "bad"&
                                              spp_zones$bins == "good"]]
)

all_frogs_zones <- all_frogs_zones[names(thresholds_out)]

####save.image####

save.image("all_frogs_enmeval_mask_biogeozones_v1_14-6-2024.RData")

