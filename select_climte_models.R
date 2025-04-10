
## Select climate models for making future projectiopns

#load bio1 files into workspace (from '/gis/download_future_projs_22-7-2022.R')
#
library(terra)
terraOptions(tempdir = "E:/backups/routinne_backup/transfer for backup/gis/raster_temp/")
library(sf)
fls_bio1 <- list.files("E:/backups/routinne_backup/transfer for backup/gis/worldclim_future_projections/bio1/wc2.1_30s/",pattern = ".tif",full.names = T)

wldclim_bio1_proj <- rast(fls_bio1)

wldclim_bio1_proj_names <- unlist(lapply(wldclim_bio1_proj@ptr$filenames(),
                                         function(x){
                                           a <- unlist(strsplit(x,"/"))
                                           a <- a[length(a)]
                                           b <- unlist(strsplit(a,"_"))
                                           b <- paste(b[3:6],collapse = "_")
                                         }))

#
wg_bnd_inter_v <- vect(wg_bnd_inter,crs = "+proj=longlat +datum=WGS84 +no_defs")
l1 <- list()

for(i in 1:nrow(wg_bnd_inter_v)){
  x <- wg_bnd_inter_v[i,]
  c <- mask(crop(wldclim_bio1_proj,x),x)
  p <- as.points(c)
  p <- data.frame(p)
  names(p) <- wldclim_bio1_proj_names
  l1[[i]] <- p
}
names(l1) <- round(wg_bnd_inter$lat,digits = 2)

l1 <- lapply(l1,
             function(x){
               names(x) <- wldclim_bio1_proj_names
               x
             })

#take median of each model within each latitude zone

l1_med <- lapply(l1,
                 function(x){
                   apply(x,2,median)
                 })

l1_med <- do.call(rbind,l1_med)

l1_med <- data.frame(lat = as.numeric(rownames(l1_med)),
                     l1_med)


bio1_extract_l <- pivot_longer(l1_med,
                               cols = !lat,
                               names_to = c("var","mod","ssp","time"),
                               names_sep = "_"
)

#Not considering "BCC.CSM2.MR" model, predictions do not appear relevant
bio1_extract_l <- bio1_extract_l[!bio1_extract_l$mod == "BCC.CSM2.MR",]




#make wide for all models
#find models that rank the lowest 3 based on median temperature in a latitude zone
bio1_extract_w <- pivot_wider(bio1_extract_l,names_from = mod,values_from = value)

t1 <- apply(bio1_extract_w[bio1_extract_w$ssp == "ssp126"&bio1_extract_w$time == "2061.2080",5:29],1,function(x)order(x,decreasing = F)[1:5])

t2 <- apply(t1,1,function(x){
  t2 <- table(x)
  t2_1 <- which.max(t2)
  as.numeric(names(t2)[t2_1])
})

names(bio1_extract_w)[t2+4]


#find models that rank the highest 3 based on median temperature in a latitude zone
t3 <- apply(bio1_extract_w[bio1_extract_w$ssp == "ssp585"&bio1_extract_w$time == "2061.2080",5:29],1,function(x)order(x,decreasing = T)[1:5])

t4 <- apply(t3,1,function(x){
  t2 <- table(x)
  t2_1 <- which.max(t2)
  as.numeric(names(t2)[t2_1])
})

names(bio1_extract_w)[t4+4]

