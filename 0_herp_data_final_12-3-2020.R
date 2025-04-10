##Final species occurrence data shared by Vijay on 10-3-2020

####@@packages@@####
library(tidyverse)
library(sf) 
library(vegan)
library(rgdal)
library(raster)
library(tmap)

rasterOptions(tmpdir = "E:/backups/routinne_backup/transfer for backup/gis/raster_temp/")

source("D:/My_Documents/OneDrive - Indian Institute of Science/My_Documents/R/ani_misc_functions.R")
####@@paths and wd@@####
# #on mac
# setwd("/Users/Murali/OneDrive - Indian Institute of Science/My_Documents/IISc_postdoc")
#on win
setwd("D:/My_Documents/OneDrive - Indian Institute of Science/My_Documents/IISc_postdoc")

load("raw_data_dumps/herp_data_final_8-5-2020.RData")

ind <- readOGR("gis/India",layer = 'India')
ind_without_islands <- readOGR("gis/India",layer = "India_without_islands")

wg_bnd <- readOGR("gis/wg_bnd","westernghats")

wg_bnd_ext <- as(extent(wg_bnd),"SpatialPolygons")
crs(wg_bnd_ext) <- crs(wg_bnd)

wg_bnd_atree <- readOGR("gis/wg_boundary_without_palghat","wg_boundary")

wg_bnd_atree_ext <- as(extent(wg_bnd_atree),Class = 'SpatialPolygons')

ext1 <- rbind(wg_bnd_ext,wg_bnd_atree_ext,makeUniqueIDs = T)

crs(wg_bnd_atree_ext) <- crs(wg_bnd_atree)


##wg massifs
massif <- read.csv("gis/wg_massif_divisions/Massifs_Bhavya_project_1-5-2020.csv")
names(massif)[2:3] <- c('ymin','ymax')

rownames(massif) <- massif$MASSIF

x_coords <- rep(c(rep(ext1@bbox["x","min"],2),
                  rep(ext1@bbox["x","max"],2),
                  ext1@bbox["x","min"]),nrow(massif))

y_coords <- apply(massif[,2:3],1,function(x) c(x[1],rep(x[2],2),rep(x[1],2)),
                  simplify = F)
y_coords <- do.call(c,y_coords)

id <- rep(massif$MASSIF,each = 5)

massif_coords <- split(data.frame(x = x_coords,y = y_coords),id)

massif_coords <- massif_coords[massif$MASSIF]


massif_poly <- SpatialPolygons(mapply(function(poly, id) {
  Polygons(list(Polygon(poly)), ID=id)
}, massif_coords, massif$MASSIF))

crs(massif_poly) <- crs(wg_bnd)

massif_poly <- SpatialPolygonsDataFrame(massif_poly,data = massif,match.ID = T)

wg_bnd_massif <- raster::intersect(wg_bnd,massif_poly)


wg_bnd_atree_massif <- raster::intersect(wg_bnd_atree,massif_poly)

writeOGR(wg_bnd_atree_massif,"gis/wg_massif_divisions","wg_atree_massif",
         driver = "ESRI Shapefile",overwrite_layer = T)


## distribution data
frogs_final <- read.csv("raw_data_dumps/wgdata_10-12-2019/CEPF_final_compilation_derived_ver21_ForKartik.csv",header = T)

header <- readr::read_lines("raw_data_dumps/wgdata_10-12-2019/CEPF_final_compilation_derived_ver21_ForKartik.csv", n_max = 1)

header <- unlist(strsplit(header,","))

names(frogs_final) <- header

summary(frogs_final)

frogs_final <- frogs_final[!is.na(frogs_final$LATITUDE),]

temp <- frogs_final$Binomial_species
while(any(spchars_end(temp))){
  temp <- clip_spchars_delim(temp," ")
}

frogs_final$Binomial_species <- temp

crs_latlong <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

frogs_final_sp <- SpatialPointsDataFrame(
  as.matrix(frogs_final[,c("LONGITUDE","LATITUDE")]),frogs_final,proj4string =  crs_latlong)


temp <- sp::over(frogs_final_sp,y =ext1)

frogs_final<- frogs_final[!is.na(temp),]

frogs_final_sp <- SpatialPointsDataFrame(
  as.matrix(frogs_final[,c("LONGITUDE","LATITUDE")]),frogs_final,proj4string =  crs_latlong)


frogs_final_grp <- frogs_final

frogs_final_grp$grp <- cut(frogs_final$LATITUDE,c(8.0,10.55,15.8,22.0), labels = c("south","central","north"))

num_grp <- tapply(frogs_final_grp$grp,frogs_final_grp$Binomial_species,
                  function(x)length(unique(x)))

frogs_final_grp$num_grp <- num_grp[match(frogs_final_grp$Binomial_species,
                                         names(num_grp))]

write.csv(frogs_final_grp,"raw_data_dumps/wgdata_10-12-2019/CEPF_final_compilation_derived_ver21_latgrp_11-6-2022.csv")

frogs_final_grp_sp <- SpatialPointsDataFrame(
  as.matrix(frogs_final_grp[,c("LONGITUDE","LATITUDE")]),frogs_final_grp,
proj4string =  crs_latlong)

inter1 <- sp::over(frogs_final_grp_sp,wg_bnd_union)

#extract dem to points
temp_env <- new.env()
load("gis/rasters_3-10-2020.RData",envir = temp_env)
dem_merge <- get("dem_merge",envir = temp_env)
rm(temp_env)
frogs_final_grp_sp$dem <- raster::extract(dem_merge,frogs_final_grp_sp)

spp_dem_sum <- data.frame(med = tapply(frogs_final_grp_sp@data$dem,frogs_final_grp_sp@data$Binomial_species,median,na.rm = T),
mean = tapply(frogs_final_grp_sp@data$dem,frogs_final_grp_sp@data$Binomial_species,mean,na.rm = T),
min = tapply(frogs_final_grp_sp@data$dem,frogs_final_grp_sp@data$Binomial_species,min,na.rm = T),
max = tapply(frogs_final_grp_sp@data$dem,frogs_final_grp_sp@data$Binomial_species,max,na.rm = T)
)

spp_dem_sum$grp <- tapply(frogs_final_grp_sp$grp,
                          frogs_final_grp_sp@data$Binomial_species,
                          function(x){
                            t <- table(x)
                            names(t)[which.max(t)]
                          })

spp_dem_sum$med_lat <- tapply(frogs_final_grp_sp$LATITUDE,
                              frogs_final_grp_sp$Binomial_species,
                              median,na.rm = T)



write.csv(frogs_final_grp_sp@data,"raw_data_dumps/wgdata_10-12-2019/CEPF_final_compilation_derived_ver21_dem_26-6-2022.csv")

write.csv(spp_dem_sum,"raw_data_dumps/wgdata_10-12-2019/CEPF_final_compilation_derived_ver21_spp_dem_sum_29-8-2022.csv")



temp_env <- new.env()
load("gis/rasters_3-10-2020.RData",envir = temp_env)
dem_merge <- get("dem_merge",envir = temp_env)

temp_env <- new.env()
load("gis/rasters_wg_3-10-2020.RData",envir = temp_env)
myexpl_uncor1 <- get("myexpl_uncor1",envir = temp_env)


grid_1deg <- readOGR("raw_data_dumps/wgdata_10-12-2019/herp_data_gis",
                     "wg_1deg_grid")

grid_1deg@data$id <- as.character(grid_1deg@data$id)


crs_latlong <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


#frogs_final_inter <- shapefiles::read.dbf("raw_data_dumps/wgdata_10-12-2019/herp_data_gis/cepf_final_point_data_20-12-2019_grid_intersect.dbf")

frogs_final_inter <- sp::over(frogs_final_sp,grid_1deg)


summary(frogs_final_inter)

#check if all latitudes are within correct grids
t1 <- frogs_final$LATITUDE>=frogs_final_inter$bottom&frogs_final$LATITUDE<=frogs_final_inter$top

all(t1,na.rm = T)#should be true, and it is

#frogs_final[is.na(t1),] #records from north east, and three incorrect longitudes

t2 <- frogs_final$LONGITUDE>=frogs_final_inter$left&frogs_final$LONGITUDE<=frogs_final_inter$right

all(t2,na.rm = T) #should be true, and it is

frogs_final[is.na(t2),] #records from north east, and three incorrect longitudes

frogs_final_inter <- data.frame(frogs_final,frogs_final_inter)

#####@@Bias Grids@@####

bias_grid_allpoints <- MASS::kde2d(frogs_final$LONGITUDE,
                                   frogs_final$LATITUDE,
                                   n = c(ncol(myexpl_uncor[[1]]),
                                         nrow(myexpl_uncor[[1]])),
                                   lims = as.vector(extent(myexpl_uncor[[1]]))
                                   )
                                   
bias_grid_allpoints <- raster(bias_grid_allpoints,)

bias_grid_allpoints <- resample(bias_grid_allpoints,myexpl_uncor[[1]])

wg_bnd_union <- rbind(as(wg_bnd_atree,'SpatialPolygons'),
                      as(wg_bnd,'SpatialPolygons'),makeUniqueIDs = T)

bias_grid_allpoints <- mask(bias_grid_allpoints,mask = wg_bnd_union)

# pts <- sample(which(!is.na(values(bias_grid_allpoints))),10000,
#               replace = F,
#               prob = values(bias_grid_allpoints)[!is.na(values(bias_grid_allpoints))])
# 
# pts_xy <- xyFromCell(bias_grid_allpoints,pts)
# 
# all_frogs_TGS_sp <- SpatialPoints(pts_xy,proj4string = crs(bias_grid_allpoints))

all_frogs_TGS <- dismo::randomPoints(bias_grid_allpoints,10000,prob = T)

all_frogs_TGS_sp <- SpatialPoints(all_frogs_TGS,
                                  proj4string = crs(frogs_final_sp)) 

all_frogs_TGS_dat <- raster::extract(myexpl_uncor,all_frogs_TGS_sp,df = T)


bias_grid_function <- function(x,long_col,lat_col,n){
  
  den <- MASS::kde2d(x[,long_col],x[,lat_col],
                     n = n,
                     lims = as.vector(extent(myexpl_uncor[[1]]))
                     )
  den_ras <- raster(den)
  
  den_ras <- resample(den_ras,myexpl_uncor[[1]])
  
  den_ras <- mask(den_ras,mask = ind_without_islands)
  
  # pts <- sample(which(!is.na(values(den_ras))),10000,replace = F,prob = values(den_ras)[!is.na(values(den_ras))])
  # 
  # pts_xy <- xyFromCell(den_ras,pts)
  # 
  # SpatialPoints(pts_xy,proj4string = crs(myexpl_uncor))
}

bias_pts_function <- function(den_ras){
  pts <- sample(which(!is.na(values(den_ras))),10000,replace = F,
                prob = values(den_ras)[!is.na(values(den_ras))])
  
  pts_xy <- xyFromCell(den_ras,pts)
  
  SpatialPoints(pts_xy,proj4string = crs(myexpl_uncor))
}


temp1 <- lapply(unique(frogs_final$Family),
                function(x)frogs_final[frogs_final$Family == x,])

names(temp1) <- unique(frogs_final$Family)

bias_grids_r <- lapply(temp1,bias_grid_function,long_col = 3,lat_col = 2,
                     n = c(ncol(myexpl_uncor1[[1]]),
                           nrow(myexpl_uncor1[[1]])
                     )
)


bias_grids <- lapply(bias_grids_r,bias_pts_function)
library(tmap)

#names(bias_grids) <- names(bias_grids_r)



grid_family_occur_plots <- grid_family_occur_plots[match(names(bias_grids),
                                                         names(grid_family_occur_plots))]



tm_shape(bias_grids[[1]])+
  tm_raster(title = "")+
  tm_shape(grid_1deg_wg)+
  tm_borders()+
  tm_layout(title = names(bias_grids)[1],
            title.position = c("RIGHT","TOP"),
            legend.position = c("RIGHT","TOP"))



bias_grid_plots <- lapply(1:length(bias_grids),
                          function(x){
                            tm_shape(bias_grids[[x]],bbox = grid_1deg_wg)+
                              tm_raster(title = "",)+
                              tm_shape(grid_1deg_wg)+
                              tm_borders()+
                              tm_layout(title = names(bias_grids)[x],
                                        title.position = c("RIGHT","TOP"),
                                        legend.position = c("RIGHT","TOP"))
                          })
tmap_arrange(bias_grid_plots)



####@@community matrices and diversity@@####
#making community matrices for each taxa 

frogs_final_mat <- with(frogs_final_inter,tapply(Binomial_species,
                                                 list(id,Binomial_species),
                                                 length))
frogs_final_mat[is.na(frogs_final_mat)] <- 0



####__estimates__####

#biodiversity estimates for all frogs

frogs_final_estimate <- t(as.matrix(estimateR(frogs_final_mat)))

grid_abund <- apply(frogs_final_mat,1,sum,na.rm = T)

#colnames(grid_abund) <- paste(colnames(grid_abund),"abund",sep = "_")


frogs_final_estimate_dat <- data.frame(grid_1deg@data[match(rownames(frogs_final_estimate),
                                                     grid_1deg@data$id),],
                                frogs_final_estimate)


frogs_final_estimate_dat$lat <- with(frogs_final_estimate_dat,(top+bottom)/2)
frogs_final_estimate_dat$long <- with(frogs_final_estimate_dat,(left+right)/2)

frogs_final_estimate_dat <- data.frame(frogs_final_estimate_dat,
                                abund = grid_abund[match(frogs_final_estimate_dat$id,
                                                 names(grid_abund))])

frogs_final_estimate_long <- frogs_final_estimate_dat %>% 
  select(id:abund) %>% 
  select(id,lat,long,everything(),-contains("se")) %>% 
  pivot_longer(cols = S.obs:abund,names_to = "estimate") 

#number of occurrencens of families in each grid
grid_family_occur <- frogs_final_inter %>%
  group_by(id,Family)%>%
  summarise(abund = length(Binomial_species)) %>% 
  pivot_wider(id_cols = id,names_from = Family,values_from = abund)

grid_1deg@data <- data.frame(grid_1deg@data,
                             frogs_final_estimate_dat[
                               match(grid_1deg@data$id,
                                     frogs_final_estimate_dat$id),-(1:5)])

grid_1deg@data <- data.frame(grid_1deg@data,
                             grid_family_occur[
                               match(grid_1deg@data$id,
                                     grid_family_occur$id),-1])

#keep cells that intersect with wg bnd

temp_intersect <- raster::intersect(x=grid_1deg,y=wg_bnd)

grid_1deg_wg <- grid_1deg[grid_1deg$id%in%temp_intersect@data$id,]

frogs_final_estimate_long_wg <- frogs_final_estimate_long[
  frogs_final_estimate_long$id%in%temp_intersect$id,]

Frogs_abund <- tm_shape(grid_1deg_wg)+
  tm_borders(col = 2)+
  tm_fill(col = "abund",title = "Occurrences")+
  tm_text("abund",size = 0.7)+
  #tm_compass(type = "arrow",position = "left")+
  tm_graticules(alpha = 0)

grid_family_occur_plots <- lapply(names(grid_1deg_wg@data)[14:21],
       function(x){
         tm_shape(grid_1deg_wg)+
           tm_borders(col = 2)+
           tm_fill(col = x,legend.show = F)+
           tm_text(x,size = 0.7)+
           #tm_compass(type = "arrow",position = "left")+
           tm_graticules(alpha = 0)+
           tm_grid(labels.show = F)+
           tm_layout(title = x,title.position = c("LEFT","TOP"))
       })

names(grid_family_occur_plots) <- names(grid_1deg_wg@data)[14:21]

tmap_arrange(grid_family_occur_plots)


abund_lat <- frogs_final_estimate_long_wg %>% 
  filter(estimate == "abund") %>% 
  ggplot()+
  aes(x=lat,y=value)+
  geom_point(size = 2)+
  labs(x="Latitude",y="Abundance")+
  theme_bw()

sprich_lat <- frogs_final_estimate_long_wg %>% 
  filter(!estimate == "abund") %>% 
  ggplot()+
  aes(x=lat,y=value)+
  geom_point(shape = 1,size = 2)+
  facet_grid(.~estimate)+
  labs(x="latitude",y="Species richness")+
  theme_bw()

sprich_lat+
  geom_point(shape = 4,size = 2,colour = "red",data = 
               frogs_final_estimate_long_wg %>% 
               filter(id%in%id[estimate == "abund"&value<11],!estimate == "abund"))
####@@RSFD and geometric constraints@@####
frogs_latrange <- frogs_final_inter %>% 
  group_by(Binomial_species) %>% 
  summarise(lat_min = min(LATITUDE),
            lat_max = max(LATITUDE),
            lat_range = lat_max - lat_min,
            lat_mid = lat_min + (lat_range/2))

summary(frogs_latrange)

ggplot(frogs_latrange)+
  geom_point(aes(x=lat_mid,y=lat_range))+
  geom_vline(xintercept = c(10.5,11),colour = "coral")+
  geom_segment(aes(x=lat_min,xend = lat_max,y=lat_range,
                   yend=lat_range))+
  scale_y_continuous(breaks = 0:12,minor_breaks = seq(0,12,0.5))+
  scale_x_continuous(breaks = 8:25,minor_breaks = seq(8,25,0.5))+
  labs(x="Latitudinal position (mid point)", y="Latitudinal extents")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


ggplot(frogs_latrange)+
  geom_histogram(aes(x=lat_range))

sites <- data.frame(min = seq(min(frogs_latrange$lat_min),
             max(frogs_latrange$lat_max),1))
sites$max = sites$min+1
sites$mid = sites$min+((sites$max - sites$min)/2)


             
calc_inter_sprich <- function(range_dat,range_min,range_max,sites){
  #get inter polated species richness based on species range extents, and
  #data frmae of intervals with min, max and midpoints for which species 
  #richness is required
  
  range <- cbind(range_dat[,range_min],range_dat[,range_max])
  
  sites1 <- apply(range, 1, function(z) {
    sites$mid[sites$min >= z[1] & sites$max <= z[2]]
  })
  
  sites1 <- do.call(c, as.list(sites1))
  data.frame(table(sites1))
}

sprich_inter <- calc_inter_sprich(frogs_latrange,2,3,sites)
sprich_inter$sites1 <- as.numeric(as.character(sprich_inter$sites1))


ggplot(sprich_inter)+
  geom_point(aes(x=sites1,y=Freq))+
  #facet_grid(.~taxa)+
  labs(x="Latitude",y="Interpolated species richness")+
  #scale_x_continuous(breaks = seq(8,19,1))+
  theme_bw()

####@@save.image@@####
save.image("raw_data_dumps/herp_data_final_8-5-2020.RData")

