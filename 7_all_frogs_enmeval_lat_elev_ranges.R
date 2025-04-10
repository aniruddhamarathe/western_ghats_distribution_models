#get latitudinal and elevational extents of species

# The main outputs from this script are:
#   1. ext_lat: latitudinal extents, range, and midpoints, of species.
#         calculated based on thresholded binary predictions. The latitudinal 
#         extents are lstitudes that cover 95% of the predicted area to avoid 
#         overestimation due to isolated pixels. 
#   2. ext_elev: elevational extents, range and midpoints. Calculatedd as 
#         maximum and minimum elevations within the latitudinal extents in 
#         'ext_lat'

temp_env <- new.env()
load("sdm_v1/frogs/all_frogs_enmeval_future_prediction_mask_biogeozones_v1_14-6-2024.RData",envir = temp_env)

bin_current_stk <- get("bin_current_stk",
                       envir = temp_env)
bin_ssp126_stk <- get("bin_ssp126_stk",envir = temp_env)
bin_ssp585_stk <- get("bin_ssp585_stk",envir = temp_env)

rm(temp_env)

wgs84_proj4 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

####@@Latitudinal extents@@####
####__current__####

##convert binary files to area 

crs(bin_current_stk) <- crs(wg_bnd_merge)

bin_current_stk_aea <- projectRaster(bin_current_stk,crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

bin_current_area <- bin_current_stk_aea*0.801520


##latitudinal  range of each species are latitudes that cover 
## 95% of the predicted area towards the south and the north
library(sf)
grd_lat_sf <- st_as_sf(grd_lat)

grd_lat_sf_aea <- st_transform(grd_lat_sf,crs(bin_current_stk_aea))

grd_lat_aea <- as_Spatial(grd_lat_sf_aea)

rm(grd_lat_sf)
rm(grd_lat_sf_aea)

area_lat <- raster::extract(bin_current_area,grd_lat_aea,fun = sum,na.rm = T,df = T)

area_lat <- area_lat[order(grd_lat_aea$lat,decreasing = T),]
x <- names(area_lat)[2]

# Calculate proportion of area within each latitudinal zone. Identify the zones
#with 95% of the area and calculate the latitudinal extents. 
ext_lat <- do.call(rbind,
                   lapply(names(area_lat)[2:ncol(area_lat)],
                          function(x){
                            y <- area_lat[,x]
                            per <- y/sum(y)
                            max <- min(which(cumsum(per)>=0.05))
                            min <- 18-min(which(cumsum(rev(per))>=0.05))
                            grd1 <- grd_lat[grd_lat$ID%in%area_lat$ID[max:min],]
                            r <- trim(crop(bin_current_stk[[x]],grd1))
                            bbox(r)[2,]
                          }))

#note: here the per is proportion of area within each latitudinal zone. 
#It is ordered NOrth to South. It is reversed to calculate the minimum latitude.
#  is reversed to calculate calculate the lati

rownames(ext_lat) <- names(area_lat)[2:ncol(area_lat)]

ext_lat <- data.frame(ext_lat)

ext_lat$range <- ext_lat$max - ext_lat$min

ext_lat$mid <- ext_lat$min + (ext_lat$range/2)

####__ssp585__####

##convert binary files to area 


crs(bin_ssp585_stk) <- crs(wg_bnd_merge)

bin_ssp585_stk_aea <- projectRaster(bin_ssp585_stk,crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

bin_ssp585_area <- bin_ssp585_stk_aea*0.801520

ssp585_area_lat <- raster::extract(bin_ssp585_area,grd_lat_aea,fun = sum,na.rm = T,df = T)

ssp585_area_lat <- ssp585_area_lat[order(grd_lat_aea$lat,decreasing = T),]

x <- names(ssp585_area_lat)[2]

ssp585_ext_lat <- do.call(rbind,
                          lapply(names(ssp585_area_lat)[2:ncol(ssp585_area_lat)],
                                 function(x){
                                   y <- ssp585_area_lat[,x]
                                   if(any(y>0)){
                                     per <- y/sum(y)
                                     max <- min(which(cumsum(per)>=0.05))
                                     min <- 18 - min(which(cumsum(rev(per))>=0.05))
                                     grd1 <- grd_lat[grd_lat$ID%in%ssp585_area_lat$ID[max:min],]
                                     r <- trim(crop(bin_ssp585_stk[[x]],grd1))
                                     bbox(r)[2,]
                                   }else{c(NA,NA)}
                                   
                                 }))

rownames(ssp585_ext_lat) <- names(ssp585_area_lat)[2:ncol(ssp585_area_lat)]


ssp585_ext_lat <- data.frame(ssp585_ext_lat)

ssp585_ext_lat$range <- ssp585_ext_lat$max - ssp585_ext_lat$min

ssp585_ext_lat$mid <- ssp585_ext_lat$min + (ssp585_ext_lat$range/2)

ext_lat$ssp585_range <- ssp585_ext_lat$range[
  match(rownames(ext_lat),rownames(ssp585_ext_lat))]

ext_lat$ssp585_mid <- ssp585_ext_lat$mid[
  match(rownames(ext_lat),rownames(ssp585_ext_lat))]


jpeg("analysis/present_future_distribution/lat_ext_21-12-2022.jpg",width = 80,
     height = 80,units = "mm",quality = 100,res = 300)

{
  layout(mat = matrix(c(1,2),ncol = 1),heights = c(1,3))
  par(mar = c(0,5,1,2),mgp = c(1.5,0.5,0))
  #plot(x = 1:10,y = seq(0,40,length.out = 10),type = "n")
  barplot(table(round(ext_lat$range)),xaxt = "n",yaxt = "n",col = NA,
          ylab = "Frequency",cex.lab = 0.6)
  axis(2,at = seq(1,13,2),labels = seq(1,13,2),cex.axis = 0.6)
  
  par(mar = c(5,5,0.2,2),mgp = c(1.6,0.5,0),bty = "L")
  boxplot(ext_lat$ssp585_range ~round(ext_lat$range),
          xlab = "Current latitudinal\nextents",
          ylab = "Future latitudinal\nextents",col = NA,
          cex.lab = 0.6,cex.axis = 0.6,ylim = c(0,12))
  points(0:12,type = "l",lwd = 1)
  
}
dev.off()


v1 <- ext_lat$ssp585_range
v1[is.na(v1)] <- 0

l1 <- lm((v1 - ext_lat$range)~scale(ext_lat$range))

l2 <- lm((v1-ext_lat$range)~0+ext_lat$range)

abline(l1)
abline(l2,col = 'red')

summary(lm((v1-ext_lat$range)~0+ext_lat$range))

plot((v1-ext_lat$range)~ext_lat$range)

jpeg("analysis/present_future_distribution/lat_mid_21-12-2022.jpg",width = 80,
     height = 80,units = "mm",quality = 100,res = 300)

{
  layout(mat = matrix(c(1,2),ncol = 1),heights = c(1,3))
  par(mar = c(0,5,1,2),mgp = c(1.5,0.5,0))
  #plot(x = 1:10,y = seq(0,40,length.out = 10),type = "n")
  barplot(table(round(ext_lat$mid)),xaxt = "n",yaxt = "n",col = NA,
          ylab = "Frequency",cex.lab = 0.6)
  axis(2,at = seq(5,30,5),labels = seq(5,30,5),cex.axis = 0.6)
  par(mar = c(5,5,0.2,2),mgp = c(1.6,0.5,0))
  boxplot(ext_lat$ssp585_mid ~round(ext_lat$mid),
          xlab = "Current latitudinal\nmidpoint",
          ylab = "Future latitudinal\nmidpoint",col = NA,cex.lab = 0.6,
          cex.axis = 0.6)
  points(9:19,type = "l",lwd = 1)
}
dev.off()


l2 <- lm((ext_lat$ssp585_mid - ext_lat$mid)~ext_lat$mid)
points((0.0996+(-0.0987*c(9:19))),type = 'l')


####__sr by latitudinal extents__####

#get latitudinal extents of each latitudinal bin in grd_lat

grd_lat_ext <- data.frame(do.call(rbind,lapply(grd_lat@polygons,
                                               function(c){
                                                 d <- c@Polygons[[1]]@coords[,'y']
                                                 c(min = min(d),max = max(d))
                                               })))

grd_lat_ext$ID <- grd_lat$ID
grd_lat_ext <- grd_lat_ext[order(grd_lat_ext$min,decreasing = T),]


x <- c(8.55,19.38)
##current
current_lat_spmat <- apply(ext_lat,1,
                           function(x){
                             v <- rep(0,nrow(grd_lat_ext))
                             v[which(grd_lat_ext$max>=x[1]&grd_lat_ext$min<=x[2])] <- 1
                             v
                           })
#rownames(current_lat_spmat) <- grd_lat$lat

sr_current_lat <- apply(current_lat_spmat,1,sum)

##ssp585
ssp585_lat_spmat <- apply(ssp585_ext_lat[complete.cases(ssp585_ext_lat),],1,
                          function(x){
                            v <- rep(0,nrow(grd_lat_ext))
                            v[which(grd_lat_ext$max>=x[1]&grd_lat_ext$min<=x[2])] <- 1
                            v
                          })
#rownames(ssp585_lat_spmat) <- grd_lat$lat

sr_ssp585_lat <- apply(ssp585_lat_spmat,1,sum)

sr_lat <- data.frame(lat = grd_lat$lat[match(grd_lat_ext$ID,grd_lat$ID)],
                     current = sr_current_lat,
                     ssp585 = sr_ssp585_lat)


plot(current~lat,data = sr_lat,type = 'n')
points(current~lat,data = sr_lat,pch = 16)
points(current~lat,data = sr_lat,type = 'l')
points(ssp585~lat,data = sr_lat,pch = 16,col = 'red')
points(ssp585~lat,data = sr_lat,type = 'l',col = 'red',lty = 3)

####regression for latitude sprich####

sr_lat$lat_scl <- scale(sr_lat$lat)

sr_lat_l <- sr_lat %>% 
  pivot_longer(c(current,ssp585),names_to = 'cat',values_to = 'richness')

sr_lat_reg1 <- glm(richness ~ lat*cat+I(lat^2)*cat,data = sr_lat_l,
                   family = poisson)

summary(sr_lat_reg1)#AIC: 218.73


sr_lat_reg2 <- glm(richness ~ lat*cat+I(lat^2),data = sr_lat_l,
                   family = poisson)

summary(sr_lat_reg2)#AIC: 216.73

sr_lat_reg3 <- glm(richness ~ cat+lat+I(lat^2),data = sr_lat_l,
                   family = poisson)

summary(sr_lat_reg3)#AIC: 214.79

sr_lat_reg4 <- glm(richness ~ lat_scl+I(lat_scl^2),data = sr_lat_l,
                   family = poisson)

summary(sr_lat_reg4)#AIC: 226.79

library(AICcmodavg)

aictab_latitude <- aictab(
  list(interaction_all = sr_lat_reg1,
       interaction_lin = sr_lat_reg2,
       no_interaction = sr_lat_reg3,
       no_change = sr_lat_reg4))




sr_lat_reg1_c <- predict(sr_lat_reg1,newdata = sr_lat_l[sr_lat_l$cat == 'current',])
sr_lat_reg1_f <- predict(sr_lat_reg1,newdata = sr_lat_l[sr_lat_l$cat == 'ssp585',])

sr_lat_reg3_c <- predict(sr_lat_reg3,newdata = sr_lat_l[sr_lat_l$cat == 'current',])
sr_lat_reg3_f <- predict(sr_lat_reg3,newdata = sr_lat_l[sr_lat_l$cat == 'ssp585',])

x11(width = 13,height = 11)
layout(matrix(1:2,nrow = 1))
plot(current~lat,data = sr_lat,xlab = 'Latitude',ylab = 'Species richness',pch = 20)
points(ssp585~lat,data = sr_lat)

points(exp(sr_lat_reg1_c)~lat,data = sr_lat,type = 'l',lty = 2,col = 3,lwd = 2)
points(exp(sr_lat_reg1_f)~lat,data = sr_lat,type = 'l',lty =3,col = 2,lwd = 2)

plot(current~lat,data = sr_lat,xlab = 'Latitude',ylab = 'Species richness',pch = 20)
points(ssp585~lat,data = sr_lat)

points(exp(sr_lat_reg3_c)~lat,data = sr_lat,type = 'l',lty = 2,col = 3,lwd = 2)
points(exp(sr_lat_reg3_f)~lat,data = sr_lat,type = 'l',lty =3,col = 2,lwd = 2)

legend(18,75,legend = c('Current','SSP5.85'),col = c(3,2),lty = c(2,3),
       lwd = c(2,3),bty = 'n')

legend(18.3,75,legend = c('',''),pch = c(20,1),bg = NULL,bty = 'n')

####@@elevational ranges and sprich@@####
####@@elevational ranges@@####
## reclass dem at 200m interval, and make polygons
temp_env <- new.env()
load("gis/rasters_3-10-2020.RData",envir = temp_env)
dem_merge <- get("dem_merge",envir = temp_env)
terrain_wg <- get("terrain_wg",temp_env)
rm(temp_env)

dem_crop <- raster::crop(dem_merge,as(extent(terrain_wg),"SpatialPolygons"))
rm(dem_merge)
rm(terrain_wg)

rcl_mat <- matrix(
  c(seq(0,2600,200),seq(200,2800,200),seq(200,2800,200)),
  ncol = 3,byrow = F)

colnames(rcl_mat) <- c("min",'max','val')

rcl_mat[1,1] <- -Inf

#dem_rcl <- raster::reclassify(dem_crop,rcl_mat,right = F)

#dem_rcl_wg <- raster::mask(dem_rcl,wg_bnd_merge)

writeRaster(dem_rcl_wg,"gis/dem_rcl_wg.tif")


#dem_rcl_wg_poly <- raster::rasterToPolygons(dem_rcl_wg,dissolve = T)
#this step is done in qgis for memory issues

library(sf)

dem_rcl_wg_poly <- st_read("gis/dem_rcl_wg/dem_rcl_wg_poly.shp")

dem_rcl_wg_poly <- as_Spatial(dem_rcl_wg_poly)

dem_rcl_wg_poly$elev <- as.numeric(dem_rcl_wg_poly$elev)

dem_rcl_wg_poly <- dem_rcl_wg_poly[order(dem_rcl_wg_poly$elev,decreasing = F),]


bin_current_area_lonlat <- projectRaster(bin_current_area,bin_current_stk)

bin_current_area_lonlat <- as(bin_current_area_lonlat,Class = "RasterStack")

####__current__####

area_zones_l <- area_lat %>%
  data.frame() %>% 
  # mutate(massif =  wg_bnd_merge_massif$MASSIF) %>% 
  pivot_longer(cols = 2:ncol(area_lat),names_to = "Binomial_species",
               values_to = "area")

##elevation range of each species within a massif is the elevations that cover 
## 95% of the predicted area

all_frogs_elev_zones <- lapply(grd_lat$ID,
                               function(x){
                                 
                                 m <- raster::crop(dem_rcl_wg_poly,grd_lat[
                                   grd_lat$ID == x,])
                                 
                                 c <- raster::crop(bin_current_area_lonlat,grd_lat[
                                   grd_lat$ID == x,])
                                 
                                 ex <- raster::extract(c,m,fun = sum,na.rm = T,df = T)
                                 
                                 ex <- ex[,!colSums(ex) == 0]
                                 
                                 d <- do.call(rbind,lapply(ex,
                                                           function(y){
                                                             per <- y/sum(y)
                                                             min <- min(which(cumsum(per)>=0.05))
                                                             max <- min(which(cumsum(rev(per))>=0.05))
                                                             
                                                             c(min = m$elev[min],max = rev(m$elev)[max])
                                                           }))
                                 d <- d[!rownames(d) == "ID",]
                                 
                                 spmat <- apply(d,1,
                                                function(x){
                                                  v <- rep(0,length(m$elev))
                                                  v[which(m$elev == x[1]):
                                                      which(m$elev == x[2])] <- 1
                                                  v
                                                })
                                 rownames(spmat) <- m$elev
                                 list(spmat = spmat,elev = d)
                               })

names(all_frogs_elev_zones) <- grd_lat$ID

#calculate species richness
all_frogs_sprich_zones <- 
  do.call("c",
          lapply(all_frogs_elev_zones,
                 function(x) rowSums(x$spmat)))

t <- do.call(rbind,strsplit(names(all_frogs_sprich_zones),"[.]"))
colnames(t) <- c("massif","elevation")

all_frogs_sprich_zones <- data.frame(t,
                                     sprich = all_frogs_sprich_zones)

all_frogs_sprich_zones$elevation <- as.numeric(all_frogs_sprich_zones$elevation)

all_frogs_sprich_zones$massif <- factor(all_frogs_sprich_zones$massif,
                                        area_lat$ID)

all_frogs_sprich_zones$elevation <- as.numeric(all_frogs_sprich_zones$elevation)


#elevational extents

all_frogs_elevrange_zones <- do.call(rbind,lapply(all_frogs_elev_zones,
                                                  function(x) data.frame(x$elev)))

all_frogs_elevrange_zones$massif <- do.call(rbind,
                                            strsplit(rownames(all_frogs_elevrange_zones),
                                                     "[.]"))[,1]
all_frogs_elevrange_zones$range <- all_frogs_elevrange_zones$max - 
  all_frogs_elevrange_zones$min

all_frogs_elevrange_zones$scnr <- "current"

all_frogs_elevrange_zones$massif <- factor(all_frogs_elevrange_zones$massif,
                                           area_lat$ID)
####ssp585####
bin_ssp585_area_lonlat <- projectRaster(bin_ssp585_area,bin_current_stk)

bin_ssp585_area_lonlat <- as(bin_ssp585_area_lonlat,Class = "RasterStack")

ssp585_area_zones_l <- ssp585_area_lat %>%
  data.frame() %>% 
  # mutate(massif =  wg_bnd_merge_massif$MASSIF) %>% 
  pivot_longer(cols = 2:ncol(ssp585_area_lat),names_to = "Binomial_species",
               values_to = "area")

##elevation range of each species within a massif is the elevations that cover 
## 95% of the predicted area

ssp585_elev_zones <- lapply(grd_lat$ID,
                            function(x){
                              
                              m <- raster::crop(dem_rcl_wg_poly,grd_lat[
                                grd_lat$ID == x,])
                              
                              c <- raster::crop(bin_ssp585_area_lonlat,grd_lat[
                                grd_lat$ID == x,])
                              
                              ex <- raster::extract(c,m,fun = sum,na.rm = T,df = T)
                              
                              ex <- ex[,!colSums(ex) == 0]
                              
                              d <- do.call(rbind,lapply(ex,
                                                        function(y){
                                                          per <- y/sum(y)
                                                          min <- min(which(cumsum(per)>=0.05))
                                                          max <- min(which(cumsum(rev(per))>=0.05))
                                                          
                                                          c(min = m$elev[min],max = rev(m$elev)[max])
                                                        }))
                              d <- d[!rownames(d) == "ID",]
                              
                              spmat <- apply(d,1,
                                             function(x){
                                               v <- rep(0,length(m$elev))
                                               v[which(m$elev == x[1]):
                                                   which(m$elev == x[2])] <- 1
                                               v
                                             })
                              rownames(spmat) <- m$elev
                              list(spmat = spmat,elev = d)
                            })

names(ssp585_elev_zones) <- grd_lat$ID


#calculate species richness
ssp585_sprich_zones <- 
  do.call("c",
          lapply(ssp585_elev_zones,
                 function(x) rowSums(x$spmat)))

t <- do.call(rbind,strsplit(names(ssp585_sprich_zones),"[.]"))
colnames(t) <- c("massif","elevation")

ssp585_sprich_zones <- data.frame(t,
                                  sprich = ssp585_sprich_zones)

ssp585_sprich_zones$elevation <- as.numeric(ssp585_sprich_zones$elevation)

ssp585_sprich_zones$massif <- factor(ssp585_sprich_zones$massif,
                                     ssp585_area_lat$ID)

all_frogs_sprich_zones$ssp585_sprich <- ssp585_sprich_zones$sprich[
  match(rownames(all_frogs_sprich_zones),rownames(ssp585_sprich_zones))
]

#elevational extents

ssp585_elevrange_zones <- do.call(rbind,lapply(ssp585_elev_zones,
                                               function(x) data.frame(x$elev)))

ssp585_elevrange_zones$massif <- do.call(rbind,
                                         strsplit(rownames(ssp585_elevrange_zones),
                                                  "[.]"))[,1]
ssp585_elevrange_zones$range <- ssp585_elevrange_zones$max - 
  ssp585_elevrange_zones$min

ssp585_elevrange_zones$scnr <- "ssp585"


elev_range_l <- rbind(all_frogs_elevrange_zones,ssp585_elevrange_zones)

elev_range_l$massif <- factor(elev_range_l$massif,area_lat$ID)

seq <- match(rownames(all_frogs_elevrange_zones),rownames(ssp585_elevrange_zones))


elev_range_w <- data.frame(
  all_frogs_elevrange_zones[,!names(all_frogs_elevrange_zones)%in%"scnr"],
  ssp585_min = ssp585_elevrange_zones$min[seq],
  ssp585_max = ssp585_elevrange_zones$max[seq],
  ssp585_range = ssp585_elevrange_zones$range[seq])

elev_range_w$spp <- do.call(rbind,
                            strsplit(rownames(elev_range_w),
                                     "[.]"))[,2]

elev_range_w$area_curr <- area_zones_l$area[
  match(rownames(elev_range_w),
        paste(area_zones_l$ID,area_zones_l$Binomial_species,sep = "."))
]

elev_range_w$area_ssp585 <- ssp585_area_zones_l$area[
  match(rownames(elev_range_w),
        paste(ssp585_area_zones_l$ID,ssp585_area_zones_l$Binomial_species,
              sep = "."))]

elev_range_w$mid = elev_range_w$min + (elev_range_w$range/2)
elev_range_w$ssp585_mid = elev_range_w$ssp585_min + (elev_range_w$ssp585_range/2)

area_per <- (elev_range_w$area_ssp585/elev_range_w$area_curr)*100


elev_range_w$area_per <- area_per


##subset elev_range_w to keep only elev zones that contain 95%of species range
#for every species

d <- elev_range_w %>% 
  select(spp,massif,area_curr) %>% 
  group_by(spp) 


l <- lapply(unique(elev_range_w$spp),
            function(x){
              d <- elev_range_w[elev_range_w$spp == x,]
              d1 <- tapply(d$area_curr,d$massif,sum)
              
              area_per <- d1/sum(d1,na.rm = T)
              
              area_per[is.na(area_per)] <- 0
              
              min <- cumsum(area_per)>=0.05
              
              max <- rev(cumsum(rev(area_per))>=0.05)
              
              min&max
            })

names(l) <- unique(elev_range_w$spp)

l1 <- do.call(c,l)

l1_names <- sapply(strsplit(names(l1),"[.]"),
                   function(x) paste(x[2],x[1],sep = "."))

elev_range_w1 <- elev_range_w[l1[match(rownames(elev_range_w),l1_names)],]

elev_range_w1$area_per <- area_per[l1[match(rownames(elev_range_w),l1_names)]]

####aggregate elevation extents####
d1 <- elev_range_w %>% 
  group_by(spp) %>% 
  summarise(min_wt_curr = weighted.mean(min,area_curr/sum(area_curr)),
            max_wt_curr = weighted.mean(max,area_curr/sum(area_curr)),
            mid_wt_curr = weighted.mean(mid,area_curr/sum(area_curr)),
            mid_wt_ssp585 = weighted.mean(ssp585_mid,
                                          area_ssp585/sum(area_ssp585,
                                                          na.rm = T),
                                          na.rm = T),
            min_wt_ssp585 = weighted.mean(ssp585_min,
                                          area_ssp585/sum(area_ssp585,na.rm = T),
                                          na.rm = T),
            max_wt_ssp585 = weighted.mean(ssp585_max,
                                          area_ssp585/sum(area_ssp585,na.rm = T),
                                          na.rm = T),
            area_curr = sum(area_curr,na.rm = T),
            area_ssp585 = sum(area_ssp585,na.rm = T)
  )

d1 <- data.frame(d1)
rownames(d1) <- d1$spp

current_elev_spmat <- apply(d1,1,
                            function(x){
                              x <- as.numeric(x[2:3])
                              v <- rep(0,length(dem_rcl_wg_poly$elev))
                              v[max(which(dem_rcl_wg_poly$elev <= x[1])):
                                  min(which(dem_rcl_wg_poly$elev >= x[2]))] <- 1
                              v
                            })
rownames(current_elev_spmat) <- dem_rcl_wg_poly$elev

sr_current_elev <- rowSums(current_elev_spmat)

ssp585_elev_spmat <- apply(d1[complete.cases(d1),],1,
                           function(x){
                             x <- as.numeric(x[6:7])
                             v <- rep(0,length(dem_rcl_wg_poly$elev))
                             v[max(which(dem_rcl_wg_poly$elev <= x[1])):
                                 min(which(dem_rcl_wg_poly$elev >= x[2]))] <- 1
                             v
                           })
rownames(ssp585_elev_spmat) <- dem_rcl_wg_poly$elev

sr_ssp585_elev <- rowSums(ssp585_elev_spmat)

sr_elev <- data.frame(elev = dem_rcl_wg_poly$elev,
                      current = sr_current_elev,
                      ssp585 = sr_ssp585_elev)




#### regression for elevation sprich####
#sr_elev$elev_scl <- scale(sr_elev$elev)

sr_elev_l <- sr_elev %>% 
  pivot_longer(c(current,ssp585),names_to = 'cat',values_to = 'richness')

sr_elev_reg1 <- glm(richness ~ elev*cat+I(elev^2)*cat,data = sr_elev_l,
                    family = poisson)

summary(sr_elev_reg1) # AIC: 186.4

sr_elev_reg1_c <- predict(sr_elev_reg1,newdata = sr_elev_l[sr_elev_l$cat == 'current',])
sr_elev_reg1_f <- predict(sr_elev_reg1,newdata = sr_elev_l[sr_elev_l$cat == 'ssp585',])


plot(current ~ elev,data = sr_elev,pch = 20,xlab = 'Elevation',
     ylab = 'Species richness' ,type = 'n',bty = 'n')
points(current ~ elev,data = sr_elev,pch = 20,xlab = 'Elevation',
       ylab = 'Species richness',frame.plot = F,xlim = c(0,2900),xaxp = c(0,2800,10))
points(sr_ssp585_elev ~ elev,data = sr_elev,pch = 0)

points(exp(sr_elev_reg1_c)~ elev,data = sr_elev,type = 'l',col = 3,lty = 2,lwd = 2)
points(exp(sr_elev_reg1_f)~ elev,data = sr_elev,type = 'l',col = 2,lty = 3,lwd = 3)

legend(1960,150,legend = c('current','SSP5.85'),col = c(3,2),lty = c(2,3),
       lwd = c(2,3),bty = 'n')
legend(1960,150,legend = c('',''),pch = c(20,0),bg = NULL,bty = 'n')

#drop interaction for quadratic term
sr_elev_reg2 <- glm(richness ~ elev*cat+I(elev^2),data = sr_elev_l,
                    family = poisson)

summary(sr_elev_reg2)# AIC: 190.41

sr_elev_reg2_c <- predict(sr_elev_reg2,newdata = sr_elev_l[sr_elev_l$cat == 'current',])
sr_elev_reg2_f <- predict(sr_elev_reg2,newdata = sr_elev_l[sr_elev_l$cat == 'ssp585',])

plot(current ~ elev,data = sr_elev,pch = 20,xlab = 'Elevation',
     ylab = 'Species richness' ,type = 'n',bty = 'n')
points(current ~ elev,data = sr_elev,pch = 20,xlab = 'Elevation',
       ylab = 'Species richness',frame.plot = F,xlim = c(0,2900),xaxp = c(0,2800,10))
points(sr_ssp585_elev ~ elev,data = sr_elev,pch = 0)

points(exp(sr_elev_reg2_c)~ elev,data = sr_elev,type = 'l',col = 3,lty = 2,lwd = 2)
points(exp(sr_elev_reg2_f)~ elev,data = sr_elev,type = 'l',col = 2,lty = 3,lwd = 3)


sr_elev_reg3 <- glm(richness ~ cat+elev+I(elev^2),data = sr_elev_l,
                    family = poisson)

summary(sr_elev_reg3)# AIC: 209.43

sr_elev_reg3_c <- predict(sr_elev_reg3,newdata = sr_elev_l[sr_elev_l$cat == 'current',])
sr_elev_reg3_f <- predict(sr_elev_reg3,newdata = sr_elev_l[sr_elev_l$cat == 'ssp585',])


plot(current ~ elev,data = sr_elev,pch = 20,xlab = 'Elevation',
     ylab = 'Species richness' ,type = 'n',bty = 'n')
points(current ~ elev,data = sr_elev,pch = 20,xlab = 'Elevation',
       ylab = 'Species richness',frame.plot = F,xlim = c(0,2900),xaxp = c(0,2800,10))
points(sr_ssp585_elev ~ elev,data = sr_elev,pch = 0)

points(exp(sr_elev_reg3_c)~ elev,data = sr_elev,type = 'l',col = 3,lty = 2,lwd = 2)
points(exp(sr_elev_reg3_f)~ elev,data = sr_elev,type = 'l',col = 2,lty = 3,lwd = 3)

sr_elev_reg4 <- glm(richness ~ elev+I(elev^2),data = sr_elev_l,
                    family = poisson)

aictab_elevation <- aictab(
  list(interaction_all = sr_elev_reg1,
       interaction_lin = sr_elev_reg2,
       no_interaction = sr_elev_reg3,
       interaction_q = glm(richness ~ elev+I(elev^2)*cat,data = sr_elev_l,
                           family = poisson),
       no_change = sr_elev_reg4))

save.image("sdm_v1/frogs/all_frogs_enmeval_lat_elev_ranges_v1.RData")