## compare present and future areas for all species

# The main outputs from this script are:
#   1. area_dat: data frame with current area masked within WG and within 
#               biogeo boudaries
#   2. total_area: data frame with current and future area masked within WG
# and plots with 2. 
#   3. Rasters of high diversity  areas (more than 50% of the regional pool)

## get current distributions and area 

temp_env <- new.env()
load("sdm_v1/frogs/all_frogs_enmeval_biased_unbiased_compare_13-6-2023.RData",envir = temp_env)
thresholds_out <- get("thresholds_out",envir = temp_env)
area_dat <- get("area_dat",envir = temp_env)
rm(temp_env)

## get future distribution and area
temp_env <- new.env()
load("sdm_v1/frogs/all_frogs_enmeval_future_predictions_mask_biogeozones_v1_14-6-2023.RData",envir = temp_env)
bin_2070_zones_area <- get("bin_2070_zones_area",envir = temp_env)
rm(temp_env)

## get masked outputs of current distribution
temp_env <- new.env()
load("all_frogs_enmeval_mask_biogeozones_v1_14-6-2024.RData",envir = temp_env)
bin_current_mask <- get("bin_current_mask",envir = temp_env)
all_frogs_zones_clip_mask <- get("all_frogs_zones_clip_mask",envir = temp_env)
all_frogs_zones_mask <- get("all_frogs_zones_mask",envir = temp_env)
rm(temp_env)


####current area####

bin_current_area <- sapply(bin_current_mask,
                           function(x){
                             
                             x <- projectRaster(x,crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
                             
                             
                             raster::cellStats(x,stat = 'sum',na.rm = T)*0.801520
                           })

bin_mask_area <- sapply(all_frogs_zones_mask,
                        function(x){
                          x <- x$out
                          x <- projectRaster(x,crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
                          
                          
                          raster::cellStats(x,stat = 'sum',na.rm = T)*0.801520
                        })

bin_clip_area <- sapply(all_frogs_zones_clip_mask,
                        function(x){
                          x <- x$out
                          x <- projectRaster(x,crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
                          
                          raster::cellStats(x,stat = 'sum',na.rm = T)*0.801520
                        })

area_dat <- data.frame(area_dat[,c(1,2,4,6)],
                       bin_current = bin_current_area,
                       bin_mask  = bin_mask_area,
                       bin_clip = bin_clip_area
)

###make data frame for total area in current and future

total_area <- data.frame(bin_2070_zones_area[,1:8],
                         area_current = area_dat$bin_current[
                           match(bin_2070_zones_area$Binomial_species,
                                 area_dat$Binomial_species)])

total_area$diff <- total_area$area_current - total_area$area
total_area$per <- (total_area$area/total_area$area_current)*100

summary(total_area$per[total_area$cat == "max"&!total_area$per>100])

# transpose for plotting
total_area_w <- total_area %>% 
  pivot_wider(id_cols = c(Family,Binomial_species,time,),id_expand = F,names_from = cat,values_from = area) 

total_area_w$area_current <- area_dat$bin_current[match(total_area_w$Binomial_species,area_dat$Binomial_species)]

total_area_w$min_per <- (total_area_w$min/total_area_w$area_current)*100
total_area_w$max_per <- (total_area_w$max/total_area_w$area_current)*100

total_area_sort <- total_area_w[order(total_area_w$area_current,decreasing = T),]

total_area_sort$cat <- cut(total_area_sort$area_current,
                           breaks = c(0,1429,5123,15571,79706),labels = 1:4)

total_area_sort$cat1 <- 
  ifelse(total_area_sort$min<total_area_sort$area_current&
           total_area_sort$max<total_area_sort$area_current,"decrease","increase")

## make plots

ggplot(total_area_sort)+
  geom_segment(aes(x = area_current,
                   xend = area_current,
                   y = min,
                   yend = max,colour = cat1))+
  geom_abline(aes(slope = 1,intercept = 0))+
  facet_wrap("cat",scales = "free")

ggplot(total_area_sort)+
  geom_segment(aes(x = area_current,
                   xend = area_current,
                   y = min,
                   yend = max),col = "blue",
               data = total_area_sort[total_area_sort$min<total_area_sort$area_current&
                                        total_area_sort$max<total_area_sort$area_current,])+
  geom_segment(aes(x = area_current,
                   xend = area_current,
                   y = min,
                   yend = max),col = "red",
               data = total_area_sort[total_area_sort$min>total_area_sort$area_current&
                                        total_area_sort$max>total_area_sort$area_current,])+
  geom_abline(aes(slope = 1,intercept = 0))+
  facet_wrap("cat")


total_area_w %>% 
  filter(max_per>min_per&min_per<100) %>% 
  mutate(diff = max-min) %>% 
  View()


plot(1:length(total_area_sort$max),
     total_area_sort$max,type = "n",ylab = "area",xlab = "rank")


with(total_area_sort,segments(1:nrow(total_area_sort),min,
                              1:nrow(total_area_sort),max,lwd = 1.2))

points(total_area_sort$area_current,cex = 0.5,pch = 16,col = "red")
points(x = which(total_area_sort$max>total_area_sort$min),
       y = total_area_sort$max[which(total_area_sort$max>total_area_sort$min)],cex = 0.5,
       pch = 16,col = "blue")

jpeg("analysis/present_future_area.jpg",width = 240,height = 130,units = "mm",quality = 100,res = 300)
par_dflt <- par()
par(mar=c(5, 5, 4,8))
plot.new()
{plot(total_area_w$area_current,total_area_w$area_current,ylim = c(0,100000),xlim = c(0,100000),type = "n",axes = F,frame.plot = T,
      xlab = expression(Current~area~(Km^2)),ylab = expression(Future~area~(Km^2)))
  axis(1,at = seq(10000,100000,10000),labels = seq(10000,100000,10000),tick = T)
  axis(2,at = seq(10000,100000,10000),labels = seq(10000,100000,10000),tick = T)
  abline(0,1,lwd = 2)
  with(total_area_w[total_area_w$min<total_area_w$area_current&
                      total_area_w$max<total_area_w$area_current,],
       segments(area_current,min,
                area_current,max,lwd = 2,col  ="blue"))
  with(total_area_w[total_area_w$min>total_area_w$area_current|
                      total_area_w$max>total_area_w$area_current,],
       segments(area_current,min,
                area_current,max,lwd = 2,lty = 4,col = "red"))
  with(total_area_w[total_area_w$max == 0,],
       points(area_current,max,pch = 4))
  par(xpd = T)
  legend(110000,90000,legend = c("Increase","Decrease"),lty = c(4,1),col = c("red","blue"),seg.len = 1)
  par(par_dflt)
}
dev.off()


####::high diversity areas::####

bin_current_stk <- stack(bin_current_mask)

bin_current_sum <- stackApply(bin_current_stk,rep(1,nlayers(bin_current_stk)),sum)

bin_current_sum <- bin_current_sum/nlayers(bin_current_stk)

bin_current_scl <- bin_current_sum/max(getValues(bin_current_sum))


writeRaster(bin_current_scl>=0.5,"E:/backups/all_frogs_enmeval_18-4-2023//high_diversity_current_18-4-2024.asc",overwrite = T)

bin_current_scl[bin_current_scl<0.5] <- NA

sum(getValues(raster::area(bin_current_scl,na.rm = T)),na.rm = T)


bin_ssp126_stk <- stack(
  lapply(thresholds_out,function(d){
    "sp126min_2061-2080_mask"
    a <- unlist(strsplit(d$bin_max@file@name,"[\\]"))
    b <- paste(a[1:5],collapse = "\\")
    
    get_pred(b,target = "ssp126min_2061-2080_mask",recursive = T)
    
  })
)

bin_ssp126_sum <- stackApply(bin_ssp126_stk,rep(1,nlayers(bin_ssp126_stk)),sum)

bin_ssp126_sum <- bin_ssp126_sum/nlayers(bin_ssp126_stk)

bin_ssp126_scl <- bin_ssp126_sum/max(getValues(bin_ssp126_sum))

writeRaster(bin_ssp126_scl>=0.5,"E:/backups/all_frogs_enmeval_18-4-2023//high_diversity_ssp126_18-4-2024.asc",overwrite = T)

bin_ssp126_scl[bin_ssp126_scl<0.5] <- NA

sum(getValues(raster::area(bin_ssp126_scl,na.rm = T)),na.rm = T)


bin_ssp585_stk <- stack(stack(
  lapply(thresholds_out,function(d){
    
    a <- unlist(strsplit(d$bin_max@file@name,"[\\]"))
    b <- paste(a[1:5],collapse = "\\")
    
    get_pred(b,target = "ssp585max_2061-2080_mask",recursive = T)
    
  })
))

bin_ssp585_sum <- stackApply(bin_ssp585_stk,rep(1,nlayers(bin_ssp585_stk)),sum)

bin_ssp585_sum <- bin_ssp126_sum/nlayers(bin_ssp585_stk)

bin_ssp585_scl <- bin_ssp585_sum/max(getValues(bin_ssp585_sum))

writeRaster(bin_ssp585_scl>=0.5,"E:/backups/all_frogs_enmeval_18-4-2023//high_diversity_ssp585_18-4-2024.asc",overwrite = T)


bin_ssp585_scl[bin_ssp585_scl<0.5] <- NA

sum(getValues(raster::area(bin_ssp585_scl,na.rm = T)),na.rm = T)

####save.image####
#save.image("sdm_v1/frogs/all_frogs_enmeval_future_prediction_mask_biogeozones_v1_14-6-2024.RData")

