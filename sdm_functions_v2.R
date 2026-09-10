## get and load files like rds
f1 <- function(x) {
  a <- load(x)
  get(a)
}

####@@ENMeval@@####

enmeval_predict <- function(mod,preds,thresh,outname){
  predict.sel<-dismo::predict(mod, preds, args="outputformat=cloglog")
  writeRaster(predict.sel, filename= paste0(spp_path,"/",outname),
              format="GTiff", overwrite=TRUE)
  
  predict.bin <- raster::reclassify(predict.sel,
                                    c(-Inf,thresh,NA,thresh,Inf,1),right = F)
  writeRaster(predict.bin,filename= paste0(spp_path,"/",outname,"_ses_bin"),
              format="GTiff", overwrite=TRUE)
}


enmeval_fit <- function(x,path,back_swd,bias){
  #x: character.species name
  #path: character.output folder
  #back_swd: either a list or data.frame if bias is T, then a named list with Family as names, and each 
  #           element a swd background points for respective family
  #bias: logical. If T the background points should be bias corrected
  spp <- x
  i <- which(maxent_swd_occ5$Binomial_species == spp)
  fam <- unique(maxent_swd_occ5$Family[i])
  path_fam <- paste0(path,"/",fam)
  path_spp <- paste0(path,"/",fam,"/",spp)
  
  occ_dat <- maxent_swd_occ5[maxent_swd_occ5$Binomial_species == spp,]
  if(bias){
    back_dat <- back_swd[[fam]]
  }else{back_dat <- back_swd}
  if(!dir.exists(path_fam)){
    dir.create(path_fam)
  }
  
  if(!dir.exists(path_spp)){
    dir.create(path_spp)
  }
  mod_out <- ENMevaluate(occs = occ_dat[,-c(1:3)],bg = back_dat[,-1],
                         tune.args = list(fc = c("L","Q","H","LQ","LH","QH","LQH"),
                                          rm = seq(0.5,5,0.5)),
                         partitions = ifelse(nrow(occ_dat)<10,"jackknife","block"),
                         algorithm = "maxent.jar",
                         doClamp = T,
                         taxon.name = spp
  )
  
  save(mod_out,file = paste0(path_spp,"/",spp,"_bioclim","_enmeval"))
  write.csv(mod_out@results,file = paste0(path_spp,"/",spp,"_bioclim","_res.csv"))
  
  
}



area_lat_elev_ext <- function(dsn_orig,dsn_dest,mod_args = NA){
  #load enmeval object with biased background
  #select model with maximum number of features, and default regularization (rm = 1)
  #predict the model and use MSS threshold to make binary
  #mask to WG 
  
  ##list files with enmeval objects
  l1 <- list.files(dsn_orig,pattern = "_enmeval$",recursive = T)
  #apply function to make predictions and save files
  lapply(l1,load_predict_mask,dsn_dest = dsn_dest,dsn_orig = dsn_orig,mod_args = mod_args)
  
  ## Calculate areas of binary predictions
  
  cmplx_withbias_curr_mask <- fetch_pred(dsn_dest,target = "_mss_mask.asc$",recursive = T)
  ## project rasters to equal area projection and calculate area
  #current
  cmplx_withbias_curr_mask_aea <- projectRaster(stack(cmplx_withbias_curr_mask),crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  
  cmplx_withbias_curr_mask_aea <- cmplx_withbias_curr_mask_aea*0.801520
  
  cmplx_withbias_curr_area <- raster::cellStats(cmplx_withbias_curr_mask_aea,
                                                stat = 'sum',na.rm = T)
  
  #ssp126
  cmplx_withbias_ssp126_mask <- fetch_pred(dsn_dest,target = "ssp126min_2061-2080_mask.asc",recursive = T)
  
  cmplx_withbias_ssp126_mask_aea <- projectRaster(stack(cmplx_withbias_ssp126_mask),crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  
  cmplx_withbias_ssp126_mask_aea <- cmplx_withbias_ssp126_mask_aea*0.801520
  
  cmplx_withbias_ssp126_area <- raster::cellStats(cmplx_withbias_ssp126_mask_aea,stat = 'sum',na.rm = T)
  
  #ssp585
  cmplx_withbias_ssp585_mask <- fetch_pred(dsn_dest,target = "ssp585max_2061-2080_mask.asc",recursive = T)
  
  cmplx_withbias_ssp585_mask_aea <- projectRaster(stack(cmplx_withbias_ssp585_mask),crs = crs("+proj=aea +lat_0=-15 +lon_0=125 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  
  cmplx_withbias_ssp585_mask_aea <- cmplx_withbias_ssp585_mask_aea*0.801520
  
  cmplx_withbias_ssp585_area <- raster::cellStats(cmplx_withbias_ssp585_mask_aea,stat = 'sum',na.rm = T)
  
  area_dat <- data.frame(current = cmplx_withbias_curr_area,
                         ssp126 = cmplx_withbias_ssp126_area,
                         ssp585 = cmplx_withbias_ssp585_area)
  
  
  ## latitude extents
  cmplx_withbias_curr_latext <- data.frame(do.call(rbind,lapply(cmplx_withbias_curr_mask,
                                                                function(b){
                                                                  tryCatch(
                                                                    {b1 <- trim(b)
                                                                    b2 <- extent(b1)
                                                                    c(lat_min = b2@ymin,lat_max = b2@ymax)},
                                                                    error = function(e){
                                                                      c(lat_min = NA,lat_max = NA)
                                                                      
                                                                    })
                                                                }))
  )
  names(cmplx_withbias_curr_latext) <- paste0(names(cmplx_withbias_curr_latext),"_curr")
  
  cmplx_withbias_ssp585_latext <- data.frame(do.call(rbind,lapply(cmplx_withbias_ssp585_mask,
                                                                  function(b){
                                                                    tryCatch(
                                                                      {b1 <- trim(b)
                                                                      b2 <- extent(b1)
                                                                      c(lat_min = b2@ymin,lat_max = b2@ymax)},
                                                                      error = function(e){
                                                                        c(lat_min = NA,lat_max = NA)
                                                                        
                                                                      })
                                                                  }))
  )
  names(cmplx_withbias_ssp585_latext) <- paste0(names(cmplx_withbias_ssp585_latext),"_ssp585")
  
  lat_ext <- data.frame(cmplx_withbias_curr_latext,
                        cmplx_withbias_ssp585_latext)
  
  ##elevational extents
  ####current####
  ## calculate current area within each latitude bin
  cmplx_withbias_curr_area_lat <- raster::extract(cmplx_withbias_curr_mask_aea,
                                                  grd_lat_aea,fun = sum,na.rm = T,df = T)
  
  cmplx_withbias_curr_area_lat <- cmplx_withbias_curr_area_lat[order(grd_lat_aea$lat,decreasing = T),]
  
  cmplx_withbias_curr_area_lat_l <- cmplx_withbias_curr_area_lat %>%
    data.frame() %>% 
    # mutate(massif =  wg_bnd_merge_massif$MASSIF) %>% 
    pivot_longer(cols = 2:ncol(.),names_to = "Binomial_species",
                 values_to = "area")
  
  cmplx_withbias_curr_elevext <- get_elev_ext(grd_lat = grd_lat,
                                              raster_stack = stack(cmplx_withbias_curr_mask),
                                              area_wt = cmplx_withbias_curr_area_lat_l,
                                              scnr = "current")
  
  names(cmplx_withbias_curr_elevext) <- paste0(names(cmplx_withbias_curr_elevext),"_curr")
  ####ssp585####
  
  ## calculate current area within each latitude bin
  cmplx_withbias_ssp585_area_lat <- raster::extract(cmplx_withbias_ssp585_mask_aea,
                                                    grd_lat_aea,fun = sum,na.rm = T,df = T)
  
  cmplx_withbias_ssp585_area_lat <- cmplx_withbias_ssp585_area_lat[order(grd_lat_aea$lat,decreasing = T),]
  
  cmplx_withbias_ssp585_area_lat_l <- cmplx_withbias_ssp585_area_lat %>%
    data.frame() %>% 
    # mutate(massif =  wg_bnd_merge_massif$MASSIF) %>% 
    pivot_longer(cols = 2:ncol(.),names_to = "Binomial_species",
                 values_to = "area")
  
  ## elevation extents
  cmplx_withbias_ssp585_elevext <- get_elev_ext(grd_lat = grd_lat,
                                                raster_stack = stack(cmplx_withbias_ssp585_mask),
                                                area_wt = cmplx_withbias_ssp585_area_lat_l,
                                                scnr = "ssp585")
  names(cmplx_withbias_ssp585_elevext) <- paste0(names(cmplx_withbias_ssp585_elevext),"_ssp585")
  
  elev_ext <- data.frame(cmplx_withbias_curr_elevext,cmplx_withbias_ssp585_elevext[,!names(cmplx_withbias_ssp585_elevext)%in%"spp_ssp585"])
  
  list(area_dat = area_dat,lat_ext = lat_ext,elev_ext = elev_ext)
}



load_predict_mask <- function(a,dsn_dest,dsn_orig,mod_args){
  #load enmeval object with biased background
  #select model with maximum number of features, and default regularization (rm = 1)
  #predict the model and use MSS threshold to make binary
  #mask to WG 
  
  #make file name, paths ect
  a1 <- unlist(strsplit(a,"/"))
  fam <- a1[1]
  spp <- a1[2]
  file <- a1[3]
  target_dir <- paste0(dsn_dest,fam,"/",spp)
  if(!file.exists(paste0(target_dir,"/",paste0(spp,"_ensemble_ssp585max_2061-2080_mask.asc")))){
  #create dir to save models
  if(!dir.exists(target_dir)){
    dir.create(target_dir,recursive = T)}
  #get model object
  a2 <- f1(paste0(dsn_orig,"/",a))
  if(is.na(mod_args)){
    sel <- a2@results %>% 
      filter(auc.train > 0.6) %>%
      filter(or.mtp.avg == min(or.mtp.avg,na.rm = T)) %>%
      filter(auc.diff.avg == min(auc.diff.avg,na.rm = T)) %>%
      filter(ncoef == min(ncoef,na.rm = T)) %>% 
      filter(rm == max(as.numeric(as.character(rm)))) %>% 
      mutate(n_f = sapply(strsplit(as.character(fc),""),length)) %>% 
      filter(n_f == min(n_f,na.rm = T)) %>% 
      data.frame()
    mod_args <- as.character(sel$tune.args)
    mod_sel_id <- which(names(a2@models)%in%sel$tune.args)
    mod_sel <- a2@models[[mod_sel_id]]
    #save selected model
    save(mod_sel,file = paste0(target_dir,"/",spp,"_",mod_args,"_sltd"))
    mod_mss <- mod_sel@results[thresh,]
  }else{
    mod_sel <- a2@models[[mod_args]]
    #save selected model
    save(mod_sel,file = paste0(target_dir,"/",spp,"_",mod_args,"_sltd"))
    mod_mss <- mod_sel@results[thresh,]
  }
  
  l1 <- mod_outputs(mod_sel)
  dir.create(paste0(target_dir,"/plots/"))
  
  ##create model plots
  jpeg(file=paste0(target_dir,"/plots/",paste(spp,"bioclim_preds","percon.jpeg",
                                            sep = "_")), height=8, width = 11, units = 'in', res=300)
  print(l1$p_con)
  dev.off()
  
  jpeg(file=paste0(target_dir,"/plots/",paste(spp,"bioclim_preds","perimp.jpeg",
                                            sep = "_")), height=8, width = 11, units = 'in', res=300)
  print(l1$p_imp)
  dev.off()
  
  jpeg(file=paste0(target_dir,"/plots/",paste(spp,"bioclim_preds","response.jpeg",
                                            sep = "_")), height=8, width = 11, units = 'in', res=300)
  dismo::response(mod_sel)
  dev.off()
  

  # create file name for output
  outname <- paste0(target_dir,"/",paste0(spp,"_",mod_args,"_bioclim_preds_mss_mask.asc"))
  # it output exists, load it as raster
  if(file.exists(outname)){
    mod_pred <- raster(outname)
  }else{# else predict model, make binary with threshold, and mask within wg
    mod_pred <- dismo::predict(mod_sel,preds,args="outputformat=cloglog")
    
    mod_pred_ses <- raster::reclassify(mod_pred,matrix(c(-Inf,mod_mss,NA,
                                                         mod_mss,Inf,1),
                                                       nrow = 2,byrow = T),right = F)
    mod_pred_ses_mask <- raster::mask(mod_pred_ses,wg_bnd_merge)
    mod_pred_mask <- raster::mask(mod_pred,wg_bnd_merge)
    
    raster::writeRaster(mod_pred_ses_mask,outname)
    raster::writeRaster(mod_pred,
                        paste0(target_dir,"/",paste(spp,mod_args,"bioclim_preds",'avg_unmasked.asc',
                                                    sep = "_")),overwrite = T)
    raster::writeRaster(mod_pred,
                paste0(target_dir,"/",paste(spp,mod_args,"bioclim_preds",'avg.asc',
                                          sep = "_")),overwrite = T)
  }
  
  ## predict model to future conditions
  #ssp126
  r1 <- dismo::predict(mod_sel,preds_ssp126_1,args="outputformat=cloglog")
  r2 <- dismo::predict(mod_sel,preds_ssp126_2,args="outputformat=cloglog")
  r3 <- dismo::predict(mod_sel,preds_ssp126_3,args="outputformat=cloglog")
  
  s <- stack(r1,r2,r3)
  r4 <- stackApply(s,c(1,1,1),median,na.rm = T)
  r5 <- raster::reclassify(r4,matrix(c(-Inf,mod_mss,NA,
                                       mod_mss,Inf,1),
                                     nrow = 2,byrow = T),right = F)
  
  outname_ssp126 <- paste0(target_dir,"/",paste0(spp,"_ensemble_ssp126min_2061-2080_avg.asc"))
  outname_ssp126_bin <- paste0(target_dir,"/",paste0(spp,"_ensemble_ssp126min_2061-2080_mask.asc"))
  
  raster::writeRaster(r4,outname_ssp126,overwrite = T)
  raster::writeRaster(r5,outname_ssp126_bin,overwrite = T)
  
  #ssp585
  r1 <- dismo::predict(mod_sel,preds_ssp585_1,args="outputformat=cloglog")
  r2 <- dismo::predict(mod_sel,preds_ssp585_2,args="outputformat=cloglog")
  r3 <- dismo::predict(mod_sel,preds_ssp585_3,args="outputformat=cloglog")
  
  s <- stack(r1,r2,r3)
  r4 <- stackApply(s,c(1,1,1),median,na.rm = T)
  r5 <- raster::reclassify(r4,matrix(c(-Inf,mod_mss,NA,
                                       mod_mss,Inf,1),
                                     nrow = 2,byrow = T),right = F)
  
  outname_ssp585 <- paste0(target_dir,"/",paste0(spp,"_ensemble_ssp585max_2061-2080_avg.asc"))
  outname_ssp585_bin <- paste0(target_dir,"/",paste0(spp,"_ensemble_ssp585max_2061-2080_mask.asc"))
  
  raster::writeRaster(r4,outname_ssp585,overwrite = T)
  raster::writeRaster(r5,outname_ssp585_bin,overwrite = T)
  }
}




##function for elevational extents of species 


get_elev_ext <- function(grd_lat,raster_stack,area_wt,scnr){
  all_frogs_elev_zones <- lapply(grd_lat$ID,
                                 function(x){
                                   c <- stack(raster::crop(raster_stack,grd_lat[
                                     grd_lat$ID == x,]))
                                   dem1 <- raster::crop(dem_wg,c) #dem_wg is loaded in the script file, it is only used here, so did not want to include it in the arguments. 
                                   d <-  do.call(rbind,lapply(c@layers,
                                                              function(c1){
                                                                
                                                                cp <- rasterToPoints(c1,spatial = T)
                                                                cp_ex <- raster::extract(dem_wg,cp)
                                                                c(min = min(cp_ex,na.rm = T),
                                                                  max = max(cp_ex,na.rm = T))
                                                              }
                                                              
                                   ))
                                   d <- data.frame(d)
                                   rownames(d) <- names(c)
                                   list(elev = d)
                                 })
  
  names(all_frogs_elev_zones) <- grd_lat$ID
  
  #elevational extents
  all_frogs_elevrange_zones <- do.call(rbind,lapply(all_frogs_elev_zones,
                                                    function(x) data.frame(x$elev)))
  all_frogs_elevrange_zones$mid <- all_frogs_elevrange_zones$min+((all_frogs_elevrange_zones$max - all_frogs_elevrange_zones$min)/2)
  all_frogs_elevrange_zones$massif <- do.call(rbind,
                                              strsplit(rownames(all_frogs_elevrange_zones),
                                                       "[.]"))[,1]
  all_frogs_elevrange_zones$spp <- do.call(rbind,
                                           strsplit(rownames(all_frogs_elevrange_zones),
                                                    "[.]"))[,2]
  all_frogs_elevrange_zones$scnr <- scnr
  
  all_frogs_elevrange_zones$massif <- factor(all_frogs_elevrange_zones$massif,
                                             grd_lat$ID)
  
  all_frogs_elevrange_zones$area <- area_wt$area[
    match(rownames(all_frogs_elevrange_zones),
          paste(area_wt$ID,area_wt$Binomial_species,sep = "."))
  ]
  
  all_frogs_elevrange_zones <- all_frogs_elevrange_zones[!(all_frogs_elevrange_zones$min == 'Inf'&
                                                             all_frogs_elevrange_zones$max == '-Inf'),]
  
  out <- all_frogs_elevrange_zones %>% 
    group_by(spp) %>% 
    summarise(min_wt = weighted.mean(min,area/sum(area)),
              max_wt = weighted.mean(max,area/sum(area)),
              mid_wt = weighted.mean(mid,area/sum(area)),
              area = sum(area,na.rm = T)) %>% 
    mutate(range = max_wt - min_wt)
  
}




mod_outputs <- function(mod_sel){
  
    ###variable contribution
    #plot variable percent contribution 
    #jpeg(file=paste0(x,"/out/",y,"/aic/",y,"_aic_varCon.jpeg"), height=8, width = 11, units = 'in', res=300)
    percon <- mod_sel@results[grep("contribution",rownames(mod_sel@results)),1]
    
    names(percon) <- do.call(rbind,strsplit(names(percon),"[.]"))[,1]
    
    d <- data.frame(var = names(percon),val = percon,stringsAsFactors = T)
    
    d <- d[order(d$val,decreasing = T),]
    
    d$var <- factor(d$var,levels = d$var)
    
    p_con <- ggplot(d)+
      aes(x = var,y = val)+
      geom_bar(stat ='identity' )+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90))+
      labs(x = "Predictor",y = "Percent contribution")
    # print(p_con)
    # dev.off()
    
    ### permutation importance
    percon <- mod_sel@results[grep("permutation.importance",rownames(mod_sel@results)),1]
    
    names(percon) <- do.call(rbind,strsplit(names(percon),"[.]"))[,1]
    
    d <- data.frame(var = names(percon),val = percon,stringsAsFactors = T)
    
    d <- d[order(d$val,decreasing = T),]
    
    d$var <- factor(d$var,levels = d$var)
    
    p_imp <- ggplot(d)+
      aes(x = var,y = val)+
      geom_bar(stat ='identity' )+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90))+
      labs(x = "Predictor",y = "permutation importance")
    
    list(p_con = p_con,p_imp = p_imp)
    
}

####@@maps and plots@@####
####biogeo zones comparison####
zones_compare_f2<- function(z){
  out1 <- zones_compare_f(z)
  grid.arrange(out1$map_masked,out1$map_clipped,out1$map_bins,nrow = 1)
  print(out1$inset_masked,vp = vp1)
  print(out1$inset_clipped,vp = vp2)
  print(out1$inset_bins,vp = vp3)
}

zones_compare_f <- function(x){
  i <- which(names(thresholds_out) == x)
  thresholds_out2 <- thresholds_out[[i]]
  raw <- thresholds_out2$bin_max
  raw <- terra::trim(rast(raw))
  raw_df <- as.data.frame(raw,xy = T,na.rm =T)
  
  all_frogs_zones_barrier2 <- all_frogs_zones_barrier[[i]]
  masked <- all_frogs_zones_barrier2$out
  masked <- terra::trim(rast(masked))
  masked_df <- as.data.frame(masked,xy = T,na.rm = T)
  masked_ext <- st_as_sfc(st_bbox(masked))
  
  all_frogs_zones_minelev2 <- all_frogs_zones_minelev[[i]]
  clipped <- all_frogs_zones_minelev2$out
  clipped <- terra::trim(rast(clipped))
  clipped_df <- as.data.frame(clipped,xy = T,na.rm = T)
  clipped_ext <- st_as_sfc(st_bbox(clipped))
  
  all_frogs_zones_bins2 <- all_frogs_zones_bins[[i]]
  bins <- all_frogs_zones_bins2$out
  bins <- terra::trim(rast(bins))
  bins_df <- as.data.frame(bins,xy = T,na.rm = T)
  bins_ext <- st_as_sfc(st_bbox(bins))
  
  pts <- thresholds_out[[x]]$pts
  
  masked_df1 <- rbind(data.frame(raw_df[,c('x','y')],var = 'out'),
                      data.frame(masked_df[,c('x','y')],var = 'in'))
  clipped_df1 <- rbind(data.frame(raw_df[,c('x','y')],var = 'out'),
                       data.frame(clipped_df[,c('x','y')],var = 'in'))
  bins_df1 <- rbind(data.frame(raw_df[,c('x','y')],var = 'out'),
                    data.frame(bins_df[,c('x','y')],var = 'in'))
  
  list(map_masked = in_out_map_fun1(masked_df1,"barrier",pts = pts),
       inset_masked = inset_fun1(masked_ext),
       map_clipped = in_out_map_fun1(clipped_df1,"minimum\nelevation",pts = pts)+
         labs(title = x),
       inset_clipped = inset_fun1(clipped_ext),
       map_bins = in_out_map_fun1(bins_df1,"latitude\nbins",pts = pts),
       inset_bins = inset_fun1(bins_ext)
  )
}

####plot biogeographical subset of sdm prediction#####
in_out_map_fun1 <- function(a,label,pts){
  ggplot()+
    geom_tile(data = a,aes(x = x, y = y, fill = var))+
    #geom_tile(data = a,aes(x = x, y = y),fill = "blue4")+
    geom_sf(data = st_as_sf(wg_bnd_merge_massif),col = "black",fill = NA)+
    geom_sf(data = pts,fill = "yellow",shape= 21,size = 1.5)+
    scale_fill_manual(values= c('out' = "coral4",'in' = "blue4"),name = label)+
    coord_sf(crs = 4326)+
    xlim(min(a$x),78.5)+
    ylim(min(a$y),max(a$y))+
    theme_classic()+
    theme(legend.position = c(0.9,0.8),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))
  
}

#### plot sdm prediction ####

bias_random_background_compare <- function(c){
  pred1 <- thresholds_out[[c]]$pred
  mask1 <- thresholds_out[[c]]$bin_max
  pts1 <- thresholds_out[[c]]$pts
  
  pred_nobias1 <- thresholds_out_nobias[[c]]$pred
  mask_nobias1 <- thresholds_out_nobias[[c]]$bin_max
  
  maps_biased <- pred_map_fun2(paste0(c,"\nBiased"),pred1,mask1,pts1)
  maps_nobiased <- pred_map_fun2(paste0(c,"\nRandom"),pred_nobias1,mask_nobias1,pts1)
  
  map_out_biased <- ggdraw(maps_biased$pred_map)+
    draw_plot(maps_biased$inset,x = 0.8,y = 0.4,
              width = w,height = h)
  map_out_nobiased <- ggdraw(maps_nobiased$pred_map)+
    draw_plot(maps_nobiased$inset,x = 0.8,y = 0.4,
              width = w,height = h)
  
  grid.arrange(map_out_biased,map_out_nobiased,nrow = 1)
}

pred_map_fun3 <- function(b,pred,mask,pts){
  require(cowplot,quietly = T)
  
  pred <- pred[[b]]
  mask <- mask[[b]]
  pts <- pts[pts$Binomial_species == b,]
  out <- pred_map_fun2(b,pred,mask,pts)
  
  ggdraw(out$pred_map)+
  draw_plot(out$inset,x = 0.7,y = 0.5,
            width = w,height = h)
}


pred_map_fun2 <- function(b,pred,mask,pts){
  
  
  pred <- raster::trim(raster::mask(pred,mask))
  pred_ext <- st_as_sfc(st_bbox(pred))
  
  #pred_df <- as.data.frame(pred,xy = T,na.rm = T)
  #names(pred_df) <- c("x","y","var")
  pred <- rast(pred)
  names(pred) <- 'var'
  list(pred_map = pred_map_fun1(pred,label = b,pts = pts),
       inset = inset_fun1(pred_ext))
}

pred_map_fun1 <- function(a,label,pts){
  ext_a <- ext(a)
  ggplot()+
    geom_spatraster(data = a,aes(fill = var),na.rm = T,maxcell = 1e+20)+
    #geom_tile(data = a,aes(x = x, y = y),fill = "blue4")+
    geom_sf(data = st_as_sf(wg_bnd_atree),col = "black",fill = NA)+
    geom_sf(data = pts,fill = "yellow",shape= 21,size = 1.5)+
    scale_fill_distiller(type = "seq",palette = "Reds",direction = 1,name = label,na.value = NA)+
    coord_sf(crs = 4326)+
    xlim(ext_a[1],78.5)+
    ylim(ext_a[3],ext_a[4])+
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title.position = "top"))+
    theme_classic()+
    theme(legend.position = c(0.9,0.85),
          legend.title = element_text(size = 9,face = "italic"),
          legend.text = element_text(size = 8))
  
}


#### gain lost area maps ####

gl_maps_fun <- function(c){
  pred1 <- thresholds_out[[c]]$pred
  mask1 <- thresholds_out[[c]]$bin_max
  pts1 <- thresholds_out[[c]]$pts
  change_ssp126_1 <- change_ssp126[[c]]
  change_ssp585_1 <- change_ssp585[[c]]
  #if(c%in%D$distribution_name){
  #  c <-D$phylogeny_name[D$distribution_name == c]
  #}
  c1 <- gsub("_","\n",c)
  curr_map <- pred_map_fun2(c1,pred1,mask1,pts1)
  gl_ssp126 <- pred_map_fun2_gl("SSP1-2.6",change_ssp126_1)
  gl_ssp585 <- pred_map_fun2_gl("SSP5-8.5",change_ssp585_1)
  
  map_curr <- ggdraw(curr_map$pred_map)#+
              #draw_plot(curr_map$inset,x = 0.85,y = 0.6,vjust = 0.5,
              #          width = w,height = h)
  map_gl_ssp126 <- ggdraw(gl_ssp126$pred_map)#+
    #draw_plot(gl_ssp126$inset,x = 0.85,y = 0.6,vjust = 0.5,
    #          width = w,height = h)
  map_gl_ssp585 <- ggdraw(gl_ssp585$pred_map)#+
    #draw_plot(gl_ssp585$inset,x = 0.85,y = 0.6,vjust = 0.5,
    #          width = w,height = h)
  grid.arrange(map_curr,map_gl_ssp126,map_gl_ssp585,nrow = 1)
}

pred_map_fun2_gl <- function(b,pred){
  
  
  pred <- raster::trim(pred)
  pred_ext <- st_as_sfc(st_bbox(pred))
  
  #pred_df <- as.data.frame(pred,xy = T,na.rm = T)
  #names(pred_df) <- c("x","y","var")
  pred <- rast(pred)
  pred <- as.factor(pred)
  names(pred) <- 'var'
  list(pred_map = pred_map_fun1_gl(pred,label = b),
       inset = inset_fun1(pred_ext))
}

pred_map_fun1_gl <- function(a,label){
  ext_a <- ext(a)
  ggplot()+
    geom_spatraster(data = a,aes(fill = var),na.rm = T,maxcell = 1e+20)+
    #geom_tile(data = a,aes(x = x, y = y),fill = "blue4")+
    geom_sf(data = st_as_sf(wg_bnd_atree),col = "black",fill = NA)+
    #geom_sf(data = pts,fill = "yellow",shape= 21,size = 1.5)+
    scale_fill_manual(values = c("-1"="red4","1"="royalblue3","2"="green4"),
                      labels = c("lost","constant","gain"),
                      name = label,na.value = "#00000000",na.translate = F)+
    coord_sf(crs = 4326)+
    xlim(ext_a[1],78.5)+
    ylim(ext_a[3],ext_a[4])+
    #guides(fill = guide_colorbar(direction = "horizontal",
     #                            title.position = "top"))+
    theme_classic()+
    theme(legend.position = c(0.9,0.85),
          legend.title = element_text(size = 9,face = "italic"),
          legend.text = element_text(size = 8))
  
}
####inset####
inset_fun1 <- function(b){
  ggplot(st_as_sf(wg_bnd_atree))+
    geom_sf()+
    geom_sf(data = b,colour = 'red',fill = NA)+
    theme_classic()+
    coord_sf(crs = 4326)+
    theme_void()+
    theme(panel.border = element_rect(colour = 'black',fill = NA,linewidth = 0.5))
}


####@@@@get pred and eval@@@@####
fetch_eval <- function(dsn){
  #take maxent eval files and load as dataframe
  
  paths_fam <- list.dirs(dsn,
                         full.names = T,recursive = F)
  
  names_fam <- list.dirs(dsn,full.names = F,recursive = F)
  
  names(paths_fam) <- names_fam
  
  # paths_spp <- unlist(lapply(paths_fam,list.dirs,full.names = T,recursive = F))
  
  eval_current <- lapply(paths_fam,
                         function(x){
                           read.csv(
                             paste0(x,"/maxentResults.csv"),
                             header = T,row.names = 1
                           )
                         })
}

get_pred <- function(path,target,crs_obj = CRS("+proj=longlat +ellps=WGS84"),recursive){
  #expects the file matching 'target' to be raster
  #default crs is CRS("+proj=longlat +ellps=WGS84")
  f <- list.files(path,pattern = target,full.names = T,recursive = recursive) 
  
  
  if(length(f)>1){
    r <- stack(f)
  }else{
    r <- raster(f)
    
  }
  crs(r) <- crs_obj
  r
  
}

fetch_pred <- function(dsn,target,recursive = F){
  #expects folder structure as 'dsn/Family/species' to use get_pred
  #with recursive = T, subfolders within '/species' will be searched,
  #returns named list with 'species' as names  
  
  paths_fam <- list.dirs(dsn,
                         full.names = T,recursive = F)
  
  # names_fam <- list.dirs(dsn,full.names = F,recursive = F)
  # 
  # names(paths_fam) <- names_fam
  
  paths_spp <- unlist(lapply(paths_fam,
                             function(x){
                               spp <- list.dirs(x,full.names = T,recursive = F)
                               if(length(c(grep("maxent.cache",spp),grep("plots",spp)))>0){
                                 spp <- spp[-c(grep("maxent.cache",spp),
                                               grep("plots",spp)
                                 )]}
                               spp
                             })
  )
  
  names_spp <- do.call(rbind,strsplit(paths_spp,"/"))
  
  names(paths_spp) <- names_spp[,ncol(names_spp)]
  
  
  pred_current <- lapply(paths_spp,
                         function(y){
                           tryCatch(
                             get_pred(y,target = target,recursive = recursive),
                             error = function(e){
                               return(NULL)}
                           )
                         })
}

dsn <- "/Volumes/ani/ani_backup/backups/plos_revision_sdm/lizards_enmeval_12-3-2026/"
target <- "_sltd$"
fun <- basename


fetch_files <- function(dsn,target,fun,...){
  #expects folder structure as 'dsn/Family/species'
  #uses 'fun' on the files that contain target
  
  paths_fam <- list.dirs(dsn,
                         full.names = T,recursive = F)
  
  # names_fam <- list.dirs(dsn,full.names = F,recursive = F)
  # 
  # names(paths_fam) <- names_fam
  
  paths_spp <- unlist(lapply(paths_fam,
                             function(x){
                               spp <- list.dirs(x,full.names = T,recursive = F)
                               if(length(c(grep("maxent.cache",spp),grep("plots",spp)))>0){
                                 spp <- spp[-c(grep("maxent.cache",spp),
                                               grep("plots",spp)
                                 )]}
                               spp
                             })
  )
  
  names_spp <- basename(paths_spp)
  
  names(paths_spp) <- names_spp
  
  paths <- lapply(paths_spp,list.files,pattern = target,full.names = T)
  pred_current <- lapply(paths,fun,...)
  names(pred_current) <- names_spp
  
  pred_current
}

######biogegraphic zones@@####
clip_sdm_wgbiogeo_mask <- function(bin,pts,zones){
  
  #zones_sf <- st_as_sf(zones)
  
  pts_zones <- unique(raster::intersect(pts,zones)@data$MASSIF)
  zones_dat <- zones[zones$MASSIF%in%pts_zones,]
  
  Nzone <- zones_dat[which.max(zones_dat$ymax),]
  Szone <- zones_dat[which.min(zones_dat$ymin),]
  
  Nzones <- zones[zones$ymax>=Nzone$ymax,]
  Szones <- zones[zones$ymin<=Szone$ymin,]
  all_zones <- zones[
    zones$ymax<Nzone$ymax&zones$ymin>Szone$ymin,
  ]
  
  
  barrier_cntr_S <- cntr_poly[
    cntr_poly$level == 
      zones$szone[zones$MASSIF == Szone$MASSIF],]
  barrier_cntr_N <-  cntr_poly[
    cntr_poly$level ==
      zones$nzone[zones$MASSIF == Nzone$MASSIF],]
  
  barriers <- rbind(barrier_cntr_N,barrier_cntr_S)
  
  barriers_in <- raster::intersect(barriers,rbind(Nzone,all_zones,Szone))
  
  
  
  if(length(unique(zones_dat$MASSIF)) == 1){
    bin_out <- mask(bin,zones_dat)
    l <- list(out = bin_out,
              barriers = cntr_single[cntr_single$id%in%barriers_in$id,])
  }else{
    
    N_mask <- mask(bin,Nzones)
    N_mask <- mask(N_mask,barrier_cntr_N[barrier_cntr_N$id%in%barriers_in$id,])
    
    S_mask <- mask(bin,Szones)
    S_mask <- mask(S_mask,barrier_cntr_S[barrier_cntr_S$id%in%barriers_in$id,])
    
    if(nrow(all_zones)>0){
      all_zones_mask <- mask(bin,all_zones)
    }
    
    # 
    #     if(exists("N_mask")){
    #       s <- stack(N_mask,all_zones_mask)
    #     }else{s <- stack(all_zones_mask)}
    # 
    #     if(exists("S_mask")){
    #       s <-addLayer(s,S_mask)
    #     }
    
    if(exists("all_zones_mask")){
      s <- stack(N_mask,all_zones_mask,S_mask)
      bin_out <- stackApply(s,indices = c(1,1,1),max,na.rm = T)
    }else{
      s <- stack(N_mask,S_mask)
      bin_out <- stackApply(s,indices = c(1,1),max,na.rm = T)
    }
    
    #   if(nlayers(s) == 3){
    # 
    # bin_out <- stackApply(s,indices = c(1,1,1),max,na.rm = T)
    #   }
    # 
    #   if(nlayers(s) == 2){
    # bin_out <- stackApply(s,indices = c(1,1),max,na.rm = T)
    # 
    #   }
    
    l <- list(out = bin_out,
              barriers = cntr_single[cntr_single$id%in%barriers_in$id,])
  }
  
}

clip_sdm_wgbiogeo_occelev <- function(bin,pts){
  #clip the sdm prediction using contour for the minimum elevation among 
  #occurrence pints
  names(pts)[1] <- 'pts_id'
  
  
  #find lowest elevation among occurrence points
  
  cntr_poly_inter <- raster::intersect(pts,cntr_poly)
  
  cntr_allpts <- tapply(cntr_poly_inter$pts_id,cntr_poly_inter$id,
                        function(x) all(pts$pts_id%in%unique(x)),simplify = T
  )
  
  if(any(cntr_allpts)){
    cntr_id_allpts <- as.numeric(names(cntr_allpts)[cntr_allpts])
    
    id_sel <- cntr_id_allpts[
      which.max(cntr_poly[cntr_id_allpts,]@data$level)]
    
    min_cntr <- cntr_poly[cntr_poly$id == id_sel,]
    
    min_cntr_mask <- mask(bin,min_cntr)
    
    list(out = min_cntr_mask ,
         barrier = cntr_single[cntr_single$id == id_sel,]) 
  }else{
    message("single contour covering all occurrences could not be identified")
    
    cntr_allpts <- tapply(cntr_poly_inter$pts_id,cntr_poly_inter$level,
                          function(x) all(pts$pts_id%in%unique(x)),simplify = T
    )
    if(any(cntr_allpts)){
      min_elev <- max(as.numeric(names(cntr_allpts)[cntr_allpts]))
      
      min_cntr <- cntr_poly[cntr_poly$id%in%
                              cntr_poly_inter$id[
                                cntr_poly_inter$level == min_elev],
      ]
      min_cntr_mask <- mask(bin,min_cntr)
      list(out = min_cntr_mask ,
           barrier = cntr_single[cntr_single$level == min_elev,])
      
    }else{
      message("Species is distributed bellow dicontinuous elevations. Therefore,
            cannot be clipped.")
      list(out = bin,barrier = cntr_single[cntr_single$level == min(cntr_single$level),
      ])
    }
  }
}

clip_sdm_wgbiogeo <- function(bin,pts,zones){
  
  #zones_sf <- st_as_sf(zones)
  
  pts_zones <- unique(raster::intersect(pts,zones)@data$MASSIF)
  zones_dat <- zones@data[zones$MASSIF%in%pts_zones,]
  Nzone <- zones_dat[which.max(zones_dat$ymax),]
  Szone <- zones_dat[which.min(zones_dat$ymin),]
  all_zones <- zones[
    zones$ymax<=Nzone$ymax&zones$ymin>=Szone$ymin,
  ]
  
  barrier_cntr_S <- cntr_poly[
    cntr_poly$level == 
      zones$szone[zones$MASSIF == Szone$MASSIF],]
  barrier_cntr_N <-  cntr_poly[
    cntr_poly$level ==
      zones$nzone[zones$MASSIF == Nzone$MASSIF],]
  
  barriers <- rbind(barrier_cntr_N,barrier_cntr_S)
  
  barriers_in <- raster::intersect(barriers,all_zones)
  
  list(out = mask(bin,all_zones),
       barriers = cntr_single[cntr_single$id%in%barriers_in$id,])
}

####@@ gain lost areas @@####
gain_lost <- function(current,future){
  #the binary raster files should have NA for not suitable
  #current: binary distribution raster for current climate
  #future: binary distribution raster for future climate
  #return: raster with -1 for lost areas, 1 for constant , and 2 for gain
  
  same <- current + future
  same[!is.na(same)] <- 1
  
  lost <- mask(current,future,inverse = T)
  lost[!is.na(lost)] <- -1
  
  gain <- mask(future,current,inverse = T)
  gain[!is.na(gain)] <- 2
  
  out <- stackApply(stack(lost,same,gain),indices = c(1,1,1),sum,na.rm = T)
  out[out == 0] <- NA
  out
  
  
}

####@@ make regression tables@@####
# a custom glance funtion for calculating pseudo r sq of glm model
glance_fun_custom <- function(mod){
  out <- broom::glance(mod)
  out %>% mutate(AICc = AICcmodavg::AICc(mod),pseudoR2 = (1-(deviance/null.deviance)))
}

format_table <- function(mod){
  mod %>% 
    tbl_regression() %>% 
    modify_table_body(
      ~ .x |>
        dplyr::filter(
          !(row_type %in% "label" & variable %in% "scnr:lat") &
            !(row_type %in% "label" & variable %in% "scnr:I(lat^2)")&
            !(row_type %in% "label" & variable %in% "scnr")
        )
    ) %>%
    modify_indent(columns = matches("label"),indent = 0L) %>% 
    add_glance_table(include = c(AICc,pseudoR2),glance_fun = glance_fun_custom)
}

# make a coef data frame for table
format_mod_coefs <- function(mod_name,mod_list){
  mod <- mod_list[[mod_name]]
  mod_s <- summary(mod)
  #add confidence intercals
  mod_confint <- confint(mod)
  colnames(mod_confint) <- c("confint_lower","confint_upper")
  
  #format confint
  mod_confint <- round(mod_confint,digits = 3)
  confint1 <- paste0("(",mod_confint[,1],",",mod_confint[,2],")")
  
  mod_coef <- data.frame(mod_s$coefficients,confint = confint1)
  names(mod_coef) <- c("Estimate","StdError","t.value","pval","confint")
  #significance codes
  mod_coef <- mod_coef %>% 
    mutate(sig_codes = case_when(
      pval>0.1 ~ "",
      pval>0.05& pval<=0.1 ~ ".",
      pval>0.01 & pval<=0.05 ~ "*",
      pval>0.001 & pval<=0.01 ~ "**",
      pval<0.001 ~ "***"
    ))
  
  #for mat a single cell in table with coef, confint, and sig code, with line break after sig code
  mod_coef_s <- data.frame(var1 =rownames(mod_coef), coef = paste0(round(mod_coef$Estimate,3),
                                                                   mod_coef$sig_codes,"<br>",
                                                                   mod_coef$confint))
  
  rownames(mod_coef_s) <- rownames(mod_coef)
  #calculate pseudo rsq and AICc
  if("glm"%in%class(mod)){
    pseudorsq <- round(1-(mod_s$deviance/mod_s$null.deviance),2)
    
    mod_coef_s <- rbind(mod_coef_s,data.frame(var1 = c("AICc","pseudo-R2"),
                                              coef = c(round(AICc(mod),2),pseudorsq)))
  }
  if(length(class(mod)) == 1){
    if(class(mod) == "lm"){
      mod_coef_s <- rbind(mod_coef_s,data.frame(var1 = c("AICc","R2"),
                                                coef = c(round(AICc(mod),2),round(mod_s$adj.r.squared,2))
      )) 
    }
  }
  
  names(mod_coef_s)[2] <- mod_name
  mod_coef_s
  
}

# format regression table in gt

#expects a data frame generated by 'format_mod_coefs

format_regression_tab_latitude <- function(dat){
  #format table with gt
  dat %>% 
    arrange(AICc)  %>%  
    gt() %>% 
    fmt_markdown() %>% 
    cols_move_to_start(columns = mod_names) %>% 
    tab_options(
      table.border.bottom.width = px(4),       # Increase thickness (e.g., 4px)
      table.border.bottom.color = "black") %>% 
    cols_label(
      mod_names = "Models",
      X.Intercept. = "Intercept",
      scnrssp585 = "SSP5-8.5",
      lat = "Latitude",
      I.lat.2. = html("Latitude<sup>2<sup>"),
      scnrssp585.lat = "SSP5-8.5:Latitude",
      scnrssp585.I.lat.2. = html("SSP5-8.5:Latitude<sup>2<sup>"))  %>%  
    tab_style(
      style = cell_borders(sides = "top",color = "black",weight = px(2),style = "solid"),
      locations = cells_body(rows = 1)
    )  %>% 
    tab_style(
      style = cell_borders(sides = "top",color = "black",weight = px(2),style = "solid"),
      locations = cells_column_labels()
    ) %>% 
    tab_style(
      style = cell_borders(sides = c("top","bottom"),color = "black",weight = px(2),style = "solid"),
      locations = cells_row_groups()
    )
  
}

format_regression_tab_elevation <- function(dat){
  dat %>%  
    arrange(AICc) %>% 
    gt() %>% 
    fmt_markdown() %>% 
    cols_move_to_start(columns = mod_names) %>% 
    tab_options(
      table.border.bottom.width = px(4),       # Increase thickness (e.g., 4px)
      table.border.bottom.color = "black") %>%
    cols_label(
      mod_names = "Models",
      X.Intercept. = "Intercept",
      scnrssp585 = "SSP5-8.5",
      elev = "Elevation",
      I.elev.2. = html("Elevation<sup>2<sup>"),
      scnrssp585.elev = "SSP5-8.5:Elevation",
      scnrssp585.I.elev.2. = html("SSP5-8.5:Elevation<sup>2<sup>"))  %>%  
    tab_style(
      style = cell_borders(sides = "top",color = "black",weight = px(2),style = "solid"),
      locations = cells_body(rows = 1)
    )  %>%  
    tab_style(
      style = cell_borders(sides = "top",color = "black",weight = px(2),style = "solid"),
      locations = cells_column_labels()
    ) %>%  
    tab_style(
      style = cell_borders(sides = c("top","bottom"),color = "black",weight = px(2),style = "solid"),
      locations = cells_row_groups())
}
