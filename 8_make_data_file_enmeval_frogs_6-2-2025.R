#Preparation of data file for present future distributions - Frogs
# on mac
setwd('/Users/shivaniagarwal/Desktop/aniruddha/My_Documents/IISc_postdoc/')

#on windows
setwd('D:/My_Documents/OneDrive - Indian Institute of Science/My_Documents/IISc_postdoc/')



library(dplyr)
temp_env <- new.env()
load("sdm/all_frogs_enmeval/all_frogs_sdm_enmevla_18-4-2023.RData",envir = temp_env)
ls(pattern = "swd", envir = temp_env)


maxent_swd_occ5 <- get("maxent_swd_occ5",envir = temp_env)
mods_name <- get('mods_name',envir = temp_env)


mod_params <- do.call(rbind,lapply(mods_name,function(a){
  a1 <- unlist(strsplit(a,"[/]"))
  a2 <- a1[2]
  a3 <- tail(unlist(strsplit(a1[3],"_")),3)[1:2]
  c(a2,a3)
  
})
)

colnames(mod_params) <- c('spp','feature_class','regularization')

mod_params[,"feature_class"] <- gsub("^fc.","",
                                     mod_params[,"feature_class"])

mod_params[,"regularization"] <- gsub("^rm.","",mod_params[,"regularization"]) 

maxent_swd <- maxent_swd_occ5

data_file <- maxent_swd %>% 
  group_by(Family,Binomial_species) %>% 
  summarise(occurrences = length(Binomial_species))

data_file <- data.frame(data_file,
                        mod_params[match(data_file$Binomial_species,
                                         mod_params[,'spp']),])


rm(temp_env)

temp_env <- new.env()
load("sdm_v1/frogs/all_frogs_enmeval_future_prediction_mask_biogeozones_v1_14-6-2024.RData",
     envir = temp_env)
total_area_w <- get("total_area_w",envir = temp_env)
spp_zones <- get("spp_zones",envir = temp_env)

data_file <- data.frame(data_file,
                        total_area_w[match(data_file$Binomial_species,
                                           total_area_w$Binomial_species),
                                     c("area_current","min","max")])
names(data_file)[5:6] <- c("area_ssp126","area_ssp585") 

rm(temp_env)

temp_env <- new.env()
load("sdm_v1/frogs/all_frogs_enmeval_lat_elev_ranges_v1.RData",
     envir =temp_env)
ext_lat <- get("ext_lat",envir = temp_env)

ssp585_ext_lat <- get("ssp585_ext_lat",envir = temp_env)

ext_lat1 <- data.frame(lat_min_curr = ext_lat$min,
                       lat_max_curr = ext_lat$max,
                       # lat_min_ssp126 = ssp126_ext_lat$min[match(rownames(ext_lat),
                       #                                           rownames(ssp126_ext_lat))],
                       # lat_max_ssp126 = ssp126_ext_lat$max[match(rownames(ext_lat),
                       #                                           rownames(ssp126_ext_lat))],
                       lat_min_ssp585 = ssp585_ext_lat$min[match(rownames(ext_lat),
                                                                 rownames(ssp585_ext_lat))],
                       lat_max_ssp585 = ssp585_ext_lat$max[match(rownames(ext_lat),
                                                                 rownames(ssp585_ext_lat))]
                       )
rownames(ext_lat1) <- rownames(ext_lat)

data_file <- data.frame(data_file,
                        ext_lat1[match(data_file$Binomial_species,rownames(ext_lat1)),])



elev_range_w <- get("elev_range_w",envir = temp_env)
all_frogs_elevrange_zones <- get("all_frogs_elevrange_zones",envir = temp_env)
# ssp126_elevrange_zones <- get("ssp126_elevrange_zones",envir = temp_env)
# 
# ssp126_area_zones_l <- get("ssp126_area_zones_l",envir = temp_env)

d1 <- get("d1",envir = temp_env)

names(d1)[2:9] <- paste("elev",names(d1)[2:9],sep = "_")

data_file <- data.frame(data_file,
                        d1[match(data_file$Binomial_species,d1$spp),!names(d1)%in%"spp"])

View(data_file)

d <- data.frame(data_file[,c("Family","Binomial_species","occurrences")],
                spp_zones[match(data_file$Binomial_species,spp_zones$spp),
                          c("bias","Nzone","Szone")],
                mod_params[match(data_file$Binomial_species,mod_params[,'spp']),
                           2:3],
                data_file[,!names(data_file)%in%c("Family","Binomial_species","occurrences")])





write.csv(d,"analysis/present_future_distribution_enmeval/data_and_scripts//data_file_enmeval_frogs_6-2-2025.csv")

save.image("analysis/present_future_distribution_enmeval/data_and_scripts/make_data_file_enmeval_frogs_6-2-2025.RData")

rm(list = ls())

