
library(tidyverse)
library(dplyr)
library(sf)
library(sp)
library(spatial)
library(rgdal)
library(rgeos)
library(raster)
library(terra)
library(rayshader)
library(lidR)
library(sp)
library(nngeo)
library(future)
library(rmapshaper)
library(concaveman)
library(parallel)
library(foreach)
library(doParallel)
library(smoothr)
library(ForestTools)
library(rmapshaper)
library(gdalUtilities)
library(exactextractr)
library(alphashape3d)
library(aRchi)
library(concaveman)

site_names = c("Harrison",
               "JordanRiver",
               "Kalamalka",
               "Skimikin",
               "TJC",
               "Whitecourt"
               )

##---------------------------------------------------------------------------##
# extract reflectance values

#cl = makeCluster(3) #not to overload your computer
#registerDoParallel(cl)

foreach(i = seq_along(site_names), 
        .combine = 'c',
        .packages = c("dplyr", "raster", "terra", "sf", "exactextractr", "lidR")
        ) %do% {

  name = site_names[i]

  dir = paste0("D://Sync//_Sites//Sx_Genecology//", name)
  # 
  # shad_mask = terra::rast(paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_mask2.tif"))
  # 
  # ms = terra::rast(paste0(dir, "\\input\\MS_ORTHO\\", name, "_MS.tif"))
  # 
  # pols_ms = sf::st_read(paste0(dir, "\\output\\HULLS\\", name, "_z50_snap.shp")) %>% 
  #   filter(!sf::st_is_empty(.)) %>% 
  #   sf::st_buffer(dist = -.05) %>%
  #   terra::vect()
  # 
  # ms_mask = terra::mask(ms, pols_ms)
  # 
  # ms_mask_trim = terra::trim(ms_mask)
  # 
  # mask_trim = terra::crop(shad_mask, ms_mask_trim)
  # 
  # # may not work with larger datasets. was crashing before being trimmed 
  # ms_mask2 = terra::mask(ms_mask_trim, mask_trim, maskvalues = 1, updatevalue = NA)
  # 
  # terra::writeRaster(ms_mask2, paste0(dir, "\\output\\RASTER\\", name, "_MS_mask.tif"),
  #                    overwrite = TRUE)

  pols = st_read(paste0(dir, "\\output\\HULLS\\", name, "_z50_snap.shp")) %>%
  filter(!st_is_empty(.)) %>%
  dplyr::filter(st_is(., c("POLYGON","MULTIPOLYGON")))

  ms_mask = rast(paste0(dir, "\\output\\RASTER\\", name, "_MS_mask.tif"))


  # ms_sum = (ms_mask[[1]] +
  #             ms_mask[[2]] +
  #             ms_mask[[3]] +
  #             ms_mask[[4]] +
  #             ms_mask[[5]] +
  #             ms_mask[[6]] +
  #             ms_mask[[7]] +
  #             ms_mask[[8]] +
  #             ms_mask[[9]] +
  #             ms_mask[[10]])

  rast_list = list(# reflectance values
      R444 = ms_mask[[1]],
      R475 = ms_mask[[2]],
      R531 = ms_mask[[3]],
      R560 = ms_mask[[4]],
      R650 = ms_mask[[5]],
      R668 = ms_mask[[6]],
      R705 = ms_mask[[7]],
      R717 = ms_mask[[8]],
      R740 = ms_mask[[9]],
      R842 = ms_mask[[10]],

      # chlorophyll
      #NDVI = (ms_mask[[10]] - ms_mask[[5]]) / (ms_mask[[10]] + ms_mask[[5]]),
      NDVI = (ms_mask[[10]] - ms_mask[[6]]) / (ms_mask[[10]] + ms_mask[[6]]),
      NIRvNDVI = ((ms_mask[[10]] - ms_mask[[6]]) / (ms_mask[[10]] + ms_mask[[6]])) * ms_mask[[10]],
      NIRvNDRE = ((ms_mask[[10]] - ms_mask[[8]]) / (ms_mask[[10]] + ms_mask[[8]])) * ms_mask[[10]],
      #NIRv3 = ((ms_mask[[10]] - ms_mask[[7]]) / (ms_mask[[10]] + ms_mask[[7]])) * ms_mask[[10]],
      NIRvCCI = (ms_mask[[3]] - ms_mask[[5]]) / (ms_mask[[3]] + ms_mask[[5]]) * ms_mask[[10]],

      NDRE1 = (ms_mask[[10]] - ms_mask[[7]]) / (ms_mask[[10]] + ms_mask[[7]]),
      NDRE2 = (ms_mask[[10]] - ms_mask[[8]]) / (ms_mask[[10]] + ms_mask[[8]]),
      NDRE3 = (ms_mask[[10]] - ms_mask[[9]]) / (ms_mask[[10]] + ms_mask[[9]]),
      EVI = (2.5 * (ms_mask[[10]] - ms_mask[[6]])) / (ms_mask[[10]] + (6 * ms_mask[[6]]) - (7.5 * ms_mask[[1]] + 1)),
      #Gcc = ms_mask[[4]] / (ms_mask[[2]] + ms_mask[[4]] + ms_mask[[6]]),
      #Gcc2 = ms_mask[[3]] / (ms_mask[[2]] + ms_mask[[3]] + ms_mask[[6]]),
      Gcc = (ms_mask[[3]] + ms_mask[[4]]) / (ms_mask[[1]] + ms_mask[[2]] + ms_mask[[3]] + ms_mask[[4]] + ms_mask[[5]] + ms_mask[[6]]),

      # carotenoids
      SIPI = (ms_mask[[10]] - ms_mask[[1]]) / (ms_mask[[10]] - ms_mask[[6]]),
      PRI = (ms_mask[[3]] - ms_mask[[4]]) / (ms_mask[[3]] + ms_mask[[4]]),
      CCI = (ms_mask[[3]] - ms_mask[[5]]) / (ms_mask[[3]] + ms_mask[[5]]),


      # Red edge
      RE_upper = (ms_mask[[9]] - ms_mask[[8]]) / 23,
      RE_lower = (ms_mask[[8]] - ms_mask[[7]]) / 12,
      RE_total = (ms_mask[[9]] - ms_mask[[7]]) / 35)

  rast_all = rast(rast_list)

  df_spectral = exact_extract(rast_all, pols, fun = "mean", append_cols = "Obs")

  saveRDS(df_spectral, paste0(dir, "\\CSV\\", name, "_spectral.rds"))
}

#future:::ClusterRegistry("stop")


#writeRaster(ndre, paste0(dir, "\\output\\RASTER\\", name, "_NDRE.tif"), overwrite = TRUE)
#writeRaster(pri, paste0(dir, "\\output\\RASTER\\", name, "_PRI.tif"), overwrite = TRUE)
#writeRaster(cci, paste0(dir, "\\output\\RASTER\\", name, "_CCI.tif"), overwrite = TRUE)
#writeRaster(re_upper, paste0(dir, "\\output\\RASTER\\", name, "_RE_upper.tif"), overwrite = TRUE)

##---------------------------------------------------------------------------##

# define the key metrics you want 

#test = readLAS("D:\\Sync\\_Sites\\Sx_Genecology\\Skimikin\\output\\CROWNS\\Skimikin_p411_r4_t2.laz")
#las = test


key_metrics = function(las) {
  
  #las_green = filter_poi(las, las@data$G > las@data$R)
  
  Z = las@data$Z
  
  chm = grid_metrics(las, func = ~max(Z), res = 0.05)
  chm_mean = grid_metrics(las, func = ~mean(Z), res = 0.1)
  chm_mean[chm_mean < (maxValue(chm_mean))] = NA 
  chm_mean_trim = raster::trim(chm_mean)
  
  apex = clip_roi(las, extent(chm_mean_trim))@data %>% 
    dplyr::filter(Z == max(.$Z)) %>% 
    dplyr::slice(1L)
  
  origin = c(apex$X, apex$Y, apex$Z)
  a_point = c(apex$X, apex$Y, 0)

  
  myangle3d = function(b1, b2, b3) {
    aRchi::angle3d(o = c(origin[1], origin[2], origin[3]),
            a = c(a_point[1], a_point[2], a_point[3]),
            b = c(b1, b2, b3))
  }
  
  ang = las@data %>%
    dplyr::select(X, Y, Z) %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(angle = myangle3d(b1 = X, b2 = Y, b3 = Z)) %>% 
    dplyr::ungroup() %>% 
    dplyr::summarise_at(dplyr::vars(angle), .funs = c(mean, sd), na.rm = TRUE)
  
  
  # las_dec = las %>% 
  #   decimate_points(homogenize(20, 10))
    
  # alphashadep3d
  a3d = cbind(las@data$X, las@data$Y, las@data$Z)

  a3d[,1] = a3d[,1] - mean(a3d[,1]) #center points around 0,0,0
  a3d[,2] = a3d[,2] - mean(a3d[,2]) #center points around 0,0,0
  
  # shape = ashape3d(x = a3d, alpha = 1)
  # plot(shape)
  
  data.frame(

    # Crown height
    Zq999 = as.numeric(quantile(Z, 0.999)),
    Zq99 = as.numeric(quantile(Z, 0.990)),
    Z_mean = mean(Z),
    
    # Crown size
    n_points = length(las@data$Z),
    area = lidR::area(las),
    #volume = raster::cellStats(chm_vol, sum),
    #n_green_points = length(las_green@data$Z),
    vol_convex = alphashape3d::volume_ashape3d(alphashape3d::ashape3d(x = a3d, alpha = Inf)),
    vol_concave = alphashape3d::volume_ashape3d(alphashape3d::ashape3d(x = a3d, alpha = 1)),
    #vol_a05= alphashape3d::volume_ashape3d(alphashape3d::ashape3d(x = a3d, alpha = .5)),
        
    apex_angle = ang$fn1,
    apex_sd = ang$fn2,
    
    # Crown complexity
    CV_Z = sd(Z) / mean(Z),
    rumple = lidR::rumple_index(chm),
    CRR = (mean(Z) - min(Z)) / (max(Z) - min(Z)))
    

}

##---------------------------------------------------------------------------##


##---------------------------------------------------------------------------##
# extract ITC structural metrics


plan(multisession, workers = 10L)

for (i in seq_along(site_names)){
  
  name = site_names[i]
  
  dir = paste0("D://Sync//_Sites//Sx_Genecology//", name)
  
  CROWNS = readLAScatalog(paste0(dir, "\\output\\CROWNS\\"))
  opt_chunk_size(CROWNS) = 0 # processing by files
  opt_chunk_buffer(CROWNS) = 0 # no buffer
  opt_wall_to_wall(CROWNS) = TRUE # disable internal checks to ensure a valid output. Free wheel mode
  
  # compute ITC structural metrics
  structural_metrics = tree_metrics(CROWNS, ~key_metrics(las))
  
  df_structural = structural_metrics@data
  saveRDS(df_structural, paste0(dir, "\\CSV\\", name, "_structural.rds"))
  } 



# pols_metrics = left_join(st_as_sf(pols), MERGED_metrics, by = c("Obs" = "treeID"), copy = TRUE)
# st_geometry(pols_metrics) = st_geometry(pols)


future:::ClusterRegistry("stop")

##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
# join the structural and spectral metrics to the census and climate data

master = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\Sx_MASTER.csv") %>% 
  dplyr::filter(Site == "Harr" |
                  Site == "Jord" |
                  Site == "Kala" |
                  Site == "Skim" |
                  Site == "Tete" |
                  Site == "Whit") %>% 
  mutate(Site = recode(Site, 
                       Harr = "Harrison",
                       Jord = "JordanRiver",
                       Kala = "Kalamalka",
                       Skim = "Skimikin",
                       Tete = "TJC",
                       Whit = "Whitecourt")) %>% 
  rename(Obs = Meas_seq) %>% 
  filter(Prov != 9999) %>% 
  group_by(Site) %>% 
  mutate(total_planted = n()) %>% 
  group_by(Site, Prov) %>% 
  mutate(n_planted = n(),
         prop = n_planted / total_planted)

surv = master %>% 
  filter(Thin_age == ".") %>% 
  group_by(Site, Prov) %>% 
  mutate(prov_n_unthinned = n()) %>% 
  group_by(Site) %>% 
  mutate(total_unthinned = n()) %>% 
  filter(HT16 != "." & HT16 != "0") %>% 
  mutate(total_alive = n()) %>% 
  summarise(survival = total_alive / total_unthinned) %>% 
  distinct()

###
area_harr = st_read("D://Sync//_Sites//Sx_Genecology//Harrison//output//HULLS//Harrison_z50_snap.shp") %>% 
  filter(Edge != 1) %>% 
  st_as_sf() %>% 
  concaveman(length_threshold = 10, concavity = 5) %>% 
  st_zm() %>% 
  st_area()    
    
pols_harr = st_read("D://Sync//_Sites//Sx_Genecology//Harrison//output//HULLS//Harrison_z50_snap.shp") %>%
  st_make_valid() %>% 
  filter(!st_is_empty(.)) %>% 
  dplyr::select(Obs, Absent, Edge, Edited) %>% 
  mutate(Dead = 0) %>% 
  dplyr::mutate(Site = "Harrison",
                plot_area = area_harr)
###
area_jord = st_read("D://Sync//_Sites//Sx_Genecology//JordanRiver//output//HULLS//JordanRiver_z50_snap.shp") %>% 
  filter(Edge != 1) %>% 
  st_as_sf() %>% 
  concaveman(length_threshold = 10, concavity = 5) %>% 
  st_zm() %>% 
  st_area()   

pols_jord = st_read("D://Sync//_Sites//Sx_Genecology//JordanRiver//output//HULLS//JordanRiver_z50_snap.shp") %>%
  st_make_valid() %>% 
  filter(!st_is_empty(.)) %>% 
  dplyr::select(Obs, Absent, Edge, Edited, Dead) %>% 
  dplyr::mutate(Site = "JordanRiver",
                plot_area = area_jord)
###
area_kal = st_read("D://Sync//_Sites//Sx_Genecology//Kalamalka//output//HULLS//Kalamalka_z50_snap.shp") %>% 
  #filter(Edge != 1) %>% 
  st_as_sf() %>% 
  concaveman(length_threshold = 10, concavity = 5) %>% 
  st_zm() %>% 
  st_area()   
      
pols_kal = st_read("D://Sync//_Sites//Sx_Genecology//Kalamalka//output//HULLS//Kalamalka_z50_snap.shp") %>%
  st_make_valid() %>% 
  filter(!st_is_empty(.)) %>% 
  dplyr::select(Obs, Absent, Edge, Edited) %>% 
  mutate(Dead = 0) %>% 
  dplyr::mutate(Site = "Kalamalka",
                plot_area = area_kal)

###
area_skim = st_read("D://Sync//_Sites//Sx_Genecology//Skimikin//output//HULLS//Skimikin_z50_snap.shp") %>% 
  filter(Edge != 1) %>% 
  st_as_sf() %>% 
  concaveman(length_threshold = 10, concavity = 5) %>% 
  st_zm() %>% 
  st_area()  

pols_skim = st_read("D://Sync//_Sites//Sx_Genecology//Skimikin//output//HULLS//Skimikin_z50_snap.shp") %>%
  st_make_valid() %>% 
  filter(!st_is_empty(.)) %>% 
  dplyr::select(Obs, Absent, Edge, Edited) %>% 
  mutate(Dead = 0) %>% 
  dplyr::mutate(Site = "Skimikin",
                plot_area = area_skim)

###
area_tjc = st_read("D://Sync//_Sites//Sx_Genecology//TJC//output//HULLS//TJC_z50_snap.shp") %>% 
  filter(Edge != 1) %>% 
  st_as_sf() %>% 
  concaveman(length_threshold = 10, concavity = 5) %>% 
  st_zm() %>% 
  st_area()  

pols_tjc = st_read("D://Sync//_Sites//Sx_Genecology//TJC//output//HULLS//TJC_z50_snap.shp") %>%
  st_make_valid() %>% 
  filter(!st_is_empty(.)) %>% 
  dplyr::select(Obs, Absent, Edge, Edited, Dead) %>% 
  dplyr::mutate(Site = "TJC",
                plot_area = area_tjc)

###
area_whit = st_read("D://Sync//_Sites//Sx_Genecology//Whitecourt//output//HULLS//Whitecourt_z50_snap.shp") %>% 
  filter(Edge != 1) %>% 
  st_as_sf() %>% 
  concaveman(length_threshold = 10, concavity = 5) %>% 
  st_zm() %>% 
  st_area()  

pols_whit = st_read("D://Sync//_Sites//Sx_Genecology//Whitecourt//output//HULLS//Whitecourt_z50_snap.shp") %>%
  st_make_valid() %>% 
  filter(!st_is_empty(.)) %>% 
  dplyr::select(Obs, Absent, Edge, Edited, Dead) %>% 
  dplyr::mutate(Site = "Whitecourt",
                plot_area = area_whit)


pols = rbind(pols_harr, pols_jord, pols_kal, pols_skim, pols_tjc, pols_whit) %>% 
  left_join(master, by = c("Obs", "Site"))

################################################################################


clim = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\output\\Sx_CC_Sites_ReferencePeriod_2005_2020_Y_S.csv") %>% 
  dplyr::select(Site, Latitude, TD, EMT, MAP, AHM) %>% 
  mutate(logMAP = log(MAP)) %>% 
  dplyr::select(-MAP) %>% 
  pivot_longer(cols = Latitude:logMAP, names_to = "clim_metric", values_to = "value") %>% 
  mutate(Seedlot = NA) %>% 
  relocate(Site, .after = last_col()) %>% 
  relocate(Seedlot, .after = last_col())


sx_pops_clim = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\output\\Sx_CC_Seedlots_Normal_1961_1990_Y_S.csv") %>% 
  dplyr::select(Seedlot, Latitude, TD, EMT, MAP, AHM) %>%  
  mutate(logMAP = log(MAP)) %>% 
  dplyr::select(-MAP) %>% 
  pivot_longer(cols = Latitude:logMAP, names_to = "clim_metric", values_to = "value") %>% 
  bind_rows(clim) %>% 
  group_by(clim_metric) %>% 
  mutate(standardized = (value - mean(value)) / sd(value))

sx_pops = sx_pops_clim %>%
  filter(!is.na(Seedlot)) %>% 
  dplyr::select(-Site) %>% 
  right_join(filter(sx_pops_clim, !is.na(Site)), by = "clim_metric") %>% 
  mutate(standardized_d = standardized.x - standardized.y) %>% 
  rename("Seedlot" = "Seedlot.x") %>% 
  dplyr::select(-standardized.y, -(value.x:value.y)) %>% 
  pivot_wider(names_from = clim_metric, values_from = standardized_d) %>% 
  mutate(Euc = sqrt(Latitude^2 + TD^2 + EMT^2 + AHM^2 + logMAP^2)) %>% 
  #dplyr::select(-(Latitude:logMAP)) %>% 
  left_join(read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\Sx_CC_Seedlots_info.csv") %>% 
              dplyr::select(-number, -donor), by = c("Seedlot" = "population")) %>% 
  # standardize the climate variables
  #mutate(Latitude_z = (Latitude - mean(.$Latitude)) / sd(.$Latitude),
  #        Elevation_z = (Elevation - mean(.$Elevation)) / sd(.$Elevation), 
  #        MAT_z = (MAT - mean(.$MAT)) / sd(.$MAT), 
  #        MCMT_z = (MCMT - mean(.$MCMT)) / sd(.$MCMT), 
  #        MWMT_z = (MWMT - mean(.$MWMT)) / sd(.$MWMT), 
  #       TD_z = (TD - mean(.$TD)) / sd(.$TD),
  #        MAP_z = (MAP - mean(.$MAP)) / sd(.$MAP), 
  #        MSP_z = (MSP - mean(.$MSP)) / sd(.$MSP), 
  #        CMD_z = (CMD - mean(.$CMD)) / sd(.$CMD), 
  #        FFP_z = (FFP - mean(.$FFP)) / sd(.$FFP), 
  #        DD5_z = (DD5 - mean(.$DD5)) / sd(.$DD5)) %>% 
  mutate(zone = case_when(region == "NM" | region == "AZ" ~ "US (south)",
                          region == "ID" | region == "MT" | region == "WA" ~ "US (north)",
                          region == "BC" ~ "British Columbia",
                          region == "AB" ~ "Alberta",
                          region == "YK" | region == "NWT" ~ "Canada (north)",
                          region == "ON" ~ "Ontario"),
         class = if_else(region == "AB", "B", class)) %>%
  left_join(read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\SxGenecology_predicted_Euc.csv") %>% 
              dplyr::select(sxProv, species, propGla, propEng, propSit),
            by = c("Seedlot" = "sxProv"))

df_spectral = bind_rows(
  read_rds("D://Sync//_Sites//Sx_Genecology//Harrison\\CSV\\Harrison_spectral.rds") %>% 
  mutate(Site = "Harrison"),
  read_rds("D://Sync//_Sites//Sx_Genecology//JordanRiver\\CSV\\JordanRiver_spectral.rds") %>% 
    mutate(Site = "JordanRiver"),
  read_rds("D://Sync//_Sites//Sx_Genecology//Kalamalka\\CSV\\Kalamalka_spectral.rds") %>% 
    mutate(Site = "Kalamalka"),
  read_rds("D://Sync//_Sites//Sx_Genecology//Skimikin\\CSV\\Skimikin_spectral.rds") %>% 
    mutate(Site = "Skimikin"),
  read_rds("D://Sync//_Sites//Sx_Genecology//TJC\\CSV\\TJC_spectral.rds") %>% 
    mutate(Site = "TJC"),
  read_rds("D://Sync//_Sites//Sx_Genecology//Whitecourt\\CSV\\Whitecourt_spectral.rds") %>% 
    mutate(Site = "Whitecourt"))

df_structural = bind_rows(
  read_rds("D://Sync//_Sites//Sx_Genecology//Harrison\\CSV\\Harrison_structural.rds") %>% 
    mutate(Site = "Harrison"),
  read_rds("D://Sync//_Sites//Sx_Genecology//JordanRiver\\CSV\\JordanRiver_structural.rds") %>% 
    mutate(Site = "JordanRiver"),
  read_rds("D://Sync//_Sites//Sx_Genecology//Kalamalka\\CSV\\Kalamalka_structural.rds") %>% 
    mutate(Site = "Kalamalka"),
  read_rds("D://Sync//_Sites//Sx_Genecology//Skimikin\\CSV\\Skimikin_structural.rds") %>% 
    mutate(Site = "Skimikin"),
  read_rds("D://Sync//_Sites//Sx_Genecology//TJC\\CSV\\TJC_structural.rds") %>% 
    mutate(Site = "TJC"),
  read_rds("D://Sync//_Sites//Sx_Genecology//Whitecourt\\CSV\\Whitecourt_structural.rds") %>% 
    mutate(Site = "Whitecourt"))

metrics = full_join(df_spectral, df_structural, by = c("Obs" = "treeID", "Site" = "Site")) %>% 
  distinct()

pols_dat = left_join(pols, metrics, by = c("Obs", "Site")) %>% 
  mutate(HT10 = as.numeric(HT10),
         HT16 = as.numeric(HT16),
         Weev16 = as.numeric(Weev16),
         DBH16 = if_else(DBH16 == 0, NA_real_, as.numeric(DBH16)),
         no_dbh = if_else(is.na(DBH16), 1, 0)) %>% 
  mutate(vol_con_DBH = vol_convex / DBH16,
         vol_con_Zq999 = vol_convex / Zq999,
         Zq999_DBH = Zq999 / DBH16,
         growth_ratio = (HT16 - HT10) / HT16,
         growth_ratio_Zq999 = ((Zq999 * 100) - HT10) / (Zq999 * 100),
         stem_vol = exp(-9.951) * ((DBH16 / 10) ^ 1.807) * ((HT16 / 100) ^ 1.080),
         hybrid_stem_vol = exp(-9.951) * ((DBH16 / 10) ^ 1.807) * ((Zq999) ^ 1.080)) %>% 
  mutate(apex_angle =  (apex_angle * 180) / pi,
         apex_sd =  (apex_sd * 180) / pi) %>% 
  left_join(sx_pops, by = c("Prov" = "Seedlot", "Site" = "Site")) %>% 
  mutate(Edge = if_else(is.na(Edge), 0, Edge),
         Dead = if_else(is.na(Dead), 0, Dead)) %>% 
  group_by(Prov, Site) %>% 
  mutate(prov_vol = sum(stem_vol, na.rm = TRUE),
         vol_per_ha = prov_vol / prop)


saveRDS(pols_dat, "D://Sync//_Sites//Sx_Genecology//all_metrics.rds")

pols_dat_write = st_cast(pols_dat, "POLYGON")

st_write(obj = pols_dat_write, 
         dsn = paste0("D:\\Sync\\_Sites\\Sx_Genecology\\pols_dat.shp"),
         #layer = outliers,
         driver = "ESRI Shapefile",
         append = FALSE)


################################################################################
# Check outliers 

pols_dat = readRDS("D://Sync//_Sites//Sx_Genecology//all_metrics.rds") %>% 
  st_as_sf()

mod = function(df) {lm(DBH16 ~ Zq999, data = df)}

pols_test = pols_dat %>% 
  drop_na(DBH16, Zq999) %>% 
  filter(Edge != 1 & Dead != 1) %>% 
  group_by(Site) %>% 
  nest() %>% 
  mutate(model = map(data, ~ lm(DBH16 ~ Zq99 + log(Zq99), data = .x))) %>% 
  mutate(resids = map(model, residuals)) %>% 
  mutate(studreds = map(model, MASS::studres)) %>% 
  mutate(glance = map(model, broom::glance)) %>% 
  unnest(c(data, resids, studreds, glance)) %>% 
  ungroup() %>% 
  mutate(studentized_dbh = resids / sigma)

# height vs DBH
p1 = ggplot(pols_test, aes(y = DBH16, x = Zq99, color = Site)) +
  geom_vline(xintercept = 1) +
  geom_point(color = "steelblue4", size = 3, alpha = .2) + 
  geom_point(data = filter(pols_test, (studentized_dbh < -4)
             | (studentized_dbh > 6 & Zq999 < 4)
             | (studentized_dbh > 3 & Zq999 < 1)),
             color = "indianred", size = 5, alpha = 1) + 
  # geom_label(data = filter(pols_test, (studentized_dbh < -3 & Z_mean < 3) | (studentized_dbh > 4 & Z_mean < 4)), aes(label = Prov)) + 
  theme_bw(base_size = 24) +
  facet_wrap(. ~ Site, nrow = 2)

p1

outlier_df = filter(pols_test, (studentized_dbh < -4)
                  | (studentized_dbh > 6 & Zq999 < 4)
                  | (studentized_dbh > 3 & Zq999 < 1)) 


ggsave(plot = p1, 
       filename = "D:\\Sync\\Figures\\_FIGURES_paper3\\figure_SI_outliers.tiff", 
       device = "tiff",
       width = 15,
       height = 10,
       units = "in",
       dpi = 600)

saveRDS(outlier_df, "D://Sync//_Sites//Sx_Genecology//outliers.rds")

# st_write(obj = outlier_df, 
#          dsn = paste0("D:\\Sync\\_Sites\\Sx_Genecology\\outliers.shp"),
#          #layer = outliers,
#          driver = "ESRI Shapefile",
#          append = FALSE)

# elimate outliers and save
outlier_join = pols_dat %>% 
  anti_join(outlier_df)

saveRDS(outlier_join, "D://Sync//_Sites//Sx_Genecology//all_metrics_no_outliers.rds")



