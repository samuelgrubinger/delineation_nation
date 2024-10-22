
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(pastecs)
library(exactextractr)

##---------------------------------------------------------------------------##

dir = "D:\\Sync\\_Sites\\Skimikin_spectral"

date_list = c(
  #"Skimikin_2020_03_22",
  "Skimikin_2020_07_02",
  "Skimikin_2020_08_07",
  "Skimikin_2020_08_25",
  "Skimikin_2021_03_31",
  "Skimikin_2021_06_29",
  "Skimikin_2021_07_29",
  "Skimikin_2021_08_14",
  "Skimikin_2023_02_24")

################################################################################
# define Olivia's function for local minimum detection

## Function to find index values where derivative switches from neg to post (aka local min)----------:
find_closest_to_zero_and_index = function(values) {
  neg_slope <- numeric()
  pos_slope <- numeric()
  index_of_closest <- numeric()
  definition = numeric()
  neg_sum = numeric()
  pos_sum = numeric()
  
  k <- 1
  
  for (i in 2:length(values)) {
    #print(i)
    if (values[i - 1] < 0 & values[i] >= 0 & i > 20) { #i>15 because local min will not be in first 15 and this stops an error occuring where a local min is found in teh first 15 values and the def_positive indexing doe snot work
      # Check if the absolute value of the previous and current values is closer to zero
      # if (abs(values[i - 1]) < abs(values[i])) {
      #   closest_to_zero[k] <- values[i - 1]
      #   index_of_closest[k] <- i - 1
      # } else {
      #   closest_to_zero[k] <- values[i]
      #   index_of_closest[k] <- i
      # }
      def_positive <- c(values[i:i+15]) # vector of next 7 gradient values
      def_pos <- def_positive[def_positive > 0] #only taking positives
      pos_sum[k] <- sum(def_pos)
      
      def_negative <- c(values[(i - 10):i-15]) # vector of next 7 gradient values # 7 was used and found first min
      def_neg <- def_negative[def_negative < 0] #only taking negatives
      neg_sum[k] <- sum(def_neg)
      
      neg_slope[k] <- values[i - 1]
      pos_slope[k] <- values[i]
      index_of_closest[k] <- i - 1
      definition[k] <- sum(abs(def_neg), abs(def_pos)) # higher the value, more pronounced the local min
      k <- k + 1
    }
  }
  
  return(data.frame(Neg_Value = neg_slope, Pos_Value = pos_slope, Index = index_of_closest, definition = definition, neg_def = neg_sum, pos_def = pos_sum))
}

# normalized difference index function
vi = function(img, k, i) {
  bk <- img[[k]]
  bi <- img[[i]]
  vi <- (bk - bi) / (bk + bi)
  return(vi)
}

# local minimum of bimodal distribution from a vector of values
find_local_min = function(x) {
  
  d = stats::density(x)
  
  dy_dt = pracma::gradient(d$y) #list of first derivatives
  
  zeros = find_closest_to_zero_and_index(dy_dt) #list of values and their indicies where slope switches from neg to pos (local min)
  zeros_filtered <- zeros[(zeros$pos_def > 0),]#filters rows with pos_def > 0, filtering out small minimums on negative slopes
  zeros_local_min <- zeros_filtered[which.max(zeros_filtered$definition), ]
  
  x_zeros <- d$x[zeros_local_min$Index] #selecting NIR (aka x_mid) values that correspond to index values where the slope switches from neg to positive (indicating a local min)
  #threshold <- max(x_zeros[x_zeros<0.6])
  threshold = x_zeros
}

################################################################################

pols_buf = vect(paste0(dir, "\\output\\HULLS\\Edits\\Skimikin_z50_updated.shp"))

################################################################################
# look at NDVI values

ndvi_all = c(#rast(paste0(dir, "\\input\\MS_ORTHO\\Skimikin_2020_03_22_MS.tif")) %>% vi(10, 6),
             rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2020_07_02_MS.tif")) %>% vi(10, 6),
             rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2020_08_07_MS.tif")) %>% vi(10, 6),
             rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2020_08_25_MS.tif")) %>% vi(10, 6),
             rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2021_03_31_MS.tif")) %>% vi(10, 6),
             rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2021_06_29_MS.tif")) %>% vi(10, 6),
             rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2021_07_29_MS.tif")) %>% vi(10, 6),
             rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2021_08_14_MS.tif")) %>% vi(10, 6),
             rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2023_02_24_MS.tif")) %>% vi(10, 6))

terra::writeRaster(ndvi_all, paste0(dir, "\\output\\RASTER\\NDVI_all.tif"),
                   overwrite = TRUE)

ndvi_all = rast(paste0(dir, "\\output\\RASTER\\NDVI_all.tif"))

# scale values to the max NDVI value of each date
#ndvi_max = minmax(ndvi_all)[2,]
# ndvi_scale = ndvi_all / ndvi_max
# 
# terra::writeRaster(ndvi_scale, paste0(dir, "\\output\\RASTER\\NDVI_scale.tif"),
#                    overwrite = TRUE)

# drop values below the .72 threshold
ndvi_clamp = ndvi_all %>% terra::clamp(lower = .75, values = FALSE)

terra::hist(ndvi_clamp, breaks = 100)

# make a mask, dropped values to 1
ndvi_mask = terra::ifel(is.na(ndvi_clamp), 1, NA)

# make another mask, for multiplying later, 1 being kept pixels
ndvi_mask2 = ifel(is.na(ndvi_mask), 1, NA)

terra::writeRaster(ndvi_mask, paste0(dir, "\\output\\RASTER\\NDVI_mask.tif"),
                   overwrite = TRUE)

# terra::hist(ndvi_clamp, breaks = 100)
# 
# ndvi_mask2 = terra::ifel(is.na(ndvi_clamp), ndvi_clamp, 1)

################################################################################

# create a shadow mask with the minimum 

nir_all = c(#rast(paste0(dir, "\\input\\MS_ORTHO\\Skimikin_2020_03_22_MS.tif"))[[10]],
            rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2020_07_02_MS.tif"))[[10]],
            rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2020_08_07_MS.tif"))[[10]],
            rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2020_08_25_MS.tif"))[[10]],
            rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2021_03_31_MS.tif"))[[10]],
            rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2021_06_29_MS.tif"))[[10]],
            rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2021_07_29_MS.tif"))[[10]],
            rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2021_08_14_MS.tif"))[[10]],
            rast(paste0(dir, "\\input\\MS_ORTHO_average\\Skimikin_2023_02_24_MS.tif"))[[10]])

nir_mask = nir_all * ndvi_mask2

hist(nir_all, breaks = 100)

#terra::writeRaster(nir_mask, paste0(dir, "\\output\\RASTER\\NIR_mask.tif"),
#                   overwrite = TRUE)

# nir_mask_scale = terra::scale(nir_mask, center = FALSE)
# 
# hist(nir_mask_scale, breaks = 100)

nir_min = app(nir_mask, "min", na.rm=TRUE)
nir_which_min = app(nir_mask, "which.min", na.rm=TRUE)

terra::writeRaster(nir_min, paste0(dir, "\\output\\RASTER\\NIR_min.tif"),
                   overwrite = TRUE)

nir_min = rast(paste0(dir, "\\output\\RASTER\\NIR_min.tif"))

hist(nir_min, breaks = 100)

################################################################################
# find the local minimum of a bimodal distribution

#NIR_na = na.omit(nir_min)

threshold = find_local_min(na.omit(as.vector(nir_min)))

# threshold = c(find_local_min(na.omit(as.vector(nir_mask[[3]]))),
#               find_local_min(na.omit(as.vector(nir_mask[[4]]))),
#               find_local_min(na.omit(as.vector(nir_mask[[5]]))))

################################################################################

shadow_mask = nir_min
shadow_mask[shadow_mask > threshold] = NA
shadow_mask[shadow_mask <= threshold] = 1

terra::writeRaster(shadow_mask, paste0(dir, "\\output\\RASTER\\NIR_min_mask.tif"),
                   overwrite = TRUE)

shadow_mask_all = nir_all
shadow_mask_all[shadow_mask_all > threshold] = NA
shadow_mask_all[shadow_mask_all <= threshold] = 1

terra::writeRaster(shadow_mask_all, paste0(dir, "\\output\\RASTER\\NIR_min_mask_all.tif"),
                   overwrite = TRUE)

shadow_mask_id = nir_which_min
shadow_mask_id[nir_min > threshold] = NA

terra::writeRaster(shadow_mask_id, paste0(dir, "\\output\\RASTER\\NIR_min_mask_id.tif"),
                   overwrite = TRUE)

# (hist = NIR_na %>% 
#     as_tibble() %>% 
#     ggplot() +
#     geom_histogram(aes(x = min), bins = 150) +
#     geom_vline(xintercept = threshold, color = "red3") +
#     labs(x = paste0(mode, " , threshold: ",round(threshold, digits = 2) ))+
#     theme_bw()+
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           plot.title = element_text(size = 13,hjust = 0.75, vjust = -28)))
# 
# 
# ggsave(hist,
#        filename = paste0(dir, "\\", date_list[x], "_NIR_shadow_hist.jpeg"),
#        device = jpeg,
#        width = 8,
#        height = 8)

# saveRDS(as_tibble(NIR), 
#         paste0(dir, "\\CSV\\NIR_hist\\", date_list[x], "_NIR.rds"))

#hist(NIR, breaks = 100) + abline(v = threshold, col='red', lwd = 3)

# shadow_mask = ms_temp[[10]]
# shadow_mask[shadow_mask > threshold] = NA
# shadow_mask[shadow_mask <= threshold] = 1
# 
# terra::writeRaster(shadow_mask, paste0(dir, "\\output\\RASTER\\", date_list[x], "\\",  date_list[x], "_NIR_shadow_mask.tif"),
#                    overwrite = TRUE)


################################################################################

#------------------------------------------------------------------------------#
# produce a hillshade mask 
# 
# library(rayshader)
# library(suncalc)
# library(exifr)
# library(lubridate)
# library(rgdal)
# library(LaplacesDemon)

# make a mean values CHM
# fill NA with a minimum filter,
# smooth aggressively
# read in norm cloud
#plan(multisession, workers = 3L)



pols = st_read(paste0(dir, "\\output\\HULLS\\Edits\\Skimikin_z50_updated.shp")) %>%
  filter(!st_is_empty(.))

pols_spat = pols %>% 
  terra::vect()

# ratios = read_rds(paste0(dir, "\\CSV\\Corrected_values//ratios.rds"))
# 
for (x in 1:8) {

  print(date_list[x])

  ms_temp = rast(paste0(dir, "\\input\\MS_ORTHO_average\\", date_list[x], "_MS.tif"))

  # NDVI mask
  NDVI_mask = ndvi_mask[[x]]
  
  # shadow mask
  #NIR_mask = rast(paste0(dir, "\\output\\RASTER\\NIR_min_mask.tif"))
  NIR_mask = shadow_mask_all[[x]]

#   terra::writeRaster(NDVI, paste0(dir,"\\output\\RASTER\\", date_list[x], "\\Updated\\",  date_list[x], "_NDVI_scale.tif"),
#                      overwrite = TRUE)
#   
#   terra::writeRaster(NDVI_mask, paste0(dir,"\\output\\RASTER\\", date_list[x], "\\Updated\\",  date_list[x], "_NDVI_mask.tif"),
#                      overwrite = TRUE)
# 
#   # make a shadow mask from NIR reflectance values
#   NIR = ms_temp[[10]] %>%
#     ifel(NDVI_mask > threshold, ., NA) %>%
#     crop(pols_spat) %>%
#     clamp(upper = 50000, values = FALSE) %>%
#     as.vector()
# 
#   #saveRDS(NIR, "D:\\Sync\\_Sites\\Skimikin_spectral\\Skimikin_2020_07_02_NIR.rds")
# 
# 
#   shadow_patches = raster::clump(raster::raster(shadow_mask), directions = 4) %>%
#     rast()
# 
#   clumps = data.frame(freq(shadow_patches))
#   # threshold 200 cm^2
#   num_pix = 0.02 / (res(shadow_patches)[1]^2)
#   flecks = clumps[clumps$count > num_pix,] #remove clump observations with frequency smaller than 9
#   flecks = as.vector(flecks$value) # record IDs from clumps which met the criteria in previous step
# 
#   new_mask = shadow_patches %in% flecks
#   new_mask[new_mask == 0] = NA
# 
#   terra::writeRaster(new_mask, paste0(dir, "\\output\\RASTER\\", date_list[x], "\\",  date_list[x], "_NIR_shadow_mask2.tif"),
#                      overwrite = TRUE)
# 
# }
# 
# 
# ##---------------------------------------------------------------------------##
# # extract reflectance values
# 
# 
# pols = st_read(paste0(dir, "\\output\\HULLS\\Edits\\Skimikin_z50_updated.shp")) %>%
#   filter(!st_is_empty(.))
# 
# pols_spat = pols %>% 
#   terra::vect()
# 
# foreach(x = 1:9, 
#         .combine = 'c',
#         .packages = c("dplyr", "raster", "terra", "sf", "exactextractr", "lidR")
# ) %do% {
# 
#   ms = rast(paste0(dir, "\\input\\MS_ORTHO\\", date_list[x], "_MS.tif"))
# 
  ms_mask = terra::mask(ms_temp, pols_spat, touches = FALSE)
# 
#   terra::writeRaster(ms_mask, paste0(dir, "\\output\\RASTER\\", date_list[x], "\\",  date_list[x],  "_MS_mask.tif"),
#                      overwrite = TRUE)

  # mask_shadow = rast(paste0(dir,"\\output\\RASTER\\", date_list[x], "\\",  date_list[x], "_NIR_shadow_mask.tif"))
  # mask_ndvi = rast(paste0(dir,"\\output\\RASTER\\", date_list[x], "\\",  date_list[x], "_NDVI_mask.tif"))
  
  mask_combined = cover(NDVI_mask, NIR_mask)

  ms_mask_trim = trim(ms_mask)
  mask_trim = terra::crop(mask_combined, ms_mask_trim)

  # may not work with larger datasets. was crashing before being trimmed
  ms_mask2 = terra::mask(ms_mask_trim, mask_trim, maskvalues = 1, updatevalue = NA)
  ms_mask2[ms_mask2 == 0] = NA

  # quant = global(ms_mask2[[1:10]], fun=quantile, probs = c(.003, .997), na.rm = TRUE)
  # 
  # # # get_mode <- function(v) {
  # # #   uniqv <- unique(na.omit(v))
  # # #   uniqv[which.max(tabulate(match(v, uniqv)))]
  # # # }
  # 
  # ms_mask3 = ms_mask2 %>%
  #   clamp(lower = c(quant$X0.3.[[1]],
  #                   quant$X0.3.[[2]],
  #                   quant$X0.3.[[3]],
  #                   quant$X0.3.[[4]],
  #                   quant$X0.3.[[5]],
  #                   quant$X0.3.[[6]],
  #                   quant$X0.3.[[7]],
  #                   quant$X0.3.[[8]],
  #                   quant$X0.3.[[9]],
  #                   quant$X0.3.[[10]]),
  #         upper = c(quant$X99.7.[[1]],
  #                   quant$X99.7.[[2]],
  #                   quant$X99.7.[[3]],
  #                   quant$X99.7.[[4]],
  #                   quant$X99.7.[[5]],
  #                   quant$X99.7.[[6]],
  #                   quant$X99.7.[[7]],
  #                   quant$X99.7.[[8]],
  #                   quant$X99.7.[[9]],
  #                   quant$X99.7.[[10]]),
  #         values = FALSE)
  # 
  # ms_mask4 = tidyterra::drop_na(ms_mask3)

  terra::writeRaster(ms_mask2, paste0(dir, "\\output\\RASTER\\", date_list[x], "\\Updated\\",  date_list[x], "_MS_use.tif"),
                     overwrite = TRUE)
  
  ms_mask = rast(paste0(dir, "\\output\\RASTER\\", date_list[x], "\\Updated\\",  date_list[x], "_MS_use.tif"))
  

  ms_count = ms_mask[[4]]
  ms_count[ms_count >= 0] = 1
  
  
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
    
    
    
    # carotenoids
    CCI = (ms_mask[[3]] - ms_mask[[5]]) / (ms_mask[[3]] + ms_mask[[5]]),   
    PRI = (ms_mask[[3]] - ms_mask[[4]]) / (ms_mask[[3]] + ms_mask[[4]]),    
    
    # greenness, blueness
    GCC = (ms_mask[[3]] + ms_mask[[4]]) / (ms_mask[[1]]  + ms_mask[[2]] + ms_mask[[3]] + ms_mask[[4]] + ms_mask[[5]] + ms_mask[[6]]),  
    BCC = (ms_mask[[1]] + ms_mask[[2]]) / (ms_mask[[1]]  + ms_mask[[2]] + ms_mask[[3]] + ms_mask[[4]] + ms_mask[[5]] + ms_mask[[6]]), 
    
    # greenness, blueness
    GCC1 = (ms_mask[[4]]) / (ms_mask[[2]] + ms_mask[[4]] + ms_mask[[6]]),  
    BCC1 = ( ms_mask[[2]]) /(ms_mask[[2]] + ms_mask[[4]] + ms_mask[[6]]),  
    
    # carotenoids, waxes
    ARI = (1 / ms_mask[[4]]) - (1 / ms_mask[[7]]),    
    EWI9 = (ms_mask[[6]] - ms_mask[[8]]) / (ms_mask[[6]] + ms_mask[[8]]),
    
    # NIR greenness
    NDVI = (ms_mask[[10]] - ms_mask[[6]]) / (ms_mask[[10]] + ms_mask[[6]]),
    mDatt = (ms_mask[[10]] - ms_mask[[8]]) / (ms_mask[[10]] - ms_mask[[6]]),
    Datt = (ms_mask[[10]] - ms_mask[[7]]) / (ms_mask[[10]] - ms_mask[[6]]),   
    
    # the red edge
    RE_total = (ms_mask[[9]] - ms_mask[[7]]) / 35, 
    RE_upper = (ms_mask[[9]] - ms_mask[[8]]) / 23,
    NDRE3 = (ms_mask[[10]] - ms_mask[[9]]) / (ms_mask[[10]] + ms_mask[[9]]))
     
    
    
    # # normalized difference red
    # NDVI2 = (ms_mask[[10]] - ms_mask[[5]]) / (ms_mask[[10]] + ms_mask[[5]]), 
    # GNDVI = (ms_mask[[10]] - ms_mask[[4]]) / (ms_mask[[10]] + ms_mask[[4]]), 
    # NIRv = ((ms_mask[[10]] - ms_mask[[6]]) / (ms_mask[[10]] + ms_mask[[6]])) * ms_mask[[10]],  
    # EVI = 2.5 * ((ms_mask[[10]] - ms_mask[[6]]) / 
    #                (ms_mask[[10]] + (6 * ms_mask[[6]]) - (7.5 * ms_mask[[1]]) + 1)), 
    # MCARI = ((ms_mask[[7]] - ms_mask[[6]]) - (.2 * (ms_mask[[7]] - ms_mask[[4]])))
    # * (ms_mask[[7]]/ms_mask[[6]]),
    # TCARI = 3 * ((ms_mask[[7]] - ms_mask[[6]]) -
    #                  (.2 * (ms_mask[[7]] - ms_mask[[4]]) * (ms_mask[[7]]/ms_mask[[6]]))),
    # # /
    # #   ((1.16 * (ms_mask[[10]] - ms_mask[[6]])) / (ms_mask[[10]] + ms_mask[[6]] + .16)),
    # reNDVI = (ms_mask[[9]] - ms_mask[[7]]) / (ms_mask[[9]] + ms_mask[[7]]), 
    # GRVI = (ms_mask[[4]] - ms_mask[[6]]) / (ms_mask[[4]] + ms_mask[[6]]),   
    # # triangles
    # #TVI1 = 0.5 * ((668-560) * (ms_mask[[7]] - ms_mask[[4]]) - (705-560) * (ms_mask[[6]] - ms_mask[[4]])),
    # #TVI = 0.5 * ((842-717) * (ms_mask[[9]] - ms_mask[[8]]) - (740-717) * (ms_mask[[10]] - ms_mask[[8]])),
    # # STVI = (0.5 * ((842-717) * (ms_mask[[9]] - ms_mask[[8]]) - (740-717) * (ms_mask[[10]] - ms_mask[[8]]))
    # #         - (0.5 * ((668-560) * (ms_mask[[7]] - ms_mask[[4]]) - (705-560) * (ms_mask[[6]] - ms_mask[[4]])))) 
    # # / (0.5 * ((842-717) * (ms_mask[[9]] - ms_mask[[8]]) - (740-717) * (ms_mask[[10]] - ms_mask[[8]]))
    # #    + (0.5 * ((668-560) * (ms_mask[[7]] - ms_mask[[4]]) - (705-560) * (ms_mask[[6]] - ms_mask[[4]])))), 
    # 
    # 
    # # normalized difference red edge 
    # NDRE1 = (ms_mask[[10]] - ms_mask[[7]]) / (ms_mask[[10]] + ms_mask[[7]]),     
    # NDRE2 = (ms_mask[[10]] - ms_mask[[8]]) / (ms_mask[[10]] + ms_mask[[8]]), 
    # 
    # MRESR = (ms_mask[[9]] - ms_mask[[1]]) / (ms_mask[[7]] + ms_mask[[1]]), 
    # NIRvNDRE = ((ms_mask[[10]] - ms_mask[[8]]) / (ms_mask[[10]] + ms_mask[[8]])) * ms_mask[[10]],  
    # CI_green = (ms_mask[[10]] / ms_mask[[4]]) - 1,
    # RG = ms_mask[[6]] / ms_mask[[4]],
    # BG1 = ms_mask[[2]] / ms_mask[[3]],
    # BG2 = ms_mask[[1]] / ms_mask[[3]],
    # BG3 = ms_mask[[2]] / ms_mask[[4]],
    # BG4 = ms_mask[[1]] / ms_mask[[4]],
    # #BB = ms_mask[[1]] / ms_mask[[2]],
    # 
    #  
    # 
    # # green chromatic coordinate
    # Gcc = ms_mask[[4]] / (ms_mask[[2]] + ms_mask[[4]] + ms_mask[[5]]), 
    # Gcc1 = ms_mask[[4]] / (ms_mask[[2]] + ms_mask[[4]] + ms_mask[[6]]),  
    # Gcc2 = ms_mask[[3]] / (ms_mask[[1]] + ms_mask[[3]] + ms_mask[[5]]),  
    # #Gcc3 = ms_mask[[3]] / (ms_mask[[2]] + ms_mask[[4]] + ms_mask[[5]]),  
    # 
    # 

    # #Bcc2 = ms_mask[[1]] / (ms_mask[[1]] + ms_mask[[3]] + ms_mask[[5]]), 
    # 
    # # carotenoids             
    # SIPI = (ms_mask[[10]] - ms_mask[[1]]) / (ms_mask[[10]] - ms_mask[[6]]),   
    # SIPI2 = (ms_mask[[10]] - ms_mask[[2]]) / (ms_mask[[10]] - ms_mask[[6]]),   
    # NIRvCCI = (ms_mask[[3]] - ms_mask[[5]]) / (ms_mask[[3]] + ms_mask[[5]]) * ms_mask[[10]],  
    # 
    # MARI = ((1 / ms_mask[[4]]) - (1 / ms_mask[[7]])) * ms_mask[[10]],
    # 
    # # waxes
    # EWI4 = ((ms_mask[[2]] * ms_mask[[2]]) - ms_mask[[6]]) / (ms_mask[[2]] - (ms_mask[[6]] * ms_mask[[6]])),
    # 
    # 
    # # Red edge  
    # ISI = ms_mask[[10]] - ms_mask[[9]],
    # RE_lower = (ms_mask[[8]] - ms_mask[[7]]) / 12,   
    # #RE_lower2 = (ms_mask[[8]] - ms_mask[[6]]) / 49,   
    # 

    # 
  
  rast_all = rast(rast_list)
  
  count = exactextractr::exact_extract(ms_count, pols, fun = "sum", append_cols = "Obs") %>% 
    dplyr::mutate(count = round(sum)) %>% 
    dplyr::select(-sum)  
  
  df_spectral = exact_extract(rast_all, pols, fun = "median", append_cols = "Obs") %>% 
    left_join(count, by = "Obs")

  
  saveRDS(df_spectral, paste0(dir, "\\CSV\\TRAITS\\Updated\\", date_list[x], "_spectral_average_median.rds"))
}

################################################################################

# all ITC spectral values
df_spectral_all = bind_rows(
  # read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2020_03_22_spectral.rds") %>%
  #   mutate(Date = ymd("2020-03-22"),
  #          timepoint = "LW"),
  # read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2020_07_02_spectral_min.rds") %>% 
  #   mutate(Date = ymd("2020-07-02"),
  #          timepoint = "ES"),
  # read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2020_08_07_spectral_min.rds") %>% 
  #   mutate(Date = ymd("2020-08-07"),
  #          timepoint = "MS"),
  # read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2020_08_25_spectral_min.rds") %>% 
  #   mutate(Date = ymd("2020-08-25"),
  #          timepoint = "LS"),
  # read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2021_03_31_spectral_min.rds") %>% 
  #   mutate(Date = ymd("2021-03-31"),
  #          timepoint = "LW"),
  # read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2021_06_29_spectral_min.rds") %>% 
  #   mutate(Date = ymd("2021-06-29"),
  #          timepoint = "ES"),
  # read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2021_07_29_spectral_min.rds") %>% 
  #   mutate(Date = ymd("2021-07-29"),
  #          timepoint = "MS"),
  # read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2021_08_14_spectral_min.rds") %>% 
  #   mutate(Date = ymd("2021-08-14"),
  #          timepoint = "LS"),
  # read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2023_02_24_spectral_min.rds") %>% 
  #   mutate(Date = ymd("2023-02-24"),
  #          timepoint = "LW")) %>% 
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2020_07_02_spectral_average_median.rds") %>%
  mutate(Date = ymd("2020-07-02"),
         timepoint = "ES"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2020_08_07_spectral_average_median.rds") %>%
  mutate(Date = ymd("2020-08-07"),
         timepoint = "MS"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2020_08_25_spectral_average_median.rds") %>%
  mutate(Date = ymd("2020-08-25"),
         timepoint = "LS"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2021_03_31_spectral_average_median.rds") %>%
  mutate(Date = ymd("2021-03-31"),
         timepoint = "LW"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2021_06_29_spectral_average_median.rds") %>%
  mutate(Date = ymd("2021-06-29"),
         timepoint = "ES"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2021_07_29_spectral_average_median.rds") %>%
  mutate(Date = ymd("2021-07-29"),
         timepoint = "MS"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2021_08_14_spectral_average_median.rds") %>%
  mutate(Date = ymd("2021-08-14"),
         timepoint = "LS"),
  read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_2023_02_24_spectral_average_median.rds") %>%
  mutate(Date = ymd("2023-02-24"),
         timepoint = "LW")) %>%
  rename_with(~ gsub("median.", "", .x, fixed = TRUE)) %>% 
  filter(count >= 10) %>% 
  pivot_longer(R444:NDRE3, names_to = "index", values_to = "value") %>% 
  dplyr::select(-timepoint, -count) %>% 
  mutate(value = if_else(is.nan(value), NA_real_, value)) %>% 
  pivot_wider(names_from = "Date", values_from = "value")
  # drop_na() %>% 
  # mutate(mean_summer = (`2020-07-02`+`2020-08-07`+
  #                         `2021-06-29`+`2021-07-29`)/4,
  #        mean_winter = (`2020-03-22`+`2021-03-31`+`2023-02-24`)/3,
  #        mean_LS = (`2020-08-25`+`2021-08-14`)/2,
  #        greenup = mean_summer - mean_winter,
  #        decline = mean_LS - mean_summer)

saveRDS(df_spectral_all, paste0(dir, "\\CSV\\TRAITS\\Updated\\Skimikin_all_spectral_average_median.rds"))

################################################################################
# 
# # all ITC spectral values
# df_spectral_all = 
#   read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2020_03_22_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2020_LW"), -Obs) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2020_07_02_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2020_ES"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2020_08_07_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2020_MS"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2020_08_25_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2020_LS"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2021_03_31_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2021_LW"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2021_06_29_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2021_ES"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2021_07_29_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2021_MS"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2021_08_14_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2021_LS"), -Obs)) %>% 
#   left_join(read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_2023_02_24_spectral.rds") %>% 
#     rename_with(~ gsub("mean.", "", .x, fixed = TRUE)) %>% 
#     rename_with(~paste0(.,"_2023_LW"), -Obs)) %>% 
#   
#   
# 
# saveRDS(df_spectral_all, paste0(dir, "\\CSV\\TRAITS\\Skimikin_all_spectral.rds"))

