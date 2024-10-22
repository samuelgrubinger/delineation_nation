
library(terra)
library(exifr)
library(tidyverse)
library(ggplot2)
library(raster)
library(lubridate)
library(exifr)
library(pracma)


dir_list = c(
  "D:\\Sync\\_Sites\\Sx_Genecology\\Harrison\\DAP\\MS",
  "D:\\Sync\\_Sites\\Sx_Genecology\\JordanRiver\\DAP\\MS",
  "D:\\Sync\\_Sites\\Sx_Genecology\\Kalamalka\\DAP\\MS",
  #"D:\\Sync\\_Sites\\Sx_Genecology\\Skimikin\\DAP\\MS",
  "D:\\Sync\\_Sites\\Sx_Genecology\\TJC\\DAP\\MS",
  "D:\\Sync\\_Sites\\Sx_Genecology\\Whitecourt\\DAP\\MS")

for (j in 4:length(dir_list)
     ){

  dir = dir_list[j]

  #dir to calibration photos
  dir_pan = paste0(dir, "\\Panels\\")
  
  #directory to folder containing masks, masks in the folder should be named identical to the images they mask with a _mask. Ie img_1, img_1_mask
  dir_mask = paste0(dir, "\\Masks\\")
  
  #Creating a list of all images in the calibration panel folder to iterate over
  img_list = list.files(dir_pan, recursive = TRUE, full.names = FALSE, pattern = ".tif")
  
  #Masking out panels in cal photos, calculating mean panel reflectance and concatenating all exif data into one dataframe
  for (i in 1:length(img_list)){
    print(i)
    
    
    #----------------------------------------------------------------
    #raster from image
    img = img_list[i]
    rast = rast(paste0(dir_pan,img))
    #plot(rast)
    
    #root name to match image to mask
    img_root <- substr(img_list[i],1,nchar(img_list[i])-4)
    img_root
    
    img_exif = read_exif(paste0(dir_pan, img))
    
    darkLevel = img_exif$BlackLevel %>% 
      str_split(" ") %>% 
      lapply(as.numeric) %>% 
      unlist() %>% 
      mean(na.rm = TRUE)
    
    cal = img_exif$RadiometricCalibration 
    a1 = cal[[1]][1] %>% as.numeric()
    a2 = cal[[1]][2] %>% as.numeric()
    a3 = cal[[1]][3] %>% as.numeric()
    
    #R = 1.0 / (1.0 + a2 * y / exposureTime - a3 * y)
    
    #distance from vignette center
    cent = img_exif$VignettingCenter
    vpoly = img_exif$VignettingPolynomial %>% 
      lapply(as.numeric) %>% 
      unlist()
    
    cent_vect = data.frame(x = cent[[1]][1],
                          y = cent[[1]][2]) %>% 
      vect(geom = c("x", "y"), crs = "epsg:26910")
      
    crs(rast) = "epsg:26910"
    # vignetting correction raster
    dist_rast = distance(rast, cent_vect)
    poly_rast = dist_rast^6 * vpoly[6] + 
      dist_rast^5 * vpoly[5] +
      dist_rast^4 * vpoly[4] +
      dist_rast^3 * vpoly[3] +
      dist_rast^2 * vpoly[2] +
      dist_rast * vpoly[1] +
      1
    
    V = 1 / poly_rast
    
    # row gradient correction
    y = rast
    values(y) = rep(seq(1, nrow(rast), 1), 
                                  each = ncol(rast))
    
    exposureTime = img_exif$ExposureTime
    gain = img_exif$ISOSpeed/100.0
    
    R = 1.0 / (1.0 + a2 * y / exposureTime - a3 * y)
    
    L = V * R * (rast - darkLevel)
    
    L[L < 0] = 0
    
    # apply the radiometric calibration - 
    # scale by the gain-exposure product and 
    #multiply with the radiometric calibration coefficient
    bitsPerPixel = img_exif$BitsPerSample
    dnMax = 2^bitsPerPixel
    radianceImage = L/(gain * exposureTime)*a1/dnMax
  #---------------------------------------------------------------
    
    #mask to raster
    mask_name = list.files(paste0(dir_mask), 
                            pattern = paste0(img_root,"+_mask.png$"))
    mask_name
    mask = rast(paste0(dir_mask, mask_name))
    mask[mask == 0] <- NA
    #plot(mask)
    
    #masking out panel
    panel = mask(radianceImage, mask)
    panel
    #plot(panel)
    
    #getting mean reflectance value of panel
    val = values(panel) #get raster values
    mean_ref = mean(val, na.rm=T)
    mean_ref
    
    #exif data of image, plus adding 3 data columns
    img_exif = img_exif %>%
      mutate(prop_scattered = ScatteredIrradiance / (ScatteredIrradiance + DirectIrradiance),
             prop_direct = DirectIrradiance / (ScatteredIrradiance + DirectIrradiance),
             mean_panel_radiance = mean_ref)
    
    #creating df with same column names as exif data
    if (i == 1){
      panel_exif_df <- as.data.frame(img_exif)
      #columns <- c(colnames(img_exif))
      #print(length(columns))
      #panel_exif_df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
      #panel_exif_df <- setNames(panel_exif_df, columns)
    } else {   #adding each new image exif to dataframe
      panel_exif_df <- merge(panel_exif_df, img_exif, by = intersect(names(panel_exif_df), names(img_exif)), all = TRUE)
      print(panel_exif_df)
    }
  }

  if (j == 1){
    all_exif_df = panel_exif_df
    #columns <- c(colnames(img_exif))
    #print(length(columns))
    #panel_exif_df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
    #panel_exif_df <- setNames(panel_exif_df, columns)
  } else {   #adding each new image exif to dataframe
    all_exif_df <- merge(all_exif_df, panel_exif_df, by = intersect(names(all_exif_df), names(panel_exif_df)), all = TRUE)
    print(all_exif_df)
  }
  }

#Plotting
#From Sam

# exif_07_20 = panel_exif_df %>% 
#   mutate(panel_DN = mean_panel_reflectance) 
# 
# panel_plot = exif_07_20
# 


#direct irradiance of panels verse time
(pan_plot = all_exif_df %>% 
  mutate(panel_irradiance = (mean_panel_radiance * pi) /
           (HorizontalIrradiance * .01)) %>% 
  group_by(ModifyDate) %>% 
  mutate(mean_panel_irr = mean(panel_irradiance)) %>% 
  ggplot(aes(x = ModifyDate, y = panel_irradiance, group = BandName)) +
  geom_point(aes(colour = BandName), size = 3) +
  geom_line(aes(colour = BandName), size = 1, alpha = .5) +
  geom_line(aes(y = mean_panel_irr), color = "black", size = 1) +
  geom_hline(yintercept = 0.538, color = "grey10") +
  theme_bw(base_size = 8) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values = c("Blue-444" = "royalblue4",
                                "Blue" = "steelblue",
                                "Green-531" = "springgreen3",
                                "Green" = "forestgreen",
                                "Red-650" = "firebrick2",
                                "Red" = "red3",
                                "Red edge-705" = "indianred",
                                "Red edge" = "lightpink3",
                                "Red edge-740" = "pink4",
                                "NIR" = "thistle4")) +
  facet_wrap(. ~ Directory, scales = "free"))

ggsave(plot = pan_plot, 
       filename = "D:\\Sync\\_Sites\\Sx_Genecology\\Sx_choose_panels.tiff", 
       device = tiff,
       width = 15,
       height = 10,
       units = "in",
       dpi = 600,
       bg = "white")

################################################################################
# look at all images to chuck bad ones

pic_path = list.files(dir, recursive = FALSE, full.names = TRUE, pattern = ".tif")

pic_exif = read_exif(pic_path)

pic_exif_2 = pic_exif %>% 
  mutate(prop_scattered = ScatteredIrradiance / HorizontalIrradiance,
         prop_scattered2 = ScatteredIrradiance / SpectralIrradiance
  )


#cleaning dataframe and creating diff_prop columns
exif_plot = pic_exif_2 %>% 
  drop_na(ScatteredIrradiance) %>% 
  mutate(Date_mod = as.factor(FileModifyDate)) %>% 
  arrange(FileModifyDate) %>% 
  group_by(Date_mod) %>% 
  mutate(prop_scat = mean(prop_scattered_i)) %>% 
  group_by(BandName) %>% 
  mutate(diff_prop1 = prop_scattered - lag(prop_scattered, k = 1, default = first(prop_scattered)),
         diff_prop2 = prop_scattered - lead(prop_scattered, k = 1, default = first(prop_scattered)),
         #        direc = as.factor(Directory))
         diff_irr1 = Irradiance / lag(Irradiance, k = 1, default = first(Irradiance)),
         diff_irr2 = Irradiance / lead(Irradiance, k = 1, default = first(Irradiance)))

exif_plot %>%
  filter(BandName == "Blue") %>% 
  ggplot(aes(x = ymd_hms(FileModifyDate), y = prop_scattered)) +
  # geom_point(aes(y = ScatteredIrradiance), color = "blue", size = 3) +
  # geom_point(aes(y = DirectIrradiance), color = "red", size = 3) +
  # geom_point(aes(y = HorizontalIrradiance), color = "purple", size = 3) +
  geom_point(aes(y = prop_scattered), color = "purple", size = 3) +
  geom_point(aes(y = prop_scattered2), color = "red", size = 3) +
  theme_bw()

see = exif_plot %>% filter(BandName == "Red" & 
                             (diff_irr2 > 1.02 | diff_irr2 < .98 |
                                diff_irr1 > 1.02 | diff_irr1 < .98))


exif_plot %>%
  filter(BandName == "Blue") %>% 
  ggplot(aes(x = ymd_hms(FileModifyDate), y = prop_scattered)) +
  geom_point(color = "indianred", size = 5) +
  geom_point(data = filter(exif_plot, BandName == "Blue" & prop_scat < .2),
             color = "purple") +
  geom_point(data = filter(exif_plot, BandName == "Blue" &
                             (abs(diff_prop1) > .01 | abs(diff_prop2) > .01)), 
             color = "steelblue3", size = 5) +
  # geom_point(data = filter(exif_plot, FileName %in% c(
  #   "IMG_0023_6.tif",
  #   "IMG_0025_6.tif",
  #   "IMG_0028_6.tif",
  #   "IMG_0664_6.tif",
  #   "IMG_0665_6.tif",
  #   # "IMG_0028_6.tif",
  #   # "IMG_0060_6.tif",
  #   "IMG_0666_6.tif")),
  #   color = "blue", size = 5, shape = 21) +
  #geom_line() +
# geom_point(data = filter(exif_plot, diff_abs2 < 1 | diff_abs < 1),
#            color = "blue") +
#geom_line(aes(y = prop_direct), size = 2) +
theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_datetime(breaks = scales::date_breaks("1 mins"), date_labels = "%H:%M:%S") 




####################################################################################################################################################  

see2 = exif_plot %>% filter(BandName == "Red" & 
                              #(diff_abs > 1 | diff_abs2 > 1))
                              (abs(diff_prop1) > .01 | 
                                 abs(diff_prop2) > .01))

see

see_blue = exif_plot %>% filter(BandName == "Blue-444" & 
                                  #(diff_abs > 1 | diff_abs2 > 1))
                                  (prop_scat < .3 |
                                     abs(diff_prop1) > .02 | 
                                     abs(diff_prop2) > .02))


################################################################################
#Sams code
# 
# pic_path_all = list.files("D:\\Sync\\_Sites\\Skimikin_spectral\\DAP\\MS\\Skimikin_2020_08_07", recursive = TRUE, full.names = TRUE)
# 
# pic_path_site_harr = list.files("E:\\_RawData\\SX_field_data\\Harrison\\August_23_2020", recursive = TRUE, full.names = TRUE)
# 
# pic_exif_skim = read_exif("D:\\Sync\\_Sites\\Skimikin_spectral\\DAP\\MS\\Skimikin_2020_08_07\\_1", recursive = TRUE) %>% 
#   mutate(prop_scattered = ScatteredIrradiance / (ScatteredIrradiance + DirectIrradiance),
#          prop_direct = DirectIrradiance / (ScatteredIrradiance + DirectIrradiance)
#   )
# 
# pic_exif_kal = read_exif('D:\\Sync\\_Sites\\Sx_Genecology\\Kalamalka\\DAP\\MS',
#                          recursive = TRUE) %>% 
#   mutate(prop_scattered = ScatteredIrradiance / (ScatteredIrradiance + DirectIrradiance),
#          prop_direct = DirectIrradiance / (ScatteredIrradiance + DirectIrradiance)
#   )
# 
pic_exif_harr = read(list.files("E:\\_RawData\\SX_field_data\\Harrison\\August_23_2020", recursive = TRUE, full.names = TRUE)) %>%
  # read_exif('D:\\Sync\\_Sites\\Sx_Genecology\\Harrison\\DAP\\MS',
  #                         recursive = TRUE) %>%
  mutate(prop_scattered = ScatteredIrradiance / (ScatteredIrradiance + DirectIrradiance),
         prop_direct = DirectIrradiance / (ScatteredIrradiance + DirectIrradiance)
  )
# 
# pic_exif_tjc = read_exif('D:\\Sync\\_Sites\\Sx_Genecology\\TJC\\DAP\\MS',
#                          recursive = TRUE) %>% 
#   mutate(prop_scattered = ScatteredIrradiance / (ScatteredIrradiance + DirectIrradiance),
#          prop_direct = DirectIrradiance / (ScatteredIrradiance + DirectIrradiance)
#   )
# 
# pic_exif_whit = read_exif('D:\\Sync\\_Sites\\Sx_Genecology\\Whitecourt\\DAP\\MS',
#                           recursive = TRUE) %>% 
#   mutate(prop_scattered = ScatteredIrradiance / (ScatteredIrradiance + DirectIrradiance),
#          prop_direct = DirectIrradiance / (ScatteredIrradiance + DirectIrradiance)
#   )
# 
# pic_exif_jr = read_exif('D:\\Sync\\_Sites\\Sx_Genecology\\JordanRiver\\DAP\\MS',
#                         recursive = TRUE) %>% 
#   mutate(prop_scattered = ScatteredIrradiance / (ScatteredIrradiance + DirectIrradiance),
#          prop_direct = DirectIrradiance / (ScatteredIrradiance + DirectIrradiance)
#   )
# 
# 
# pic_harr_panels = read_exif("D:\\Sync\\_Sites\\Sx_Genecology\\Harrison\\DAP\\MS\\Panels",
#                             recursive = TRUE) %>% 
#   mutate(prop_scattered = ScatteredIrradiance / (ScatteredIrradiance + DirectIrradiance),
#          prop_direct = DirectIrradiance / (ScatteredIrradiance + DirectIrradiance))
# 
# 
# # pic_exif %>% 
# #   ggplot(aes(x = FileModifyDate, y = DirectIrradiance, color = BandName)) +
# #   geom_point() +
# #   geom_line() +
# #   theme_bw(base_size = 18)

exif_plot  = exif_03_20 %>% 
  drop_na(ScatteredIrradiance) %>% 
  mutate(Date_mod = as.factor(FileModifyDate)) %>% 
  arrange(FileModifyDate) %>% 
  group_by(Date_mod) %>% 
  mutate(prop_scat = mean(prop_scattered)) %>% 
  group_by(BandName) %>% 
  mutate(diff_prop1 = prop_scattered - lag(prop_scattered, k = 1, default = first(prop_scattered)),
         diff_prop2 = prop_scattered - lead(prop_scattered, k = 1, default = first(prop_scattered)),
         direc = as.factor(Directory))


exif_plot %>%
  filter(BandName == "NIR") %>% 
  # filter(FileName %in% c("IMG_0022_6.tif",
  #                                 "IMG_0023_6.tif",
  #                                 "IMG_0024_6.tif",
  #                                 "IMG_0025_6.tif",
  #                                 "IMG_0026_6.tif",
  #                                 "IMG_0027_6.tif",
  #                                 "IMG_0028_6.tif",
  #                                 "IMG_0029_6.tif")) %>% 
  ggplot(aes(x = ymd_hms(FileModifyDate), y = DirectIrradiance)) +
  geom_point(color = "indianred", size = 5) +
  geom_point(data = filter(exif_plot, prop_scat < .3),
             color = "purple") +
  geom_point(data = filter(exif_plot, 
                           abs(diff_prop1) > .02 | abs(diff_prop2) > .02), 
             color = "steelblue3", size = 5) +
  geom_point(data = filter(exif_plot, FileName %in% c(
    "IMG_0023_6.tif",
    "IMG_0025_6.tif",
    "IMG_0028_6.tif",
    "IMG_0664_6.tif",
    "IMG_0665_6.tif",
    # "IMG_0028_6.tif",
    # "IMG_0060_6.tif",
    "IMG_0666_6.tif")),
    color = "blue", size = 5, shape = 21) +
  #geom_line() +
  # geom_point(data = filter(exif_plot, diff_abs2 < 1 | diff_abs < 1),
  #            color = "blue") +
  #geom_line(aes(y = prop_direct), size = 2) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_datetime(breaks = scales::date_breaks("1 mins"), date_labels = "%H:%M:%S") 
+
  facet_wrap(. ~ direc, scales = "free_x")

see = exif_plot %>% filter(BandName == "Red" & 
                             #(diff_abs > 1 | diff_abs2 > 1))
                             (prop_scat < .3 |
                                abs(diff_prop1) > .02 | 
                                abs(diff_prop2) > .02))

see_blue = exif_plot %>% filter(BandName == "Blue-444" & 
                                  #(diff_abs > 1 | diff_abs2 > 1))
                                  (prop_scat < .3 |
                                     abs(diff_prop1) > .02 | 
                                     abs(diff_prop2) > .02))