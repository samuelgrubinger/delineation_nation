
library(tidyverse)
library(sf)
library(sp)
library(spatial)
library(rgdal)
library(rgeos)
library(raster)
library(terra)
library(lidR)
library(sp)
library(nngeo)
library(future)
library(rmapshaper)
library(concaveman)
library(parallel)
library(foreach)
library(smoothr)
library(ForestTools)
library(rmapshaper)
library(gdalUtilities)
library(exactextractr)

site = "Harrison"
name = "Harrison"

dir = "D:/Sync/_Sites/Sx_Genecology/Harrison"

# Create output file directories if they don't exist

if (!dir.exists(paste0(dir, "\\output\\CROWNS"))) {
  dir.create(paste0(dir, "\\output\\CROWNS"))
}

if (!dir.exists(paste0(dir, "\\output\\GROUND"))) {
  dir.create(paste0(dir, "\\output\\GROUND"))
}

if (!dir.exists(paste0(dir, "\\output\\HULLS"))) {
  dir.create(paste0(dir, "\\output\\HULLS"))
}

if (!dir.exists(paste0(dir, "\\output\\MERGED"))) {
  dir.create(paste0(dir, "\\output\\MERGED"))
}

if (!dir.exists(paste0(dir, "\\output\\NORM"))) {
  dir.create(paste0(dir, "\\output\\NORM"))
}

if (!dir.exists(paste0(dir, "\\output\\NORM_clean"))) {
  dir.create(paste0(dir, "\\output\\NORM_clean"))
}

if (!dir.exists(paste0(dir, "\\output\\RASTER"))) {
  dir.create(paste0(dir, "\\output\\RASTER"))
}

if (!dir.exists(paste0(dir, "\\output\\SEGMENTED"))) {
  dir.create(paste0(dir, "\\output\\SEGMENTED"))
}

########################

plan(multisession, workers = 6L)

# 20x20m tiles output from Metashape
CTG = readLAScatalog(folder = paste0(dir, "\\input\\TILES\\"))
opt_chunk_buffer(CTG) = 15
opt_laz_compression(CTG) = FALSE
# selecting only confidence above 2 will help with belowground noise
# taking a random 5% of the point cloud speeds things up
opt_filter(CTG) = "-keep_random_fraction 0.03 -keep_attribute_above 0 1"
opt_output_files(CTG) = paste0(dir, "\\output\\GROUND\\{*}_GROUND")
opt_progress(CTG) = TRUE


# classify ground from a point cloud
# these parameters are for dense DAP data on relatively flat terrain 


GROUND = classify_ground(CTG, algorithm = csf(
  sloop_smooth = FALSE,
  class_threshold = 0.07,
  # larger resolution = coarser DTM
  cloth_resolution = .7,
  rigidness = 2L,
  # even for non flat sites, this should be 3L in my experience
  iterations = 500L,
  time_step = 0.65)) 
# higher time step = more ground covered

# Read in the ground tiles
GROUND = readLAScatalog(folder = paste0(dir, "\\output\\GROUND\\"))
opt_chunk_buffer(GROUND) = 20
opt_laz_compression(CTG) = TRUE
opt_output_files(GROUND) = ""
opt_progress(GROUND) = TRUE

# Create a DTM from the ground points
# smooth the terrain model a bit
DTM = grid_terrain(GROUND, res = .1, tin(), full_raster = FALSE) %>% 
  focal(w = matrix(1, 25, 25),
        fun = mean,
        na.rm = TRUE,
        pad = TRUE)


# save the DTM 
writeRaster(DTM, paste0(dir, "\\output\\RASTER\\", name, "_DTM.tif"), overwrite = TRUE) 
#DTM = raster(paste0(dir, "\\output\\RASTER\\", name, "_DTM.tif"))

# read the tiles again, but with high-res parameters
CTG = readLAScatalog(folder = paste0(dir, "\\input\\TILES\\"))
# small buffer; we're just normalizing
opt_chunk_buffer(CTG) = .5
opt_laz_compression(CTG) = TRUE
# 1cm voxelization
opt_filter(CTG) = "-thin_with_voxel 0.02"
opt_output_files(CTG) = paste0(dir, "\\output\\NORM\\{*}_NORM")
opt_progress(CTG) = TRUE

# normalize the heights
NORM = normalize_height(CTG, DTM, na.rm = TRUE)

# now read them in again and filter out points below -0.5
NORM = readLAScatalog(paste0(dir, "\\output\\NORM\\"))
# no buffer, just filtering
opt_chunk_buffer(NORM) = 0
opt_laz_compression(NORM) = TRUE
opt_filter(NORM) = "-drop_z_below -.25"
opt_output_files(NORM) = paste0(dir, "\\output\\NORM_clean\\{*}")
opt_progress(NORM) = TRUE

# now write over the NORM files but filter points less than -.5
NORM_clean = catalog_retile(NORM)

# read in norm cloud
NORM = readLAScatalog(paste0(dir, "\\output\\NORM_clean\\"))
opt_chunk_buffer(NORM) = .5
opt_laz_compression(NORM) = TRUE
opt_output_files(NORM) = ""
opt_progress(NORM) = TRUE

# create a CHM for tree crown delineation
# start with a simple 2cm max grid
CHM_max = grid_metrics(NORM,
                       res = .03,
                       func = max(Z))

writeRaster(CHM_max, paste0(dir, "\\output\\RASTER\\", name, "_CHM_max.tif"), overwrite = TRUE)
CHM_max = raster(paste0(dir, "\\output\\RASTER\\", name, "_CHM_max.tif"))

# smooth the max raster with a gaussian blur
CHM_max_smooth = CHM_max %>% 
  focal(w = focalWeight(CHM_max, .03, "Gauss"), na.rm = TRUE, NAonly = FALSE)

# set any NAs from the original CHM to 0
# these are areas of no points and likely gaps between crowns
CHM_max_smooth[is.na(CHM_max)] = 0

writeRaster(CHM_max_smooth, paste0(dir, "\\output\\RASTER\\", name, "_CHM_max_smooth.tif"), overwrite = TRUE)

future:::ClusterRegistry("stop")


##################################################
# read in a shapefile of gridded points located within 
# 40 cm of each tree apex,
# with census information appended
# buffer trees must have points with Row and Col == 0
# rowid will provide a unique id for joining later


grid_sf = st_read(paste0(dir, "\\GIS\\", site, "_grid_manual.shp")) %>% 
  st_set_crs(26910) %>% 
  mutate(Absent = replace_na(Absent, 0)) %>% 
  dplyr::filter(Absent != 1) %>% 
  rowid_to_column()

# load the smoothed CHM
CHM_max_smooth = raster(paste0(dir, "\\output\\RASTER\\", name, "_CHM_max_smooth.tif"))

buffers = st_buffer(grid_sf, dist = .4) %>% 
  as_Spatial() %>% 
  gdalUtilities::gRasterize(CHM_max_smooth, 
                            field = "Obs")


chm_mask = CHM_max_smooth
chm_mask[is.na(buffers)] = NA

# we will search for local maxima within a buffer of these points,
# and "snap" to the correct nearby maximum for each point
tops = find_trees(chm_mask, algorithm = lmf(ws = .4,
                                            hmin = 0,
                                            shape = "circular")) %>% 
  remove.duplicates() %>% 
  st_as_sf() %>% 
  st_set_crs(26910)

# here we snap the census data grid to the true local maxima
joined1 = st_join(tops, grid_sf, join = st_nearest_feature) %>%
  group_by(rowid) %>% 
  filter(Z == max(Z))


# save the new treetops if you want
st_write(obj = joined1,
         dsn = paste0(dir, "\\output\\HULLS\\", name, "_treetops.shp"),
         append = FALSE)

# read in the treetops as an sp object
ttops_sp = st_read(paste0(dir, "\\output\\HULLS\\", name, "_treetops.shp")) %>% 
  as_Spatial()

# apply the marker-controlled watershed algorithm,
# faster when OSGeo64 path used
ws = ForestTools::mcws(CHM = CHM_max_smooth,
                       treetops = ttops_sp,
                       OSGeoPath = "C:\\OSGeo4W64",
                       format = "raster",
                       minHeight = 0.2)

ws_max = raster::subs(x = ws, 
                      y = as.data.frame(zonal(CHM_max_smooth, 
                                              ws, 
                                              max, 
                                              na.rm = TRUE)))

chm_mask = CHM_max_smooth
chm_mask[chm_mask < (0.5 * ws_max)] = NA


ws2 = ForestTools::mcws(CHM = chm_mask,
                        treetops = ttops_sp,
                        OSGeoPath = "C:\\OSGeo4W64",
                        format = "polygons",
                        minHeight = 0.2)

ws2_smooth = st_as_sf(ws2) %>% 
  fill_holes(threshold = units::set_units(2000, cm^2)) 

ws3 = ws2_smooth %>%
  dplyr::filter(Obs != 0)

st_write(obj = ws3, 
         dsn = paste0(dir, "\\output\\HULLS\\", name, "_z50_new.shp"),
         layer = name,
         driver = "ESRI Shapefile",
         append = FALSE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# produce a hillshade mask 

library(rayshader)
library(suncalc)
library(exifr)
library(lubridate)
library(rgdal)


# make a mean values CHM
# fill NA with a minimum filter,
# smooth aggressively
# read in norm cloud
plan(multisession, workers = 3L)

# read in norm cloud
CTG = readLAScatalog(paste0(dir, "\\input\\TILES_untrimmed\\"))
opt_chunk_buffer(CTG) = .5
opt_filter(CTG) = "-thin_with_voxel 0.02"
opt_output_files(CTG) = ""
opt_progress(CTG) = TRUE

# create a DSM for tree crown delineation
# start with a simple 2cm max grid
DSM_max = pixel_metrics(CTG,
                        res = .04,
                        func = ~max(Z))

writeRaster(DSM_max, paste0(dir, "\\output\\RASTER\\", name, "_DSM_max.tif"), overwrite = TRUE)

future:::ClusterRegistry("stop")


################################################################################


pols = st_read(paste0(dir, "\\output\\HULLS\\", name, "_z50_snap.shp")) 

pols_bound = pols %>% 
  st_centroid() %>% 
  st_union() %>% 
  st_convex_hull() %>%
  st_buffer(12) %>% 
  st_sf() %>% 
  vect()

DSM_max = rast(paste0(dir, "\\output\\RASTER\\", name, "_DSM_max.tif")) %>% 
  crop(pols_bound, mask = TRUE) 

ms_temp = rast(paste0(dir, "\\input\\MS_ORTHO\\", name, "_MS.tif"))

DSM_max_fill = DSM_max %>% 
  terra::focal(w = 3, fun = "min", na.policy = "only", na.rm = TRUE) %>% 
  terra::focal(w = 3, fun = "median", na.policy = "only", na.rm = TRUE) %>% 
  terra::focal(w = 5, fun = "min", na.policy = "only", na.rm = TRUE) %>% 
  terra::focal(w = 5, fun = "median", na.policy = "only", na.rm = TRUE) %>% 
  ifel(terra::terrain(., "TRI") > 1,
       terra::focal(., w = 3, fun = "median", na.policy = "omit", na.rm = TRUE),
       .) %>% 
  ifel(terra::terrain(., "TRI") > .5,
       terra::focal(., w = 3, fun = "median", na.policy = "omit", na.rm = TRUE),
       .)

terra::writeRaster(DSM_max_fill, paste0(dir, "\\output\\RASTER\\", name, "_DSM_max_fill.tif"),
                   overwrite = TRUE)


DSM_for_shadow = rast(paste0(dir, "\\output\\RASTER\\", name, "_DSM_max_fill.tif"))

z_scale = terra::res(DSM_for_shadow)[1]

# read in the metadata from a directory of micasense images


pic_path = list.files(path = c(paste0(dir, "\\DAP\\MS\\Blue\\"), paste0(dir, "\\DAP\\MS\\Red\\")), recursive = TRUE, full.names = TRUE)


pic_exif = read_exif(pic_path)

pic_exif = pic_exif %>% 
  mutate(prop_scattered = ScatteredIrradiance / (ScatteredIrradiance + DirectIrradiance),
         prop_direct = DirectIrradiance / (ScatteredIrradiance + DirectIrradiance),
         prop_scattered_new = ScatteredIrradiance / DirectIrradiance
  )

ggplot(pic_exif, aes(x = TimeStamp, y = prop_scattered, color = BandName)) +
  geom_point(aes(y = prop_scattered_new), color = "steelblue", size = 3) +
  geom_point(aes(y = prop_scattered)) +
  #geom_point(aes(y = DirectIrradiance), color = "red3", size = 3) +
  #geom_line(aes(y = prop_direct), size = 2) +
  theme_bw(base_size = 18)

scattered = mean(pic_exif$prop_scattered, na.rm = TRUE)
direct = mean(pic_exif$prop_direct, na.rm = TRUE)

pic_date_first = min(pic_exif$ModifyDate, na.rm = TRUE) %>% 
  ymd_hms()
pic_date_last = max(pic_exif$ModifyDate, na.rm = TRUE) %>% 
  ymd_hms()

pic_lat = mean(pic_exif$GPSLatitude, na.rm = TRUE)
pic_lon = mean(pic_exif$GPSLongitude, na.rm = TRUE)

sun_angle = suncalc::getSunlightPosition(date = c(pic_date_first, pic_date_last),
                                         lat = pic_lat,
                                         lon = pic_lon) %>% 
  mutate(alt_deg = altitude * 180 / pi,
         az_deg = 180 + azimuth * 180 / pi)

# elevation matrix 
elmat = DSM_for_shadow %>% raster() %>% raster_to_matrix()


elmat2 = elmat %>% 
  ray_shade(sunangle = mean(sun_angle$az_deg),
            anglebreaks = c(sun_angle$alt_deg[1], sun_angle$alt_deg[2]),
            zscale = z_scale,
            multicore = FALSE)

DSM_shadow_direct = DSM_for_shadow %>% 
  rast() %>% 
  setValues(elmat2) %>% 
  flip(direction = "horizontal")

terra::writeRaster(DSM_shadow_direct, paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_direct.tif"),
                   overwrite = TRUE)

elmat3 = elmat %>% 
  ambient_shade(anglebreaks = 90 * cospi(seq(5, 85, by = 10)/180),
                maxsearch = 120,
                subreaks = 8,
                multicore = TRUE,
                zscale = z_scale,
                progbar = TRUE)

DSM_shadow_scattered = DSM_for_shadow %>% 
  rast() %>% 
  setValues(elmat3) %>% 
  flip(direction = "horizontal")

terra::writeRaster(DSM_shadow_scattered, paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_scattered.tif"),
                   overwrite = TRUE)


stack = c(DSM_shadow_direct, DSM_shadow_scattered)

DSM_shadow_weighted = terra::weighted.mean(stack,
                                           c(direct, scattered)) %>% 
  resample(ms_temp, "bilinear")


terra::writeRaster(DSM_shadow_weighted, paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_weighted.tif"),
                   overwrite = TRUE)

DSM_shadow_weighted = rast(paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_weighted.tif"))

DSM_shadow_canopy = DSM_shadow_weighted

threshold = 1/3

DSM_shadow_mask = DSM_shadow_weighted
DSM_shadow_mask[DSM_shadow_mask > threshold] = NA
DSM_shadow_mask[DSM_shadow_mask <= threshold] = 1

terra::writeRaster(DSM_shadow_mask, paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_mask.tif"),
                   overwrite = TRUE)


shadow_patches = raster::clump(raster::raster(DSM_shadow_mask), directions = 4) %>% 
  rast()

clumps = data.frame(freq(shadow_patches))
# threshold 200 cm^2
num_pix = 0.02 / (res(shadow_patches)[1]^2)
flecks = clumps[clumps$count > num_pix,] #remove clump observations with frequency smaller than 9
flecks = as.vector(flecks$value) # record IDs from clumps which met the criteria in previous step

new_mask = shadow_patches %in% flecks
new_mask[new_mask == 0] = NA

terra::writeRaster(new_mask, paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_mask2.tif"),
                   overwrite = TRUE)




#-------------------------------------------------------------------------------
# manual editing of crown polygons here
#-------------------------------------------------------------------------------

library(lwgeom)

CHM_max = terra::rast(paste0(dir, "\\output\\RASTER\\", name, "_CHM_max.tif"))



pols = st_read(paste0(dir, "\\output\\HULLS\\", name, "_z50_manual_2.shp")) %>% 
  st_snap_to_grid(res(CHM_max)) %>% 
  filter(!st_is_empty(.)) %>% 
  #st_make_valid() %>% 
  st_cast("POLYGON", do_split = FALSE)

# check to make sure that there are no multiple geometries
# this df should come out empty 
pols_dup = pols %>% 
  group_by(Seedlot, Tree, Rep) %>% 
  filter(n()>1)

pols_dup = pols %>% 
  group_by(Obs) %>% 
  filter(n()>1)

st_write(obj = pols, 
         dsn = paste0(dir, "\\output\\HULLS\\", name, "_z50_snap.shp"),
         layer = name,
         driver = "ESRI Shapefile",
         append = FALSE)

pols = st_read(paste0(dir, "\\output\\HULLS\\", name, "_z50_snap.shp")) %>% 
  filter(!st_is_empty(.)) %>% 
  st_make_valid() %>% 
  st_cast("POLYGON")

################################################################################

plan(multisession, workers = 8L)

# read in norm cloud
NORM = readLAScatalog(paste0(dir, "\\output\\NORM_clean\\"))
opt_chunk_buffer(NORM) = 5
#opt_filter(NORM) = "-thin_with_voxel 0.04"
opt_laz_compression(NORM) = TRUE
opt_output_files(NORM) = paste0(dir, "\\output\\SEGMENTED\\{*}_SEGMENTED")
opt_progress(NORM) = TRUE

#chunk = "D:\\Sync\\_Sites\\Sx_Genecology\\Skimikin\\output\\NORM_clean\\Skimikin-0-0_NORM.laz"
zq = 0.5

polys_to_las = function(chunk, zq = 0.5, polygons = pols) {
  
  las = readLAS(chunk)                  
  if (is.empty(las)) {
    return(NULL) }
  
  #polys_crop = raster::crop(polys, extent(las))
  las2 = merge_spatial(las, polygons, "Obs")
  
  las_df = las2@data %>%
    dplyr::group_by(Obs) %>%
    dplyr::mutate(Zq999 = quantile(Z, 0.999)) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(Obs = if_else(Z > Zq999 * zq, as.numeric(Obs), 0))
  
  las3 = las2 %>% 
    add_lasattribute(las_df$Obs, name = "treeID", desc = "Observation") %>% 
    filter_poi(treeID > 0) %>% 
    filter_poi(buffer == 0)
  
  if (is.empty(las3)) {
    return(NULL)
  } else {
    return(las3)
  }
  
}

# new chunks have only segmented portions of tree crowns
SEGMENTED = catalog_apply(NORM, polys_to_las)

SEGMENTED = readLAScatalog(paste0(dir, "\\output\\SEGMENTED\\"))
opt_chunk_buffer(SEGMENTED) = 0
opt_chunk_size(SEGMENTED) = 10000
opt_laz_compression(SEGMENTED) = TRUE
opt_progress(SEGMENTED) = TRUE
opt_output_files(SEGMENTED) = paste0(dir, "\\output\\MERGED\\", name, "_HULLS_merged")

# merge all the segmented trees into a single point cloud 
MERGED = catalog_retile(SEGMENTED)

# save each individual tree as a point cloud
opt_output_files(SEGMENTED) = paste0(dir, "\\output\\CROWNS\\", name, "_p{Seedlot}_r{Rep}_t{Tree}")
CROWNS = clip_roi(SEGMENTED, pols) 

# Now we have to make sure that each ITC cloud has only one treeID value;
# it seems that clip_roi does not assure this

# write a function to drop remnant treeID values
clean_crowns = function(chunk) {
  
  las = readLAS(chunk)                  
  if (is.empty(las)) return(NULL) 
  
  las = readLAS(chunk)                  
  if (is.empty(las)) return(NULL) 
  
  treeID_true = as.numeric(names(sort(table(las@data$treeID), decreasing = TRUE))[1])
  
  las2 = filter_poi(las, treeID == treeID_true)
}

CROWNS = readLAScatalog(paste0(dir, "\\output\\CROWNS\\"))
opt_chunk_size(CROWNS) = 0 # processing by files
opt_laz_compression(CROWNS) = TRUE
opt_chunk_buffer(CROWNS) = 0 # no buffer
opt_wall_to_wall(CROWNS) = TRUE # disable internal checks to ensure a valid output. Free wheel mode
opt_output_files(CROWNS) = paste0(dir, "\\output\\CROWNS\\{*}")

CROWNS_clean = catalog_apply(CROWNS, clean_crowns)

future:::ClusterRegistry("stop")



##---------------------------------------------------------------------------##

# mask the MS layer
# and MS polygons

ms = rast(paste0(dir, "\\input\\MS_ORTHO\\", name, "_MS.tif"))

pols_ms = st_read(paste0(dir, "\\output\\HULLS\\", name, "_z50_snap.shp")) %>% 
  filter(!st_is_empty(.)) %>% 
  st_buffer(dist = -.05) %>%
  vect()


# 
# st_write(obj = pols_ms, 
#          dsn = paste0(dir, "\\output\\HULLS\\", name, "_z50_ms.shp"),
#          layer = name,
#          driver = "ESRI Shapefile",
#          append = FALSE)


mask = rast(paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_mask2.tif"))

ms_mask = terra::mask(ms, pols_ms) %>% 
  terra::mask(mask, maskvalues = 1, updatevalue = NA)

terra::writeRaster(ms_mask, paste0(dir, "\\output\\RASTER\\", name, "_MS_mask.tif"),
                   overwrite = TRUE)

##---------------------------------------------------------------------------##
