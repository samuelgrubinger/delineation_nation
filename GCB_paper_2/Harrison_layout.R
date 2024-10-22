
# This script takes (1) a csv representing census data from 
# common garden experimental trials, and (2) shapefiles 
# respresenting the location of (1st row, 1st column)
# and (1st row, last column) and produces a shapefile with a point
# for each tree, at the correct spacing and generally correct location.
# This shapefile can then be manually adjusted in a GIS. 
# User must create buffer trees with fields == 0

library(tidyverse)
library(sf)
library(spatial)
library(rgdal)
library(raster)
library(lidR)
library(sp)
library(nngeo)
library(future)

#opts_knit$set(root.dir = normalizePath("D:/Sync/_PhD/Analyses/Skimikin_TEST/"))
dir = "D:\\Sync\\_Sites\\Sx_Genecology\\Harrison"

name = "Harrison"

csv_path = paste0(dir, "\\CSV\\", name, ".csv")


# set up here for units in meters, UTM

# give some inputs for rough geolocation
# a shapefile point of the location of the first tree (1st col, 1st row)
# and last tree in the same column (1st col, last row)
# This script assumes ROWS to be numbered east-west
# COLUMNS are numbered north-south

first_first = st_read(paste0(dir, "\\GIS\\", name, "_first_first.shp")) %>% 
  st_geometry()

first_last = st_read(paste0(dir, "\\GIS\\", name, "_first_last.shp")) %>% 
  st_geometry()

# column spacing in m 
c_space = 2

# row spacing in m 
r_space = 1

# read in the CSV
# eliminate unplanted, dead, and thinned trees
# notation may vary by trial 
# 8888 - thinned
# 9999 - unplanted
# 0 - dead
# . - dead

csv = read.csv(csv_path) %>% 
  dplyr::select(Obs, Site, Col, Row, Rep, Blk, Seedlot, 
         Tree, Ht10, Ht16, DBH16, Weevil16, Absent) 

# produce a grid using row and column numbers, to be rotated
grid = csv %>% 
  mutate(
    x_coord = st_coordinates(first_first)[1] + (Col * c_space) - c_space,
    y_coord = st_coordinates(first_first)[2] + (Row * r_space) - c_space)

# to sf object
grid_sf = st_as_sf(grid, coords = c("x_coord", "y_coord"), crs = 26910)


# get the coordinates of the tree in the 1st col, last row
first_last_sf = grid_sf %>% 
  filter(Row == 1 & Col == max(.$Col)) %>% 
  st_geometry()

# define a FUNCTION which rotates coords for a given angle
rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)

# get the angle of rotation between the true (first, last) and
# the current position of (first, last)

a = st_coordinates(first_last - first_first)
b = st_coordinates(first_last_sf - first_first)

theta = acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )

# rotate all points 
grid_rot = (st_geometry(grid_sf) - first_first) * rot((theta)) + first_first

# in order to merge the new coordinates and the census info,
# copy the sf object created from the master csv
grid_rot_sf = grid_sf
# then replace its geometry with the rotated coords
st_geometry(grid_rot_sf) = grid_rot
# and set its CRS to the correct one
st_crs(grid_rot_sf) = 26910

grid_rot_sf_clean = grid_rot_sf %>% 
  filter(Ht16 != 8888 &
           Ht16 != 88888 &
           Ht16 != 9999 &
           Ht16 != 0 &
           Ht16 != ".") %>%
  mutate(Weevil16 = replace_na(Weevil16, 0),
         Absent = replace_na(Absent, 0))




# save the rotated grid as a shapefile
st_write(obj = grid_rot_sf_clean,
         dsn = paste0(dir, "\\GIS\\", name, "_grid.shp"),
         append = FALSE,
         crs = 26910)

# now edit manually and rename site_grid_manual.shp

# for the missing Column 26


grid_rot_sf_clean_26 = grid_rot_sf %>% 
  filter(Col == 26 &
           Ht16 != 8888 &
           Ht16 != 88888 &
           Ht16 != 9999 &
           Ht16 != 0) %>%
  mutate(Weevil16 = replace_na(Weevil16, 0),
         Absent = replace_na(Absent, 0))


# save the rotated grid as a shapefile
st_write(obj = grid_rot_sf_clean_26,
         dsn = paste0(dir, "\\GIS\\", name, "_26.shp"),
         append = FALSE,
         crs = 26910)


