
library(zoo)
library(cowplot)
library(ggpubr)
library(raster)
library(terra)
library(sf)
library(RStoolbox)
library(tidyterra)
library(tidyverse)
library(grid)


###########################################
rect = rectGrob(
  x = unit(.1, "in"),
  y = unit(1, "npc") - unit(.1, "in"),
  width = unit(.5, "in"),
  height = unit(.5, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

font_s = 22

lab_a = textGrob(
  label = "(a)",
  x = unit(.15, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))


lab_b = textGrob(
  label = "(b)",
  x = unit(.15, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_c = textGrob(
  label = "(c)",
  x = unit(.15, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_d = textGrob(
  label = "(d)",
  x = unit(.15, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_e = textGrob(
  label = "(e)",
  x = unit(.15, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_f = textGrob(
  label = "(f)",
  x = unit(.2, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))


lab_g = textGrob(
  label = "(g)",
  x = unit(.15, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_h = textGrob(
  label = "(h)",
  x = unit(.15, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_i = textGrob(
  label = "(i)",
  x = unit(.2, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_j = textGrob(
  label = "(j)",
  x = unit(.2, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_k = textGrob(
  label = "(k)",
  x = unit(.16, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_l = textGrob(
  label = "(l)",
  x = unit(.2, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_m = textGrob(
  label = "(m)",
  x = unit(.12, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_n = textGrob(
  label = "(n)",
  x = unit(.15, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_o = textGrob(
  label = "(o)",
  x = unit(.15, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_p = textGrob(
  label = "(p)",
  x = unit(.15, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

#######################################

site = "Skimikin"
name = "Skimikin"


timepoint_cols = c("late winter" = "darkslategray3",
                   "early summer" = "darkolivegreen4",
                   "mid summer" = "goldenrod",
                   "late summer" =  "goldenrod4")

pols_dat = read_rds("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\traits//pols_dat_skim.rds") %>% 
  st_as_sf() %>% 
  #filter(Edge != 1 & Dead != 1) %>% 
  mutate(Tree = as.numeric(Tree)) %>% 
  mutate(Plot = "Corrected")

pols_raw = st_read("D:\\Sync\\_Sites\\Skimikin_spectral\\GIS\\Skimikin_z50_spectral.shp") %>% 
  mutate(Plot = "2020")

pols_updated = st_read("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\HULLS\\Edits\\Skimikin_z50_updated.shp") %>% 
  mutate(Plot = "Updated")

pols_plot = bind_rows(dplyr::select(pols_raw, Plot, geometry),
                      dplyr::select(pols_updated, Plot, geometry)) %>% 
  arrange(desc(Plot)) %>% 
  mutate(plot = ordered(Plot, levels = c("2020", "Updated")))


tree_pol = pols_dat %>% 
  filter(Prov == 382,
         Rep == 6,
         Tree == 2)

ext = st_centroid(tree_pol) %>% 
  ext() + c(7.5, 1.5, 1.5, 7.5)

ext2 = st_centroid(tree_pol) %>% 
  ext() + c(1.5, 1.5, 1.15, 1.85)

trees_buf = pols_dat %>% 
  st_crop(ext)

trees_pols = pols_plot %>% 
  st_cast("LINESTRING") %>% 
  st_crop(ext)

trees_pols2 = pols_plot %>% 
  st_cast("LINESTRING") %>% 
  st_crop(ext2)

segs = st_read("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\HULLS\\SEGS_step6_c5.shp") %>% 
  st_crop(ext2) %>% 
  mutate(label = "Supercells")

segs_keep = st_read("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\HULLS\\SEGS_KEEP.shp") %>% 
  st_crop(ext2)

tree_chm = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\RASTER\\Crowns\\fill_2021_CHM_max.tif") %>% 
  crop(ext2)

tree_NDVI = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\RASTER\\NDVI_all.tif")[[6]] %>% 
  crop(ext)

tree_NDVI_mask = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\RASTER\\NDVI_mask.tif")[[6]] %>% 
  crop(ext)

tree_NIR_mask = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\RASTER\\NIR_min_mask_all.tif")[[6]] %>%
  crop(ext)

tree_ms = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\input\\MS_ORTHO_average\\Skimikin_2021_07_29_MS.tif") %>% 
  crop(ext)

tree_ms2 = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\input\\MS_ORTHO_average\\Skimikin_2021_07_29_MS.tif") %>% 
  crop(ext2)

tree_NIR = tree_ms[[10]]

tree_ms_mask = tree_ms %>% 
  crop(st_buffer(trees_buf, -.05)) %>% 
  extend(ext)

################################################################################

# CHM
p_1 = ggR(raster(crop(tree_chm, ext2)),
           geom_raster = TRUE) +
   scale_fill_continuous(limits = c(0, 5.45),
                         na.value = "black",
                         high = "white",
                         low = "black",                   
                         breaks = seq(1, 5, 2),
                         name = " CHM (m)",
                         guide = FALSE) +
   theme_void() +
   ggnewscale::new_scale_fill() +
   geom_sf(data = segs, aes(color = label), fill = NA, linewidth = .5) +
  scale_color_manual(values = c("Supercells" = "orange")) +
   scale_x_continuous(expand = c(0, 0)) + 
   scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.185, .98),
        legend.justification = c(0, 1),
        panel.background = element_rect(fill = "grey50"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key.width = unit(.3, "in"),
        legend.key.height = unit(.3, "in"),
        legend.text = element_text(family = "Cambria", size = 14),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(6, 6, 6, 6))

(p1 = ggdraw(p_1) +
  draw_grob(rect) +
  draw_grob(lab_a))

# CHM
p_2 = ggR(raster(tree_chm),
          geom_raster = TRUE) +
  scale_fill_continuous(limits = c(0, 5.45),
                        na.value = "black",
                        high = "white",
                        low = "black",                   
                        breaks = seq(1, 5, 2),
                        name = " CHM (m)",
                        guide = guide_colorbar(frame.colour = "black",
                                               ticks = FALSE,
                                               title.position = "top",
                                               title.hjust = 0.5)) +
  theme_void() +
  ggnewscale::new_scale_fill() +
    geom_sf(data = segs_keep, aes(color = plot), fill = NA, linewidth = .5, color = "orange") +
   geom_sf(data = filter(trees_pols2, plot == "2020"), aes(color = plot), fill = NA, linewidth = .6) +
   scale_color_manual(values = c("2020" = "red3", "Updated" = "orange"),
                      labels = c("2020" = "2020  ", "Updated" = "Updated"),
                      guide = FALSE) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.18, .98),
        legend.justification = c(0, 1),
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key.width = unit(.355, "in"),
        legend.key.height = unit(.1, "in"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.background = element_rect(color = "black", fill = "white", size = .25),
        legend.margin = margin(4, 13, 4, 13))

(p2 = ggdraw(p_2) +
    draw_grob(rect) +
    draw_grob(lab_b))

# RGB
p_3 = ggRGB(tree_ms2, r = 6, g = 4, b = 2,
             scale = 255, 
             stretch = "lin",
             quantiles = c(0, 1),
            geom_raster = TRUE) +
   geom_sf(data = trees_pols2, aes(color = plot), fill = NA, linewidth = .6) +
   scale_color_manual(values = c("2020" = "red3", "Updated" = "orange"),
                      labels = c("2020" = "2020  ", "Updated" = "Updated")) +
  theme_void() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.185, .98),
        legend.justification = c(0, 1),
        panel.background = element_rect(fill = "grey50"),
        #legend.position = "none",
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key = element_rect(fill = "white", color = "white"),
        legend.key.width = unit(.3, "in"),
        legend.key.height = unit(.3, "in"),
        legend.text = element_text(family = "Cambria", size = 14),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(6, 6, 6, 6))

(p3 = ggdraw(p_3) +
    draw_grob(rect) +
    draw_grob(lab_c))

# NDVI
p_4 = ggR(raster(tree_NDVI),
           geom_raster = TRUE) +
    scale_fill_continuous(limits = c(minmax(tree_NDVI)[[1]], minmax(tree_NDVI)[[2]]),
                          na.value = "grey50",
                          high = "white",
                          low = "black",                   
                          breaks = c(.3, .6, .9),
                          name = "NDVI",
                          guide = guide_colorbar(frame.colour = "black",
                                                 ticks = FALSE,
                                                 title.position = "top",
                                                 title.hjust = 0.5)) +
    theme_void() +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.18, .98),
        legend.justification = c(0, 1),
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key.width = unit(.355, "in"),
        legend.key.height = unit(.1, "in"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.background = element_rect(color = "black", fill = "white", size = .25),
        legend.margin = margin(4, 13, 4, 13))

(p4 = ggdraw(p_4) +
    draw_grob(rect) +
    draw_grob(lab_e))

# NDVI mask
p_5 = ggR(raster(tree_NDVI),
           geom_raster = TRUE) +
  scale_fill_continuous(limits = c(minmax(tree_NDVI)[[1]], minmax(tree_NDVI)[[2]]),
                        na.value = "#571010",
                        high = "white",
                        low = "black",                   
                        breaks = c(.4, .65, .9),
                        name = "NDVI",
                        guide = guide_colorbar(frame.colour = "black",
                                               ticks = FALSE,
                                               title.position = "top",
                                               title.hjust = 0.5)) +
    ggnewscale::new_scale_fill() +
    ggR(tree_NDVI_mask, ggLayer=T, geom_raster=T) + 
    scale_fill_gradient(high = "#571010", low = "#571010", na.value = NA, guide = "none") +
    theme_void() +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          legend.position = "none",#c(.24, .973),
          legend.justification = c(0, 1),
          legend.direction = "horizontal",
          legend.key.width = unit(.35, "in"),
          legend.key.height = unit(.09, "in"),
          legend.text = element_text(family = "Cambria", size = 15, margin = margin(-3, 3, 3, 3)),
          legend.title = element_text(family = "Cambria", size = 17),
          legend.background = element_rect(color = "black", fill = "white", size = .5),
          legend.margin = margin(3, 14, 3, 12))

(p5 = ggdraw(p_5) +
    draw_grob(rect) +
    draw_grob(lab_f))


# NIR
p_6 = ggR(raster(tree_NIR),
           geom_raster = TRUE) +
    scale_fill_continuous(limits = c(minmax(tree_NIR)[[1]], minmax(tree_NIR)[[2]]),
                          na.value = "grey40",
                          high = "white",
                          low = "black",                   
                          breaks = c(.1, .3, .5),
                          labels = c("10", "30", "50"),
                          name = "NIR reflectance (%)",
                          guide = guide_colorbar(frame.colour = "black",
                                                 ticks = FALSE,
                                                 title.position = "top",
                                                 title.hjust = 0.5)) +
    theme_void() +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.18, .98),
        legend.justification = c(0, 1),
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key.width = unit(.355, "in"),
        legend.key.height = unit(.1, "in"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.background = element_rect(color = "black", fill = "white", size = .25),
        legend.margin = margin(4, 13, 4, 13))

(p6 = ggdraw(p_6) +
    draw_grob(rect) +
    draw_grob(lab_g))

# NIR mask
p_7 = ggR(raster(tree_NIR),
           geom_raster = TRUE) +
  scale_fill_continuous(limits = c(minmax(tree_NIR)[[1]], minmax(tree_NIR)[[2]]),
                        na.value = "#571010",
                        high = "white",
                        low = "black",                   
                        breaks = c(.10, .30, .50),
                        name = "NIR reflectance",
                        guide = guide_colorbar(frame.colour = "black",
                                               ticks = FALSE,
                                               title.position = "top",
                                               title.hjust = 0.5)) +
    ggnewscale::new_scale_fill() +
    ggR(tree_NIR_mask, ggLayer=T, geom_raster=T) + 
    scale_fill_gradient(high = "#571010", low = "#571010", na.value = NA, guide = "none") +
    theme_void() +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          legend.position = "none", #c(.24, .973),
          legend.justification = c(0, 1),
          legend.direction = "horizontal",
          legend.key.width = unit(.35, "in"),
          legend.key.height = unit(.09, "in"),
          legend.text = element_text(family = "Cambria", size = 15, margin = margin(-3, 3, 3, 3)),
          legend.title = element_text(family = "Cambria", size = 17),
          legend.background = element_rect(color = "black", fill = "white", size = .5),
          legend.margin = margin(3, 14, 3, 12))

(p7 = ggdraw(p_7) +
    draw_grob(rect) +
    draw_grob(lab_h))

# RGB
p_8 = ggRGB(brick(tree_ms), r = 6, g = 4, b = 2,
             scale = 255, 
             stretch = "lin",
             quantiles = c(0, .999),
            geom_raster = TRUE) +
    ggnewscale::new_scale_fill() +
    # ggR(tree_NIR_mask, ggLayer=T, geom_raster=T) + 
    # scale_fill_gradient(high = "#0E3057", low = "#1E3057", na.value = NA, guide = "none") +
    # ggnewscale::new_scale_fill() +
    # ggR(tree_NDVI_mask, ggLayer=T, geom_raster=T) + 
    # scale_fill_gradient(high = "#1E3057", low = "#1E3057", na.value = NA, guide = "none") +
    geom_sf(data = filter(trees_pols, plot == "Updated"), aes(color = plot), fill = NA, size = .8) +
    scale_color_manual(values = c("2020" = "red3", "Updated" = "orange"),
                       labels = c("2020" = "2020  ", "Updated" = "Updated")) +
    theme_void() +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.background = element_rect(fill = "grey50"),
          legend.position = "none",
          #legend.position = c(.18, .97),
          legend.justification = c(0, 1),
          legend.title.align = 1,
          legend.direction = "horizontal",
          legend.key.width = unit(.3, "in"),
          legend.key.height = unit(.3, "in"),
          legend.text = element_text(family = "Cambria", size = 14),
          legend.title = element_blank(),
          legend.background = element_rect(color = "black", fill = "white", size = .5),
          legend.margin = margin(6, 6, 6, 6))

(p8 = ggdraw(p_8) +
    draw_grob(rect) +
    draw_grob(lab_d))

################################################################################

rect_date = rectGrob(
  y = unit(1, "npc") - unit(.1, "in"),
  x = unit(1, "npc") - unit(.1, "in"),
  width = unit(1.2, "in"),
  height = unit(.4, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

lab_2020_07_02 = textGrob(
  label = "2020-07-02",
  y = unit(1, "npc") - unit(.2, "in"),
  x = unit(1, "npc") - unit(.17, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = 15))

lab_2020_08_07 = textGrob(
  label = "2020-08-07",
  y = unit(1, "npc") - unit(.2, "in"),
  x = unit(1, "npc") - unit(.17, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = 15))

lab_2020_08_25 = textGrob(
  label = "2020-08-25",
  y = unit(1, "npc") - unit(.2, "in"),
  x = unit(1, "npc") - unit(.17, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = 15))

lab_2021_03_31 = textGrob(
  label = "2021-03-31",
  y = unit(1, "npc") - unit(.2, "in"),
  x = unit(1, "npc") - unit(.17, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = 15))

lab_2021_06_29 = textGrob(
  label = "2021-06-29",
  y = unit(1, "npc") - unit(.2, "in"),
  x = unit(1, "npc") - unit(.17, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = 15))

lab_2021_07_29 = textGrob(
  label = "2021-07-29",
  y = unit(1, "npc") - unit(.2, "in"),
  x = unit(1, "npc") - unit(.17, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = 15))

lab_2021_08_14 = textGrob(
  label = "2021-08-14",
  y = unit(1, "npc") - unit(.2, "in"),
  x = unit(1, "npc") - unit(.17, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = 15))

lab_2023_02_24 = textGrob(
  label = "2023-02-24",
  y = unit(1, "npc") - unit(.2, "in"),
  x = unit(1, "npc") - unit(.17, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = 15))

# Index NIRvCCI

#(ms[[3]] - ms[[5]]) / (ms[[3]] + ms[[5]]) * ms[[10]])

pri = function(img, x = 3, y = 5) {
  bx = img[[x]]
  by = img[[y]]
  #bz = img[[z]]
  vi = ((bx - by) / (bx + by))
  return(vi)
}

gcc = function(img, b1 = 1, b2 = 2, b3 = 3, b4 = 4, b5 = 5, b6 = 6) {
  b1 = img[[b1]]
  b2 = img[[b2]]
  b3 = img[[b3]]
  b4 = img[[b4]]
  b5 = img[[b5]]
  b6 = img[[b6]]
  vi = ((b3 + b4) / (b1 + b2 + b3 + b4 + b5 + b6))
  return(vi)
}

re_total = function(img, x = 9, y = 7) {
  bx = img[[x]]
  by = img[[y]]
  vi = (bx - by)
  return(vi)
}

ari = function(img, x = 4, y = 7) {
  bx = img[[x]]
  by = img[[y]]
  vi = (1 / bx) - (1 / by)
  return(vi)
}

ES_20 = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\RASTER\\Skimikin_2020_07_02\\Updated\\Skimikin_2020_07_02_MS_use.tif") %>% 
  crop(ext) %>% 
  #nirvcci(3, 4, 10)
  pri()
  
MS_20 = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\RASTER\\Skimikin_2020_08_07\\Updated\\Skimikin_2020_08_07_MS_use.tif") %>% 
  crop(ext)%>% 
  #nirvcci(3, 4, 10)
  pri()

LS_20 = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\RASTER\\Skimikin_2020_08_25\\Updated\\Skimikin_2020_08_25_MS_use.tif") %>% 
  crop(ext)%>% 
  #nirvcci(3, 4, 10)
  pri()

LW_21 = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\RASTER\\Skimikin_2021_03_31\\Updated\\Skimikin_2021_03_31_MS_use.tif") %>% 
  crop(ext)%>% 
  #nirvcci(3, 4, 10)
  pri()

ES_21 = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\RASTER\\Skimikin_2021_06_29\\Updated\\Skimikin_2021_06_29_MS_use.tif") %>% 
  crop(ext)%>% 
  #nirvcci(3, 4, 10)
  pri()

MS_21 = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\RASTER\\Skimikin_2021_07_29\\Updated\\Skimikin_2021_07_29_MS_use.tif") %>% 
  crop(ext)%>% 
  #nirvcci(3, 4, 10)
  pri()

LS_21 = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\RASTER\\Skimikin_2021_08_14\\Updated\\Skimikin_2021_08_14_MS_use.tif") %>% 
  crop(ext)%>% 
  #nirvcci(3, 4, 10)
  pri()

LW_23 = rast("D:\\Sync\\_Sites\\Skimikin_spectral\\output\\RASTER\\Skimikin_2023_02_24\\Updated\\Skimikin_2023_02_24_MS_use.tif") %>% 
  crop(ext)%>% 
  #nirvcci(3, 4, 10)
  pri()

VI_stack = c(ES_20, MS_20, LS_20, LW_21,
             ES_21, MS_21, LS_21, LW_23) 

quants = global(VI_stack, fun = quantile, probs = c(.03, .5, .97), na.rm = TRUE) %>% 
  summarise(X3 = min(X3.),
            x97 = max(X97.),
            Xmid = x97 - ((x97 - X3) * .5))

VI_clamp = clamp(VI_stack, upper = quants[[2]], lower = quants[[1]])
###

vi_lims = c(quants[[1]] - .01, quants[[2]] + .01)

breaks = c(round(quants[[1]], 2), round(quants[[3]], 2), round(quants[[2]], 2))

# CHM

p_9 = ggR(raster(VI_clamp[[1]]),
           geom_raster = TRUE) +
    theme_void() +
    scale_fill_gradient2(limits = vi_lims,
                         na.value = "white",
                         high = "#0080D4",
                         mid = "#EFC139",
                         low = "#AF0842",
                         midpoint = quants[[3]],
                         breaks = breaks,
                         #name = expression(" NDRE"["740"]),
                         name = "CCI",
                         guide = guide_colorbar(frame.colour = "black",
                                                ticks = FALSE,
                                                title.position = "top",
                                                title.hjust = 0.5)) +
    theme_void(base_size = 14) +
    # ggnewscale::new_scale_fill() +
    # geom_sf(data = filter(trees_pols, plot == "Updated"), aes(color = plot), fill = NA, size = .9) +
    # scale_color_manual(values = c("2020" = "red3", "Updated" = "orange"),
    #                    labels = c("2020" = "2020  ", "Updated" = "Updated")) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          legend.position = c(.185, .98),
          legend.justification = c(0, 1),
          legend.title.align = 1,
          legend.direction = "horizontal",
          legend.key.width = unit(.25, "in"),
          legend.key.height = unit(.1, "in"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.background = element_rect(color = "black", fill = "white", size = .25),
          legend.margin = margin(4, 13, 4, 13))

(p9 = ggdraw(p_9) +
    draw_grob(rect) +
    draw_grob(lab_i) +
    draw_grob(rect_date) +
    draw_grob(lab_2020_07_02))


p_10 = ggR(raster(VI_clamp[[2]]),
           geom_raster = TRUE) +
    theme_void() +
  scale_fill_gradient2(limits = vi_lims,
                       na.value = "white",
                       high = "#0080D4",
                       mid = "#EFC139",
                       low = "#AF0842",
                       midpoint = quants[[3]],
                       breaks = breaks,
                       #name = expression(" NDRE"["740"]),
                       name = " PRI ",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks = FALSE)) +
    theme_void(base_size = 14) +
    # ggnewscale::new_scale_fill() +
    # geom_sf(data = filter(trees_pols, plot == "Updated"), aes(color = plot), fill = NA, size = .9) +
    # scale_color_manual(values = c("2020" = "red3", "Updated" = "orange"),
    #                    labels = c("2020" = "2020  ", "Updated" = "Updated")) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          legend.position = "none")

(p10 = ggdraw(p_10) +
    draw_grob(rect) +
    draw_grob(lab_j) +
    draw_grob(rect_date) +
    draw_grob(lab_2020_08_07))


p_11 = ggR(raster(VI_clamp[[3]]),
            geom_raster = TRUE) +
    theme_void() +
  scale_fill_gradient2(limits = vi_lims,
                       na.value = "white",
                       high = "#0080D4",
                       mid = "#EFC139",
                       low = "#AF0842",
                       midpoint = quants[[3]],
                       breaks = breaks,
                       #name = expression(" NDRE"["740"]),
                       name = " PRI ",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks = FALSE)) +
    theme_void(base_size = 14) +
    # ggnewscale::new_scale_fill() +
    # geom_sf(data = filter(trees_pols, plot == "Updated"), aes(color = plot), fill = NA, size = .9) +
    # scale_color_manual(values = c("2020" = "red3", "Updated" = "orange"),
    #                    labels = c("2020" = "2020  ", "Updated" = "Updated")) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          legend.position = "none")

(p11 = ggdraw(p_11) +
    draw_grob(rect) +
    draw_grob(lab_k) +
    draw_grob(rect_date) +
    draw_grob(lab_2020_08_25))



p_12 = ggR(raster(VI_clamp[[4]]),
            geom_raster = TRUE) +
    theme_void() +
  scale_fill_gradient2(limits = vi_lims,
                       na.value = "white",
                       high = "#0080D4",
                       mid = "#EFC139",
                       low = "#AF0842",
                       midpoint = quants[[3]],
                       breaks = breaks,
                       #name = expression(" NDRE"["740"]),
                       name = " PRI ",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks = FALSE)) +
    theme_void(base_size = 14) +
    # ggnewscale::new_scale_fill() +
    # geom_sf(data = filter(trees_pols, plot == "Updated"), aes(color = plot), fill = NA, size = .9) +
    # scale_color_manual(values = c("2020" = "red3", "Updated" = "orange"),
    #                    labels = c("2020" = "2020  ", "Updated" = "Updated")) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          legend.position = "none")

(p12 = ggdraw(p_12) +
    draw_grob(rect) +
    draw_grob(lab_l) +
    draw_grob(rect_date) +
    draw_grob(lab_2021_03_31))



p_13 = ggR(raster(VI_clamp[[5]]),
            geom_raster = TRUE) +
    theme_void() +
  scale_fill_gradient2(limits = vi_lims,
                       na.value = "white",
                       high = "#0080D4",
                       mid = "#EFC139",
                       low = "#AF0842",
                       midpoint = quants[[3]],
                       breaks = breaks,
                       #name = expression(" NDRE"["740"]),
                       name = " PRI ",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks = FALSE)) +
    theme_void(base_size = 14) +
    # ggnewscale::new_scale_fill() +
    # geom_sf(data = filter(trees_pols, plot == "Updated"), aes(color = plot), fill = NA, size = .9) +
    # scale_color_manual(values = c("2020" = "red3", "Updated" = "orange"),
    #                    labels = c("2020" = "2020  ", "Updated" = "Updated")) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          legend.position = "none")

(p13 = ggdraw(p_13) +
    draw_grob(rect) +
    draw_grob(lab_m) +
    draw_grob(rect_date) +
    draw_grob(lab_2021_06_29))


p_14 = ggR(raster(VI_clamp[[6]]),
            geom_raster = TRUE) +
    theme_void() +
  scale_fill_gradient2(limits = vi_lims,
                       na.value = "white",
                       high = "#0080D4",
                       mid = "#EFC139",
                       low = "#AF0842",
                       midpoint = quants[[3]],
                       breaks = breaks,
                       #name = expression(" NDRE"["740"]),
                       name = " PRI ",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks = FALSE)) +
    theme_void(base_size = 14) +
    # ggnewscale::new_scale_fill() +
    # geom_sf(data = filter(trees_pols, plot == "Updated"), aes(color = plot), fill = NA, size = .9) +
    # scale_color_manual(values = c("2020" = "red3", "Updated" = "orange"),
    #                    labels = c("2020" = "2020  ", "Updated" = "Updated")) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          legend.position = "none")

(p14 = ggdraw(p_14) +
    draw_grob(rect) +
    draw_grob(lab_n) +
    draw_grob(rect_date) +
    draw_grob(lab_2021_07_29))


p_15 = ggR(raster(VI_clamp[[7]]),
            geom_raster = TRUE) +
    theme_void() +
  scale_fill_gradient2(limits = vi_lims,
                       na.value = "white",
                       high = "#0080D4",
                       mid = "#EFC139",
                       low = "#AF0842",
                       midpoint = quants[[3]],
                       breaks = breaks,
                       #name = expression(" NDRE"["740"]),
                       name = " PRI ",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks = FALSE)) +
    theme_void(base_size = 14) +
    # ggnewscale::new_scale_fill() +
    # geom_sf(data = filter(trees_pols, plot == "Updated"), aes(color = plot), fill = NA, size = .9) +
    # scale_color_manual(values = c("2020" = "red3", "Updated" = "orange"),
    #                    labels = c("2020" = "2020  ", "Updated" = "Updated")) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          legend.position = "none")

(p15 = ggdraw(p_15) +
    draw_grob(rect) +
    draw_grob(lab_o) +
    draw_grob(rect_date) +
    draw_grob(lab_2021_08_14))



p_16 = ggR(raster(VI_clamp[[8]]),
            geom_raster = TRUE) +
    theme_void() +
  scale_fill_gradient2(limits = vi_lims,
                       na.value = "white",
                       high = "#0080D4",
                       mid = "#EFC139",
                       low = "#AF0842",
                       midpoint = quants[[3]],
                       breaks = breaks,
                       #name = expression(" NDRE"["740"]),
                       name = " PRI ",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks = FALSE)) +
    theme_void(base_size = 14) +
    # ggnewscale::new_scale_fill() +
    # geom_sf(data = filter(trees_pols, plot == "Updated"), aes(color = plot), fill = NA, size = .9) +
    # scale_color_manual(values = c("2020" = "red3", "Updated" = "orange"),
    #                    labels = c("2020" = "2020  ", "Updated" = "Updated")) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          legend.position = "none")

(p16 = ggdraw(p_16) +
    draw_grob(rect) +
    draw_grob(lab_p) +
    draw_grob(rect_date) +
    draw_grob(lab_2023_02_24))


################################################################################
g2 = ggarrange(NULL, p4, NULL, p5, NULL, p6, NULL, p7, NULL,
               nrow = 1, ncol = 9,
               widths = c(.02, 1, .02, 1, .02, 1, .02, 1, .02))

g1 = ggarrange(NULL, p1, NULL, p2, NULL, p3, NULL, p8, NULL,
               nrow = 1, ncol = 9,
               widths = c(.02, 1, .02, 1, .02, 1, .02, 1, .02))

g3 = ggarrange(NULL, p9, NULL, p10, NULL, p11, NULL, p12, NULL,
               nrow = 1, ncol = 9,
               widths = c(.02, 1, .02, 1, .02, 1, .02, 1, .02))

g4 = ggarrange(NULL, p13, NULL, p14, NULL, p15, NULL, p16, NULL,
               nrow = 1, ncol = 9,
               widths = c(.02, 1, .02, 1, .02, 1, .02, 1, .02))

g3 = ggarrange(NULL, g1, NULL, g2, NULL, g3, NULL, g4, NULL, nrow = 9, ncol = 1, 
               heights = c(.008, 1, .008, 1, .008, 1, .008, 1, .008))

ggsave(plot = g3, 
       filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\figure_2.jpeg", 
       device = jpeg,
       width = 15,
       height = 15,
       units = "in",
       dpi = 300,
       bg = "white")

################################################################################

