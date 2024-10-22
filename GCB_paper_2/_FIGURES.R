

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(cowplot)
library(ggpubr)
library(lemon)
library(raster)
library(stars)
library(elevatr)
library(ggspatial)
library(ggrepel)
library(RStoolbox)
library(ggsn)
library(terra)
library(extrafont)

windowsFonts(Times = windowsFont("Times New Roman"))
windowsFonts(Cambria = windowsFont("Cambria"))


all_plot = read_rds("D:\\Sync\\_Sites\\Sx_Genecology\\all_plot_pc.rds")

sx_pops = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\output\\Sx_CC_Seedlots_Normal_1961_1990_Y_S.csv") %>% 
  dplyr::select(Seedlot, Longitude, Latitude, TD, EMT, MAP, AHM) %>%  
  mutate(logMAP = log(MAP)) %>% 
  dplyr::select(-MAP) %>% 
  left_join(read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\Sx_CC_Seedlots_info.csv") %>% 
              dplyr::select(population, region, class),
            by = c("Seedlot" = "population")) %>% 
  left_join(read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\SxGenecology_predicted_Euc.csv") %>% 
              dplyr::select(sxProv, species, propGla, propEng, propSit),
            by = c("Seedlot" = "sxProv")) %>% 
  dplyr::rename("Class" = "class") %>% 
  mutate(Class = if_else(region == "AB", "B", Class)) %>% 
  mutate(Longitude = as.numeric(Longitude) * -1) %>% 
  mutate(improved = if_else(Class == "B", "Wildstand", "Orchard")) %>% 
  mutate(species = if_else(Class != "B", "Orchard", species)) %>% 
  left_join(dplyr::select(all_plot, -Site, -geometry))


    
    
sx_pops$Class = factor(sx_pops$Class, levels = c("B", "A", "A+"))

pops2 = sx_pops %>% 
  mutate(Longitude = Longitude * -1) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

pops2$Class = ordered(pops2$Class, levels = c("B", "A", "A+"))



sites = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\output\\Sx_CC_Sites_ReferencePeriod_2005_2020_Y_S.csv") %>% 
  dplyr::select(Site, Longitude, Latitude, TD, EMT, MAP, AHM) %>%  
  mutate(logMAP = log(MAP)) %>% 
  dplyr::select(-MAP) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  mutate(hjust1 = case_when(Site == "JordanRiver" ~ -.1,
                            Site == "Harrison" ~ -.1,
                            Site == "Skimikin" ~ 1.1,
                            Site == "Kalamalka" ~ -.1,
                            Site == "TJC" ~ 1.05,
                            Site == "Whitecourt" ~ -.1)) %>% 
  mutate(hjust2 = case_when(Site == "JordanRiver" ~ 0,
                           Site == "Harrison" ~ 1.1,
                           Site == "Skimikin" ~ 1.1,
                           Site == "Kalamalka" ~ -.1,
                           Site == "TJC" ~ -.05,
                           Site == "Whitecourt" ~ 0)) %>% 
  mutate(vjust1 = case_when(Site == "JordanRiver" ~ -.1,
                            Site == "Harrison" ~ 1.1,
                            Site == "Skimikin" ~ -.1,
                            Site == "Kalamalka" ~ -.1,
                            Site == "TJC" ~ -.3,
                            Site == "Whitecourt" ~ 1.1)) %>% 
  mutate(vjust2 = case_when(Site == "JordanRiver" ~ -1.1,
                            Site == "Harrison" ~ 1,
                            Site == "Skimikin" ~ -.2,
                            Site == "Kalamalka" ~ 1,
                            Site == "TJC" ~ -.1,
                            Site == "Whitecourt" ~ -.6))


sites$Site = recode_factor(sites$Site, "JordanRiver" = "Jordan River",
                              "Harrison" = "Harrison",
                              "Skimikin" = "Skimikin",
                              "Kalamalka" = "Kalamalka",
                              "TJC" = "Tête Jaune Cache",
                              "Whitecourt" = "Whitecourt")

sites2 = sites %>% 
  st_join(dplyr::select(all_plot, -geometry))

#saveRDS(sites, "E:\\Retreat\\sites.rds")

#################################

states = st_read("D:\\Sync\\_Sites\\Sx_Genecology\\_rnaturalearth\\ne_50m_admin_1_states_provinces.shp") %>% 
  filter(admin %in% c("United States of America", "Canada"))
countries = ne_countries(scale = 50, returnclass = "sf") %>% 
  filter(brk_name %in% c("Canada", "United States"))
coastlines = ne_coastline(scale = 50, returnclass = "sf")
ocean = st_read("D:\\Sync\\_Sites\\Sx_Genecology\\_rnaturalearth\\ne_50m_ocean.shp")

elevation_data1 = get_elev_raster(locations = countries, z = 4) %>% 
  terra::crop((st_bbox(pops2) + c(xmin = -7, ymin = -2, xmax = -13, ymax = 1))) %>% 
  projectRaster(crs = 3005)

rast_df1 = as.data.frame(elevation_data1, xy = TRUE)
names(rast_df1)[3] <- "elev"

saveRDS(rast_df1, "D:\\Sync\\_Sites\\Sx_Genecology\\rast_df1")

elevation_data2 = get_elev_raster(locations = countries, z = 5, clip = "locations") %>% 
  crop((st_bbox(sites) + c(xmin = -2, ymin = -3, xmax = 2, ymax = 2))) %>% 
  projectRaster(crs = 3005)

rast_df2 = as.data.frame(elevation_data2, xy = TRUE)
names(rast_df2)[3] <- "elev"


class_cols =  c("Wildstand" = #"#0098C7"
                 # "#255698",
                "#1752A2",
               "Orchard" = #"#815090"
                 "#0A8C35"
                 )
#class_alpha = c("B" = "steelblue1", "A" = "royalblue1", "A+" = "royalblue4")

###########################################
rect = rectGrob(
  x = unit(.15, "in"),
  y = unit(1, "npc") - unit(.15, "in"),
  width = unit(.6, "in"),
  height = unit(.6, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

lab_a = textGrob(
  label = "a",
  x = unit(.32, "in"),
  y = unit(1, "npc") - unit(.24, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))


lab_b = textGrob(
  label = "b",
  x = unit(.32, "in"),
  y = unit(1, "npc") - unit(.25, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

lab_c = textGrob(
  label = "c",
  x = unit(.32, "in"),
  y = unit(1, "npc") - unit(.24, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

lab_d = textGrob(
  label = "d",
  x = unit(.32, "in"),
  y = unit(1, "npc") - unit(.25, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

lab_e = textGrob(
  label = "e",
  x = unit(.32, "in"),
  y = unit(1, "npc") - unit(.24, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

lab_f = textGrob(
  label = "f",
  x = unit(.33, "in"),
  y = unit(1, "npc") - unit(.25, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))


lab_g = textGrob(
  label = "g",
  x = unit(.32, "in"),
  y = unit(1, "npc") - unit(.19, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

lab_h = textGrob(
  label = "h",
  x = unit(.32, "in"),
  y = unit(1, "npc") - unit(.25, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

lab_i = textGrob(
  label = "i",
  x = unit(.35, "in"),
  y = unit(1, "npc") - unit(.24, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

lab_j = textGrob(
  label = "j",
  x = unit(.35, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

lab_k = textGrob(
  label = "k",
  x = unit(.32, "in"),
  y = unit(1, "npc") - unit(.25, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

lab_l = textGrob(
  label = "l",
  x = unit(.36, "in"),
  y = unit(1, "npc") - unit(.25, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

lab_m = textGrob(
  label = "m",
  x = unit(.24, "in"),
  y = unit(1, "npc") - unit(.24, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

lab_n = textGrob(
  label = "n",
  x = unit(.32, "in"),
  y = unit(1, "npc") - unit(.24, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

#######################################

# sx range
p_1 = ggplot() +
  geom_raster(data = rast_df1,  aes(x = x, y = y, fill = elev)) +
  scale_fill_gradient(high = "white", low = "grey20") +
  geom_sf(data = countries, fill = NA, color = "grey20") +
  geom_sf(data = states, fill = NA, color = "grey20") +
  geom_sf(data = ocean, fill = "white", color = "grey20") + 
 # geom_sf(data = coastlines, color = "grey20", size = 1) +
  geom_sf(data = pops2, aes(color = improved), shape = 16, size = 2.2, alpha = .8) +
  geom_sf(data = sites, size = 3.1, shape = 21, fill = "red2") + 
  annotate("rect",
           xmin = 1010000,
           ymin = 300000,
           xmax = 1790000,
           ymax = 1140000,
           fill = NA,
           colour = "black",
           size = 0.6) +
  annotate("text", x = 1900000, y = 1490000, label = "C    A    N    A    D    A",
           size = 10, color = "grey25", family = "Cambria", angle = 13) +
  annotate("text", x = 2170000, y = -284500, label = "U  N  I  T  E  D    S  T  A  T  E  S\no f\nA  M  E  R  I  C  A",
           size = 7.5, color = "grey25", family = "Cambria", angle = 11) +
  annotate("text", x = 650000, y = -13000, label = "P   A   C   I   F   I   C\n\nO   C   E   A   N",
           size = 5.5, color = "grey25", family = "Cambria", angle = 0) +
  coord_sf(xlim = c(225000, 2800000), ylim = c(-768333, 2150000),
           crs = st_crs(3005)) +
  scale_color_manual(values = class_cols) +
  theme_nothing(14) +
  theme(panel.border = element_rect(fill = NA, color = "grey25", size = 1),
        legend.position = "none",
        plot.background = element_rect(fill = NA, color = "black", size = 1),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#p_1 

p1 = ggdraw(p_1) +
  draw_grob(rect) +
  draw_grob(lab_b)
  

e = ext(1723000, 1760000,
        340000, 368000)

el_crop = crop(rast(elevation_data2), e)

rast_df_crop = as.data.frame(el_crop, xy = TRUE)
names(rast_df_crop)[3] <- "elev"

# sx test
p_2 = ggplot() +
        geom_raster(data = rast_df2,  aes(x = x, y = y, fill = elev)) +
        scale_fill_gradient(high = "white", low = "grey55") +
        geom_sf(data = countries, fill = NA, color = "grey20") +
        geom_sf(data = states, fill = NA, color = "grey20") +
        geom_sf(data = ocean, fill = "white", color = "grey20") + 
        geom_sf(data = sites, size = 4.8, shape = 21, fill = "red2") + 
        geom_sf_label(data = sites, aes(label = Site, vjust = -.4, hjust = .4, 
                                        angle = 0), size = 7, color = NA, fill = "white", alpha = .6, fontface = "bold",
                      family = "Cambria") +
        geom_sf_text(data = sites, aes(label = Site, vjust = -.9, hjust = .4, 
                                       angle = 0), size = 7, color = "black", fontface = "bold",
                     family = "Cambria") +
        coord_sf(xlim = c(1080000, 1740000), ylim = c(352000, 1100000),
                 crs = st_crs(3005)) +
        scalebar(sites,
                 st.bottom = FALSE,
                 anchor = c(x = 1720000, y = 325000),
                 st.size = 5,
                 st.dist = 3200,
                 dist = 100,
                 dist_unit = "km",
                 height = 2000,
                 family = "Cambria",
                 transform = FALSE) +
        geom_raster(data = rast_df_crop,  aes(x = x, y = y, fill = elev)) +
        annotate("text", x = 1742000, y = 355000, label = " km",
           size = 5, family = "Cambria") +
        annotate("text", x = 1600000, y = 1000000, label = "A  L  B  E  R  T  A",
                 size = 9, color = "grey25", family = "Cambria", angle = 0) +
        annotate("text", x = 1280000, y = 800000, label = "B  R  I  T  I  S  H\nC  O  L  U  M  B  I  A",
                 size = 9, color = "grey25", family = "Cambria", angle = 0) +
        theme_nothing(14) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1),
              legend.position = "none",
              plot.background = element_rect(fill = NA, color = "black", size = 1),
              axis.text = element_blank(),
              axis.ticks = element_blank())

p2 = ggdraw(p_2) +
  draw_grob(rect) +
  draw_grob(lab_a)
#p2

g1 = ggarrange(NULL, p2, NULL, p1, NULL, ncol = 1,  heights = c(.005, 1, 0.01, 1, .005), align = "h")
#g1

# #------------------------------------------------------------------------------#
# 
# var1 = "AHM"
# var2 = "TD"
# 
# p_3 = ggplot() +
#   geom_point(data = pops2, aes_string(x = var1, y = var2, color = "improved"), size = 4.5, alpha = .6) +
#   geom_point(data = sites, aes_string(x = var1, y = var2), size = 5, shape = 21, fill = "red2") + 
#   geom_text(data = sites, 
#                   aes_string(x = var1, y = var2, label = "Site", hjust = "hjust1", vjust = "vjust1"), 
#                   size = 6.5, color = "black", fontface = "bold",
#             family = "Cambria") +
#   theme_bw(20) +
#   labs(x = "Annual heat:moisture index (°C/m)", y = "Temperature differential (°C)") +
#   scale_color_manual(values = class_cols, name = "Seedlot Class") +
#   # annotate(geom = "text", label = "MOIST & \nMARITIME", size = 5,
#   #          x = 5,
#   #          y = 10,
#   #          fontface = "bold",
#   #          alpha = .5) +
#   # annotate(geom = "text", label = "ARID & \nCONTINENTAL", size = 5,
#   #          x = 44,
#   #          y = 41,
#   #          fontface = "bold",
#   #          alpha = .5) +
#   #geom_text(data = ab, geom = "text", aes(x = -5, y = 3200, label = letter), size = 9, fontface = "bold") +
#   #facet_rep_grid(. ~ tree, repeat.tick.labels = TRUE) +
#   theme(strip.text = element_blank(),
#         plot.background = element_rect(fill = NA, color = "black", size = 1),
#         legend.text = element_text(family = "Cambria"),
#         legend.title = element_text(family = "Cambria"),
#         axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0), family = "Cambria"),
#         axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0), family = "Cambria"),
#         legend.position = "none")
# 
# p3 = ggdraw(p_3) +
#   draw_grob(rect) +
#   draw_grob(lab_c)
# 
# pops2_new = mutate(pops2, MAP = exp(logMAP))
# sites_new = mutate(sites, MAP = exp(logMAP))
# 
# var1 = "EMT"
# var2 = "MAP"
# 
# p_4 = ggplot() +
#   geom_point(data = pops2_new, aes_string(x = var1, y = var2, color = "improved"), size = 4.5, alpha = .6) +
#   geom_point(data = sites_new, aes_string(x = var1, y = var2), size = 5, shape = 21, fill = "red2") + 
#   geom_text(data = sites_new, 
#                   aes_string(x = var1, y = var2, label = "Site", hjust = "hjust2", vjust = "vjust2"), 
#                   size = 6.5, color = "black", fontface = "bold",
#             family = "Cambria") +
#   theme_bw(20) +
#   labs(x = "Extreme minimum temperature (°C)", 
#        y = "Mean annual precipitation (mm)") +
#   scale_color_manual(values = class_cols, name = "Class") +
#   coord_trans(y = "log10") +
#   scale_y_continuous(breaks= c(250, 500, 1000, 2000)) +
#   expand_limits(y = c(250, 2400)) +
#   # annotate(geom = "text", label = "DRY, WITH \nHARSH WINTERS", size = 5,
#   #          x = -55,
#   #          y = 5.5,
#   #          fontface = "bold",
#   #          alpha = .5) +
#   # annotate(geom = "text", label = "WET, WITH \nMILD WINTERS", size = 5,
#   #          x = -17,
#   #          y = 8.5,
#   #          fontface = "bold",
#   #          alpha = .5) +
#   # annotate("segment", x = 21, xend = 21, y = 100, yend = 100, arrow=arrow(angle=30, type = "closed")) +
#   #geom_text(data = ab, geom = "text", aes(x = -5, y = 3200, label = letter), size = 9, fontface = "bold") +
#   #facet_rep_grid(. ~ tree, repeat.tick.labels = TRUE) +
#   theme(strip.text = element_blank(),
#         plot.background = element_rect(color = "black", fill = NA, size = 1),
#         legend.position = c(0.82, 0.12),
#         legend.background = element_rect(fill = "white", color = "black"),
#         legend.text = element_text(family = "Cambria"),
#         legend.title = element_text(family = "Cambria"),
#         axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0), family = "Cambria"),
#         axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0), family = "Cambria"),
#         panel.grid.minor.y = element_blank()) +
#   guides(color = guide_legend(title = "Seed class"))
# 
# p4 = ggdraw(p_4) +
#   draw_grob(rect) +
#   draw_grob(lab_d)
# 
# p4

#------------------------------------------------------------------------------#


p_3 = ggplot() +
  geom_point(data = pops2, aes(x = PC1, y = PC2, color = improved), size = 4.5, alpha = .6) +
  geom_point(data = sites2, aes(x = PC1, y = PC2), size = 5, shape = 21, fill = "red2") + 
  geom_text(data = sites2, 
            aes(x = PC1, y = PC2,, label = Site.x, hjust = hjust1, vjust = vjust1), 
            size = 6.5, color = "black", fontface = "bold",
            family = "Cambria") +
  theme_bw(20) +
    annotate(geom = "text", label = "C O L D →",
             x = max(all_plot$PC1),
             y = min(all_plot$PC2),
             fontface = "bold",
             family = "Cambria",
             size = 7,
             hjust = 1,
             vjust = 2.5,
             alpha = .5) +
  annotate(geom = "text", label = "← W A R M",
           x = min(all_plot$PC1),
           y = min(all_plot$PC2),
           fontface = "bold",
           family = "Cambria",
           size = 7,
           hjust = .1,
           vjust = 2.5,
           alpha = .5) +
  annotate(geom = "text", label = "M O I S T →",
           x = min(all_plot$PC1),
           y = max(all_plot$PC2),
           fontface = "bold",
           family = "Cambria",
           size = 7,
           hjust = 1,
           vjust = -1,
           alpha = .5,
           angle = 90) +
  annotate(geom = "text", label = "← A R I D",
           x = min(all_plot$PC1),
           y = min(all_plot$PC2),
           fontface = "bold",
           family = "Cambria",
           size = 7,
           hjust = .2,
           vjust = -1,
           alpha = .5,
           angle = 90) +
  labs(x = "PC1 (42.4%)", y = "PC2 (28.2%)") +
  scale_color_manual(values = class_cols, name = "Seedlot Class") +
  scale_y_continuous(expand = expansion(mult = c(.12, .07))) +
  scale_x_continuous(expand = expansion(mult = c(.12, .07))) +
  theme(strip.text = element_blank(),
      plot.background = element_rect(fill = NA, color = "black", size = 1),
      legend.text = element_text(family = "Cambria"),
      legend.title = element_text(family = "Cambria"),
      axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0), family = "Cambria"),
      axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0), family = "Cambria"),
      legend.position = "none")

p3 = ggdraw(p_3) +
  draw_grob(rect) +
  draw_grob(lab_c)


p_4 = ggplot() +
  geom_point(data = pops2, aes(x = PC2, y = PC3, color = improved), size = 4.5, alpha = .6) +
  geom_point(data = sites2, aes(x = PC2, y = PC3), size = 5, shape = 21, fill = "red2") + 
  geom_text(data = sites2, 
            aes(x = PC2, y = PC3,, label = Site.x, hjust = hjust2, vjust = vjust2
                ), 
            size = 6.5, color = "black", fontface = "bold",
            family = "Cambria") +
  theme_bw(20) +
  annotate(geom = "text", label = "M O I S T →",
           x = max(all_plot$PC2),
           y = min(all_plot$PC3),
           fontface = "bold",
           family = "Cambria",
           size = 7,
           hjust = 1,
           vjust = 2.5,
           alpha = .5) +
  annotate(geom = "text", label = "← A R I D",
           x = min(all_plot$PC2),
           y = min(all_plot$PC3),
           fontface = "bold",
           family = "Cambria",
           size = 7,
           hjust = .1,
           vjust = 2.5,
           alpha = .5) +
  annotate(geom = "text", label = "M O N T A N E →",
           x = min(all_plot$PC2),
           y = max(all_plot$PC3),
           fontface = "bold",
           family = "Cambria",
           size = 7,
           hjust = 1,
           vjust = -1,
           alpha = .5,
           angle = 90) +
  annotate(geom = "text", label = "← L O W L A N D",
           x = min(all_plot$PC2),
           y = min(all_plot$PC3),
           fontface = "bold",
           family = "Cambria",
           size = 7,
           hjust = .1,
           vjust = -1,
           alpha = .5,
           angle = 90) +
  labs(x = "PC2 (28.2%)", y = "PC3 (14.0%)") +
  scale_y_continuous(expand = expansion(mult = c(.12, .07))) +
  scale_x_continuous(expand = expansion(mult = c(.12, .07))) +
  scale_color_manual(values = class_cols, name = "Seedlot Class") +
  theme(strip.text = element_blank(),
        plot.background = element_rect(fill = NA, color = "black", size = 1),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(family = "Cambria"),
        legend.title = element_text(family = "Cambria"),
        axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0), family = "Cambria"),
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0), family = "Cambria"),
        legend.position = c(0.82, 0.88),) +
  guides(color = guide_legend(title = "Seed class"))

p4 = ggdraw(p_4) +
  draw_grob(rect) +
  draw_grob(lab_d)

p4


g2 = ggarrange(NULL, p3, NULL, p4, NULL, ncol = 1, heights = c(.005, 1, .01, 1, .005), align = "h")
g2


g3 = ggarrange(g1, NULL, g2, ncol = 3, heights = c(1, 1, 1), widths = c(1, .008, 1), align = "h")
#g3


ggsave(plot = g3, 
       filename = "D:\\Sync\\Figures\\_FIGURES_paper3\\figure_1.tiff", 
       device = tiff,
       width = 13.5,
       height = 15.3,
       units = "in",
       dpi = 600,
       bg = "white")


ggsave(plot = g3, 
       filename = "D:\\Sync\\Figures\\_FIGURES_paper3\\figure_1.jpeg", 
       device = "jpeg",
       width = 13.5,
       height = 15.3,
       units = "in",
       dpi = 300,
       bg = "white")


################################################################################
################################################################################
################################################################################
# TRANSFER FUNCTIONS

pols_dat = read_rds("D://Sync//_Sites//Sx_Genecology//all_metrics.rds") %>% 
  filter(Edge != 1 & Dead != 1) %>% 
  mutate(Tree = as.numeric(Tree))

by_pop = pols_dat %>%  
  as.data.frame() %>% 
  dplyr::select(-geometry) %>% 
  #filter(class == "A+") %>% 
  group_by(Prov, Site) %>% 
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))


ggplot(by_pop, aes(x = Euc, y = Zq99)) + 
  #geom_abline(slope = 1, intercept = 0, size = 1.5, alpha = .5) +
  geom_point(color = "steelblue4", alpha = .6, size = 3) + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "red", size = 2, se = FALSE) +
  theme_bw(base_size = 30) + 
  #coord_equal() +
  facet_wrap(. ~ Site)



################################################################################
################################################################################
################################################################################
# METRICS

library(rgl)
library(threed)
library(magick)
library(alphashape3d)
library(lidR)
library(RStoolbox)
library(aRchi)

site = "Whitecourt"
name = "Whitecourt"

dir = "D://Sync//_Sites//Sx_Genecology//Whitecourt"

ang_func = function(las) {
    
    #las_green = filter_poi(las, las@data$G > las@data$R)
    
    Z = las@data$Z
    
    chm = grid_metrics(las, func = ~max(Z), res = 0.05)
    chm_mean = grid_metrics(las, func = ~mean(Z), res = 0.05)
    chm_mean[chm_mean < (maxValue(chm_mean))] = NA 
    chm_mean_trim = raster::trim(chm_mean)
    
    apex = clip_roi(las, extent(chm_mean_trim))@data %>% 
      filter(Z == max(.$Z)) %>% 
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
      mutate(apex_Y = apex$Y,
             apex_Z = apex$Z) %>% 
      rowwise() %>% 
      mutate(angle = myangle3d(b1 = X, b2 = Y, b3 = Z))
}

key_metrics = function(las) {
  
  #las_green = filter_poi(las, las@data$G > las@data$R)
  
  Z = las@data$Z
  
  chm = grid_metrics(las, func = ~max(Z), res = 0.05)
  chm_mean = grid_metrics(las, func = ~mean(Z), res = 0.2)
  chm_mean[chm_mean < (maxValue(chm_mean))] = NA 
  chm_mean_trim = raster::trim(chm_mean)
  
  apex = clip_roi(las, extent(chm_mean_trim))@data %>% 
    filter(Z == max(.$Z))
  
  origin = c(apex$X, apex$Y, apex$Z)
  a_point = c(apex$X, apex$Y, 0)
  
  
  myangle3d = function(b1, b2, b3) {
    aRchi::angle3d(o = c(origin[1], origin[2], origin[3]),
                   a = c(a_point[1], a_point[2], a_point[3]),
                   b = c(b1, b2, b3))
  }
  
  ang = las@data %>%
    dplyr::select(X, Y, Z) %>%
    rowwise() %>% 
    mutate(angle = myangle3d(b1 = X, b2 = Y, b3 = Z)) %>% 
    ungroup() %>% 
    summarise_at(vars(angle), .funs = c(mean, sd), na.rm = TRUE)
  
  
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

################################################################################

seedlot_1 = 375
rep_1 = 1
tree_1 = 1

seedlot_2 = 335
rep_2 = 1
tree_2 = 1

seedlot_3 = 309
rep_3 = 1
tree_3 = 1

##############
lab_size = 5
lab_y_1 = 1
lab_y_2 = .7

##############

las1 = lidR::readLAS(paste0(dir, "\\output\\CROWNS\\", site, "_p", seedlot_1, "_r", rep_1, "_t", tree_1, ".laz"))
las2 = lidR::readLAS(paste0(dir, "\\output\\CROWNS\\", site, "_p", seedlot_2, "_r", rep_2, "_t", tree_2, ".laz"))
las3 = lidR::readLAS(paste0(dir, "\\output\\CROWNS\\", site, "_p", seedlot_3, "_r", rep_3, "_t", tree_3, ".laz"))

dat1 = las1@data %>% 
  dplyr::mutate(RGB = rgb(R, G, B, maxColorValue = 65535))
dat2 = las2@data %>% 
  dplyr::mutate(RGB = rgb(R, G, B, maxColorValue = 65535))
dat3 = las3@data %>% 
  dplyr::mutate(RGB = rgb(R, G, B, maxColorValue = 65535))

met1 = key_metrics(las1) %>% 
  mutate(Zq99_lab = paste0(format(round(Zq99, 2), nsmall = 2), " m"),
         Zq999_lab = paste0(format(round(Zq999, 2), nsmall = 2), " m"),
         apex_angle_lab = paste0(format(round((apex_angle * 180 / pi), 2), nsmall = 2), "°"),
         apex_sd_lab = paste0(format(round((apex_sd * 180 / pi), 2), nsmall = 2), "°"),
         vol_convex_lab = format(round(vol_convex, 2), nsmall = 2))
met2 = key_metrics(las2) %>% 
  mutate(Zq99_lab = paste0(format(round(Zq99, 2), nsmall = 2), " m"),
         Zq999_lab = paste0(format(round(Zq999, 2), nsmall = 2), " m"),
         apex_angle_lab = paste0(format(round((apex_angle * 180 / pi), 2), nsmall = 2), "°"),
         apex_sd_lab = paste0(format(round((apex_sd * 180 / pi), 2), nsmall = 2), "°"),
         vol_convex_lab = format(round(vol_convex, 2), nsmall = 2))
met3 = key_metrics(las3) %>% 
  mutate(Zq99_lab = paste0(format(round(Zq99, 2), nsmall = 2), " m"),
         Zq999_lab = paste0(format(round(Zq999, 2), nsmall = 2), " m"),
         apex_angle_lab = paste0(format(round((apex_angle * 180 / pi), 2), nsmall = 2), "°"),
         apex_sd_lab = paste0(format(round((apex_sd * 180 / pi), 2), nsmall = 2), "°"),
         vol_convex_lab = format(round(vol_convex, 2), nsmall = 2))

p_1 = ggplot() +
    geom_point(data = dat1, size = .8, aes(x = Y, y = Z, color = RGB)) +
    geom_point(data = dat2, size = .8, aes(x = Y, y = Z, color = RGB)) +
    geom_point(data = dat3, size = .8, aes(x = Y, y = Z, color = RGB)) +
    scale_color_identity() +
  geom_linerange(color = "red3", size = .9,
                 aes(y = met1$Zq999),
                 xmin = min(dat1$Y, na.rm = TRUE),
                 xmax = max(dat1$Y, na.rm = TRUE)) +
  geom_linerange(color = "red3", size = .9,
                 aes(y = met2$Zq999),
                 xmin = min(dat2$Y, na.rm = TRUE),
                 xmax = max(dat2$Y, na.rm = TRUE)) +
  geom_linerange(color = "red3", size = .9,
                 aes(y = met3$Zq999),
                 xmin = min(dat3$Y, na.rm = TRUE),
                 xmax = max(dat3$Y, na.rm = TRUE)) +
    ##
    geom_linerange(color = "red3", size = .9,
                   aes(y = met1$Zq99),
                   xmin = min(dat1$Y, na.rm = TRUE),
                   xmax = max(dat1$Y, na.rm = TRUE)) +
    geom_linerange(color = "red3", size = .9,
                   aes(y = met2$Zq99),
                   xmin = min(dat2$Y, na.rm = TRUE),
                   xmax = max(dat2$Y, na.rm = TRUE)) +
    geom_linerange(color = "red3", size = .9,
                   aes(y = met3$Zq99),
                   xmin = min(dat3$Y, na.rm = TRUE),
                   xmax = max(dat3$Y, na.rm = TRUE)) +
    annotate(geom = "text",
             color = "black",
             label = deparse(bquote(Z[Q999] == .(met1$Zq999_lab))),
             parse = TRUE,
             size = lab_size,
             x = mean(dat1$Y, na.rm = TRUE),
             y = lab_y_1,
             hjust = .5,
             vjust = 1.2,
             fontface = "bold",
             family = "Cambria") +
    annotate(geom = "text",
             color = "black",
             label = deparse(bquote(Z[Q999] == .(met2$Zq999_lab))),
             parse = TRUE,
             size = lab_size,
             x = mean(dat2$Y, na.rm = TRUE) - .06,
             y = lab_y_1,
             hjust = .5,
             vjust = 1.2,
             fontface = "bold",
             family = "Cambria") +
    annotate(geom = "text",
             color = "black",
             label = deparse(bquote(Z[Q999] == .(met3$Zq999_lab))),
             parse = TRUE,
             size = lab_size,
             x = mean(dat3$Y, na.rm = TRUE) - .2,
             y = lab_y_1,
             hjust = .5,
             vjust = 1.2,
             fontface = "bold",
             family = "Cambria") +
  ###
  annotate(geom = "text",
           color = "black",
           label = deparse(bquote(Z[Q99] == .(met1$Zq99_lab))),
           parse = TRUE,
           size = lab_size,
           x = mean(dat1$Y, na.rm = TRUE),
           y = lab_y_2,
           hjust = .5,
           vjust = 1.2,
           fontface = "bold",
           family = "Cambria") +
  annotate(geom = "text",
           color = "black",
           label = deparse(bquote(Z[Q99] == .(met2$Zq99_lab))),
           parse = TRUE,
           size = lab_size,
           x = mean(dat2$Y, na.rm = TRUE) -.06,
           y = lab_y_2,
           hjust = .5,
           vjust = 1.2,
           fontface = "bold",
           family = "Cambria") +
  annotate(geom = "text",
           color = "black",
           label = deparse(bquote(Z[Q99] == .(met3$Zq99_lab))),
           parse = TRUE,
           size = lab_size,
           x = mean(dat3$Y, na.rm = TRUE) - .2,
           y = lab_y_2,
           hjust = .5,
           vjust = 1.2,
           fontface = "bold",
           family = "Cambria") +
  # LABEL
  # annotate(geom = "label",
  #          color = "black",
  #          label = "a",
  #          size = 18,
  #          x = min(met1),
  #          y = Inf,
  #          fontface = "bold") +
    #scale_alpha(range = c(0, .5)) +
    #geom_line() +
  expand_limits(y = 0) + 
  coord_equal() +
  #theme_bw(base_size = 18) +
  theme_void() +
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        legend.position = "none",
        aspect.ratio = 1)

p1 = ggdraw(p_1) +
  draw_grob(rect) +
  draw_grob(lab_a)

################################################################################


ang1 = ang_func(decimate_points(las1, random_per_voxel(.2, 2)))
ang2 = ang_func(decimate_points(las2, random_per_voxel(.2, 2)))
ang3 = ang_func(decimate_points(las3, random_per_voxel(.2, 2)))

p_2 = ggplot() +
    geom_point(data = dat1, size = 1, aes(x = Y, y = Z), color = "grey50", alpha = 0) +
    geom_point(data = dat2, size = 1, aes(x = Y, y = Z), color = "grey50", alpha = 0) +
    geom_point(data = dat3, size = 1, aes(x = Y, y = Z), color = "grey50", alpha = 0) +
    geom_segment(data = ang1, aes(x = apex_Y, y = apex_Z, xend = Y, yend = Z), color = "red4", alpha = .2, size = .6) +
    geom_segment(data = ang2, aes(x = apex_Y, y = apex_Z, xend = Y, yend = Z), color = "red4", alpha = .2, size = .6) +
    geom_segment(data = ang3, aes(x = apex_Y, y = apex_Z, xend = Y, yend = Z), color = "red4", alpha = .2, size = .6) +
   annotate(geom = "text",
            color = "black",
            label = deparse(bquote(~theta[apex] == .(met1$apex_angle_lab))),
            parse = TRUE,
            size = lab_size,
            x = mean(dat1$Y, na.rm = TRUE),
            y = lab_y_1,
            hjust = .5,
            vjust = 1.2,
            fontface = "bold",
            family = "Cambria") +
   annotate(geom = "text",
            color = "black",
            label = deparse(bquote(~theta[apex] == .(met2$apex_angle_lab))),
            parse = TRUE,
            size = lab_size,
            x = mean(dat2$Y, na.rm = TRUE) - .06,
            y = lab_y_1,
            hjust = .5,
            vjust = 1.2,
            fontface = "bold",
            family = "Cambria") +
   annotate(geom = "text",
            color = "black",
            label = deparse(bquote(~theta[apex] == .(met3$apex_angle_lab))),
            parse = TRUE,
            size = lab_size,
            x = mean(dat3$Y, na.rm = TRUE) - .2,
            y = lab_y_1,
            hjust = .5,
            vjust = 1.2,
            fontface = "bold",
            family = "Cambria") +
  ###
  annotate(geom = "text",
           color = "black",
           label = deparse(bquote(~sigma[apex] == .(met1$apex_sd_lab))),
           parse = TRUE,
           size = lab_size,
           x = mean(dat1$Y, na.rm = TRUE),
           y = lab_y_2,
           hjust = .5,
           vjust = 1.2,
           fontface = "bold",
           family = "Cambria") +
  annotate(geom = "text",
           color = "black",
           label = deparse(bquote(~sigma[apex] == .(met2$apex_sd_lab))),
           parse = TRUE,
           size = lab_size,
           x = mean(dat2$Y, na.rm = TRUE) - .06,
           y = lab_y_2,
           hjust = .5,
           vjust = 1.2,
           fontface = "bold",
           family = "Cambria") +
  annotate(geom = "text",
           color = "black",
           label = deparse(bquote(~sigma[apex] == .(met3$apex_sd_lab))),
           parse = TRUE,
           size = lab_size,
           x = mean(dat3$Y, na.rm = TRUE) - .2,
           y = lab_y_2,
           hjust = .5,
           vjust = 1.2,
           fontface = "bold",
           family = "Cambria") +
  expand_limits(y = 0) + 
  coord_equal() +
  #theme_bw(base_size = 18) +
  theme_void() +
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        legend.position = "none",
        aspect.ratio = 1)


p2 = ggdraw(p_2) +
  draw_grob(rect) +
  draw_grob(lab_b)

#p2
################################################################################

# alphashadep3d
a3d_1 = cbind(las1@data$X, las1@data$Y, las1@data$Z)
meany_1 = mean(a3d_1[,2])
a3d_1[,1] = a3d_1[,1] - mean(a3d_1[,1]) #center points around 0,0,0
a3d_1[,2] = a3d_1[,2] - mean(a3d_1[,2]) #center points around 0,0,0

a3d_2 = cbind(las2@data$X, las2@data$Y, las2@data$Z)
meany_2 = mean(a3d_2[,2])
a3d_2[,1] = a3d_2[,1] - mean(a3d_2[,1]) #center points around 0,0,0
a3d_2[,2] = a3d_2[,2] - mean(a3d_2[,2]) #center points around 0,0,0

a3d_3 = cbind(las3@data$X, las3@data$Y, las3@data$Z)
meany_3 = mean(a3d_3[,2])
a3d_3[,1] = a3d_3[,1] - mean(a3d_3[,1]) #center points around 0,0,0
a3d_3[,2] = a3d_3[,2] - mean(a3d_3[,2]) #center points around 0,0,0


shape1 = ashape3d(x = a3d_1, alpha = 1)
shape2 = ashape3d(x = a3d_2, alpha = 1)
shape3 = ashape3d(x = a3d_3, alpha = 1)
#plot(shape)

mesh1 = rgl::as.mesh3d(shape1) %>% 
  as_data_frame() %>% 
  mutate(alpha_ = percent_rank(desc(fcx)),
         y_new = y + meany_1)

mesh2 = rgl::as.mesh3d(shape2) %>% 
  as_data_frame() %>% 
  mutate(alpha_ = percent_rank(desc(fcx)),
         y_new = y + meany_2)

mesh3 = rgl::as.mesh3d(shape3) %>% 
  as_data_frame() %>% 
  mutate(alpha_ = percent_rank(desc(fcx)),
         y_new = y + meany_3)

p_3 = ggplot(mesh1, aes(x = y_new, y = z, group = zorder, alpha = alpha_)) +
  geom_polygon(data = mesh1, fill = "steelblue4") +
  geom_polygon(data = mesh2, fill = "steelblue4") +
  geom_polygon(data = mesh3, fill = "steelblue4") +
  scale_alpha(range = c(.2, .8)) +
  geom_line(data = mesh1) +
  geom_line(data = mesh2) +
  geom_line(data = mesh3) +
  annotate(geom = "text",
           color = "black",
           label = deparse(bquote(V["convex"] == .(met1$vol_convex_lab)~m^3)),
           parse = TRUE,
           size = lab_size,
           x = mean(dat1$Y, na.rm = TRUE),
           y = 1,
           hjust = .5,
           vjust = 1.2,
           fontface = "bold",
           family = "Cambria") +
  annotate(geom = "text",
           color = "black",
           label = deparse(bquote(V["convex"] == .(met2$vol_convex_lab)~m^3)),
           parse = TRUE,
           size = lab_size,
           x = mean(dat2$Y, na.rm = TRUE) - .06,
           y = 1,
           hjust = .5,
           vjust = 1.2,
           fontface = "bold",
           family = "Cambria") +
  annotate(geom = "text",
           color = "black",
           label = deparse(bquote(V["convex"] == .(met3$vol_convex_lab)~m^3)),
           parse = TRUE,
           size = lab_size,
           x = mean(dat3$Y, na.rm = TRUE) - .2,
           y = 1,
           hjust = .5,
           vjust = 1.2,
           fontface = "bold",
           family = "Cambria") +
  expand_limits(y = 0) + 
  coord_equal() +
  #theme_bw(base_size = 18) +
  theme_void() +
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        legend.position = "none",
        aspect.ratio = 1)

p3 = ggdraw(p_3) +
  draw_grob(rect) +
  draw_grob(lab_c)

################################################################################


pols_dat = read_rds("D://Sync//_Sites//Sx_Genecology//all_metrics.rds") %>% 
  filter(Edge != 1 & Dead != 1) %>% 
  mutate(Tree = as.numeric(Tree)) %>% 
  mutate(Plot = "Corrected")

pols_raw = st_read("D:\\Sync\\_Sites\\Sx_Genecology\\Whitecourt\\output\\HULLS\\Whitecourt_z50.shp") %>% 
  mutate(Plot = "Watershed")

pols_plot = bind_rows(dplyr::select(pols_dat, Plot, geometry),
                      dplyr::select(pols_raw, Plot, geometry)) %>% 
  arrange(desc(Plot)) %>% 
  mutate(plot = ordered(Plot, levels = c("Watershed", "Corrected")))

tree_pol = pols_dat %>% 
  filter(Site == "Whitecourt",
         Prov == 387,
         Rep == 4,
         Tree == 1)

ext = st_centroid(tree_pol) %>% 
  ext() + c(5.4, 5.5, 5, 5.9)

trees_buf = pols_dat %>% 
  st_crop(ext + 0.05)

trees_pols = pols_plot %>% 
  st_cast("LINESTRING") %>% 
  st_crop(ext)

tree_chm = rast(paste0(dir, "\\output\\RASTER\\", name, "_CHM_max.tif")) %>% 
  crop(ext)
quants = global(tree_chm, fun = quantile, probs = c(0, .97), na.rm = TRUE)
tree_chm_stretch = clamp(tree_chm, upper = quants[[2]], lower = quants[[1]])

tree_shadow = rast(paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_weighted.tif")) %>% 
  crop(ext)

tree_shadow_mask = rast(paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_mask2.tif")) %>% 
  crop(ext)

tree_rgb = rast(paste0(dir, "\\input\\ORTHO\\", name, ".tif")) %>% 
  crop(ext)


tree_ms = rast(paste0(dir, "\\input\\MS_ORTHO\\", name, "_MS.tif")) %>% 
  crop(ext)

tree_ms_mask = rast(paste0(dir, "\\output\\RASTER\\", name, "_MS_mask.tif")) %>% 
  crop(st_buffer(trees_buf, -.05)) %>% 
  extend(ext)

#cci = (tree_ms[[3]] - tree_ms[[5]]) / (tree_ms[[3]] + tree_ms[[5]])
cci_tree = (tree_ms_mask[[3]] - tree_ms_mask[[5]]) / (tree_ms_mask[[3]] + tree_ms_mask[[5]])
# perform a .01 - .99 stretch
quants = global(cci_tree, fun = quantile, probs = c(.01, .99), na.rm = TRUE)
cci_tree_stretch = clamp(cci_tree, upper = quants[[2]], lower = quants[[1]])
cci_midpoint = quants[[2]] - ((quants[[2]] - quants[[1]]) / 2)

#ndre = (tree_ms[[10]] - tree_ms[[8]]) / (tree_ms[[10]] + tree_ms[[8]])
ndre_tree = (tree_ms_mask[[10]] - tree_ms_mask[[8]]) / (tree_ms_mask[[10]] + tree_ms_mask[[8]])
# perform a .01 - .99 stretch
quants = global(ndre_tree, fun = quantile, probs = c(.01, .99), na.rm = TRUE)
ndre_tree_stretch = clamp(ndre_tree, upper = quants[[2]], lower = quants[[1]])
ndre_midpoint = quants[[2]] - ((quants[[2]] - quants[[1]]) / 2)

#ndre = (tree_ms[[10]] - tree_ms[[8]]) / (tree_ms[[10]] + tree_ms[[8]])
ndvi_tree = (tree_ms_mask[[10]] - tree_ms_mask[[6]]) / (tree_ms_mask[[10]] + tree_ms_mask[[6]])
# perform a .01 - .99 stretch
quants = global(ndvi_tree, fun = quantile, probs = c(.01, .99), na.rm = TRUE)
ndvi_tree_stretch = clamp(ndvi_tree, upper = quants[[2]], lower = quants[[1]])
ndvi_midpoint = quants[[2]] - ((quants[[2]] - quants[[1]]) / 2)

#re_upper = (tree_ms[[9]] - tree_ms[[8]]) / 23
re_upper_tree = (tree_ms_mask[[9]] - tree_ms_mask[[8]]) / 23
# perform a .01 - .99 stretch
quants = global(re_upper_tree, fun = quantile, probs = c(.01, .99), na.rm = TRUE)
re_upper_tree_stretch = clamp(re_upper_tree, upper = quants[[2]], lower = quants[[1]])
re_upper_midpoint = quants[[2]] - ((quants[[2]] - quants[[1]]) / 2)

#re_upper = (tree_ms[[9]] - tree_ms[[8]]) / 23
evi_tree = (2.5 * (tree_ms_mask[[10]] - tree_ms_mask[[6]])) / (tree_ms_mask[[10]] + (6 * tree_ms_mask[[6]]) - ( 6 * tree_ms_mask[[1]] + 1))
# perform a .01 - .99 stretch
quants = global(evi_tree, fun = quantile, probs = c(.01, .99), na.rm = TRUE)
evi_tree_stretch = clamp(evi_tree, upper = quants[[2]], lower = quants[[1]])
evi_midpoint = quants[[2]] - ((quants[[2]] - quants[[1]]) / 2)


p_4 = ggRGB(brick(tree_rgb), r = 1, g = 2, b = 3,
      scale = 255, 
      stretch = "lin",
      quantiles = c(0, 1)) +
  geom_sf(data = filter(trees_pols, Plot == "Corrected"), color = "red3", fill = NA, size = .6) +
  # scale_color_manual(values = c("Watershed" = "orange", "Corrected" = "red"),
  #                    labels = c("Watershed" = "Watershed  ", "Corrected" = "Corrected")) +
  theme_void() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.18, .97),
        legend.justification = c(0, 1),
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key.width = unit(.3, "in"),
        legend.key.height = unit(.3, "in"),
        legend.text = element_text(family = "Cambria", size = 14),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(6, 6, 6, 6))

p4 = ggdraw(p_4) +
  draw_grob(rect) +
  draw_grob(lab_d)

p4

p_5 = ggR(raster(ndre_tree_stretch),
          geom_raster = TRUE) +
  #geom_sf(data = trees_pols, color = "red", fill = NA, size = 1) +
  scale_fill_gradient2(na.value = "white",
                       high = "#0080D4",
                       mid = "#EFC139",
                       low = "#AF0842",
                       midpoint = ndre_midpoint + .01,
                       breaks = c(.3,
                                  .4,
                                  #0.2,
                                  .5),
                       name = expression(" NDRE"["705"]),
                       #name = " NDRE ",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks = FALSE)) +
  theme_void(base_size = 14) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.18, .97),
        legend.justification = c(0, 1),
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key.width = unit(.21, "in"),
        legend.key.height = unit(.1, "in"),
        legend.text = element_text(family = "Cambria"),
        legend.title = element_text(family = "Cambria"),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(6, 4, 6, 6))

p5 = ggdraw(p_5) +
  draw_grob(rect) +
  draw_grob(lab_e)

p5

c("Wildstand" = class_cols[1], "Orchard" = class_cols[2])

p_6 = ggR(raster(cci_tree_stretch),
         geom_raster = TRUE) +
  #geom_sf(data = trees_pols, color = "red", fill = NA, size = 1) +
  scale_fill_gradient2(na.value = "white",
                       high = "#0080D4",
                       mid = "#EFC139",
                       low = "#AF0842",
                       midpoint = cci_midpoint,
                       breaks = c(#0.0,
                                  0.1,
                                  0.2),
                       # name = expression("NIR"["V"^"CCI"]),
                       name = "CCI  ",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks = FALSE)) +
  theme_void(base_size = 14) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.18, .97),
        legend.justification = c(0, 1),
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key.width = unit(.24, "in"),
        legend.key.height = unit(.1, "in"),
        legend.text = element_text(family = "Cambria", size = 11),
        legend.title = element_text(family = "Cambria", size = 14),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(6, 6, 6, 6))

p6 = ggdraw(p_6) +
  draw_grob(rect) +
  draw_grob(lab_f)

p6

################################################################################
g1 = ggarrange(NULL, p1, NULL, p2, NULL, p3, NULL,
               nrow = 1, ncol = 7,
               widths = c(.004, 1, .015, 1, .015, 1, .004))

g2 = ggarrange(NULL, p4, NULL, p5, NULL, p6, NULL,
               nrow = 1, ncol = 7,
               widths = c(.004, 1, .015, 1, .015, 1, .004))

g3 = ggarrange(NULL, g1, NULL, g2, NULL, nrow = 5, ncol = 1, 
               heights = c(.001, 1, .01, 1, .001))

ggsave(plot = g3, 
       filename = "D:\\Sync\\Figures\\_FIGURES_GCB\\figure_4.tiff", 
       device = tiff,
       width = 15,
       height = 10,
       units = "in",
       dpi = 300,
       bg = "white")

ggsave(plot = g3, 
       filename = "D:\\Sync\\Figures\\_FIGURES_GCB\\figure_4.jpeg", 
       device = jpeg,
       width = 15,
       height = 10,
       units = "in",
       dpi = 300,
       bg = "white")


# gx = ggarrange(p4, NULL, p6, nrow = 1, ncol = 3, 
#                widths = c(1, .02, 1))
# 
# ggsave(plot = gx, 
#        filename = "D:\\Sync\\Figures\\3MT_2022.tiff", 
#        device = tiff,
#        width = 12,
#        height = 6,
#        units = "in",
#        dpi = 600,
#        bg = "white")


################################################################################
#
# CROWNS / SHADOWS
#

site = "Skimikin"
name = "Skimikin"

dir = "D://Sync//_Sites//Sx_Genecology//Skimikin"

pols_dat = read_rds("D://Sync//_Sites//Sx_Genecology//all_metrics.rds") %>% 
  filter(Edge != 1 & Dead != 1) %>% 
  mutate(Tree = as.numeric(Tree)) %>% 
  mutate(Plot = "Corrected")

pols_raw = st_read("D:\\Sync\\_Sites\\Sx_Genecology\\Skimikin\\output\\HULLS\\Skimikin_z50.shp") %>% 
  mutate(Plot = "Watershed")

pols_plot = bind_rows(dplyr::select(pols_dat, Plot, geometry),
                      dplyr::select(pols_raw, Plot, geometry)) %>% 
  arrange(desc(Plot)) %>% 
  mutate(plot = ordered(Plot, levels = c("Watershed", "Corrected")))

seeds = st_read("D:\\Sync\\_Sites\\Sx_Genecology\\Skimikin\\GIS\\Skimikin_grid_manual.shp") %>% 
  mutate(Plot = "Seed points")

ttops = st_read("D:\\Sync\\_Sites\\Sx_Genecology\\Skimikin\\output\\HULLS\\Skimikin_treetops.shp") %>% 
  mutate(Plot = "Maxima")



tree_pol = pols_dat %>% 
  filter(Site == site,
         Prov == 406,
         Rep == 1,
         Tree == 4)

ext = st_centroid(tree_pol) %>% 
  ext() + c(3.25, 3.35, 4.7, 5.85)

trees_buf = pols_dat %>% 
  st_crop(ext + 0.05)

ws_tops = bind_rows(dplyr::select(seeds, Plot, geometry),
                    dplyr::select(ttops, Plot, geometry)) %>% 
  arrange(desc(Plot)) %>% 
  mutate(plot = ordered(Plot, levels = c("Seed points", "Maxima"))) %>% 
  st_crop(ext + 0.05)

bufs = seeds %>% 
  st_crop(ext) %>% 
  st_buffer(dist = .4) %>% 
  st_crop(ext)

trees_pols = pols_plot %>% 
  st_cast("LINESTRING") %>% 
  st_crop(ext)

tree_chm = rast(paste0(dir, "\\output\\RASTER\\", name, "_CHM_max.tif")) %>% 
  crop(ext)

#plot(tree_chm)

# quants = global(tree_chm, fun = quantile, probs = c(.2, 99.5), na.rm = TRUE)
# tree_chm_stretch = clamp(tree_chm, upper = quants[[2]], lower = quants[[1]])

tree_chm_smooth = rast(paste0(dir, "\\output\\RASTER\\", name, "_CHM_max_smooth.tif")) %>% 
  crop(ext)

tree_chm_for_shadow = rast(paste0(dir, "\\output\\RASTER\\", name, "_DSM_max_fill.tif")) %>% 
  crop(ext)

tree_direct = rast(paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_direct.tif")) %>% 
  crop(ext)

tree_scattered = rast(paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_scattered.tif")) %>% 
  crop(ext)

tree_shadow = rast(paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_weighted.tif")) %>% 
  crop(ext)

tree_shadow_mask = rast(paste0(dir, "\\output\\RASTER\\", name, "_DSM_shadow_mask2.tif")) %>% 
  crop(ext)

tree_rgb = rast(paste0(dir, "\\input\\ORTHO\\", name, ".tif")) %>% 
  crop(ext)


tree_ms = rast(paste0(dir, "\\input\\MS_ORTHO\\", name, "_MS.tif")) %>% 
  crop(ext)

tree_ms_mask = rast(paste0(dir, "\\output\\RASTER\\", name, "_MS_mask.tif")) %>% 
  crop(st_buffer(trees_buf, -.05)) %>% 
  extend(ext)

###

p_1 = ggR(raster(tree_chm),
          geom_raster = TRUE) +
  geom_sf(data = filter(ws_tops, plot == "Seed points"), 
          shape = 21, fill = "orange", size = 3, stroke = 1) +
  scale_fill_continuous(limits = c(0, 7.25),
                        na.value = "black",
                        high = "white",
                        low = "black",                   
                        breaks = seq(1, 7, 2),
                        name = "CHM (m)",
                        guide = guide_colorbar(frame.colour = "black",
                                               ticks = FALSE,
                                               title.position="top",
                                               title.hjust = 0.5)) +
  theme_void() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.24, .973),
        legend.justification = c(0, 1),
        legend.title.align = 0,
        legend.direction = "horizontal",
        legend.key.width = unit(.3, "in"),
        legend.key.height = unit(.09, "in"),
        legend.text = element_text(family = "Cambria", size = 15, margin = margin(-3, 3, 3, 3)),
        legend.title = element_text(family = "Cambria", size = 17),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(5, 6, 3, 6)) 

p1 = ggdraw(p_1) +
  draw_grob(rect) +
  draw_grob(lab_a)

p1
###

###

p_2 = ggR(raster(tree_chm_smooth),
          geom_raster = TRUE) +
  scale_fill_continuous(limits = c(0, 7.25),
                        na.value = "black",
                        high = "white",
                        low = "black",                   
                        breaks = seq(1, 7, 2),
                        name = " CHM (m)",
                        guide = FALSE) +
  theme_void() +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values = c("Seed points" = "orange", "Maxima" = "red3")) +
    geom_sf(data = bufs, 
          shape = 21, fill = NA, color = "red3", size = .9) +
  geom_sf(data = ws_tops, 
          shape = 21, aes(fill = Plot), size = 3, stroke = 1) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.24, .973),
        legend.justification = c(0, 1),
        legend.title.align = 0,
        legend.direction = "vertical",
        legend.key.width = unit(.005, "in"),
        legend.key.height = unit(.28, "in"),
        legend.text = element_text(family = "Cambria", size = 17),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(1, 1, 5, 5)) 

p2 = ggdraw(p_2) +
  draw_grob(rect) +
  draw_grob(lab_b)

p2
###

###

p_3 = ggR(raster(tree_chm_smooth),
          geom_raster = TRUE) +
  scale_fill_continuous(limits = c(0, 7.25),
                        na.value = "black",
                        high = "white",
                        low = "black",                   
                        breaks = seq(1, 7, 2),
                        name = " CHM (m)",
                        guide = "none") +
  theme_void(base_size = 14) +
  ggnewscale::new_scale_fill() +
  scale_fill_manual(values = c("Seed points" = "orange", "Local maxima" = "red3")) +
  geom_sf(data = filter(ws_tops, Plot == "Local maxima"), 
          shape = 21, fill = "red3", size = 3, stroke = 1.2) +
  geom_sf(data = filter(trees_pols, Plot == "Watershed"), 
          color = "orange", size = .9) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

p3 = ggdraw(p_3) +
  draw_grob(rect) +
  draw_grob(lab_c)

p3
###

p_4 = ggRGB(brick(tree_rgb), r = 1, g = 2, b = 3,
            scale = 255, 
            stretch = "lin",
            quantiles = c(0, 1)) +
  geom_sf(data = trees_pols, aes(color = plot), fill = NA, size = .9) +
  scale_color_manual(values = c("Watershed" = "orange", "Corrected" = "red3"),
                     labels = c("Watershed" = "Watershed  ", "Corrected" = "Corrected")) +
  theme_void() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.24, .973),
        legend.justification = c(0, 1),
        legend.title.align = 0,
        legend.direction = "vertical",
        legend.key.width = unit(.28, "in"),
        legend.key.height = unit(.26, "in"),
        legend.text = element_text(family = "Cambria", size = 17),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(1, -1, 5, 5))  +
  guides(color = guide_legend(byrow = TRUE))

p4 = ggdraw(p_4) +
  draw_grob(rect) +
  draw_grob(lab_d)

p4

p_5 = ggR(raster(tree_chm),
          geom_raster = TRUE) +
  geom_sf(data = filter(trees_pols, plot == "Corrected"), 
          color = "red3", fill = NA, size = .9) +
  scale_fill_continuous(limits = c(0, 7.25),
                        na.value = "black",
                        high = "white",
                        low = "black",                   
                        breaks = seq(1, 7, 2),
                        name = " CHM (m)",
                        guide = "none") +
  theme_void(base_size = 14) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

p5 = ggdraw(p_5) +
  draw_grob(rect) +
  draw_grob(lab_e)

p5
###

# which tree?
site = "Skimikin"
prov = 383
rep = 8
tree = 2

tree_pol_for_laz = pols_dat %>% 
  filter(Site == site,
         Prov == prov,
         Rep == rep,
         Tree == tree)


# read in the points for the tree only 
NORM = readLAScatalog(paste0(dir, "\\output\\NORM_clean\\"))
#opt_filter(NORM) = paste("-keep_xyz", ext[1] - 45, ext[3] - 55, "0", ext[2] + 35, ext[4] + 25, "30")
tree_laz = clip_roi(NORM, tree_pol_for_laz)

tree_laz_clip = readLAS(paste0(dir, "\\output\\CROWNS\\", 
                               site, "_p", prov, "_r", rep, "_t", tree, ".laz"))

tree_met = key_metrics(tree_laz) %>% 
  mutate(Zq999_lab = paste0(format(round(Zq999, 2), nsmall = 2), " m"),
         p50 = Zq999 * .5,
         p50_lab = "50% threshold")

tree_laz_dat = tree_laz@data %>% 
  dplyr::mutate(RGB = rgb(R, G, B, maxColorValue = 90000))

tree_laz_clip_dat = tree_laz_clip@data %>% 
  dplyr::mutate(RGB = rgb(R, G, B, maxColorValue = 90000))


# get ratio for bounding box
x_dist = ext[2] - ext[1]
y_dist = ext[4] - ext[3]
x_y = y_dist / x_dist

tree_x_min = min(tree_laz_dat$X) - .7
tree_x_max = max(tree_laz_dat$X) + .7
tree_y_max = (tree_x_max - tree_x_min) * x_y


###

p_6 = ggplot() +
  geom_point(data = tree_laz_dat, size = 1.2, aes(x = X, y = Z, color = RGB)) +
  scale_color_identity() +
  scale_x_continuous(limits = c(tree_x_min, tree_x_max), expand = c(.03, 0)) +
  scale_y_continuous(limits = c(0, tree_y_max), expand = c(.03, 0)) +
  coord_equal() +
  theme_void() +
  geom_linerange(color = "red3", size = 1,
                 aes(y = tree_met$Zq999),
                 xmin = tree_x_min + .6,
                 xmax = tree_x_max - .75) +
  geom_linerange(color = "red3", size = 1,
                 aes(y = tree_met$p50),
                 xmin = tree_x_min + .6,
                 xmax = tree_x_max - .75) +
  geom_linerange(color = "grey35", size = 1.2, alpha = .7,
                 aes(y = 0),
                 xmin = tree_x_min + .6,
                 xmax = tree_x_max - .7) +
  annotate(geom = "text",
           color = "black",
           label = deparse(~bold(Z["Q999"])),
           parse = TRUE,
           x = tree_x_max - .7,
           y = tree_met$Zq999,
           hjust = 0,
           vjust = .5,
           fontface = "bold",
           size = 7,
           family = "Cambria") +
  annotate(geom = "text",
           color = "black",
           label = "50%",
           parse = FALSE,
           x = tree_x_max - .6,
           y = tree_met$p50,
           hjust = 0,
           vjust = .5,
           fontface = "bold",
           size = 7,
           family = "Cambria") +
  theme(panel.background = element_rect(color = "black", fill = "white"),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        legend.position = "none")

p6 = ggdraw(p_6) +
  draw_grob(rect) +
  draw_grob(lab_f)

p6

###

p_7 = ggplot() +
  geom_point(data = tree_laz_clip_dat, size = 1.2, aes(x = X, y = Z, color = RGB)) +
  scale_color_identity() +
  scale_x_continuous(limits = c(tree_x_min, tree_x_max), expand = c(.03, 0)) +
  scale_y_continuous(limits = c(0, tree_y_max), expand = c(.03, 0)) +
  coord_equal() +
  theme_void() +
  geom_linerange(color = "grey35", size = 1.2, alpha = .7,
                 aes(y = 0),
                 xmin = tree_x_min + .7,
                 xmax = tree_x_max - .7) +
  theme(panel.background = element_rect(color = "black", fill = "white"),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 1))

p7 = ggdraw(p_7) +
  draw_grob(rect) +
  draw_grob(lab_g)

p7

###

# p_8 = ggR(raster(tree_chm_for_shadow),
#           geom_raster = TRUE) +
#   scale_fill_continuous(limits = c(0, 7.25),
#                         na.value = "black",
#                         high = "white",
#                         low = "black",                   
#                         breaks = seq(1, 7, 2),
#                         name = " CHM (m)",
#                         guide = "none") +
#   theme_void(base_size = 14) +
#   scale_x_continuous(expand = c(0, 0)) + 
#   scale_y_continuous(expand = c(0, 0)) +
#   theme(panel.border = element_rect(color = "black", fill = NA),
#         legend.position = c(.1, .9885),
#         legend.justification = c(0, 1),
#         legend.title.align = 1,
#         legend.direction = "horizontal",
#         legend.key.width = unit(.395, "in"),
#         legend.key.height = unit(.088, "in"),
#         legend.text = element_text(family = "Cambria"),
#         legend.title = element_blank(),
#         legend.background = element_rect(color = "black", fill = "white", size = .5),
#         legend.margin = margin(6, 6, 6, 6))
# 
# p8 = ggdraw(p_8) +
#   draw_grob(rect) +
#   draw_grob(lab_h)
# 
# p8
###

p_8 = ggR(raster(tree_direct),
          geom_raster = TRUE) +
  geom_sf(data = filter(trees_pols, plot == "Corrected"), 
          color = "red3", fill = NA, size = .9) +
  scale_fill_continuous(na.value = "black",
                        high = "white",
                        low = "black",                   
                        breaks = seq(1, 7, 2),
                        name = " CHM (m)",
                        guide = "none") +
  theme_void(base_size = 14) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.18, .97),
        legend.justification = c(0, 1),
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key.width = unit(.395, "in"),
        legend.key.height = unit(.088, "in"),
        legend.text = element_text(family = "Cambria"),
        legend.title = element_text(family = "Cambria"),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(6, 6, 6, 6))

p8 = ggdraw(p_8) +
  draw_grob(rect) +
  draw_grob(lab_h)

p8
###

p_9 = ggR(raster(tree_scattered),
          geom_raster = TRUE) +
  geom_sf(data = filter(trees_pols, plot == "Corrected"), 
          color = "red3", fill = NA, size = .9) +
  scale_fill_continuous(na.value = "black",
                        high = "white",
                        low = "black",                   
                        breaks = seq(1, 7, 2),
                        name = " CHM (m)",
                        guide = "none") +
  theme_void(base_size = 14) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.18, .97),
        legend.justification = c(0, 1),
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key.width = unit(.395, "in"),
        legend.key.height = unit(.088, "in"),
        legend.text = element_text(family = "Cambria"),
        legend.title = element_text(family = "Cambria"),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(6, 6, 6, 6))

p9 = ggdraw(p_9) +
  draw_grob(rect) +
  draw_grob(lab_i)

p9

###

tree_shadow2 = tree_shadow * 100

p_10 = ggR(raster(tree_shadow2),
          geom_raster = TRUE) +
  geom_sf(data = filter(trees_pols, plot == "Corrected"), 
          color = "red3", fill = NA, size = .8) +
  scale_fill_continuous(limits = c(10, 100),
                        na.value = "black",
                        high = "white",
                        low = "black",                   
                        breaks = c(15, 55, 95),
                        name = "Illumination (%)",
                        guide = guide_colorbar(frame.colour = "black",
                                               ticks = FALSE,
                                               title.position = "top",
                                               title.hjust = 0.5)) +
  theme_void() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(.24, .973),
        legend.justification = c(0, 1),
        legend.direction = "horizontal",
        legend.key.width = unit(.35, "in"),
        legend.key.height = unit(.09, "in"),
        legend.text = element_text(family = "Cambria", size = 15, margin = margin(-3, 3, 3, 3)),
        legend.title = element_text(family = "Cambria", size = 17),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(3, -3, 3, 5)) 

p10 = ggdraw(p_10) +
  draw_grob(rect) +
  draw_grob(lab_j)

p10
###

p_11 = ggR(raster(tree_shadow2),
           geom_raster = TRUE) +
  scale_fill_continuous(limits = c(10, 100),
                        na.value = "black",
                        high = "white",
                        low = "black",  
                        guide = "none") +
  ggnewscale::new_scale_fill() +
  ggR(tree_shadow_mask, ggLayer=T, geom_raster=T) + 
  scale_fill_gradient(high = "black", low = "black", na.value = NA, guide = "none") +
  geom_sf(data = filter(trees_pols, plot == "Corrected"), 
          color = "red3", fill = NA, size = .8) +
  theme_void() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) 

p11 = ggdraw(p_11) +
  draw_grob(rect) +
  draw_grob(lab_k)

p11
###

p_12 = ggRGB(brick(tree_ms), r = 6, g = 4, b = 2,
            scale = 255, 
            stretch = "lin",
            quantiles = c(0, 1)) +
  geom_sf(data = filter(trees_pols, plot == "Corrected"), 
          color = "red3", fill = NA, size = .8) +
  theme_void() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "grey50"),
        legend.position = c(.18, .97),
        legend.justification = c(0, 1),
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key.width = unit(.3, "in"),
        legend.key.height = unit(.3, "in"),
        legend.text = element_text(family = "Cambria", size = 14),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(6, 6, 6, 6))

p12 = ggdraw(p_12) +
  draw_grob(rect) +
  draw_grob(lab_l)

p12

###

p_13 = ggRGB(brick(tree_ms), r = 6, g = 4, b = 2,
             scale = 255, 
             stretch = "lin",
             quantiles = c(0, 1)) +
  ggR(tree_shadow_mask, ggLayer=T, geom_raster=T) + 
  scale_fill_gradient(high = "black", low = "black", na.value = NA, guide = "none") +
  geom_sf(data = filter(trees_pols, plot == "Corrected"), 
          color = "red3", fill = NA, size = .8) +
  theme_void() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "grey50"),
        legend.position = c(.18, .97),
        legend.justification = c(0, 1),
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key.width = unit(.3, "in"),
        legend.key.height = unit(.3, "in"),
        legend.text = element_text(family = "Cambria", size = 14),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(6, 6, 6, 6))

p13 = ggdraw(p_13) +
  draw_grob(rect) +
  draw_grob(lab_m)

p13

###

p_14 = ggRGB(brick(tree_ms_mask), r = 6, g = 4, b = 2,
             scale = 255, 
             stretch = "lin",
             quantiles = c(0, 1)) +
  geom_sf(data = filter(trees_pols, plot == "Corrected"), 
          color = "red3", fill = NA, size = .8) +
  theme_void() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "white"),
        legend.position = c(.18, .97),
        legend.justification = c(0, 1),
        legend.title.align = 1,
        legend.direction = "horizontal",
        legend.key.width = unit(.3, "in"),
        legend.key.height = unit(.3, "in"),
        legend.text = element_text(family = "Cambria", size = 14),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", size = .5),
        legend.margin = margin(6, 6, 6, 6))

p14 = ggdraw(p_14) +
  draw_grob(rect) +
  draw_grob(lab_n)

p14


################################################################################
g1 = ggarrange(NULL, p1, NULL, p2, NULL, p3, NULL, p4, NULL, p5, NULL, p6, NULL, p7, NULL,
               nrow = 1, ncol = 15,
               widths = c(.02, 1, .02, 1, .02, 1, .02, 1, .02, 1, .02, 1, .02, 1, .02))

g2 = ggarrange(NULL, p8, NULL, p9, NULL, p10, NULL, p11, NULL, p12, NULL, p13, NULL, p14, NULL,
               nrow = 1, ncol = 15,
               widths = c(.02, 1, .02, 1, .02, 1, .02, 1, .02, 1, .02, 1, .02, 1, .02))

g3 = ggarrange(NULL, g1, NULL, g2, NULL, nrow = 5, ncol = 1, 
               heights = c(.008, 1, .008, 1, .008))

ggsave(plot = g3, 
       filename = "D:\\Sync\\Figures\\_FIGURES_paper3\\figure_3.tiff", 
       device = tiff,
       width = 25,
       height = 11.3,
       units = "in",
       dpi = 300,
       bg = "white")

ggsave(plot = g3, 
       filename = "D:\\Sync\\Figures\\_FIGURES_paper3\\figure_3.jpeg", 
       device = jpeg,
       width = 25,
       height = 11.3,
       units = "in",
       dpi = 300,
       bg = "white")


################################################################################
# Create a text grob for the label square

library(jpeg)
library(magick)
library(grid)

rect = rectGrob(
  x = unit(.15, "in"),
  y = unit(1, "npc") - unit(.15, "in"),
  width = unit(.6, "in"),
  height = unit(.6, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

lab_a = textGrob(
  label = "a",
  x = unit(.32, "in"),
  y = unit(1, "npc") - unit(.24, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))


lab_b = textGrob(
  label = "b",
  x = unit(.32, "in"),
  y = unit(1, "npc") - unit(.25, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

lab_c = textGrob(
  label = "c",
  x = unit(.32, "in"),
  y = unit(1, "npc") - unit(.24, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 40,
            fontfamily = "Cambria"))

image1 = jpeg::readJPEG("D:\\Sync\\Figures\\Drone_field_pics\\Whitecourt3.jpg")

image2 = jpeg::readJPEG("D:\\Sync\\Figures\\Drone_field_pics\\M200.jpg")

p1 = ggdraw() +
  draw_image(image1) +
  theme(panel.background = element_rect(color = NA, fill = "black", size = .5),
        panel.border = element_rect(color = "black", fill = NA, size = .5)) +
  draw_grob(rect) +
  draw_grob(lab_a)


p2 = ggdraw() +
  draw_image(image2) +
  theme(panel.background = element_rect(color = NA, fill = "black", size = .4),
        panel.border = element_rect(color = "black", fill = NA, size = .4),
        legend.position = "none") +
  draw_grob(rect)+
  draw_grob(lab_b)

#p2  

bands = read.csv("D:\\Sync\\Figures\\bands.csv") %>% 
  mutate(rig_num = if_else(rig == "blue", 1.3, 2.3))



p_3 = ggplot(bands, aes(y = center, x = rig)) +
  coord_flip(ylim=c(430, 870), clip = "off") +
  geom_crossbar(aes(fill = Band, color = Band,
                    ymin = (center - (width / 2)), ymax = (center + (width / 2))),
                width = .5) +
  geom_errorbar(aes(ymin = (center - (width / 2)), ymax = (center + (width / 2)), x = rig_num),
                width = .05) +
  scale_fill_manual(values = c("blue 444" = "royalblue4",
                               "blue 475" = "steelblue",
                               "green 531" = "springgreen3",
                               "green 560" = "forestgreen",
                               "red 650" = "firebrick2",
                               "red 668" = "red3",
                               "red edge 705" = "indianred",
                               "red edge 717" = "lightpink3",
                               "red edge 740" = "pink4",
                               "NIR 842" = "thistle4")) +
  scale_color_manual(values = c("blue 444" = "royalblue4",
                                "blue 475" = "steelblue",
                                "green 531" = "springgreen3",
                                "green 560" = "forestgreen",
                                "red 650" = "firebrick2",
                                "red 668" = "red3",
                                "red edge 705" = "indianred",
                                "red edge 717" = "lightpink3",
                                "red edge 740" = "pink4",
                                "NIR 842" = "thistle4")) +
  geom_text(aes(label = center), vjust = 6.2, size = 4, fontface = "bold", family = "Cambria") +
  geom_text(aes(label = width), vjust = -6.5, size = 4, fontface = "bold", family = "Cambria") +
  # annotate(geom = "label", label = "c", y = 450, x = 2.5,
  #          hjust = 4,
  #          size = 20, family = "Cambria") +
  lims(y = c(450, 850)) +
  labs(x = "", y = "Wavelength (nm)") +
  scale_y_continuous(breaks = seq(400, 850, 50), limits = c(410, 890)) +
  scale_x_discrete(labels = c("red" = "RedEdge\nMX Red", "blue" = "RedEdge\nMX Blue")) +
  theme_bw() +
  theme(text=element_text(family="Cambria", size=14, face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(size = 14, color = "black", face = "plain"),
        axis.text.y = element_text(margin = margin(r = 10, l = -10), vjust = 0.5, hjust = 0.5, color = "black", size = 16),
        plot.background = element_rect(color = "black", fill = NA, size = 1.1),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        plot.margin = margin(.4, .4, .4, .4, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = .5))


p3 = ggdraw(p_3) +
  draw_grob(rect)+
  draw_grob(lab_c)


g1 = ggarrange(NULL, p2, NULL, p3, NULL, ncol = 1,
               heights = c(.004, 1, .012, 1, .004))

g2 = ggarrange(NULL, p1, NULL, ncol = 1,
               heights = c(0.0025, 1, 0.0025))

g3 = ggarrange(NULL, g2, NULL, g1, NULL, ncol = 5,
               widths = c(.002, 1, .01, 1, .002))

#g2

ggsave(plot = g3, 
       filename = "D:\\Sync\\Figures\\_FIGURES_paper3\\figure_2.tiff", 
       device = tiff,
       width = 15,
       height = 10,
       units = "in",
       dpi = 600,
       bg = "white")

ggsave(plot = g3, 
       filename = "D:\\Sync\\Figures\\_FIGURES_paper3\\figure_2.jpeg", 
       device = jpeg,
       width = 15,
       height = 10,
       units = "in",
       dpi = 600,
       bg = "white")


################################################################################


ggplot(las@data, aes(x = X, y = Z)) +
  geom_point() +
  #scale_alpha(range = c(0, .5)) +
  #geom_line() +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "none")

las_green = las %>% 
  filter_poi(G > R)

dat = las@data %>% 
  dplyr::mutate(RGB = rgb(R, G, B, maxColorValue = 65535))

ggplot(dat, aes(X, Z, color = RGB)) + 
  geom_point(size = 1) + 
  scale_color_identity() +
  coord_equal() +
  scale_y_continuous(limits = c(0, 6.25),
                     breaks = seq(0, 6, by = 1),
                     expand = c(0, .2)) +
  scale_x_continuous(limits = c((mean(dat$X) - 1.9), (mean(dat$X) + 1.9))) +
  geom_hline(yintercept = Zq999, color = "red", size = 1) +
  geom_hline(yintercept = Zq99, color = "red", size = 1) +
  geom_hline(yintercept = Z_mean, color = "red", size = 1) +
  annotate(geom = "text",
           color = "white",
           label = paste0(expression("Z[Q999]")," == ", round(Zq999, digits = 2)),
           parse = TRUE,
           size = 8,
           x = mean(dat$X) + 1.5,
           y = Zq999,
           hjust = .5,
           vjust = -.25,
           fontface = "bold",) +
  annotate(geom = "text",
           color = "white",
           label = paste0(expression("Z[Q99]")," == ", round(Zq99, digits = 2)),
           parse = TRUE,
           size = 8,
           x = mean(dat$X) + 1.5,
           y = Zq99,
           hjust = .5,
           vjust = 1.2,
           fontface = "bold",) +
  annotate(geom = "text",
           color = "white",
           label = paste0(expression("mu[Z]")," == ", round(Z_mean, digits = 2)),
           parse = TRUE,
           size = 8,
           x = mean(dat$X) + 1.5,
           y = Z_mean,
           hjust = .5,
           vjust = 1.2,
           fontface = "bold",) +
  theme_bw(base_size = 30) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "grey20"),
        panel.grid.major.y = element_line(color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())



