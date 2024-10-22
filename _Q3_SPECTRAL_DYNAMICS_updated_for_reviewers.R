
#save.image(file = "D:\\Sync\\_Sites\\Skimikin_spectral\\R_scripts\\_Q3_SPECTRAL_DYNAMICS.RData")
#load("D:\\Sync\\_Sites\\Skimikin_spectral\\R_scripts\\_Q3_SPECTRAL_DYNAMICS.RData")

library(sf)
library(terra)
library(grid)
library(ggh4x)
library(lmerTest)
library(ggpubr)
library(cowplot)
library(corrr)
library(emmeans)
library(modelr)
library(ggnewscale)
library(cluster)
library(factoextra)

library(tidyverse)

###########################################
# ggplot aesthetic variables

rect = rectGrob(
  x = unit(.1, "in"),
  y = unit(1, "npc") - unit(.1, "in"),
  width = unit(.5, "in"),
  height = unit(.5, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

font_s = 22


timepoint_cols = c("late winter" = "darkslategray3",
                   "early summer" = "darkolivegreen4",
                   "mid summer" = "goldenrod",
                   "late summer" =  "goldenrod4")


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

# read in polygons, join to seedlot info, hybrid status, field assessment, save.
################################################################################


pols = st_read(paste0(dir, "\\output\\HULLS\\Edits\\Skimikin_z50_updated.shp")) %>%
  filter(!st_is_empty(.)) %>% 
  dplyr::mutate(DBH18 = if_else(Obs == 1912, "70", DBH18)) %>% # seems to be missing a zero 
  left_join(read.csv("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\Sx_MASTER_ver_21_ADDITIONAL_TRAITS.csv") %>% 
              filter(Site == "Skim") %>% 
              dplyr::select(Meas_seq, Blk, Weev16),
            by = c("Obs" = "Meas_seq"))

seedlot_info = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\Sx_CC_Seedlots_info.csv") %>% 
  dplyr::select(-number, -donor)

hybrid_prop = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\SxGenecology_predicted_Euc.csv") %>% 
  dplyr::select(sxProv, species, propGla, propEng, propSit)

pols_dat_skim = pols %>% 
  dplyr::rename("Seedlot" = "Prov") %>% 
  left_join(seedlot_info,
            by = c("Seedlot" = "population")) %>% 
  dplyr::mutate(zone = case_when(region == "NM" | region == "AZ" ~ "US (south)",
                          region == "ID" | region == "MT" | region == "WA" ~ "US (north)",
                          region == "BC" ~ "British Columbia",
                          region == "AB" ~ "Alberta",
                          region == "YK" | region == "NWT" ~ "Canada (north)",
                          region == "ON" ~ "Ontario"),
                HT16 = as.numeric(HT16),
                DBH16 = as.numeric(DBH16),
                DBH18 = as.numeric(DBH18)) %>%
  left_join(hybrid_prop,
            by = c("Seedlot" = "sxProv")) %>% 
  dplyr::rename("Prov" = "Seedlot") %>% 
  dplyr::mutate(Tree = as.numeric(Tree)) %>% 
  filter(class == "B" & HT16 >= 130 & DBH18 > 0) %>% 
  mutate(Class = class) %>% 
  #dplyr::mutate(Class = if_else(region %in% c("BC", "ON") & class != "B", "Orchard", "Wildstand")) %>% 
  as.data.frame() %>% 
  st_drop_geometry()


saveRDS(pols_dat_skim, paste0(dir, "\\CSV\\traits//pols_dat_skim.rds"))

################################################################################
# means, range, standrad deviation of tree heights and DBH

tree_size = pols_dat_skim %>% 
  dplyr::select(Obs, HT16, DBH16, DBH18, Prov) %>% 
  distinct()

tree_size %>% 
  ungroup() %>% 
  summarise_all("max")

################################################################################

ratios = readRDS(paste0(dir, "\\CSV\\Corrected_values//ratios.rds")) %>% 
  filter(!grepl("_80m", Date)) %>% 
  dplyr::mutate(Date = as_date(Date_mean))


df_spectral_all = readRDS(paste0(dir, "\\CSV\\TRAITS\\Updated\\Skimikin_all_spectral_average_median.rds")) %>% 
  #dplyr::select(-`2020-03-22`) %>% 
  pivot_longer(`2020-07-02`:`2023-02-24`, names_to = "trait", values_to = "value") %>% 
  #dplyr::select(-(mean_summer:decline)) %>% 
  left_join(pols_dat_skim %>% dplyr::select(Obs, Blk, Rep, Prov, Class, HT16:DBH18), by = "Obs") %>% 
  dplyr::mutate(timepoint = case_when(substr(trait,7,7) %in% c("2", "3") ~ "late winter",
                                      trait %in% c("2020-07-02", "2021-06-29") ~ "early summer",
                                      trait %in% c("2020-08-07", "2021-07-29") ~ "mid summer",
                                      trait %in% c("2020-08-25", "2021-08-14") ~ "late summer")) %>% 
  dplyr::rename(Date = trait) %>% 
  dplyr::mutate(value = if_else(index %in% c("RE_upper", "RE_lower", "RE_total"), value * 100, value),
                index2 = recode(index,
                "RE_upper" = "RE['upper']",
                "NDRE3" = "NDRE['740']",
                "EWI9" = "EWI['9']",
                "Datt" = "mDatt"),
                index2 = factor(index2,
                levels = c("CCI", "PRI", "GCC", "BCC", "ARI", "EWI['9']",
                           "mDatt", "RE['upper']", "NDVI", "NDRE['740']")))

################################################################################

keep_list = c(
  "CCI",
  "PRI",
  "ARI",
  "GCC",
  #"EWI9",
  "mDatt",
  "NDVI",
  "RE_upper"
    )

# keep_list = c(
#   "CCI", "BCC", "ARI", 
#   "EWI9", "mDatt", "NDRE3"
# )

R_list = c("R444",
           "R475",
           "R531",
           "R560",
           "R650",
           "R668",
           "R705",
           "R717",
           "R740",
           "R842")

field_list = c("HT16", "DBH16", "DBH18")

(timepoint_boxplots = df_spectral_all %>%  
       filter(index %in% keep_list) %>% 
       #dplyr::mutate(Class = factor(Class, levels = c("Wildstand", "Orchard"))) %>% 
       ggplot() +
       geom_boxplot(#outlier.alpha = 0, 
                    width = .5, alpha = .7,
                    position = position_dodge(width = .75),
                    linewidth = .3,
                    aes(fill = timepoint,
                        x = factor(Date),
                        y = value)) +
       theme_bw(base_size = 15) +
       #scale_y_discrete(limits = rev(levels(blups_r2_plot$clim))) +
       scale_fill_manual(values = timepoint_cols) +
       scale_y_continuous(expand = expansion(mult = c(.09, .09))) +
       theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, margin = margin(4, 0, -8, 0)),
             panel.grid.major.x = element_blank(),
             axis.title = element_blank(),
             strip.background = element_blank(),
             axis.ticks.x = element_blank(),
             aspect.ratio = 1.3,
             legend.position = "none",
             strip.placement = "outside",
             axis.text.y = element_text(size = 9.5)) +
       #coord_cartesian(expand = FALSE) +
       facet_wrap2(. ~ index2,
                   labeller = label_parsed,
                   nrow = 2,
                   strip.position = "left",
                   scales = "free_y") 
  +
       facetted_pos_scales(
         y = list(
           index2 == "CCI" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                breaks = c(.05, .15, .25)),
           index2 == "GCC" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                breaks = seq(.42, .5, .04)),
           index2 == "PRI" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                breaks = c(-.03, -.09, -.15)),
           index2 == "BCC" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                breaks = c(.24, .27, .3)),
           index2 == "ARI" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                breaks = c(3, 6, 9),
                                                labels = c("3.0", "6.0", "9.0")),
           index2 == "EWI['9']" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                     breaks = c(-.52, -.60, -.68)),
           index2 == "mDatt" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                  breaks = c(.84, .86, .88)),
           index2 == "RE['upper']" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                        breaks = c(.45, .60, .75)),
           index2 == "NDVI" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                 breaks = c(.82, .85, .88)),
           index2 == "NDRE['740']" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                        breaks = c(.12, .15, .18))
         )
       )
  )
  


(timepoint_boxplots = df_spectral_all %>%  
    filter(index %in% R_list) %>% 
    #dplyr::mutate(Class = factor(Class, levels = c("Wildstand", "Orchard"))) %>% 
    ggplot() +
    geom_boxplot(#outlier.alpha = 0, 
                 width = .5, alpha = .7,
                 position = position_dodge(width = .75),
                 linewidth = .3,
                 aes(fill = timepoint,
                     x = factor(Date),
                     y = value)) +
    # geom_point(position = position_dodge(width=0.75), 
    #            shape = 21, stroke = NA, size = 3, alpha = .3,
    #            aes(group=Date, fill = timepoint,
    #                x = factor(Date),
    #                y = value)) +
    # geom_text(data = filter(blups_p_VI, p.comp < .1), aes(x = x_loc, y = max_plot, label = p_plot),
    #           size = 6, vjust = .2) +
    # geom_text(data = filter(blups_p_VI, p_plot == "NS"), aes(x = x_loc, y = max_plot, label = p_plot),
    #           size = 3.5, vjust = -.4) +
    # geom_linerange(data = blups_p_VI, aes(xmin = x_start, xmax = x_end, y = max_plot),
    #                color = "black") +
    # # geom_text(data = blups_p_VI, aes(x = x_loc, y = max_comp),
    # #           label = "—",
    # #           size = 6,
    # #           vjust = -.2) +
    theme_bw(base_size = 15) +
    #scale_y_discrete(limits = rev(levels(blups_r2_plot$clim))) +
    scale_fill_manual(values = timepoint_cols) +
    scale_y_continuous(expand = expansion(mult = c(.09, .09))) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, margin = margin(4, 0, -8, 0)),
          panel.grid.major.x = element_blank(),
          axis.title = element_blank(),
          strip.background = element_blank(),
          axis.ticks.x = element_blank(),
          aspect.ratio = 1.3,
          legend.position = "none",
          strip.placement = "outside",
          axis.text.y = element_text(size = 9.5)) +
    #coord_cartesian(expand = FALSE) +
    facet_wrap2(. ~ index,
                labeller = label_parsed,
                nrow = 2,
                strip.position = "left",
                scales = "free_y") 
  # +
  #      facetted_pos_scales(
  #        y = list(
  #          index2 == "CCI" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
  #                                               breaks = c(.05, .15, .25)),
  #          index2 == "GCC" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
  #                                               breaks = seq(.42, .5, .04)),
  #          index2 == "PRI" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
  #                                               breaks = c(-.03, -.09, -.15)),
  #          index2 == "BCC" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
  #                                               breaks = c(.24, .27, .3)),
  #          index2 == "ARI" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
  #                                               breaks = c(3, 6, 9),
  #                                               labels = c("3.0", "6.0", "9.0")),
  #          index2 == "EWI['9']" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
  #                                                    breaks = c(-.52, -.60, -.68)),
  #          index2 == "mDatt" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
  #                                                 breaks = c(.84, .86, .88)),
  #          index2 == "RE['upper']" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
  #                                                       breaks = c(.45, .60, .75)),
  #          index2 == "NDVI" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
  #                                                breaks = c(.82, .85, .88)),
  #          index2 == "NDRE['740']" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
  #                                                       breaks = c(.12, .15, .18))
  #        )
  #      )
)


################################################################################
# variance components for spectral indices
mm_dat_init = df_spectral_all %>% 
  #filter(trait != "2020-03-22") %>% 
 dplyr::rename("var_value" = "value") %>% 
  dplyr::mutate(Blk = factor(Blk),
         Rep = factor(Rep),
         Prov = factor(Prov),
         Class = factor(Class)) %>% 
  dplyr::mutate(timepoint = case_when(substr(Date,7,7) %in% c("2", "3") ~ "late winter",
                               Date %in% c("2020-07-02", "2021-06-29") ~ "early summer",
                               Date %in% c("2020-08-07", "2021-07-29") ~ "mid summer",
                               Date %in% c("2020-08-25", "2021-08-14") ~ "late summer"),
         Season = if_else(substr(Date,7,7) %in% c("2", "3"), "Winter", "Summer"),
         #Date = as_date(trait),
         Year = if_else(Date < "2021-06-01", "2020", "2021")) %>% 
  dplyr::group_by(index, Obs) %>% 
  # drop any trees with incomplete data
  dplyr::mutate(mean_obs = mean(var_value, na.rm = FALSE))

# number of trees dropped because complete spectra not available
mm_dat_init %>% 
  filter(is.na(mean_obs)) %>% 
  ungroup() %>% 
  dplyr::distinct(Obs) %>% 
  tally()
  
# 220 mosaic median
# 250 average mean 15 pixels
# 219 average mean 5 pixels
# 235 average mean 10 pixels
# 270 average mean 30 pixels
  
mm_dat = mm_dat_init %>% 
  drop_na(mean_obs) %>% # drop incomplete datasets
  drop_na(Prov) %>% 
  dplyr::select(-mean_obs)

mm_dat_count = mm_dat %>% 
  ungroup() %>% 
  distinct(Obs, Prov) %>% 
  group_by(Prov) %>% 
  add_count()

# min 3, max 16
  
#df = mm_dat %>% filter(index == "EWI9")

blup_mod = function(df) {
  # fit a random effects model
  mm = lmer(var_value ~ 
              (1|Prov) +
              (1|Rep) +
              (1|Blk:Rep),
            data = df, REML = TRUE)
  }

blup_mod_all = function(df) {
  mm_all = lmer(var_value ~ 
                  (1|timepoint:Prov) + 
                  (1|timepoint:Rep) + 
                  (1|timepoint:Blk:Rep) +
                  (1|timepoint) +
                  (1|timepoint:Date),
            data = df, REML = TRUE)
  }

df_mm = mm_dat %>%
  filter(is.finite(var_value)) %>% 
  group_by(index, Date, timepoint) %>% 
  nest() %>% 
  dplyr::mutate(mm = map(data, blup_mod)) %>%
  dplyr::mutate(var_comp = map(mm, VarCorr)) %>% 
  dplyr::mutate(var_df = map(var_comp, as.data.frame)) %>% 
  unnest(c(var_df)) %>% 
  dplyr::select(-var_comp, -mm, -data) 

df_mm_field = mm_dat %>%
  ungroup() %>% 
  dplyr::select(Obs, Blk:Prov, HT16:DBH18) %>% 
  pivot_longer(HT16:DBH18, names_to = "index", values_to = "var_value") %>% 
  group_by(index) %>% 
  nest() %>% 
  dplyr::mutate(mm = map(data, blup_mod)) %>%
  dplyr::mutate(var_comp = map(mm, VarCorr)) %>% 
  dplyr::mutate(var_df = map(var_comp, as.data.frame)) %>% 
  unnest(c(var_df)) %>% 
  dplyr::select(-var_comp, -mm, -data) 

df_mm_all = mm_dat %>%
  filter(is.finite(var_value)) %>% 
  group_by(index) %>% 
  nest() %>% 
  dplyr::mutate(mm = map(data, blup_mod_all)) %>%
  dplyr::mutate(var_comp = map(mm, VarCorr)) %>% 
  dplyr::mutate(var_df = map(var_comp, as.data.frame)) %>% 
  unnest(c(var_df)) %>% 
  dplyr::select(-var_comp, -mm, -data)

saveRDS(df_mm, "D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\BLUPS\\Updated\\df_mm.rds")
saveRDS(df_mm_field, "D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\BLUPS\\Updated\\df_mm_field.rds")
saveRDS(df_mm_all, "D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\BLUPS\\Updated\\df_mm_all.rds")

df_mm = readRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\BLUPS\\Updated\\df_mm.rds")
df_mm_field = readRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\BLUPS\\Updated\\df_mm_field.rds")
df_mm_all = readRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\BLUPS\\Updated\\df_mm_all.rds") 

components = df_mm %>% 
  group_by(index, timepoint, Date) %>% 
  dplyr::mutate(var_sum = sum(vcov),
         vcov_percent = round((vcov/var_sum) * 100, 1)) %>% 
  dplyr::select(index, Date, grp, vcov_percent) %>% 
  pivot_wider(names_from = grp, values_from = vcov_percent) %>%
  bind_rows(
    df_mm_field %>% 
      group_by(index) %>% 
      dplyr::mutate(var_sum = sum(vcov),
                    vcov_percent = round((vcov/var_sum) * 100, 1)) %>% 
      dplyr::select(index, grp, vcov_percent) %>% 
      pivot_wider(names_from = grp, values_from = vcov_percent)
  ) %>% 
  dplyr::mutate(type = case_when(index %in% R_list ~ "band",
                                 index %in% keep_list ~"VI",
                                 index %in% c("HT16", "DBH16", "DBH18") ~ "field"))


components_all = df_mm_all %>% 
  group_by(index) %>% 
    dplyr::mutate(var_sum = sum(vcov),
           vcov_percent = round((vcov/var_sum) * 100, 1)) %>% 
    dplyr::select(index, grp, vcov_percent) %>% 
    pivot_wider(names_from = grp, values_from = vcov_percent) %>% 
  dplyr::mutate(type = if_else(index %in% c("R444","R475", "R531","R560","R650",               
                                     "R668","R705", "R717","R740","R842"),
                        "band", "VI")) %>% 
  ungroup() %>% 
  dplyr::mutate(order = if_else(type == "VI",
                         dense_rank(`timepoint:Prov` + `timepoint`),
                         dense_rank(as.numeric(substr(index, 2, 4)))))

vpop = components %>% 
  dplyr::mutate(vpop = Prov / (Prov + Residual)) %>% 
  dplyr::mutate(trait = if_else(!is.na(Date), Date, timepoint))

vpop_all = components_all %>% 
  dplyr::mutate(vpop = `timepoint:Prov` / (`timepoint:Prov` + Residual))

vpop_plot = vpop %>%
  filter(type == "band") %>% 
  dplyr::mutate(R_num = as.numeric(substr(index, 2,4)),
         Date = as_date(trait)) %>% 
  full_join(ratios %>% 
              dplyr::select(Date, mean_scattered), by = "Date") %>% 
  group_by(Date) %>% 
  dplyr::mutate(vpop_mean = mean(vpop),
         vpop_rel = vpop - vpop_mean)

vpop_plot_all = vpop_all %>%
  filter(type == "band") %>% 
  dplyr::mutate(R_num = as.numeric(substr(index, 2,4))) %>% 
  ungroup() %>% 
  dplyr::mutate(vpop_mean = mean(vpop),
         vpop_rel = vpop - vpop_mean)



(p2 = vpop_plot %>% 
  ggplot(aes(x = R_num, y = vpop_rel)) +
  geom_line(aes(color = timepoint, group = trait), linewidth = 1.1) +
    geom_line(data = filter(vpop_plot_all, type == "band"), 
              linewidth = 1.1, linetype = 5, color = "grey30") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5),
        legend.position = "none",
        aspect.ratio = 1) +
  scale_x_continuous(breaks = as.numeric(levels(factor(vpop_plot$R_num))),
                     minor_breaks = NULL) +
  scale_color_manual(values = timepoint_cols) +
  labs(x = "Band central wavelength (nm)",
       y = "Vpop, de-meaned"))

# vpop_plot %>% 
#   ggplot(aes(x = timepoint, y = vpop_rel, color = index, group = year)) +
#   geom_line(linewidth = 1.1) +
#   theme_bw(base_size = 16) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_color_manual(values = timepoint_cols) +
#   #geom_text(aes(label = year)) +
#   labs(x = "Band central wavelength (nm)",
#        y = "Proportion of variance due to population differences (Vpop)")



# vpop_plot %>%
#   group_by(timepoint, index) %>% 
#   dplyr::mutate(vpop_mean_nm = mean(vpop_rel),
#          timepoint = factor(timepoint, levels = c("late winter",
#                                                   "early summer",
#                                                   "mid summer",
#                                                   "late summer"))) %>% 
#   #filter(index == "R740") %>% 
#   ggplot(aes(x = Date, y = vpop_rel, group = index)) +
#   geom_line(aes(color = index), alpha = .5) +
#   geom_point(aes(color = index), size = 1.5) +
#   theme_bw() +
#   labs(y = "Proportion of variance due to population differences (Vpop)",
#        x = "Proportion of scattered irradiance") +
#   theme(aspect.ratio = 1)

# (p4 = ggarrange(p1, p2,
#                     ncol = 1, align = "h"))
# 
# (p5 = ggarrange(p3, p4, widths = c(1, .75)))
# 
# ggsave(plot = p5, 
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\figure_2.jpeg", 
#        device = jpeg,
#        width = 15,
#        height = 8,
#        units = "in",
#        dpi = 300,
#        bg = "white")

################################################################################

index_cols = c("Carotenoid" = "steelblue4",
                   "Greenness" = "forestgreen",
                   "Red Edge" =  "orchid4")

comp_cols = c("Residual" = "grey50", 
             "timepoint:Blk:Rep" = "chocolate3", 
             "timepoint:Rep" = "chocolate4",
             "timepoint:Date" = "orange4", 
             "timepoint:Prov" = "steelblue4", 
             #"timepoint:Class" = "steelblue", 
             "timepoint" = "red4")

################################################################################

# 
# (compare = mm_dat %>%
#     filter(!(index %in% c("R444",
#                           "R475",
#                           "R531",
#                           "R560",               
#                           "R650",               
#                           "R668",               
#                           "R705",               
#                           "R717",               
#                           "R740",               
#                           "R842")) 
#            #& Date != "2020-03-22"
#            ) %>% 
#     ungroup() %>% 
#    dplyr::rename("value" = "var_value") %>% 
#     # dplyr::mutate(order = dense_rank(`Date:timepoint`)) %>% 
#     # pivot_longer(`Blk:Rep`:Residual, names_to = "component") %>% 
#     ggplot(aes(x = Prov,
#                y = value)) +
#     geom_point(size = 1, alpha = .2, aes(color = factor(Date))) +
#     stat_smooth(aes(group = Date)) +
#     theme_bw(base_size = 10) +
#     theme(axis.text.x = #element_text(angle = 90)) +
#             element_blank()) +
#     facet_grid(index ~ Season,
#                #ncol = 3,
#                scales = "free"))

################################################################################


index_labs = c(RE_upper = expression(RE['upper']),
               NDRE3 = expression(NDRE['740']),
               EWI9 = expression(EWI['9']),
               Datt = "mDatt")

(p1 = vpop_plot %>% 
  ggplot(aes(x = R_num, y = vpop)) +
   coord_cartesian(expand = FALSE) +
  geom_vline(aes(xintercept = R_num)) +
    geom_line(aes(color = timepoint, group = trait), linewidth = 1.5) +
  geom_line(data = filter(vpop_plot_all, type == "band"), 
            linewidth = 1.5, color = "grey30") +
  theme_bw(base_size = 16) +
   scale_x_continuous(breaks = as.numeric(levels(factor(vpop_plot$R_num))),
                      minor_breaks = NULL,
                      expand = expansion(mult = c(.01, .01))) +
   scale_y_continuous(expand = expansion(mult = c(.01, .01)),
                      limits = c(0, .55),
                      breaks = seq(0, .5, .1),
                      labels = c("0", "10", "20", "30", "40", " 50")) +
  scale_color_manual(values = timepoint_cols) +
  labs(x = "Band central wavelength (nm)",
       y = expression(Population~differentiation~(V[pop]~", %"))) +
   theme_bw(base_size = 16) +
   #scale_fill_manual(values = timepoint_cols) +
   theme(axis.text.x = element_text(angle = 90, vjust = .5, size = 14),
         axis.text.y = element_text(size = 14),
         axis.title.x = element_blank(),
         axis.title.y = element_text(size = 18),
         axis.ticks = element_blank(),
         plot.margin = unit(c(.5,.5,.5,.5), 'lines'),
         aspect.ratio = 1,
         legend.position = "none"))

# (p1 = vpop %>% 
#     left_join(components_all %>% 
#                 mutate(vpop_all = (`timepoint:Prov`) / (`timepoint:Prov` + Residual)) %>% 
#                 dplyr::select(index, vpop_all)) %>% 
#     filter(type == "band" & index %in% R_list) %>% 
#     ungroup() %>% 
#     dplyr::mutate(Timepoint = factor(timepoint,
#                                      levels = c("early summer", "mid summer", "late summer", "late winter")),
#                   time_num = as.numeric(dense_rank(Timepoint))) %>% 
#     dplyr::arrange(time_num) %>% 
#     ggplot(aes(x = reorder(index, vpop_all),
#                y = vpop, 
#                #group = Date,
#                fill = Timepoint,
#                group = Date)) +
#     coord_cartesian(expand = FALSE) +
#     geom_bar(stat = "identity", position = "dodge") +
#     # geom_line(aes(y = vpop_all),
#     #           linewidth = 1.5, color = "grey30") +
#     theme_bw(base_size = 16) +
#     scale_fill_manual(values = timepoint_cols) +
#     scale_x_discrete(labels = index_labs) +
#     scale_y_continuous(expand = expansion(mult = c(.01, .01)),
#                        limits = c(0, .55),
#                        breaks = seq(0, .5, .1),
#                        labels = c("0", "10", "20", "30", "40", " 50")) +
#     labs(x = "Vegetation index",
#          y = expression(Population~differentiation~(V[pop]~", %"))) +
#     theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 14),
#           axis.text.y = element_text(size = 14),
#           panel.grid.major.x = element_blank(),
#           #axis.title. = element_blank(),
#           axis.ticks = element_blank(),
#           legend.position = "right",#c(.2,.85),
#           legend.background = element_rect(color = "white"),
#           plot.margin = unit(c(1,1,1,1), 'lines'),
#           aspect.ratio = 1) 
#   # +
#   #   facet_wrap(. ~ Date,
#   #              nrow = 2)
# )

(p6 = vpop %>%
    left_join(components_all %>%
                mutate(vpop_all = (`timepoint:Prov`) / (`timepoint:Prov` + Residual)) %>%
                dplyr::select(index, vpop_all)) %>%
    filter(index %in% keep_list | index %in% R_list | index %in% field_list) %>%
    mutate(timepoint = if_else(index %in% field_list, "mid summer", timepoint)) %>%
    ungroup() %>%
    dplyr::mutate(order = if_else(type == "VI" | type == "field",
                                  dense_rank(vpop_all),
                                  dense_rank(as.numeric(substr(index, 2, 4))))) %>%
    dplyr::mutate(Timepoint = factor(timepoint,
                                     levels = c("early summer", "mid summer", "late summer", "late winter")),
                  time_num = as.numeric(dense_rank(Timepoint))) %>%
    dplyr::arrange(time_num) %>%
    ggplot(aes(x = reorder(index, order),
               y = vpop,
               #group = Date,
               fill = Timepoint,
               group = Date)) +
    coord_cartesian(expand = FALSE) +
    geom_bar(stat = "identity", position = "dodge") +
    # geom_line(aes(y = vpop_all),
    #           linewidth = 1.5, color = "grey30") +
    theme_bw(base_size = 16) +
    scale_fill_manual(values = timepoint_cols) +
    scale_x_discrete(labels = index_labs) +
    scale_y_continuous(expand = expansion(mult = c(.01, .01)),
                       limits = c(0, .55),
                       breaks = seq(0, .5, .1),
                       labels = c("0", "10", "20", "30", "40", " 50")) +
    labs(x = "Vegetation index",
         y = expression(Population~differentiation~(V[pop]~", %"))) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing = unit(2, "lines"),
          legend.position = c(.14,.82),
          legend.background = element_rect(color = "black", linewidth = .4),
          plot.margin = unit(c(1,1,1,1), 'lines'))
  +
    facet_grid(. ~ factor(type, levels = c("band", "VI", "field"),
                          labels = c("Spectral band reflectance", "Vegetation index", "Field-assessed trait")),
               #nrow = 1,
               scales = "free_x",
               space = "free")
  )

(p6 = vpop %>%
    left_join(components_all %>%
                mutate(vpop_all = (`timepoint:Prov`) / (`timepoint:Prov` + Residual)) %>%
                dplyr::select(index, vpop_all)) %>%
    filter(index %in% keep_list | index %in% R_list | index %in% field_list) %>%
    mutate(timepoint = if_else(index %in% field_list, "mid summer", timepoint)) %>%
    ungroup() %>%
    dplyr::mutate(order = if_else(type == "VI" | type == "field",
                                  dense_rank(vpop_all),
                                  dense_rank(as.numeric(substr(index, 2, 4))))) %>%
    dplyr::mutate(Timepoint = factor(timepoint,
                                     levels = c("early summer", "mid summer", "late summer", "late winter")),
                  time_num = as.numeric(dense_rank(Timepoint))) %>%
    dplyr::arrange(time_num) %>%
    ggplot(aes(x = reorder(index, order),
               y = vpop,
               #group = Date,
               fill = Timepoint,
               group = Date)) +
    coord_cartesian(expand = FALSE) +
    geom_bar(stat = "identity", position = "dodge") +
    stat_summary(fun="mean", geom="segment",  color = "grey20", linewidth = 1.5, alpha = .7,
                 mapping=aes(y = vpop_all, 
                             xend=..x.. - 0.5, 
                             yend=..y..)) +
    stat_summary(fun="mean", geom="segment", color = "grey20", linewidth = 1.5, alpha = .7,
                 mapping=aes(y = vpop_all, 
                             xend=..x.. + 0.5, 
                             yend=..y..)) +
    # geom_line(aes(y = vpop_all),
    #           linewidth = 1.5, color = "grey30") +
    theme_bw(base_size = 16) +
    scale_fill_manual(values = timepoint_cols) +
    scale_x_discrete(labels = index_labs) +
    scale_y_continuous(expand = expansion(mult = c(.01, .01)),
                       limits = c(0, .55),
                       breaks = seq(0, .5, .1),
                       labels = c("0", "10", "20", "30", "40", " 50")) +
    labs(x = "Vegetation index",
         y = expression(Population~differentiation~(V[pop]~", %"))) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing = unit(2, "lines"),
          legend.position = c(.14,.82),
          legend.background = element_rect(color = "black", linewidth = .4),
          plot.margin = unit(c(1,1,1,1), 'lines'))
  +
    facet_grid(. ~ factor(type, levels = c("band", "VI", "field"),
                          labels = c("Spectral band reflectance", "Vegetation index", "Field-assessed trait")),
               #nrow = 1,
               scales = "free_x",
               space = "free")
  )

(p_field = vpop %>% 
    left_join(components_all %>% 
                mutate(vpop_all = (`timepoint:Prov`) / (`timepoint:Prov` + Residual)) %>% 
                dplyr::select(index, vpop_all)) %>% 
    filter(type == "field") %>% 
    ungroup() %>% 
    dplyr::arrange(vpop) %>% 
    ggplot(aes(x = reorder(index, vpop_all),
               y = vpop)) +
    coord_cartesian(expand = FALSE) +
    geom_bar(stat = "identity", position = "dodge",
             fill = "indianred4") +
    theme_bw(base_size = 16) +
    scale_x_discrete(labels = index_labs) +
    scale_y_continuous(expand = expansion(mult = c(.01, .01)),
                       limits = c(0, .55),
                       breaks = seq(0, .5, .1),
                       labels = c("0", "10", "20", "30", "40", " 50")) +
    labs(x = "Field measurements") +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 14),
          axis.text.y = element_blank(),#element_text(size = 14),
          panel.grid.major.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          legend.justification = "right",
          plot.margin = unit(c(1,1,1,1), 'lines'),
          aspect.ratio =7/3) 
  # +
  #   facet_wrap(. ~ Date,
  #              nrow = 2)
)

(p_3_ = ggarrange(p6, p_field,
                  ncol = 2,
                  heights = c(1,1),
                  widths = c(8,3.5)))

# ggsave(plot = p_3_, 
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated\\figure_3.jpeg", 
#        device = jpeg,
#        width = 11,
#        height = 7,
#        units = "in",
#        dpi = 300,
#        bg = "white")


(p6_legend = as_ggplot(get_legend(p6)))


(p4 = components_all %>% 
    filter(type == "band") %>% 
    ungroup() %>% 
    pivot_longer(`timepoint:Blk:Rep`:Residual, names_to = "component") %>% 
    dplyr::mutate(component = factor(component, 
                              levels = c("Residual", "timepoint:Blk:Rep", "timepoint:Rep", "timepoint:Date", "timepoint:Prov", "timepoint")),
           R_num = as.numeric(substr(index, 2, 4))) %>%
    ggplot(aes(x = R_num,
               y = value,
               fill = component)) +
    coord_cartesian(expand = FALSE) +
    geom_density(stat = "identity", position = "stack",
                 color = NA) +
    geom_vline(aes(xintercept = R_num)) +
    theme_bw(base_size = 16) +
    scale_fill_manual(values = comp_cols) +
    scale_x_continuous(breaks = as.numeric(levels(factor(vpop_plot$R_num))),
                       minor_breaks = NULL,
                       expand = expansion(mult = c(.01, .01))) +
    scale_y_continuous(expand = expansion(mult = c(.01, .01)),
                       labels = c("0", "25", "50", "75", " 100")) +
    labs(y = "Variance component (%)",
         x = "Band central wavelength (nm)") +
    scale_fill_manual(values = comp_cols) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 16),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 18),
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(1,1,1,1), 'lines'),
          aspect.ratio = 1))

(p5 = components_all %>% 
    filter(index %in% keep_list | index %in% R_list) %>% 
    ungroup() %>% 
    dplyr::mutate(order = if_else(type == "VI",
                           dense_rank(`timepoint` + `timepoint:Prov`),
                           dense_rank(as.numeric(substr(index, 2, 4))))) %>%
    pivot_longer(`timepoint:Blk:Rep`:Residual, names_to = "component") %>% 
    dplyr::mutate(Effect = factor(component, 
                                  levels = c("Residual", "timepoint:Blk:Rep", "timepoint:Rep", "timepoint:Date", "timepoint:Prov", "timepoint"))) %>% 
    ggplot(aes(x = reorder(index, order),
               y = value,
               fill = Effect)) +
    coord_cartesian(expand = FALSE) +
    geom_bar(stat = "identity") +
    theme_bw(base_size = 16) +
    labs(x = "Vegetation Index", y = "Percent Variance") +
    scale_fill_manual(values = comp_cols,
                      labels = c("timepoint:Blk:Rep" = "Block in Rep",
                                 "timepoint:Rep" =  "Rep",
                                 "timepoint:Date" = "Date",
                                "timepoint:Prov" = "Population",
                                             #"timepoint:Class" = "steelblue", 
                                             "timepoint" = "Timepoint")) +
    scale_x_discrete(labels = index_labs) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          panel.spacing = unit(2, "lines"),
          legend.box.margin = margin(l = -12, r = 40),
          strip.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(r = 8)),
          axis.ticks = element_blank(),
          legend.background = element_rect(color = "black", linewidth = .4),
          legend.justification = "right") +
    facet_grid(. ~ factor(type, levels = c("band", "VI"),
                          labels = c("Spectral band reflectance", "Vegetation index")),
               #nrow = 1,
               scales = "free_x",
               space = "free"))

(p5_legend = as_ggplot(get_legend(p5)))


# (p_1 = ggarrange(p1 + theme(legend.position = "none",
#                             plot.margin = unit(c(.5,0,.5,0), 'lines')),
#                  p6 + theme(legend.position = "none",
#                             plot.margin = unit(c(.5,0,.5,0), 'lines')), 
#                  p6_legend,
#                  p4 + theme(legend.position = "none",
#                             plot.margin = unit(c(.5,0,.5,0), 'lines')), 
#                  p5 + theme(legend.position = "none",
#                             plot.margin = unit(c(.5,0,.5,0), 'lines')), 
#                  p5_legend,
#                  ncol = 3,
#                  nrow = 2,
#                  widths = c(1, 1, 1/4, 1, 1,  1/4),
#                 heights = c(1, 1, 1, 1, 1, 1),
#                  align = "hv"))

(p_1 = ggarrange(p6 + theme(#legend.position = "none",
                            plot.margin = unit(c(.5,.5,.5,.5), 'lines')), 
                 #p6_legend,
                 p5 + theme(#legend.position = "none",
                            plot.margin = unit(c(.5,.5,.5,.5), 'lines')), 
                 #p5_legend,
                 ncol = 1,
                 nrow = 2,
                 widths = c(1, 1),
                 # heights = c(1, 1, 1, 1, 1, 1),
                 align = "hv"))

rect = rectGrob(
  x = unit(1.08, "in"),
  y = unit(.97, "npc") - unit(.2, "in"),
  width = unit(.5, "in"),
  height = unit(.5, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

font_s = 21.5

lab_a = textGrob(
  label = "(a)",
  x = unit(1.15, "in"),
  y = unit(.97, "npc") - unit(.3, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

rect2 = rectGrob(
  x = unit(7.25, "in"),
  y = unit(1, "npc") - unit(.2, "in"),
  width = unit(.5, "in"),
  height = unit(.5, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))


lab_b = textGrob(
  label = "(b)",
  x = unit(7.32, "in"),
  y = unit(1, "npc") - unit(.3, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

rect3 = rectGrob(
  x = unit(1.08, "in"),
  y = unit(.47, "npc") - unit(.2, "in"),
  width = unit(.5, "in"),
  height = unit(.5, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

lab_c = textGrob(
  label = "(b)",
  x = unit(1.15, "in"),
  y = unit(.47, "npc") - unit(.3, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))


rect4 = rectGrob(
  x = unit(7.25, "in"),
  y = unit(.5, "npc") - unit(.2, "in"),
  width = unit(.5, "in"),
  height = unit(.5, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

lab_d = textGrob(
  label = "(d)",
  x = unit(7.32, "in"),
  y = unit(.499, "npc") - unit(.3, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))



(p_3 = ggdraw(p_1) +
    draw_grob(rect) +
    #draw_grob(rect2) +
    draw_grob(rect3) +
    # draw_grob(rect4) +
    draw_grob(lab_a) +
    #draw_grob(lab_b) #+
     draw_grob(lab_c)# +
    # draw_grob(lab_d)
  )

  

ggsave(plot = p_3, 
       filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\figure_3.jpeg", 
       device = jpeg,
       width = 14,
       height = 12,
       units = "in",
       dpi = 300,
       bg = "white")



(p3 = vpop_plot %>% 
    #filter(index == "R740") %>% 
    ggplot(aes(x = mean_scattered, y = vpop)) +
    geom_line(aes(group = Date), alpha = .5) +
    geom_point(aes(group = trait, color = timepoint), size = 1.5, show.legend = FALSE) +
    stat_summary(aes(fill = timepoint), fun.y=mean, shape = 21, size = 3, geom ="point") +
    geom_smooth(method = "lm", se = FALSE, color = "black", 
                linewidth = 1, alpha = .3, linetype = 3) +
    scale_color_manual(values = timepoint_cols) +
    scale_fill_manual(values = timepoint_cols) +
    theme_bw(base_size = 20) +
    labs(y = expression(Population~differentiation~(V[pop]~", %")),
         x = "Proportion of scattered irradiance") +
    theme(aspect.ratio = 1,
          legend.position = c(.2, .85),
          legend.background = element_rect(color = "black"),
          legend.title = element_blank()))

(p_6 = ggarrange(p2, p3,
                 ncol = 2, align = "h"))

# ggsave(plot = p_6, 
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\SI\\figure_s1.jpeg", 
#        device = jpeg,
#        width = 15,
#        height = 8,
#        units = "in",
#        dpi = 300,
#        bg = "white")

ggsave(plot = p3, 
       filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\SI\\figure_s1.jpeg", 
       device = jpeg,
       width = 9,
       height = 8,
       units = "in",
       dpi = 300,
       bg = "white")


################################################################################

# (compare = mm_dat %>%
#     filter(!(index %in% c("R444",
#                           "R475",
#                           "R531",
#                           "R560",               
#                           "R650",               
#                           "R668",               
#                           "R705",               
#                           "R717",               
#                           "R740",               
#                           "R842"))) %>% 
#     ungroup() %>% 
#    dplyr::rename("value" = "var_value") %>% 
#     # dplyr::mutate(order = dense_rank(`Date:timepoint`)) %>% 
#     # pivot_longer(`Blk:Rep`:Residual, names_to = "component") %>% 
#     ggplot(aes(x = Date,
#                y = value)) +
#     geom_point(size = 2.5, alpha = .05, aes(color = timepoint)) +
#     theme_bw(base_size = 10) +
#     theme(axis.text.x = #element_text(angle = 90)) +
#             element_blank()) +
#     facet_grid(index ~ timepoint,
#                #ncol = 3,
#                scales = "free"))
# 
# (compare = mm_dat %>%
#     filter(!(index %in% c("R444",
#                           "R475",
#                           "R531",
#                           "R560",               
#                           "R650",               
#                           "R668",               
#                           "R705",               
#                           "R717",               
#                           "R740",               
#                           "R842")) &
#              Date != "2020-03-22") %>% 
#     ungroup() %>% 
#    dplyr::rename("value" = "var_value") %>% 
#     # dplyr::mutate(order = dense_rank(`Date:timepoint`)) %>% 
#     # pivot_longer(`Blk:Rep`:Residual, names_to = "component") %>% 
#     ggplot(aes(x = reorder_within(Prov, value, timepoint),
#                y = value)) +
#     geom_point(size = 1, alpha = .2, aes(color = Date)) +
#     theme_bw(base_size = 10) +
#     theme(axis.text.x = #element_text(angle = 90)) +
#             element_blank()) +
#     facet_grid(index ~ timepoint,
#                #ncol = 3,
#                scales = "free"))
# 
# (vpop_plot_R = components_VI %>%
#     filter((index %in% c("R444",
#                           "R475",
#                           "R531",
#                           "R560",               
#                           "R650",               
#                           "R668",               
#                           "R705",               
#                           "R717",               
#                           "R740",               
#                           "R842"))) %>% 
#     ungroup() %>% 
#     dplyr::mutate(order = dense_rank(`Date:timepoint`)) %>% 
#     pivot_longer(`Blk:Rep`:Residual, names_to = "component") %>% 
#     ggplot(aes(x = reorder(index, vpop),
#                y = value,
#                fill = component)) +
#     geom_bar( stat = "identity") +
#     theme_bw(base_size = 16) +
#     theme(axis.text.x = element_text(angle = 90)))
# 
# (p_VI = vpop_plot_VI %>% 
#   ggplot(aes(x = reorder(index, desc(vpop)), y = vpop, group = trait)) +
#   geom_line(aes(color = trait), linewidth = 1.1) +
#   theme_bw(base_size = 16) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   labs(x = "Vegetation index",
#        y = "Proportion of variance due to population differences (Vpop)") +
#   facet_wrap(. ~ model,
#              ncol = 5, 
#              scales = "free_x"))
# 
# ggsave(plot = p_VI, 
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\indices_vpop.jpeg", 
#        device = jpeg,
#        width = 10,
#        height = 8,
#        units = "in",
#        dpi = 300,
#        bg = "white")
# 
# (p_VI = vpop_plot_VI %>% 
#     ggplot(aes(x = reorder(index, desc(vpop)), y = vpop, group = trait)) +
#     geom_line(aes(color = trait), linewidth = 1.1) +
#     theme_bw(base_size = 16) +
#     theme(axis.text.x = element_text(angle = 90)) +
#     labs(x = "Vegetation index",
#          y = "Proportion of variance due to population differences (Vpop)"))
# 
# (p_VI = vpop_plot_VI %>% 
#     filter(trait %in% c("decline", "greenup")) %>% 
#     ggplot(aes(x = reorder_within(index, by = desc(vpop), within = trait), y = vpop, group = trait)) +
#     geom_point(aes(color = trait)) +
#     theme_bw(base_size = 16) +
#     theme(axis.text.x = element_text(angle = 90)) +
#     labs(x = "Vegetation index",
#          y = "Proportion of variance due to population differences (Vpop)") +
#     facet_wrap(. ~ trait,
#                ncol = 5, 
#                scales = "free_x"))
# 
# ggsave(plot = p_VI, 
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\figure_3.jpeg", 
#        device = jpeg,
#        width = 15,
#        height = 8,
#        units = "in",
#        dpi = 300,
#        bg = "white")

# (table_re_upper = components %>% 
#   bind_rows(vpop) %>% 
#   dplyr::mutate(across(2:7, round, 2)) %>% 
#   flextable::as_flextable())
# 
# flextable::save_as_docx("NDRE" = table_ndre1, "CCI" = table_cci, "REupper" = table_re_upper,
#                         path = "D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\prelim.docx")


################################################################################

sx_pops_clim = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\output\\Sx_CC_Seedlots_Normal_1961_1990_Y_S.csv") %>% 
  dplyr::rename("Prov" = "Seedlot") %>% 
  dplyr::select(-X, -Longitude, -MAR) %>% 
  dplyr::mutate(MSP_log = log(MSP),
         PAS_log = log(PAS))
  

(pops_clim_use = sx_pops_clim %>% 
  dplyr::select("Latitude", 
                "Elevation", 
                
                "bFFP", # spring frost risk
                "eFFP", # fall frost risk
                
                "MCMT", # cold winters
                "MWMT", # hot summers
                
                "MSP_log", # summer moisture
                "PAS_log", # spring moisture / winter protection
                
                "CMD", # drought
                "AHM", # aridity
           ) )

(clim_hist = pops_clim_use %>% 
  rownames_to_column() %>% 
  pivot_longer(Latitude:AHM) %>% 
  dplyr::mutate(name = factor(name, levels = c(
    "Latitude", "Elevation", "bFFP", "eFFP", "MCMT", "MWMT",
    "MSP_log", "PAS_log", "CMD", "AHM"))) %>%  
  ggplot(aes(x = value)) +
  geom_histogram(bins = 25) +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title = element_blank()) +
  facet_wrap(. ~ name, 
             ncol = 5,
             scales = "free"))

(clim_corr = pops_clim_use %>% 
  correlate() %>% 
  pivot_longer(-term) %>% 
    dplyr::mutate(term = factor(term, levels = rev(c(
      "Latitude", "Elevation", "bFFP", "eFFP", "MCMT", "MWMT",
      "MSP_log", "PAS_log", "CMD", "AHM"))),
      name = factor(name, levels = rev(c(
      "Latitude", "Elevation", "bFFP", "eFFP", "MCMT", "MWMT",
      "MSP_log", "PAS_log", "CMD", "AHM")))) %>% 
  dplyr::mutate(across(where(is.numeric), round, 2)) %>% 
    filter(as.integer(term) > as.integer(name)) %>%
  ggplot(aes(x = name, y = term, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), size = 2) +
  scale_fill_gradient2(low = "steelblue4", mid = "white", high = "orchid4") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90),
        aspect.ratio = 1))

#hist(sx_pops_clim$PAS, breaks = 100)
#hist(clim_cor$PAS, breaks = 100)#

(p_clim = ggarrange(clim_corr, clim_hist,
                ncol = 2, align = "h", widths = c(.4, .6)))


ggsave(plot = p_clim, 
       filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\SI\\clim_cor.jpeg", 
       device = jpeg,
       width = 12,
       height = 5,
       units = "in",
       dpi = 300,
       bg = "white")



################################################################################

class_cols =  c("Wildstand" = "#1752A2",
                "Orchard" = "#0A8C35"
)

# blups for spectral indices
# mm_dat = df_spectral_all %>% 
#  dplyr::rename("var_value" = "value") %>% 
#   dplyr::mutate(Blk = factor(Blk),
#          Rep = factor(Rep),
#          Prov = factor(Prov),
#          Class = factor(Class),
#          Season = facdtor(timepoint))

df = mm_dat %>%
  left_join(sx_pops_clim %>%
              dplyr::mutate(Prov = factor(Prov)) %>%
              dplyr::select("Prov", "Latitude",
                            "Elevation",

                            "bFFP", # spring frost risk
                            "eFFP", # fall frost risk

                            "MCMT", # cold winters
                            "MWMT", # hot summers

                            "MSP_log", # summer moisture
                            "PAS_log", # spring moisture / winter protection

                            "CMD", # drought
                            "AHM"),# aridity),
            by = "Prov") %>%
  pivot_longer(Latitude:AHM, names_to = "clim", values_to = "clim_value") %>%
  filter(index == "CCI" & clim == "MWMT")

library(emmeans)
emm_options(pbkrtest.limit = 3000)

blup_mod = function(df) {
  
  # fit a random effects model
  mm = lmer(var_value ~ 
              (1|timepoint:Prov) + 
              (1|timepoint:Rep) +
              (1|timepoint:Rep:Blk) +
              #(1|timepoint:Date) + 
              timepoint,
            data = df, REML = TRUE)

  
  
  # sum = summary(mm)$coefficients %>% 
  #   as_tibble(rownames = "coefficient") %>% 
  #  dplyr::rename("t_value_class" = `t value`,
  #          "p_value_class" = `Pr(>|t|)`) %>% 
  #   dplyr::select(coefficient, t_value_class, p_value_class) %>% 
  #   dplyr::mutate(coefficient = gsub("timepoint", "", coefficient)) %>% 
  #   separate(coefficient, ":", into = c("timepoint", "drop")) %>% 
  #     drop_na(drop) %>% 
  #   dplyr::select(-drop)
  
  vals = df %>% 
    #drop_na(Class) %>% 
    dplyr::select(Prov, timepoint) %>% 
    distinct() %>% 
    # Prov in Class
    left_join(ranef(mm)$`timepoint:Prov` %>% 
                rownames_to_column() %>% 
                separate(rowname, sep = ":", into = c("timepoint", "Prov")) %>% 
               dplyr::rename(BLUP = `(Intercept)`),
              by = c("timepoint", "Prov")) %>% 
    dplyr::mutate(INTERCEPT = fixef(mm)[[1]],
           TIMEPOINT = case_when(timepoint == "early summer" ~ 0,
                               timepoint == "late summer" ~ fixef(mm)[[2]],
                               timepoint == "late winter" ~ fixef(mm)[[3]],
                               timepoint == "mid summer" ~ fixef(mm)[[4]]))

  
  sum = emmeans(mm, pairwise ~ timepoint, lmer.df = "asymp", adjust = "none")$contrasts %>% 
    as_tibble() %>% 
    dplyr::mutate(contrast = case_when(
      contrast == "early summer - late summer" ~ "ES_LS_p",
      contrast == "early summer - late winter" ~ "ES_LW_p",
      contrast == "early summer - mid summer" ~ "ES_MS_p",
      contrast == "late summer - late winter" ~ "LS_LW_p",
      contrast == "late summer - mid summer" ~ "LS_MS_p",
      contrast == "late winter - mid summer" ~ "LW_MS_p")) %>% 
    dplyr::select(contrast, p.value) %>% 
    pivot_wider(names_from = contrast, values_from = p.value) %>% 
    slice(rep(1:n(), each = nrow(vals)))
  
  vals2 = bind_cols(vals, sum)
    
}

df_mm_blups = mm_dat %>% 
  dplyr::select(Obs, Prov, Blk, Rep, index, var_value, timepoint) %>% 
  dplyr::group_by(index, Obs) %>% 
  # drop any trees with incomplete data
  dplyr::mutate(mean_obs = mean(var_value, na.rm = FALSE)) %>% 
  drop_na(mean_obs) %>% 
  dplyr::select(-mean_obs) %>% 
  dplyr::group_by(index, Obs, Prov, Blk, Rep, timepoint) %>% 
  dplyr::summarise(var_value = mean(var_value)) %>% 
  drop_na() %>% 
  dplyr::filter(is.finite(var_value)) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(index, Prov, timepoint) %>% 
  dplyr::add_count() %>% 
  group_by(index) %>% 
  nest() %>% 
  dplyr::mutate(mm = map(data, blup_mod)) %>%
  unnest(c(mm)) %>%
  dplyr::select(-data)

# df_mm_blups = mm_dat %>% 
#   filter(is.finite(var_value)) %>% 
#   drop_na(var_value) %>% 
#   group_by(index) %>% 
#   nest() %>% 
#   dplyr::mutate(mm = map(data, blup_mod)) %>%
#   unnest(c(mm)) %>%
#   dplyr::select(-data)

saveRDS(df_mm_blups, "D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\BLUPS\\Updated\\blups_spectral_index_avg.rds")

df_mm_blups = readRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\BLUPS\\Updated\\blups_spectral_index_avg.rds")


mm_blups_plot = df_mm_blups %>% 
  dplyr::mutate(ESTIMATE = BLUP + INTERCEPT + TIMEPOINT)

mm_blups_change = mm_blups_plot %>% 
  dplyr::select(index:timepoint, ESTIMATE) %>% 
  distinct() %>% 
  #dplyr::mutate(timepoint = str_replace_all(timepoint, " ", "_")) %>% 
  pivot_wider(names_from = timepoint, values_from = ESTIMATE) %>% 
  dplyr::mutate(`green-up` = `mid summer` - `late winter`,
         decline = `late summer` - `early summer`) %>% 
  pivot_longer(`early summer`:decline, names_to = "timepoint", values_to = "timepoint_value") 
  # %>% 
  # left_join(mm_blups_plot %>% 
  #             dplyr::select(index:timepoint) %>% 
  #             distinct() %>% 
  #             #dplyr::mutate(timepoint = str_replace_all(timepoint, " ", "_")) %>% 
  #             pivot_wider(names_from = timepoint, values_from = CLASS) %>% 
  #             dplyr::mutate(greenup = `mid summer` - `late winter`,
  #                    decline = `late summer` - `mid summer`) %>% 
  #             pivot_longer(`mid summer`:decline, names_to = "timepoint", values_to = "class_value"))


mm_blups_clim = mm_blups_change %>%
  left_join(sx_pops_clim %>% 
              dplyr::mutate(Prov = factor(Prov)) %>% 
              dplyr::select("Prov", "Latitude", 
                            "Elevation", 
                            
                            "bFFP", # spring frost risk
                            "eFFP", # fall frost risk
                            
                            "MCMT", # cold winters
                            "MWMT", # hot summers
                            
                            "MSP_log", # summer moisture
                            "PAS_log", # spring moisture / winter protection
                            
                            "CMD", # drought
                            "AHM"),# aridity), 
            by = "Prov") %>% 
  pivot_longer(Latitude:AHM, names_to = "clim", values_to = "clim_value")
  
  
mod_fit = function(df) {
    
    # fit a random effects model
    mm = lm(var_value ~ clim_value, data = df)
    
    }
                
blups_model = mm_blups_clim %>% 
  #filter(index %in% keep_list) %>% 
    group_by(timepoint, clim, index) %>% 
    nest() %>% 
    dplyr::mutate(linear = map(data, ~lm(timepoint_value ~ clim_value, data = .)),
           quadratic = map(data, ~lm(timepoint_value ~ clim_value + I(clim_value^2), data = .))) %>% 
  pivot_longer(cols = c(quadratic, linear), names_to = "type", values_to = "model") %>% 
  dplyr::mutate(clim = factor(clim, levels = c("Latitude", 
                "Elevation", 
                
                "MCMT", # cold winters
                "MWMT", # hot summers
                
                "bFFP", # spring frost risk
                "eFFP", # fall frost risk

                "MSP_log", # summer moisture
                "PAS_log", # spring moisture / winter protection
                
                "CMD", # drought
                "AHM")))
  
blups_pred = blups_model %>% 
    dplyr::mutate(preds = map2(data, model, add_predictions)) %>%
    dplyr::mutate(preds = map(preds, dplyr::select, pred)) %>% 
  unnest(c(data, preds)) %>% 
  dplyr::select(-model)

################################################################################

(VI_list = vpop %>% 
   as_tibble() %>% 
   # left_join(components_all %>% 
   #             dplyr::mutate(Prov_all = `timepoint:Prov`,
   #                           timepoint_all = timepoint) %>% 
   #             dplyr::select(index, Prov_all, timepoint_all, order)) %>% 
   left_join(vpop_all %>%
               dplyr::mutate(Prov_all = `timepoint:Prov`,
                             timepoint_all = timepoint,
                             vpop_all = vpop) %>%
               dplyr::select(index, Prov_all, timepoint_all, vpop_all)) %>%
   filter(type == "VI") %>% 
   group_by(index) %>% 
   dplyr::mutate(vpop_mean = mean(vpop),
                 vpop_max = max(vpop),
                 vpop_min = min(vpop)) %>% 
   ungroup() %>% 
   dplyr::mutate(order = dense_rank(vpop_all)) %>% 
   filter(type == "VI" #& vpop_mean > .2
   ))

# 
# (test = vpop_all %>%
#     dplyr::filter(index %in% VI_list$index) %>% 
#     dplyr::mutate(timepoint_plot = (timepoint / (timepoint + Residual))) %>% 
#     ggplot(aes(x = timepoint_plot, 
#                y = vpop)) +
#     geom_point() +
#     geom_text(aes(label = index)) +
#     theme_bw())
# 

#keep_list = c("CCI", "PRI", "MCARI")

# tree level cor
(mm_dat_VI = mm_dat %>% 
    filter(index %in% VI_list$index) %>% 
    dplyr::select(Obs, index, Date, timepoint, var_value)  %>% 
    distinct() %>% 
    filter(is.finite(var_value)) %>% 
    drop_na(var_value) %>% 
    pivot_wider(names_from = index, values_from = var_value))

(mm_dat_VI = mm_blups_clim %>%
    # left_join(VI_list %>%
    #             dplyr::select(index, order)) %>%
    # filter(index %in% VI_list$index
    #        # &
    #        # timepoint %in% c("mid summer", "late summer", "late winter")
    #        ) %>%
    dplyr::select(Prov, index, timepoint, timepoint_value)  %>%
    distinct() %>%
    mutate(ESTIMATE = timepoint_value) %>%
    dplyr::select(-timepoint_value) %>%
    filter(is.finite(ESTIMATE)) %>%
    drop_na(ESTIMATE) %>%
    pivot_wider(names_from = index, values_from = ESTIMATE))

mm_dat_VI %>% 
  ggplot(aes(x = mDatt, y = NDRE3, color = factor(timepoint))) +
  geom_point(size = 4) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  facet_wrap(. ~ timepoint,
             scales = "free",
             nrow = 1)


# df_cluster %>% 
#   ggplot(aes(x = `R668_mid summer`, y = `R668_late winter`)) +
#   geom_point() +
#   theme_bw() +
#   theme(aspect.ratio = 1)

# (mm_dat_VI = mm_dat %>%
#     # left_join(VI_list %>%
#     #             dplyr::select(index, order)) %>%
#     dplyr::filter(index %in% VI_list$index) %>%
#     dplyr::distinct(Prov, index, var_value, Date)  %>%
#     #distinct(Prov, Date, var_value) %>%
#     dplyr::filter(is.finite(var_value)) %>%
#     tidyr::drop_na(var_value) %>%
#     tidyr::pivot_wider(id_cols = c(Prov, Date), id_expand = TRUE, names_from = index, values_from = var_value))

cor_all = mm_dat_VI %>% 
  #dplyr::select(-Prov) %>% 
  correlate() %>% 
  pivot_longer(-term) %>% 
  dplyr::mutate(across(where(is.numeric), round, 3)) %>% 
  left_join(VI_list %>% 
              dplyr::rename("order_term" = "order") %>% 
              dplyr::select(index, order_term),
            by = c("term" = "index")) %>% 
  left_join(VI_list%>% 
              dplyr::rename("order_name" = "order") %>% 
              dplyr::select(index, order_name),
            by = c("name" = "index"))

cor_early_summer = mm_dat_VI %>%
  filter(timepoint == "early summer") %>%
  dplyr::select(-Prov) %>%
  correlate() %>%
  pivot_longer(-term) %>%
  dplyr::mutate(across(where(is.numeric), round, 3)) %>%
  left_join(VI_list %>%
              dplyr::rename("order_term" = "order") %>%
              dplyr::select(index, order_term),
            by = c("term" = "index")) %>%
  left_join(VI_list%>%
              dplyr::rename("order_name" = "order") %>%
              dplyr::select(index, order_name),
            by = c("name" = "index")) %>%
  dplyr::mutate(timepoint = "early summer")

cor_late_winter = mm_dat_VI %>%
  filter(timepoint == "late winter") %>%
  dplyr::select(-Prov) %>%
  correlate() %>%
  pivot_longer(-term) %>%
  dplyr::mutate(across(where(is.numeric), round, 3)) %>%
  left_join(VI_list %>%
              dplyr::rename("order_term" = "order") %>%
              dplyr::select(index, order_term),
            by = c("term" = "index")) %>%
  left_join(VI_list%>%
              dplyr::rename("order_name" = "order") %>%
              dplyr::select(index, order_name),
            by = c("name" = "index")) %>%
  dplyr::mutate(timepoint = "late winter")

cor_late_summer = mm_dat_VI %>%
  filter(timepoint == "late summer") %>%
  dplyr::select(-Prov) %>%
  correlate() %>%
  pivot_longer(-term) %>%
  dplyr::mutate(across(where(is.numeric), round, 3)) %>%
  left_join(VI_list %>%
              dplyr::rename("order_term" = "order") %>%
              dplyr::select(index, order_term),
            by = c("term" = "index")) %>%
  left_join(VI_list%>%
              dplyr::rename("order_name" = "order") %>%
              dplyr::select(index, order_name),
            by = c("name" = "index")) %>%
  dplyr::mutate(timepoint = "late summer")

cor_mid_summer = mm_dat_VI %>%
  filter(timepoint == "mid summer") %>%
  dplyr::select(-Prov) %>%
  correlate() %>%
  pivot_longer(-term) %>%
  dplyr::mutate(across(where(is.numeric), round, 3)) %>%
  left_join(VI_list %>%
              dplyr::rename("order_term" = "order") %>%
              dplyr::select(index, order_term),
            by = c("term" = "index")) %>%
  left_join(VI_list%>%
              dplyr::rename("order_name" = "order") %>%
              dplyr::select(index, order_name),
            by = c("name" = "index")) %>%
  dplyr::mutate(timepoint = "mid summer")

cor_greenup = mm_dat_VI %>%
  filter(timepoint == "green-up") %>%
  dplyr::select(-Prov) %>%
  correlate() %>%
  pivot_longer(-term) %>%
  dplyr::mutate(across(where(is.numeric), round, 3)) %>%
  left_join(VI_list %>%
              dplyr::rename("order_term" = "order") %>%
              dplyr::select(index, order_term),
            by = c("term" = "index")) %>%
  left_join(VI_list%>%
              dplyr::rename("order_name" = "order") %>%
              dplyr::select(index, order_name),
            by = c("name" = "index")) %>%
  dplyr::mutate(timepoint = "green-up")

cor_decline = mm_dat_VI %>%
  filter(timepoint == "decline") %>%
  dplyr::select(-Prov) %>%
  correlate() %>%
  pivot_longer(-term) %>%
  dplyr::mutate(across(where(is.numeric), round, 3)) %>%
  left_join(VI_list %>%
              dplyr::rename("order_term" = "order") %>%
              dplyr::select(index, order_term),
            by = c("term" = "index")) %>%
  left_join(VI_list%>%
              dplyr::rename("order_name" = "order") %>%
              dplyr::select(index, order_name),
            by = c("name" = "index")) %>%
  dplyr::mutate(timepoint = "decline")

cor_all = bind_rows(cor_early_summer, cor_mid_summer, cor_late_summer, cor_late_winter, cor_decline, cor_greenup)

(cor_all %>% 
    #filter(as.integer(term) > as.integer(name)) %>%
    ggplot(aes(x = reorder(name, order_name), y = reorder(term, order_term), fill = value)) +
    geom_tile() +
    geom_text(aes(label = value), size = 2) +
    scale_fill_gradient2(low = "steelblue4", mid = "white", high = "orchid4") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90),
          aspect.ratio = 1,
          strip.background = element_blank()) +
    facet_grid(. ~ timepoint))

# # correlation threshold of .9
# keep_list = factor(VI_list$index) %>%
#   factor() %>%
#   levels() %>%
#   .[!(. %in% c("NIRv", "EVI", "BG1"
# 
#                # "MRESR", "SIPI2", "STVI", "CI_green", "NDRE1"
#                #"MARI", "MCARI", "MRESR", "CI_green", "Gcc2"
#   ))]
# 


keep_list = 
  # factor(VI_list$index) %>%
  # factor() %>%
  # levels() %>%
  # .[. %in% 
      c(
    "CCI", "PRI", 
    "GCC", "ARI",
     "mDatt", "RE_upper", "NDVI"#, "EWI9"
    )
#]

(VI_cor = cor_all %>% 
    filter(term %in% keep_list & name %in% keep_list) %>%
    dplyr::mutate(timepoint2 = factor(timepoint, c("early summer", "mid summer", "late summer", "late winter",
                                                    "green-up", "decline")),
                  term = recode(term,
                                  "RE_upper" = "RE['upper']",
                                  "NDRE3" = "NDRE['740']",
                                  "EWI9" = "EWI['9']",
                                "Datt" = "mDatt"),
                  name = recode(name,
                                "RE_upper" = "RE['upper']",
                                "NDRE3" = "NDRE['740']",
                                "EWI9" = "EWI['9']",
                                "Datt" = "mDatt"),
                  term = factor(term,
                                  levels = c("CCI", "PRI", "GCC", "BCC", "ARI", "EWI['9']",
                                             "mDatt", "RE['upper']", "NDVI", "NDRE['740']")),
                  name = factor(name,
                                levels = c("CCI", "PRI", "GCC", "BCC", "ARI", "EWI['9']",
                                           "mDatt", "RE['upper']", "NDVI", "NDRE['740']"))) %>% 
    ggplot(aes(x = name, y = term, fill = value)) +
    geom_tile() +
    geom_text(aes(label = gsub("0\\.", "\\.", round(value, 2))), size = 4) +
    scale_fill_gradient2(low = "steelblue4", mid = "white", high = "red4") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 17) +
    coord_cartesian(expand = FALSE) +
    theme(axis.text.x = element_text(angle = 90, size = 11),
          axis.text.y = element_text(size = 11),
          axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          aspect.ratio = 1,
          legend.position = "none",
          legend.direction = "horizontal",
          legend.key.width = unit(.4, "in"),) +
    scale_x_discrete(labels=scales::parse_format()) +
    scale_y_discrete(labels=scales::parse_format()) +
    guides(fill = guide_colourbar(title.position="top", 
                                  title = "correlation",
                                  title.hjust = .5,
                                  )) +
    facet_wrap2(. ~ timepoint2, axes = "all"))


(test = vpop_all %>%
    dplyr::filter(index %in% keep_list) %>%
    dplyr::mutate(timepoint_plot = (timepoint / (timepoint + Residual))) %>%
    ggplot(aes(x = timepoint_plot,
               y = vpop)) +
    geom_point() +
    geom_text(aes(label = index)) +
    theme_bw())

ggsave(plot = VI_cor, 
       filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\SI\\VI_cor.jpeg", 
       device = jpeg,
       width = 12,
       height = 8.5,
       units = "in",
       dpi = 300,
       bg = "white")


################################################################################
# perform post-hoc test comparing means of timepoints

# class_cols =  c("Wildstand" = "#1752A2",
#                 "Orchard" = "#0A8C35")


blups_class = df_mm_blups %>% 
  distinct(index, Prov, timepoint, BLUP, INTERCEPT, TIMEPOINT, 
           ES_LS_p, ES_LW_p, ES_MS_p, LS_LW_p, LS_MS_p, LW_MS_p)  %>% 
  left_join(components_all %>% 
              dplyr::select(index, order)) %>% 
  dplyr::mutate(ESTIMATE = BLUP + INTERCEPT + TIMEPOINT,
                ESTIMATE = if_else(index == "RE_upper", ESTIMATE * 100, ESTIMATE),
                index2 = recode(index,
                                "RE_upper" = "RE['upper']",
                                "NDRE3" = "NDRE['740']",
                                "EWI9" = "EWI['9']",
                                "Datt" = "mDatt"),
         index2 = factor(index2,
                        levels = c("CCI", "PRI", "GCC", "BCC", "ARI", "EWI['9']",
                                   "mDatt", "RE['upper']", "NDVI", "NDRE['740']")))

blups_class_values = blups_class %>% 
  group_by(index) %>% 
  dplyr::mutate(range_comp = abs(max(ESTIMATE) - min(ESTIMATE))) %>% 
  group_by(index, timepoint) %>% 
  dplyr::mutate(max_val = max(ESTIMATE)) %>% 
  distinct(index, index2, Prov, order, timepoint, ES_LS_p, ES_LW_p, ES_MS_p, LS_LW_p, LS_MS_p, LW_MS_p, max_val, range_comp)  %>% 
  left_join(components_all %>% 
              dplyr::select(index, order))

blups_p = blups_class_values %>% 
  pivot_wider(names_from = timepoint, values_from = max_val) %>% 
  dplyr::select(-Prov) %>% 
  distinct() %>% 
  pivot_longer(ES_LS_p:LW_MS_p, names_to = "comparison", values_to = "p.comp") %>% 
  dplyr::group_by(index) %>% 
  dplyr::mutate(max_comp = case_when(comparison == "ES_LS_p" ~ max(c(`early summer`, `mid summer`, `late summer`)),
                                     comparison == "ES_LW_p" ~ max(c(`early summer`, `mid summer`, `late summer`, `late winter`)),
                                     comparison == "ES_MS_p" ~ max(c(`early summer`, `mid summer`)),
                                     comparison == "LS_LW_p" ~ max(c(`late summer`, `late winter`)),
                                     comparison == "LS_MS_p" ~ max(c(`late summer`, `mid summer`)),
                                     comparison == "LW_MS_p" ~ max(c(`early summer`, `mid summer`, `late summer`, `late winter`))),
                max_plot = case_when(comparison == "ES_LS_p" ~ max_comp + (range_comp * .12),
                                     comparison == "ES_LW_p" ~ max_comp + (range_comp * .24),
                                     comparison == "ES_MS_p" ~ max_comp + (range_comp * .06),
                                     comparison == "LS_LW_p" ~ max_comp + (range_comp * .06),
                                     comparison == "LS_MS_p" ~ max_comp + (range_comp * .06),
                                     comparison == "LW_MS_p" ~ max_comp + (range_comp * .18)),
                x_loc = case_when(comparison == "ES_LS_p" ~ 2,
                                  comparison == "ES_LW_p" ~ 2.5,
                                  comparison == "ES_MS_p" ~ 1.5,
                                  comparison == "LS_LW_p" ~ 3.5,
                                  comparison == "LS_MS_p" ~ 2.5,
                                  comparison == "LW_MS_p" ~ 3),
                x_start = case_when(comparison == "ES_LS_p" ~ 1.1,
                                    comparison == "ES_LW_p" ~ 1.1,
                                    comparison == "ES_MS_p" ~ 1.1,
                                    comparison == "LS_LW_p" ~ 3.1,
                                    comparison == "LS_MS_p" ~ 2.1,
                                    comparison == "LW_MS_p" ~ 2.1),
                x_end = case_when(comparison == "ES_LS_p" ~ 2.9,
                                  comparison == "ES_LW_p" ~ 3.9,
                                  comparison == "ES_MS_p" ~ 1.9,
                                  comparison == "LS_LW_p" ~ 3.9,
                                  comparison == "LS_MS_p" ~ 2.9,
                                  comparison == "LW_MS_p" ~ 3.9)) %>% 
  dplyr::mutate(p_plot = case_when(p.comp < 0.001 ~ '***',
                                   p.comp >=  0.001 & p.comp < 0.01 ~ '**',
                                   p.comp >=  0.01 & p.comp < 0.05 ~ '*',
                                   p.comp >=  0.05 & p.comp < 0.1 ~ '•',
                                   p.comp >  0.1 ~ 'NS'
  )) 
# %>% 
#   filter(p.comp <= .1)


R_list = c("R444",
           "R475",
           "R531",
           "R560",
           "R650",
           "R668",
           "R705",
           "R717",
           "R740",
           "R842")

blups_p_VI = filter(blups_p, index %in% keep_list)
blups_p_R = filter(blups_p, index %in% R_list)


(plot_class = blups_class %>%  
    filter(index %in% keep_list) %>% 
    #dplyr::mutate(Class = factor(Class, levels = c("Wildstand", "Orchard"))) %>% 
    ggplot() +
    geom_boxplot(outlier.alpha = 0, width = .5, alpha = .7,
                 position = position_dodge(width = .75),
                 linewidth = .3,
                 aes(fill = timepoint,
                     x = factor(timepoint, levels = 
                                  c("early summer", "mid summer", "late summer", "late winter")),
                     y = ESTIMATE)) +
    geom_point(position = position_dodge(width=0.75), 
               shape = 21, stroke = NA, size = 3, alpha = .3,
               aes(group=timepoint, fill = timepoint,
                   x = factor(timepoint, levels = 
                                c("early summer", "mid summer", "late summer", "late winter")),
                   y = ESTIMATE)) +
    geom_text(data = filter(blups_p_VI, p.comp < .1), aes(x = x_loc, y = max_plot, label = p_plot),
              size = 5, vjust = .2) +
    geom_text(data = filter(blups_p_VI, p_plot == "NS"), aes(x = x_loc, y = max_plot, label = p_plot),
              size = 2.8, vjust = -.4) +
    geom_linerange(data = blups_p_VI, aes(xmin = x_start, xmax = x_end, y = max_plot),
                   color = "black") +
    # geom_text(data = blups_p_VI, aes(x = x_loc, y = max_comp),
    #           label = "—",
    #           size = 6,
    #           vjust = -.2) +
    theme_bw(base_size = 15) +
    #scale_y_discrete(limits = rev(levels(blups_r2_plot$clim))) +
    scale_fill_manual(values = timepoint_cols) +
    scale_y_continuous(expand = expansion(mult = c(.09, .09))) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, margin = margin(4, 0, -8, 0)),
          panel.grid.major.x = element_blank(),
          axis.title = element_blank(),
          strip.background = element_blank(),
          axis.ticks.x = element_blank(),
          aspect.ratio = 1.3,
          legend.position = "none",
          strip.placement = "outside",
          axis.text.y = element_text(size = 9.5)) +
    #coord_cartesian(expand = FALSE) +
    facet_wrap2(. ~ index2,
                labeller = label_parsed,
                nrow = 2,
                strip.position = "left",
                scales = "free_y") +
  facetted_pos_scales(
    y = list(
      index2 == "CCI" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                   breaks = c(.05, .15, .25)),
      index2 == "GCC" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                           breaks = seq(.42, .5, .04)),
      index2 == "PRI" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                           breaks = c(-.03, -.09, -.15)),
      index2 == "BCC" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                           breaks = c(.24, .27, .3)),
      index2 == "ARI" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                           breaks = c(3,5,7),
                                           labels = c("3.0", "5.0", "7.0")),
      index2 == "EWI['9']" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                           breaks = c(-.52, -.60, -.68)),
      index2 == "mDatt" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                breaks = c(.62, .70, .78)),
      index2 == "RE['upper']" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                breaks = c(50, 65, 80)),
      index2 == "NDVI" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                   breaks = c(.82, .85, .88)),
      index2 == "NDRE['740']" ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                                   breaks = c(.12, .15, .18))
    )
  ))


ggplot2::ggsave(plot = plot_class,
                filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\figure_4.jpeg",
                device = jpeg,
                width = 14,
                height = 8.5,
                units = "in",
                dpi = 300,
                bg = "white")

(plot_class_band = blups_class %>% 
    filter(index %in% R_list) %>% 
    ggplot() +
    geom_boxplot(outlier.alpha = 0, width = .5, alpha = .7,
                 position = position_dodge(width = .75),
                 linewidth = .3,
                 aes(fill = timepoint,
                     x = factor(timepoint, levels = 
                                  c("early summer", "mid summer", "late summer", "late winter")),
                     y = ESTIMATE)) +
    geom_point(position = position_dodge(width=0.75), 
               shape = 21, stroke = NA, size = 3, alpha = .3,
               aes(group=timepoint, fill = timepoint,
                   x = factor(timepoint, levels = 
                                c("early summer", "mid summer", "late summer", "late winter")),
                   y = ESTIMATE)) +
    geom_text(data = filter(blups_p_R, p.comp < .1), aes(x = x_loc, y = max_plot, label = p_plot),
              size = 5, vjust = .2) +
    geom_text(data = filter(blups_p_R, p_plot == "NS"), aes(x = x_loc, y = max_plot, label = p_plot),
              size = 2.8, vjust = -.4) +
    geom_linerange(data = blups_p_R, aes(xmin = x_start, xmax = x_end, y = max_plot),
                   color = "black") +
    # geom_text(data = blups_p_VI, aes(x = x_loc, y = max_comp),
    #           label = "—",
    #           size = 6,
    #           vjust = -.2) +
    theme_bw(base_size = 15) +
    #scale_y_discrete(limits = rev(levels(blups_r2_plot$clim))) +
    scale_fill_manual(values = timepoint_cols) +
    scale_y_continuous(expand = expansion(mult = c(.09, .09))) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, margin = margin(4, 0, -8, 0)),
          panel.grid.major.x = element_blank(),
          axis.title = element_blank(),
          strip.background = element_blank(),
          axis.ticks.x = element_blank(),
          aspect.ratio = 1.3,
          legend.position = "none",
          strip.placement = "outside",
          axis.text.y = element_text(size = 9.5)) +
    #coord_cartesian(expand = FALSE) +
    facet_wrap2(. ~ index,
                nrow = 2,
                strip.position = "left",
                scales = "free_y") +
  facetted_pos_scales(
    y = list(
      !(index %in% c("R717", "R740", "R842")) ~ scale_y_continuous(expand = expansion(mult = c(.09, .09)),
                                           breaks = seq(0, 1, .01)),
      index == "R717" ~ scale_y_continuous(expand = expansion(mult = c(.1, .1)),
                                           breaks = seq(0, 1, .02)),
      index == "R742" ~ scale_y_continuous(expand = expansion(mult = c(.1, .1)),
                                           breaks = seq(0, 1, .05)),
      index == "R842" ~ scale_y_continuous(expand = expansion(mult = c(.1, .1)),
                                           breaks = seq(0, 1, .1))
    )
  ))


ggplot2::ggsave(plot = plot_class_band,
                filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\SI\\bands_boxplot.jpeg",
                device = jpeg,
                width = 16,
                height = 8.5,
                units = "in",
                dpi = 300,
                bg = "white")

(plot_band = blups_class %>% 
    filter(index %in% R_list) %>% 
    dplyr::mutate(#Class = factor(Class, levels = c("Wildstand", "Orchard")),
                  EST = INTERCEPT + TIMEPOINT) %>% 
    distinct(index, EST, timepoint) %>% 
    ggplot(aes(color = factor(timepoint, levels = 
                                c("early summer", "mid summer", "late summer", "late winter")),
               fill = factor(timepoint, levels = 
                               c("early summer", "mid summer", "late summer", "late winter")),
               x = index,
               y = EST
               #group = factor(Class, levels = c("Wildstand", "Orchard"))
    )) +
    geom_line(aes(group = timepoint), linewidth = 1) +
    geom_point(size = 5, alpha = .7) +
    theme_bw(base_size = 22) +
    labs(x = "Spectral band", y = "Reflectance") +
    #scale_y_discrete(limits = rev(levels(blups_r2_plot$clim))) +
    scale_fill_manual(values = timepoint_cols) +
    scale_color_manual(values = timepoint_cols) +
    #scale_y_continuous(expand = expansion(mult = c(.08, .18))) +
    theme(aspect.ratio = 1,
          legend.position = c(.24, .79),
          legend.text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = .5),
          axis.title.x = element_blank(),
          legend.background = element_rect(color = "black", linewidth = .5),
          legend.title = element_blank()))


ggplot2::ggsave(plot = plot_band,
                filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\SI\\bands_timepoint.jpeg",
                device = jpeg,
                width = 10,
                height = 10,
                units = "in",
                dpi = 300,
                bg = "white")



################################################################################


blups_r2 = blups_model %>%
    dplyr::mutate(glance = map(model, broom::glance),
           tidy = map(model, broom::tidy),
           slope = tidy %>% map_dbl(function(x) x$estimate[2]),
           ) %>% 
    unnest(c(data, glance, slope)) %>% 
  dplyr::select(-model, -tidy) %>% 
  dplyr::select(index:clim, type, adj.r.squared:slope) %>% 
  distinct() %>% 
  dplyr::group_by(clim, timepoint) %>% 
  nest() %>% 
  dplyr::mutate(p.adjusted = map(data, function(df) p.adjust(df$p.value, method = "bonferroni"))) %>% 
  #dplyr::mutate(q.value = map(data, function(x) qvalue::qvalue(x$p.value)$qvalues)) %>% 
  unnest(c(data, p.adjusted))

# 2 model fits * 10 clim indices * 6 timepoints
p_threshold = .005 / (2 * 10 * 6) 


blups_r2_plot = blups_r2 %>% 
  left_join(VI_list %>% 
              dplyr::select(index, order)) %>% 
  # group_by(index, clim, type) %>% 
  # dplyr::mutate(max_r2 = max(adj.r.squared),
  #        mean_r2 = mean(adj.r.squared),
  #        rel_r2 = adj.r.squared - mean_r2) %>% 
  group_by(index, clim, timepoint) %>% 
  dplyr::mutate(slope_lin = if_else(type == "linear", slope, NA_real_),
                r2_lin = if_else(type == "linear", adj.r.squared, NA_real_),
                slope_lin_all = mean(slope_lin, na.rm = TRUE),
                r2_lin_all = mean(r2_lin, na.rm = TRUE)) %>% 
  group_by(index, clim, type, timepoint) %>% 
  dplyr::mutate(r2_sign = if_else(slope_lin_all > 0, r2_lin_all, r2_lin_all * -1)) %>% 
  distinct() %>% 
  dplyr::group_by(timepoint, clim, index) %>% 
  # dplyr::mutate(count = if_else(p.value < p_threshold, 1, 0),
  #               sum = sum(count),
  #               # keep quadratic if only quadratic is signif at alpha level
  #               keep_flag = if_else((type == "linear" & sum != 1) |
  #                                     (type != "linear" & sum == 1),
  #                                   1, 0))
  dplyr::mutate(AIC_min = min(AIC),
                count = if_else(AIC < (AIC_min + 10), 1, 0),
                sum = sum(count),
                # keep quadratic if only quadratic is signif at alpha level
                keep_flag = if_else((type == "linear" & sum != 1) |
                                      (type != "linear" & sum == 1),
                                    1, 0)) %>% 
  ungroup() %>% 
  mutate(max_r2 = max(adj.r.squared),
         alpha = adj.r.squared / max_r2,
         alpha2 = 1 / (1 + exp(-10*(alpha - .5))) )




# (test = blups_r2 %>%
#   dplyr::filter(type == "quadratic") %>%
#   dplyr::group_by(index) %>%
#   dplyr::mutate(max_r2 = max(adj.r.squared)) %>%
#   dplyr::filter(adj.r.squared == max_r2) %>%
#   dplyr::full_join(VI_list, by = "index") %>%
#   ggplot(aes(x = adj.r.squared, y = vpop_all)) +
#   geom_point() +
#   geom_label(aes(label = index)) +
#   theme_bw())

scale_this = function(x){
  (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
  #(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

scale_this_robust = function(x){
  (x - median(x, na.rm=TRUE)) / stats::IQR(x, na.rm=TRUE)
}

scale_this_new = function(x){
  (x - median(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
}

mm_blups_change %>%
  filter(index == "NDVI" & timepoint == "late summer") %>%
  mutate(new = scale_this_new(timepoint_value),
         this = scale_this(timepoint_value)) %>%
  ggplot() +
  geom_histogram(aes(x = new), bins = 20, color = "black")

blups_summary = blups_pred %>% 
  left_join(blups_r2_plot) %>% 
  group_by(index, clim, type, timepoint) %>% 
  dplyr::mutate(pred_new = scale_this(pred),
         timepoint_value_new = scale_this(timepoint_value))


################################################################################

R_list = c("R444",
           "R475",
            "R531",
            "R560",
            "R650",
            "R668",
            "R705",
            "R717",
            "R740",
            "R842")

R_cols = c("R444" = "royalblue4",
            "R475" = "steelblue",
            "R531" = "springgreen3",
            "R560" = "forestgreen",
            "R650" = "firebrick2",
            "R668" = "red3",
            "R705" = "indianred",
            "R717" = "lightpink3",
            "R740" = "pink4",
            "R842" = "thistle4")

 

# 
# (blups_summary %>% 
#     dplyr::mutate(timepoint = factor(timepoint, levels = 
#                                 c("mid_summer", "late_summer", "late_winter", "decline", "greenup"))) %>% 
#     filter(index %in% R_list
#            #& clim == "MWMT"
#            & !(timepoint %in% c("decline", "greenup"))) %>% 
#     #filter(index == "Gcc1" & type == "quadratic"
#     #& clim %in% c("MWMT", "AHM", "bFFP")
#     #) %>% 
#     ggplot(aes(x = clim_value, y = pred)) +
#     #geom_point(aes(y = timepoint_value, color = index), size = 2, alpha = .05) +
#     geom_smooth(aes(color = timepoint, linewidth = adj.r.squared), se = FALSE) +
#     #scale_color_manual(values = timepoint_cols) +
#     scale_alpha_continuous() +
#     theme_bw() +
#     theme(aspect.ratio = 1) +
#     # scale_y_continuous(expand = expansion(mult = c(.08, .24))) +
#     # scale_x_continuous(expand = expansion(mult = c(.08, .24))) +
#     facet_grid2(index ~ clim, scales = "free") +
#     facetted_pos_scales(
#       x = list(
#         clim == "Latitude" ~ scale_x_reverse(),
#         clim == "Elevation" ~ scale_x_reverse(),
#         clim == "eFFP" ~ scale_x_reverse(),
#         clim == "CMD" ~ scale_x_reverse(),
#         clim == "AHM" ~ scale_x_reverse()
#       )))

# 
# (band_scaled_summary = blups_summary %>% 
#     dplyr::mutate(timepoint = factor(timepoint, levels = 
#                                 c("mid_summer", "late_summer", "late_winter", "decline", "greenup"))) %>% 
#     filter(index %in% R_list
#            #& clim == "MWMT"
#            #& !(timepoint %in% c("decline", "greenup"))
#            ) %>% 
#     #filter(index == "Gcc1" & type == "quadratic"
#     #& clim %in% c("MWMT", "AHM", "bFFP")
#     #) %>% 
#     ggplot(aes(x = clim_value, y = pred_new, 
#                color = index, alpha = adj.r.squared, group = adj.r.squared)) +
#     geom_line(stat = "smooth", linewidth = 1.1) +
#     theme_bw(base_size = 12) +
#     theme(aspect.ratio = 1,
#           panel.grid = element_blank(),
#           axis.text = element_text(size = 7)) +
#     scale_alpha_continuous(range = c(0, 1)) +
#     scale_color_manual(values = R_cols) +
#     theme_bw() +
#     theme(aspect.ratio = 1) +
#     # scale_y_continuous(expand = expansion(mult = c(.08, .24))) +
#     # scale_x_continuous(expand = expansion(mult = c(.08, .24))) +
#     facet_grid(timepoint ~ clim, scales = "free") +
#     facetted_pos_scales(
#       x = list(
#         clim == "Latitude" ~ scale_x_reverse(),
#         clim == "Elevation" ~ scale_x_reverse(),
#         clim == "eFFP" ~ scale_x_reverse(),
#         clim == "CMD" ~ scale_x_reverse(),
#         clim == "AHM" ~ scale_x_reverse()
#       )))
# 
# ggsave(plot = band_scaled_summary,
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\band_scaled_summary.jpeg",
#        device = jpeg,
#        width = 15,
#        height = 8,
#        units = "in",
#        dpi = 300,
#        bg = "white")

dat_plot = blups_summary %>% 
  dplyr::mutate(timepoint = factor(timepoint, levels = 
                                     c("early summer", "mid summer", "late summer", "late winter", "decline", "green-up"))) %>% 
  dplyr::filter(index %in% R_list
         & clim == "Elevation"
  )  %>% 
  dplyr::group_by(index, clim, timepoint) %>% 
  dplyr::mutate(max_val = max(timepoint_value),
                max_r2 = max(adj.r.squared),
                r_lab = round(adj.r.squared, 2))


(r2_band_plot = dat_plot %>% 
    ggplot(aes(x = clim_value, y = pred, alpha = adj.r.squared, group = adj.r.squared)) +
    geom_point(data = filter(dat_plot, type == "linear"), 
               aes(y = timepoint_value, color = r2_sign), size = 2, alpha = .6) +
    geom_line(data = filter(dat_plot, p.value < p_threshold
                            & keep_flag == 1),
              stat = "smooth", linewidth = .7, alpha = .8) +
    geom_text(data = filter(dat_plot, timepoint_value == max_val 
                            & keep_flag == 1
                            &  p.value < p_threshold
                            & r2_sign > 0), 
              aes(x = -Inf, y = max_val, 
                  label = paste0(expression(R^2), "==", r_lab)),
              size = 3,
              hjust = -.15,
              vjust = -.25,
              alpha = 1,
              parse = TRUE) +
    geom_text(data = filter(dat_plot, timepoint_value == max_val 
                            & keep_flag == 1
                            & p.value < p_threshold
                            & r2_sign < 0), 
              aes(x = Inf, y = max_val, 
                  label = paste0(expression(R^2), "==", r_lab)),
              size = 3,
              hjust = 1.15,
              vjust = -.25,
              alpha = 1,
              parse = TRUE) +
    scale_color_gradient2(low = "steelblue4", mid = "grey70", high = "red4",
                          limits = c(-.38, .38),
                          breaks = seq(-.4, .4, .2)) +
    theme_bw(base_size = 12) +
    labs(y = "Reflectance", x = "Elevation (m)") +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          strip.text = element_text(size = 11.5),
          axis.text = element_text(size = 8),
          axis.title.y = element_text(face = "bold", margin = margin(8, 10, 8, -2)),
          axis.title.x = element_text(face = "bold", margin = margin(10, -2, 8, 8)),
          legend.position = "none",
          strip.background = element_blank(),
          strip.placement = "outside") +
    scale_alpha_continuous(range = c(.1, 1)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01),
                       expand = expansion(mult = c(.12, .26))) +
    scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                       breaks = seq(0, 3000, 1500)) +
    facet_grid2(timepoint ~ index, scales = "free", independent = "y",
                switch = "y"))


ggsave(plot = r2_band_plot,
       filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\SI\\_SI_r2_band_Elevation.jpeg",
       device = jpeg,
       width = 15,
       height = 6.5,
       units = "in",
       dpi = 300,
       bg = "white")


################################################################################

# summary figure 
dat_plot = blups_summary %>% 
  dplyr::mutate(timepoint = factor(timepoint, levels = 
                              c("early summer", "mid summer", "late summer", "late winter", "decline", "green-up"))) %>% 
  filter(index == "NDVI"
  ) %>% 
  group_by(index, clim, timepoint) %>% 
  dplyr::mutate(max_val = max(timepoint_value),
         max_r2 = max(adj.r.squared),
         r_lab = round(adj.r.squared, 2)) %>% 
  dplyr::mutate(clim = factor(clim,
                levels = c( "Latitude", "Elevation", 
                            "MWMT", "MCMT",  
                            "bFFP", "eFFP",
                            "MSP_log", "PAS_log", 
                            "CMD", "AHM"))) %>% 
  distinct() %>% 
  dplyr::mutate(clim = dplyr::recode(clim, 
                                     Latitude = "Latitude (°N)",
                                     Elevation = "Elevation (m)",
                                     MWMT = "MWMT (°C)",
                                     MCMT = "MCMT (°C)",
                                     MSP_log = "MSP (mm)",
                                     PAS_log = "PAS (mm)"))


(r2_index_plot = dat_plot %>% 
    ggplot(aes(x = clim_value, y = pred, alpha = adj.r.squared, group = adj.r.squared)) +
    geom_point(data = filter(dat_plot, keep_flag == 1),#filter(dat_plot, type == "linear"), 
               aes(y = timepoint_value, color = adj.r.squared), size = 2, alpha = .6) +
    geom_line(data = filter(dat_plot, p.value < p_threshold
                            & keep_flag == 1),
              stat = "smooth", linewidth = .7, alpha = .8) +
    geom_text(data = filter(dat_plot, timepoint_value == max_val 
                            & keep_flag == 1
                            &  p.value < p_threshold
                            & r2_sign > 0), 
              aes(x = -Inf, y = max_val, 
                  label = paste0(expression(R^2), "==", r_lab)),
              size = 3.5,
              hjust = -.15,
              vjust = -.25,
              alpha = 1,
              parse = TRUE) +
    geom_text(data = filter(dat_plot, timepoint_value == max_val 
                            & keep_flag == 1
                            &  p.value < p_threshold
                            & r2_sign < 0), 
              aes(x = Inf, y = max_val, 
                  label = paste0(expression(R^2), "==", r_lab)),
              size = 3.5,
              hjust = 1.15,
              vjust = -.25,
              alpha = 1,
              parse = TRUE) +
    # scale_color_gradient2(low = "steelblue4", mid = "grey70", high = "red4",
    #                      limits = c(-.64, .64),
    #                      breaks = seq(-.4, .4, .2)) +
    scale_color_gradient(low = "grey70", high = "red4") +
    theme_bw(base_size = 12) +
    labs(x = expression(bold(~~~~~~~~~~~~~~~~~~~~GEOGRAPHY
                             ~~~~~~~~~~~~~~~~~~~~~~~~~
                               ~~~~~~~~~~~TEMPERATURE
                             ~~~~~~~~~~~~~~~~~~~~~~~~~~
                               ~~~~GROWING~SEASON
                             ~~~~~~~~~~~~~~~~~~~~~~~~~
                               ~~~~~~PRECIPITATION
                             ~~~~~~~~~~~~~~~~~~~~~~~~~~
                               ~~~~~~~~~~~~~~~ARIDITY)),
         y = expression(bold("mDatt"))) +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          strip.text = element_text(size = 11.5),
          axis.text = element_text(size = 9),
          axis.title.y = element_text(face = "bold"),
          axis.title.x = element_text(hjust = 0, face = "bold"),
          legend.position = "none",
          strip.background = element_blank(),
          strip.placement = "outside") +
    scale_alpha_continuous(range = c(.1, 1)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01),
      # breaks = seq(-10, 10, 5),
      # labels = seq(-10, 10, 5),
      expand = expansion(mult = c(.12, .26))) +
    scale_x_continuous(labels = scales::label_number(accuracy = 1),
                       expand = expansion(mult = c(.12, .12))) +
    facet_grid2(timepoint ~ clim, scales = "free", independent = "x",
                switch = "both") +
    facetted_pos_scales(
      y = list(timepoint == "decline" ~ scale_y_continuous(labels = scales::label_number(accuracy = 0.001),
                                                           expand = expansion(mult = c(.12, .26)))),
      x = list(
        clim == "Elevation (m)" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                                     breaks = seq(0, 3000, 1500)),
        clim == "MWMT (°C)" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                                 breaks = seq(12, 18, 3)),
        clim == "MCMT (°C)" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                                 breaks = seq(-30, 0, 10)),
        clim == "bFFP" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                            breaks = c(121, 152, 182),
                                            labels = c("May 1", "Jun 1", "Jul 1")),
        clim == "eFFP" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                            breaks = c(244, 274, 305),
                                            labels = c("Sep 1", "Oct 1", "Nov 1")),
        clim == "MSP (mm)" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                                labels = c(125, 250, 500),
                                                breaks = c(log(125), log(250), log(500))),
        clim == "PAS (mm)" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                                labels = c(125, 250, 500, 1000),
                                                breaks = c(log(125), log(250), log(500), log(1000))),
        clim == "CMD" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                                     breaks = seq(100, 500, 200)),
        clim == "AHM" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                           breaks = seq(10, 40, 10))
      )
    ))



ggsave(plot = r2_index_plot,
       filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\SI\\_SI_r2_NDVI.jpeg",
       device = jpeg,
       width = 15,
       height = 8.5,
       units = "in",
       dpi = 300,
       bg = "white")


####
#keep quadratic if quad is highly significant but linear isn't.
####

blups_r2_keep = blups_r2_plot %>% 
  dplyr::filter((index %in% keep_list
                | index %in% R_list)
                & p.value < p_threshold 
                & keep_flag == 1) %>% 
  dplyr::mutate(category = if_else(timepoint %in% c("early summer", "mid summer", "late summer"), "seasonal", timepoint),
                kind = if_else(index %in% R_list, "band", "index")) %>% 
  dplyr::group_by(index, timepoint) %>% 
  dplyr::mutate(r2_sum = sum(adj.r.squared)) %>% 
  dplyr::group_by(clim, timepoint, kind) %>%
  dplyr::mutate(r2_max = max(adj.r.squared),
                max_flag = if_else(adj.r.squared == r2_max, 1, 0)) %>%
  dplyr::group_by(index, timepoint) %>%
  dplyr::mutate(count = n()) 

# test = blups_r2_plot %>% 
#   dplyr::filter((index %in% keep_list
#                  #| index %in% R_list
#                  )
#                 & p.value < p_threshold 
#                 & keep_flag == 1) %>% 
#   dplyr::group_by(clim) %>%
#   dplyr::mutate(r2_max = max(adj.r.squared),
#                 max_flag = if_else(adj.r.squared == r2_max, 1, 0)) %>%
#   dplyr::filter(adj.r.squared == r2_max) %>% 
#   dplyr::mutate(index = factor(index,
#                                levels = c("CCI", "PRI", "GCC", "BCC", "ARI", "EWI9",
#                                           "NDVI", "mDatt", "RE_upper", 'NDRE3')))


################################################################################


(heatmap_band = blups_r2_plot %>%
    dplyr::filter(index %in% R_list) %>%
    ggplot(aes(x = factor(index),
               y = clim,
               group = index,
               fill = adj.r.squared,
               color = adj.r.squared)) +
    #coord_cartesian(expand = FALSE) +
    geom_tile(data = filter(blups_r2_plot, index %in% R_list), fill = "white", color = NA, show.legend = FALSE) +
    geom_tile(data = filter(blups_r2_keep, index %in% R_list
                            & type == "linear")) +
    # geom_point(data = filter(blups_r2_plot, index %in% keep_list & type == "quadratic" & p.value < .001),
    #            shape = 21, size = 6, stroke = NA, fill = "white") +
    geom_point(data = filter(blups_r2_keep, index %in% R_list
                             & type == "quadratic"),
               shape = 16, size = 6, stroke = NA, fill = "white",
    ) +
    # geom_tile(data = filter(blups_r2_keep, index %in% keep_list & max_flag == 1),
    #           color = "black",
    #           fill = NA,
    #           linewidth = 1,
    #           width=0.9,
    #           height=0.9,
    #           alpha = 1) +
    theme_bw(base_size = 12.5) +
    scale_y_discrete(limits = rev(levels(blups_r2_plot$clim)),
                     labels = c("MSP_log" = "MSP",
                                "PAS_log" = "PAS")) +
    scale_x_discrete(labels = index_labs) +
    labs(y = "GEOGRAPHY\n\n\nTEMPERATURE\n\n\nGROWING SEASON\n\n\nPRECIPITATION\n\n\nARIDITY",
         x = "Micasense spectral band",
         title = "SEASONAL                                                                                                                              DERIVED                                                 ") +
    # scale_fill_gradient2(low = "steelblue4", mid = "grey50", high = "red4",
    #                      limits = c(-.64, .64)) +
    scale_fill_gradient(name = "      linear  ",#expression(bold(R["adj."]^2)),
                        guide = guide_legend(#reverse = TRUE, 
                          title.position = "left", 
                          label.position = "bottom",
                          label = TRUE,
                          label.vjust = 1,
                          nrow = 1,
                          title.vjust = .8,
                          override.aes = list(linetype = 0)),
                        breaks = seq(.15, .704, .1), labels = c(expression(atop(" ", atop(".15", " "))), 
                                                               expression(atop(" ", atop(".25", " "))), 
                                                               expression(atop(" ", atop(".35", " "))), 
                                                               expression(atop(" ", atop(".45", bold(R["adj."]^2)))), 
                                                               expression(atop(" ", atop(".55", " "))), 
                                                               expression(atop(" ", atop(".65", " ")))), 
                        low = "grey95", high = "red4",
                        limits = c(.15, .65)) +
    scale_color_gradient(name = " quadratic ", guide = guide_legend(#reverse = TRUE, 
      title.position = "left", 
      label.position = "top",
      label = FALSE,
      nrow = 1,
      order = 1,
      # label.theme = element_text(size = 1, color = "white"),
      override.aes = list(fill="white", linetype = 0)),
      breaks = seq(.15, .704, .1), labels = c(".15", ".25", ".35", ".45", ".55", ".65"),
      low = "grey95", high = "red4",
      limits = c(.15, .704)) +
    #annotate("text", label = expression(bold(R["adj."]^2)), x = -.2, y)
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 12.5),
          axis.text.y = element_text(size = 12.5),
          panel.grid.major = element_blank(),
          plot.title = element_text(angle = 0, 
                                    vjust = .5, 
                                    hjust = 1, 
                                    size = 11,
                                    face = "bold",
                                    margin = margin(0, 0, 15, 0)),
          axis.title.x = element_text(angle = 0, 
                                      vjust = .0, 
                                      hjust = .5, 
                                      size = 15,
                                      margin = margin(10, 0, 0, 0)),
          axis.title.y = element_text(angle = 0, 
                                      vjust = .5, 
                                      hjust = 1, 
                                      size = 11,
                                      face = "bold",
                                      margin = margin(0, 10, 0, 0)),
          strip.background = element_blank(),
          strip.text = element_text(size = 15),
          axis.ticks = element_blank(),
          aspect.ratio = 1,
          legend.position = c(-.1, -.27),
          legend.box = "vertical",
          legend.spacing.y = unit(-.22, "cm"),
          legend.spacing.x = unit(-.01, "cm"),
          legend.text = element_text(size = 14)) +
    coord_cartesian(expand = FALSE) +
    facet_grid(. ~ factor(timepoint, levels =
                            c("early summer", "mid summer", "late summer", "late winter", "green-up", "decline"))))

ggsave(plot = heatmap_band,
       filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated\\SI\\heatmap_band.jpeg",
       device = jpeg,
       width = 15,
       height = 4.5,
       units = "in",
       dpi = 300,
       bg = "white")


index_labs = c(RE_upper = expression(RE['upper']),
               NDRE3 = expression(NDRE['740']),
               EWI9 = expression(EWI['9']),
               mDatt = "mDatt")

(heatmap_index = blups_r2_plot %>%
    dplyr::filter(index %in% keep_list) %>%
    ggplot(aes(x = factor(index,
                          levels = c("CCI", "PRI", "GCC", "BCC", "ARI", "EWI9",
                                     "mDatt", "RE_upper", "NDVI", 'NDRE3')),
               y = clim,
               group = index,
               fill = adj.r.squared,
               color = adj.r.squared)) +
    #coord_cartesian(expand = FALSE) +
    geom_tile(data = filter(blups_r2_plot, index %in% keep_list), fill = "white", color = NA, show.legend = FALSE) +
    geom_tile(data = filter(blups_r2_keep, index %in% keep_list
                            & type == "linear")) +
    # geom_point(data = filter(blups_r2_plot, index %in% keep_list & type == "quadratic" & p.value < .001),
    #            shape = 21, size = 6, stroke = NA, fill = "white") +
    geom_point(data = filter(blups_r2_keep, index %in% keep_list
                             & type == "quadratic"),
               shape = 16, size = 6, stroke = NA, fill = "white",
               ) +
    # geom_tile(data = filter(blups_r2_keep, index %in% keep_list & max_flag == 1),
    #           color = "black",
    #           fill = NA,
    #           linewidth = 1,
    #           width=0.9,
    #           height=0.9,
    #           alpha = 1) +
    theme_bw(base_size = 12.5) +
    scale_y_discrete(limits = rev(levels(blups_r2_plot$clim)),
                     labels = c("MSP_log" = "MSP",
                                "PAS_log" = "PAS")) +
    scale_x_discrete(labels = index_labs) +
    labs(y = "GEOGRAPHY\n\n\nTEMPERATURE\n\n\nGROWING SEASON\n\n\nPRECIPITATION\n\n\nARIDITY",
         x = "Vegetation index",
         title = "                      SEASONAL                                                                                                              DERIVED                                    ") +
    # scale_fill_gradient2(low = "steelblue4", mid = "grey50", high = "red4",
    #                      limits = c(-.64, .64)) +
    scale_fill_gradient(name = "      linear  ",#expression(bold(R["adj."]^2)),
                        guide = guide_legend(#reverse = TRUE, 
                                             title.position = "left", 
                                             label.position = "bottom",
                                             label = TRUE,
                                             label.vjust = 1,
                                             nrow = 1,
                                             title.vjust = .8,
                                             override.aes = list(linetype = 0)),
                        breaks = seq(.15, .704, .1), labels = c(expression(atop(" ", atop(".15", " "))), 
                                                               expression(atop(" ", atop(".25", " "))), 
                                                               expression(atop(" ", atop(".35", " "))), 
                                                               expression(atop(" ", atop(".45", bold(R["adj."]^2)))), 
                                                               expression(atop(" ", atop(".55", " "))), 
                                                               expression(atop(" ", atop(".65", " ")))), 
                        low = "grey95", high = "red4",
                        limits = c(.15, .704)) +
    scale_color_gradient(name = " quadratic ", guide = guide_legend(#reverse = TRUE, 
                                                                  title.position = "left", 
                                                                  label.position = "top",
                                                                  label = FALSE,
                                                                  nrow = 1,
                                                                  order = 1,
                                                                  # label.theme = element_text(size = 1, color = "white"),
                                                                  override.aes = list(fill="white", linetype = 0)),
                         breaks = seq(.15, .65, .1), labels = c(".15", ".25", ".35", ".45", ".55", ".65"),
                         low = "grey95", high = "red4",
                         limits = c(.15, .65)) +
    #annotate("text", label = expression(bold(R["adj."]^2)), x = -.2, y)
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 12.5),
          axis.text.y = element_text(size = 12.5),
          panel.grid.major = element_blank(),
          plot.title = element_text(angle = 0, 
                                    vjust = .5, 
                                    hjust = 1, 
                                    size = 12,
                                    face = "bold",
                                    margin = margin(0, 0, 15, 0)),
          axis.title.x = element_text(angle = 0, 
                                      vjust = .0, 
                                      hjust = .5, 
                                      size = 15,
                                      margin = margin(10, 0, 15, 0)),
          axis.title.y = element_text(angle = 0, 
                                      vjust = .5, 
                                      hjust = 1, 
                                      size = 12,
                                      face = "bold",
                                      margin = margin(0, 10, 0, 0)),
          strip.background = element_blank(),
          strip.text = element_text(size = 15),
          axis.ticks = element_blank(),
          aspect.ratio = 10/7,
          legend.background = element_rect(fill = NA),
          legend.position = c(-.13, -.25),
          legend.box = "vertical",
          legend.key.spacing = unit(.05, "cm"),
          legend.spacing.y = unit(-.25, "cm"),
          legend.spacing.x = unit(-.02, "cm"),
          legend.text = element_text(size = 14)) +
    coord_cartesian(expand = FALSE) +
    facet_grid(. ~ factor(timepoint, levels =
                            c("early summer", "mid summer", "late summer", "late winter", "green-up", "decline"))))


ggsave(plot = heatmap_index,
       filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\Figure_5.jpeg",
       device = jpeg,
       width = 15,
       height = 4.8,
       units = "in",
       dpi = 300,
       bg = "white")

# 
# (heatmap_dummy = blups_r2_plot %>% 
#     ggplot() +
#     geom_point(data = filter(blups_r2_plot,
#                              index %in% keep_list & r2_sign < 0),
#                aes(x = r2_sign, y = r2_sign, fill = adj.r.squared)) +
#     scale_fill_gradient(high = "red4", low = "white",
#                         breaks = seq(.1, .4, .1),
#                         labels = c("", "", "", ""),
#                         guide = guide_colorbar(#title.vjust = 0,
#                           title = NULL,
#                           title.position = "bottom")) +
#     theme(legend.title = element_blank(),
#           legend.spacing = unit(.01, 'cm'),
#           legend.title.align=0.5,
#           legend.text = element_blank()) +
#   new_scale_fill() +
#   geom_point(data = filter(blups_r2_plot,
#                            index %in% keep_list & r2_sign > 0),
#              aes(x = r2_sign, y = r2_sign, fill = adj.r.squared)) +
#     scale_fill_gradient(high = "steelblue4", low = "white",
#                         guide = guide_colorbar(title = expression(R^2),
#                                                title.position = "top")) +
#     
#     theme(legend.box = "vertical",
#           legend.position = "bottom",
#           legend.title = element_text(color = "black"),
#           legend.title.align=0.5,
#           legend.text = element_text(margin = margin(t = 0, b = 0, r = 0, l = 0),
#                                      vjust = 0)))
#

# 
# 
# (heatmap = ggarrange(heatmap_band + theme(legend.position = "none"),
#                      heatmap_index + theme(legend.position = "none"),
#                      widths = c(1, 1),
#                      ncol = 1))
# 
# 
# (heatmap_leg = ggarrange(dummy_plot, heatmap, 
#                         widths = c(.112, 1),
#                         ncol = 2))
# 
# ggplot2::ggsave(plot = heatmap_leg,
#                 filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated\\figure_5.jpeg",
#                 device = jpeg,
#                 width = 15,
#                 height = 7,
#                 units = 'in',
#                 dpi = 300,
#                 bg = 'white')
# 

################################################################################
# K means clustering

  
# # if no dignif diff between timepoints, keep one with highest r2
# timepoint_drop_list = blups_class_values %>% 
#   dplyr::select(-Prov, -order, -max_val) %>% 
#   dplyr::distinct() %>% 
#   dplyr::left_join(blups_r2_keep, by = c("index", "timepoint")) %>% 
#   drop_na(clim) %>% 
#   slice(which.max(adj.r.squared)) %>% 
#   dplyr::select(index:MS_LW_p, adj.r.squared) %>% 
#   pivot_wider(names_from = timepoint, values_from = adj.r.squared) %>% 
#   dplyr::mutate(MS = case_when(LS_MS_p > .05 & `mid summer` < `late summer` ~ 1,
#                                MS_LW_p > .05 & `mid summer` < `late winter` ~ 1),
#                 LS = case_when(LS_MS_p > .05 & `late summer` < `mid summer` ~ 1,
#                                LS_LW_p > .05 & `late summer` < `late winter`  ~ 1),
#                 LW = case_when(MS_LW_p > .05 & `late winter` < `mid summer` ~ 1,
#                                LS_LW_p > .05 & `late winter` < `late summer`  ~ 1)) %>% 
#   pivot_longer(`late summer`:`mid summer`, names_to = "timepoint", values_to = "adj.r.squared") %>% 
#   drop_na(adj.r.squared) %>% 
#   filter((MS == 1 & timepoint == "mid summer") |
#            (LS == 1 & timepoint == "late summer") |
#            (LW == 1 & timepoint == "late winter")) %>% 
#   unite("trait", c(index, timepoint))

(blups_r2_keep %>% 
  dplyr::filter(kind == "index"))$adj.r.squared %>% 
  hist(breaks = 100)

index_keep_list = blups_r2_keep %>% 
  #filter(max_flag == 1) %>% 
  #dplyr::filter(p.value < .0000005) %>% 
  #dplyr::filter(max_flag == 1) %>% 
  #                 adj.r.squared > .2) %>%
  unite("trait", c(index, timepoint))

# class_keep_list = blups_class %>% 
#   filter(p_value_class < .1) %>% 
#   unite("trait", c(index, timepoint))

df_cluster = mm_blups_change %>% 
  dplyr::mutate(ESTIMATE = timepoint_value) %>% 
  distinct(index, Prov, timepoint, ESTIMATE) %>% 
  dplyr::mutate(Prov = as.character(Prov)) %>% 
  saveRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\BLUPS\\Updated\\df_cluster.rds")

df_cluster = readRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\BLUPS\\Updated\\df_cluster.rds") %>% 
  mutate(cat = if_else(timepoint %in% c("early summer", "mid summer", "mid summer", "late summer", "late winter"), "seasonal", timepoint)) %>% 
  group_by(index, timepoint) %>% 
  dplyr::filter(index %in% keep_list) %>% 
  dplyr::mutate(across(ESTIMATE, ~ scale_this_new(.))) %>% 
                         #scale(.)[,1])) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-cat) %>% 
  #dplyr::mutate(ESTIMATE = scale(ESTIMATE[,1])) %>% 
  #filter(timepoint %in% c("decline", "greenup")) %>% 
  unite("trait", c(index, timepoint)) %>% 
  filter(trait %in% index_keep_list$trait
  #        #& !(trait %in% timepoint_drop_list$trait)
         ) %>%
  pivot_wider(names_from = trait, values_from = ESTIMATE) %>% 
  drop_na()
  

set.seed(123)

# #gap_stat = clusGap(dplyr::select(df_cluster, -Prov), FUN = kmeans, iter.max = 200, nstart = 100, K.max = 15, B = 200)
# #gap_stat = clusGap(dplyr::select(df_cluster, -Prov), FUN = pam, K.max = 15, B = 500)
# gap_stat = clusGap(dplyr::select(df_cluster, -Prov), FUN = hcut, K.max = 15, B = 500)
# 
# spec_wss = fviz_nbclust(dplyr::select(df_cluster, -Prov), 
#                         hcut,
#                         #pam, 
#                         #kmeans, iter.max = 200,
#                         nstart = 100, method = "wss", k.max = 15) +
#   #scale_y_log10() +
#     geom_vline(xintercept = 8, linewidth = 5, color = "red4", alpha = .2) +
#   theme(panel.grid.major = element_line(),
#         aspect.ratio = 1,
#         plot.title = element_blank())
# 
# spec_sil = fviz_nbclust(dplyr::select(df_cluster, -Prov),
#                         hcut,
#                         #pam, 
#                         #kmeans, iter.max = 200, 
#                         nstart = 100, method = "silhouette", k.max = 15) +
#     geom_vline(xintercept = 8, linewidth = 5, color = "red4", alpha = .2) +
#     theme(panel.grid.major = element_line(),
#         aspect.ratio = 1,
#         plot.title = element_blank())
# 
# spec_sil$layers[[3]]$aes_params$colour <- NA
# 
# spec_gap = fviz_gap_stat(gap_stat) +
#     geom_vline(xintercept = 8, linewidth = 5, color = "red4", alpha = .2) +
#   theme(panel.grid.major = element_line(),
#         aspect.ratio = 1,
#         plot.title = element_blank())
# 
# spec_gap$layers[[4]]$aes_params$colour <- NA
# 
# ################################################################################
# 
# ################################################################################
# 
# (n_clust = ggarrange(spec_wss, spec_sil, spec_gap,
#                     #clim_wss, clim_sil, clim_gap,
#                     ncol = 3,
#                     nrow = 1))
# 
# ggplot2::ggsave(plot = n_clust,
#                 filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated\\SI\\nclust.jpeg",
#                 device = jpeg,
#                 width = 14,
#                 height = 5,
#                 units = 'in',
#                 dpi = 300,
#                 bg = 'white')

################################################################################
# PCA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggbiplot)

df_pca = df_cluster %>% 
  dplyr::select(-Prov) %>% 
  rename_with(~ gsub(" ", "_", .x, fixed = TRUE)) %>% 
  rename_with(~ gsub("-", "_", .x, fixed = TRUE))

pca_spec = prcomp(~ ., df_pca, scale = FALSE)
pca_spec$x[,1:3]

pca_plot = bind_cols(df_cluster, pca_spec$x[,1:4]) %>% 
  # left_join(pops_plot %>% dplyr::select(Prov, cluster),
  #           by = "Prov") %>%
  distinct()


################################################################################
# plot relationships between PCs and climate

mm_pca = pca_plot %>% 
  dplyr::select(Prov, PC1, PC2, PC3) %>% 
  dplyr::mutate(PC1 = PC1 * -1,
                PC2 = PC2 * -1) %>% 
  dplyr::left_join(blups_summary %>% 
                     ungroup() %>% 
                     dplyr::select(Prov, clim, clim_value)) %>% 
  pivot_longer(PC1:PC3, names_to = "PC", values_to = "PC_value") %>% 
  dplyr::distinct()

pca_model = mm_pca %>% 
  #filter(index %in% keep_list) %>% 
  group_by(PC, clim) %>% 
  nest() %>% 
  dplyr::mutate(linear = map(data, ~lm(PC_value ~ clim_value, data = .)),
                quadratic = map(data, ~lm(PC_value ~ clim_value + I(clim_value^2), data = .))) %>% 
  pivot_longer(cols = c(quadratic, linear), names_to = "type", values_to = "model") %>% 
  dplyr::mutate(clim = factor(clim, levels = c("Latitude", 
                                               "Elevation", 
                                               
                                               "MCMT", # cold winters
                                               "MWMT", # hot summers
                                               
                                               "bFFP", # spring frost risk
                                               "eFFP", # fall frost risk
                                               
                                               "MSP_log", # summer moisture
                                               "PAS_log", # spring moisture / winter protection
                                               
                                               "CMD", # drought
                                               "AHM")))


pca_r2 = pca_model %>%
  dplyr::mutate(glance = map(model, broom::glance),
                tidy = map(model, broom::tidy),
                slope = tidy %>% map_dbl(function(x) x$estimate[2]),
  ) %>% 
  unnest(c(data, glance, slope)) %>% 
  dplyr::select(-model, -tidy) %>% 
  dplyr::select(PC:clim, type, adj.r.squared:slope) %>% 
  distinct() %>% 
  dplyr::group_by(clim, PC) %>% 
  nest() %>% 
  dplyr::mutate(p.adjusted = map(data, function(df) p.adjust(df$p.value, method = "bonferroni"))) %>% 
  #dplyr::mutate(q.value = map(data, function(x) qvalue::qvalue(x$p.value)$qvalues)) %>% 
  unnest(c(data, p.adjusted))


pca_r2_plot = pca_r2 %>% 
  # left_join(VI_list %>% 
  #             dplyr::select(index, order)) %>% 
  group_by(PC, clim) %>% 
  dplyr::mutate(r2_lin_slope = if_else(type == "linear", slope, NA_real_),
                r2_slope = mean(r2_lin_slope, na.rm = TRUE)) %>% 
  group_by(PC, clim, type) %>% 
  dplyr::mutate(r2_sign = if_else(r2_slope > 0, adj.r.squared, adj.r.squared * -1)) %>% 
  distinct() %>% 
  dplyr::group_by(PC, clim) %>% 
  dplyr::mutate(AIC_min = min(AIC),
                count = if_else(AIC < (AIC_min + 40), 1, 0),
                sum = sum(count),
                # keep quadratic if only quadratic is signif at alpha level
                keep_flag = if_else((type == "linear" & sum != 1) |
                                      (type != "linear" & sum == 1),
                                    1, 0))

pca_pred = pca_model %>% 
  dplyr::mutate(preds = map2(data, model, add_predictions)) %>%
  dplyr::mutate(preds = map(preds, dplyr::select, pred)) %>% 
  unnest(c(data, preds)) %>% 
  dplyr::select(-model) %>% 
  dplyr::distinct()

pca_summary = pca_pred %>% 
  dplyr::left_join(pca_r2_plot) %>% 
  dplyr::group_by(clim, type, PC)


# summary figure 
dat_plot_pca = pca_summary %>% 
  filter(PC != "PC4") %>% 
  group_by(PC) %>% 
  dplyr::mutate(max_val = max(PC_value)) %>% 
  group_by(PC, clim, type) %>% 
  dplyr::mutate(max_r2 = max(adj.r.squared),
                r_lab = as.character(round(adj.r.squared, 2)),
                r_lab = str_remove(r_lab, "^0+"),
                r_lab = str_pad(r_lab, 3, "right", "0")) %>% 
  dplyr::mutate(clim = factor(clim,
                              levels = c( "Latitude", "Elevation", 
                                          "MWMT", "MCMT",  
                                          "bFFP", "eFFP",
                                          "MSP_log", "PAS_log", 
                                          "CMD", "AHM"))) %>% 
  distinct() %>% 
  dplyr::mutate(clim = dplyr::recode(clim, 
                                     Latitude = "Latitude (°N)",
                                     Elevation = "Elevation (m)",
                                     MWMT = "MWMT (°C)",
                                     MCMT = "MCMT (°C)",
                                     MSP_log = "MSP (mm)",
                                     PAS_log = "PAS (mm)"),
                # PC = dplyr::recode(PC,
                #                    "PC1" = "PC1 — greenness",
                #                    "PC2" = "PC2 — green-up",
                #                    "PC3" = "PC3 — red edge stress"))
                PC = dplyr::recode(PC,
                                   "PC1" = "PC1\ngreenness",
                                   "PC2" = "PC2\ngreen-up",
                                   "PC3" = "PC3\nred edge stress"))

# test = blups_summary %>% 
#   filter(index == "mDatt" &
#            clim == "AHM" &
#            timepoint == "green-up"
#   )
# 
# lin = lm(timepoint_value ~ clim_value, data = test)
# 
# quad = lm(timepoint_value ~ clim_value + I(clim_value^2), data = test)

PC_labs = c("PC1" = "PC1 — greenness",
            "PC2" = "PC2 — green-up",
            "PC3" = "PC3 — red edge stress")

# clim_labs = dplyr::recode(dat_plot_pca$clim, "MSP_log" = "log(MSP)",
#               "PAS_log" = "log(PAS)")

# 10 clim vars * 3 PCs
p_threshold_pcs = .005 / (10)

(r2_pca_plot = dat_plot_pca %>% 
    ggplot(aes(x = clim_value, y = pred, group = adj.r.squared)) +
    geom_point(data = filter(dat_plot_pca, type == "linear"), 
               aes(y = PC_value, color = adj.r.squared), size = 2.5, alpha = .55) +
    geom_line(data = filter(dat_plot_pca, p.value < p_threshold_pcs
                            & keep_flag == 1
                            ),
              stat = "smooth", linewidth = .7, alpha = .8) +
    geom_text(data = filter(dat_plot_pca, PC_value == max_val 
                            & keep_flag == 1
                            &  p.value < p_threshold_pcs
                            & r2_sign > 0), 
              aes(x = -Inf, y = max_val, 
                  label = paste0(expression(R^2))),
              size = 4,
              hjust = -.3,
              vjust = -.28,
              alpha = 1,
              parse = TRUE) +
    geom_text(data = filter(dat_plot_pca, PC_value == max_val 
                            & keep_flag == 1
                            &  p.value < p_threshold_pcs
                            & r2_sign > 0), 
              aes(x = -Inf, y = max_val, 
                  label = paste0("=", r_lab)),
              size = 4,
              hjust = -.8,
              vjust = -.4,
              alpha = 1,
              parse = FALSE) +
    geom_text(data = filter(dat_plot_pca, PC_value == max_val 
                            & keep_flag == 1
                            &  p.value < p_threshold_pcs
                            & r2_sign < 0), 
              aes(x = Inf, y = max_val, 
                  label = paste0(expression(R^2))),
              size = 4,
              hjust = 3.05,
              vjust = -.28,
              alpha = 1,
              parse = TRUE) +
    geom_text(data = filter(dat_plot_pca, PC_value == max_val 
                            & keep_flag == 1
                            &  p.value < p_threshold_pcs
                            & r2_sign < 0), 
              aes(x = Inf, y = max_val, 
                  label = paste0("=", r_lab)),
              size = 4,
              hjust = 1.15,
              vjust = -.4,
              alpha = 1,
              parse = FALSE) +
    scale_color_gradient(low = "grey90", high = "red4",
                          limits = c(-.1, .605),
                          breaks = seq(-.4, .4, .2)) +
    theme_bw(base_size = 12) +
    labs(x = expression(bold(~~~~~~~~~~~~~~~~~~GEOGRAPHY
                             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               ~~~~~~~~~~~TEMPERATURE
                             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               ~~~~~GROWING~SEASON
                             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               ~~~~~~PRECIPITATION
                             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               ~~~~~~~~~~~~~~~ARIDITY))) +
    theme(aspect.ratio = 1.1,
          panel.grid = element_blank(),
          strip.text = element_text(size = 11.5),
          axis.text = element_text(size = 9),
          axis.title.y = element_blank(),
          axis.title.x = element_text(hjust = 0, face = "bold"),
          legend.position = "none",
          strip.background = element_blank(),
          strip.placement = "outside") +
    scale_alpha_continuous(range = c(.1, 1)) +
    scale_y_continuous(#labels = scales::label_number(accuracy = 0.01,
      breaks = seq(-2, 2, 1),
      labels = seq(-2, 2, 1),
      expand = expansion(mult = c(.14, .25))) +
    scale_x_continuous(labels = scales::label_number(accuracy = 1),
                       expand = expansion(mult = c(.12, .12))) +
    facet_grid2(PC ~ clim, scales = "free", independent = "x",
                switch = "both") +
    facetted_pos_scales(
      x = list(
        clim == "Elevation (m)" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                                     breaks = seq(0, 3000, 1500)),
        clim == "MWMT (°C)" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                                 breaks = seq(12, 18, 3)),
        clim == "MCMT (°C)" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                                 breaks = seq(-30, 0, 10)),
        clim == "bFFP" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                            breaks = c(121, 152, 182),
                                            labels = c("May 1", "Jun 1", "Jul 1")),
        clim == "eFFP" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                            breaks = c(244, 274, 305),
                                            labels = c("Sep 1", "Oct 1", "Nov 1")),
        clim == "MSP (mm)" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                                labels = c(125, 250, 500),
                                                breaks = c(log(125), log(250), log(500))),
        clim == "PAS (mm)" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                                labels = c(125, 250, 500, 1000),
                                                breaks = c(log(125), log(250), log(500), log(1000))),
        clim == "CMD" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                           breaks = seq(100, 500, 200)),
        clim == "AHM" ~ scale_x_continuous(expand = expansion(mult = c(.14, .14)),
                                           breaks = seq(10, 40, 10))
      )
      # ,
      # y = list(PC == "PC1 — greenness" ~ scale_y_continuous(expand = expansion(mult = c(.14, .225)),
      #                                                       breaks = c(-10, 0, 10)),
      #          PC == "PC2 — green-up" ~ scale_y_continuous(expand = expansion(mult = c(.14, .23)),
      #                                                      breaks = c(-5, 0, 5)),
      #          PC == "PC3 — red edge" ~ scale_y_continuous(expand = expansion(mult = c(.14, .15)),
      #                                                      breaks = c(-4, 0, 4))
      #)
    ))



ggsave(plot = r2_pca_plot,
       filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\figure_6_PCs.jpeg",
       device = jpeg,
       width = 15,
       height = 6,
       units = "in",
       dpi = 300,
       bg = "white")


################################################################################

test = df_cluster
  # pca_plot %>%
  # dplyr::select(PC1:PC3,#contains(keep_list)
  #               Prov)
#%>% 
#   remove_rownames() %>% 
#   column_to_rownames(var="Prov")

#cluster_7 = kmeans(dplyr::select(df_cluster, -Prov), 9, iter.max = 200, nstart = 100)
#cluster_7 = cluster::pam(dplyr::select(df_cluster, -Prov), k = 9, nstart = 100)
#cluster_7 = factoextra::hkmeans(dplyr::select(df_cluster, -Prov), 6, iter.max = 200)
cluster_7 = hclust(dist(test %>%
                          #dplyr::select(-(PC1:PC4)) %>% 
                          remove_rownames() %>%
                          column_to_rownames(var="Prov")),
                   method = "ward.D2")

fviz_dend(cluster_7, k = 8)


loadings = pca_spec$rotation %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  dplyr::select(rowname:PC3) %>% 
  #pivot_longer(cols = PC1:PC3, names_to = "PC", values_to = "loading") %>% 
  dplyr::mutate(PC1_abs = abs(PC1),
                PC2_abs = abs(PC2),
                PC3_abs = abs(PC3))

# fviz_cluster(cluster_7, data = dplyr::select(df_cluster, -Prov), axes = c(1, 2))+
#   #scale_colour_manual(values = c("steelblue", "springgreen3", "firebrick2")) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme_bw()
# 
# fviz_cluster(cluster_7, data = dplyr::select(df_cluster, -Prov), axes = c(2, 3))+
#   #scale_colour_manual(values = c("steelblue", "springgreen3", "firebrick2")) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme_bw()

# fviz_cluster(cluster_7, data = dplyr::select(df_cluster, -Prov), axes = c(1, 3))+
#   #scale_colour_manual(values = c("steelblue", "springgreen3", "firebrick2")) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme_bw()

################################################################################

library(rnaturalearth)
library(rnaturalearthdata)

#clust_cols = c("#7E1A0A", "#BE5A2A", "#CA9A28", "#7E9037", "#1A7D62", "#3AA2BC", "#5C79C4", "#9A86B5",  "#9A3675", "#BE5A2A" )
clust_cols = c("1" = "#9E2A0A", "2" = "#cd6f37", "3" = "#CA9A28", 
               "4" = "#7e9633", "5" = "#4d8d57", "6" = "#24759f",
                "7" = "#9080b0", "8" = "#b06a70")#, "9" = "#8A2665")

cols_species = c("#7F8F9F", "#556555", "#111121")

states = st_read("D:\\Sync\\_Sites\\Sx_Genecology\\_rnaturalearth\\ne_50m_admin_1_states_provinces.shp") %>% 
  filter(admin %in% c("United States of America", "Canada"))

countries = ne_countries(scale = 50, returnclass = "sf") %>% 
  filter(brk_name %in% c("Canada", "United States"))

ocean = st_read("D:\\Sync\\_Sites\\Sx_Genecology\\_rnaturalearth\\ne_50m_ocean.shp")

rast_df1 = read_rds("D:\\Sync\\_Sites\\Sx_Genecology\\rast_df1")

###

cluster_plot_start = df_cluster %>% 
  dplyr::mutate(clust = cutree(cluster_7, k = 8))
  #dplyr::mutate(clust = cluster_7$cluster)

# number clusters by elevation

(cluster_plot = sx_pops_clim %>% 
  dplyr::select("Prov",
                "Latitude") %>% 
  dplyr::mutate(Prov = as.character(Prov)) %>% 
  left_join(cluster_plot_start %>% 
              dplyr::select(Prov, clust)) %>% 
  drop_na(clust) %>% 
  group_by(clust) %>% 
  dplyr::mutate(mean_lat = mean(Latitude)) %>% 
  ungroup() %>% 
  dplyr::mutate(cluster = dense_rank(mean_lat),
         cluster = factor(cluster)) %>% 
  distinct())
  
pops_plot = cluster_plot %>% 
  left_join(read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\SxGenecology_predicted_Euc.csv") %>% 
              mutate(Prov = as.character(sxProv)) %>% 
              dplyr::select(Prov, lat, long) %>% 
              distinct()) %>% 
  sf::st_as_sf(coords = c("long", "lat"), crs = 4326)

#autoplot(pca_clim, colour = 'name', label = TRUE)
(bp1 = ggbiplot::ggbiplot(pca_spec, choices = c(1,2), groups = pca_plot$cluster,
                varname.size = 3, labels.size = 3,
                #var.axes = FALSE, 
                color = "black", labels = df_cluster$Prov
                ) +
    theme_bw(base_size = 20) +
    theme(aspect.ratio = 1))

(bp2 = ggbiplot::ggbiplot(pca_spec, choices = c(2, 3), groups = pca_plot$cluster,
                varname.size = 3, labels.size = 3,
                #var.axes = FALSE, 
                color = "black", labels = df_cluster$Prov) +
    theme_bw(base_size = 20) +
    theme(aspect.ratio = 1))

### facets

pops_plot_facet = pops_plot %>% 
  dplyr::mutate(group_plot = case_when(cluster %in% c(1, 2) ~ "low",
                                cluster %in% c(3, 4, 5) ~ "mid",
                                cluster %in% c(6, 7, 8, 9, 10) ~ "high"),
           group_plot = factor(group_plot, levels = c("low", "mid", "high"))) %>% 
  distinct()

# # sx range
# cluster_map = ggplot() +
#   geom_raster(data = rast_df1,  aes(x = x, y = y, fill = elev)) +
#   scale_fill_gradient(high = "white", low = "grey20", guide = "none") +
#   geom_sf(data = countries, fill = NA, color = "grey20") +
#   geom_sf(data = states, fill = NA) +
#   geom_sf(data = ocean, fill = "white", color = "grey20") + 
#   geom_sf(data = pops_plot, aes(color = cluster), shape = 19, size = 2) +
#   scale_color_manual(values = clust_cols) +
#   coord_sf(xlim = c(225000, 2800000), ylim = c(-768333, 2150000),
#            crs = st_crs(3005)) +
#   theme(panel.border = element_rect(fill = NA, color = "grey25", size = .5),
#         legend.position = "none",#c(0.94, 0.5),
#         legend.background = element_rect(fill = "white", color = "black"),
#         legend.text = element_text(size = 16),
#         legend.title = element_text(size = 24),
#         legend.margin = margin(8, 8, 14, 8),
#         panel.background = element_rect(fill = "white", color = "black", linewidth = .5),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank() +
#           facet_wrap(. ~ factor(group_plot), ncol = 3))


# ggplot2::ggsave(plot = cluster_map,
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\cluster_map_9.jpeg",
#        device = jpeg,
#        width = 12.5,
#        height = 14,
#        units = 'in',
#        dpi = 300,
#        bg = 'white')

cluster_centroids = pops_plot_facet %>% 
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::summarise(lat = mean(lat),
            lon = mean(lon),
            group_plot = group_plot) %>%
  dplyr::mutate(
                lon = if_else(cluster %in% c(6), lon - 0, lon),
                lon = if_else(cluster %in% c(5), lon - 2.6, lon),
                lon = if_else(cluster %in% c(3), lon - 1.3, lon),
                lon = if_else(cluster %in% c(1), lon - .2, lon),
                lon = if_else(cluster %in% c(2), lon + .1, lon),
                lon = if_else(cluster %in% c(4), lon - 1, lon),
                lon = if_else(cluster %in% c(7), lon + .5, lon),
    #             lat = if_else(cluster %in% c(4, 7), lat - .5, lat),
                lat = if_else(cluster %in% c(3, 6, 4), lat - 1.3, lat),
                lat = if_else(cluster %in% c(2), lat + 1.9, lat),
                lat = if_else(cluster %in% c(8), lat - 1, lat),
                lat = if_else(cluster %in% c(5), lat + 1.3, lat)
                ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  distinct()


# sx range
cluster_map_facet = ggplot() +
    geom_raster(data = rast_df1,  aes(x = x, y = y, fill = elev)) +
    scale_fill_gradient(high = "white", low = "grey20", guide = "none") +
    geom_sf(data = countries, fill = NA, color = "grey20") +
    geom_sf(data = states, fill = NA) +
    geom_sf(data = ocean, fill = "white", color = "grey20") + 
    geom_sf(data = pops_plot_facet, aes(color = factor(cluster)), shape = 16, size = 1.8) +
  geom_sf(data = cluster_centroids, aes(color = factor(cluster)), shape = 19, size = 9) +
  annotate("text", x = 1900000, y = 1490000, label = "C    A    N    A    D    A",
           size = 5.5, color = "grey36", angle = 13) +
  annotate("text", x = 2170000, y = -354500, label = "U  N  I  T  E  D    S  T  A  T  E  S\no f\nA  M  E  R  I  C  A",
           size = 4.1, color = "grey36", angle = 11) +
  annotate("text", x = 650000, y = -13000, label = "P   A   C   I   F   I   C\n\nO   C   E   A   N",
           size = 3, color = "grey36", angle = 0) +
  annotate("text", vjust = 1,  x = 1966000, y = 494000, label = "   R  O  C  K  Y           M  O  U  N  T  A  I  N  S",
           size = 2.7, color = "grey36", angle = -45) +
    annotate("text", vjust = 1,  x = 890000, y = 1100000, label = " C  O  A  S  T         M  O  U  N  T  A  I  N  S    ",
             size = 2.7, color = "grey36", angle = -61) +
  scale_color_manual(values = clust_cols) +
  geom_sf_text(data = cluster_centroids, aes(label = cluster), color = "black", fontface = "bold", size = 4.5) +
    coord_sf(xlim = c(225000, 2800000), ylim = c(-768333, 2150000),
             crs = st_crs(3005)) +
    theme(panel.border = element_rect(fill = NA, color = "grey25", size = .5),
          legend.position = "none",#c(0.94, 0.5),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 24),
          legend.margin = margin(8, 8, 14, 8),
          panel.background = element_rect(fill = "white", color = "black", linewidth = .5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_blank()) +
    facet_wrap(. ~ factor(group_plot), ncol = 3)

#cluster_map_facet

# ggplot2::ggsave(plot = cluster_map_facet,
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\cluster_map_facet_9.jpeg",
#        device = jpeg,
#        width = 12.5,
#        height = 14,
#        units = 'in',
#        dpi = 300,
#        bg = 'white')


pops_plot_clim = sx_pops_clim %>% 
    dplyr::select("Prov",
                  "Latitude", 
                  "Elevation",
                  
                  "MCMT", # cold winters
                  "MWMT", # hot summers
                  
                  "bFFP", # spring frost risk
                  "eFFP", # fall frost risk
                  
                  "MSP_log", # summer moisture
                  "PAS_log", # spring moisture / winter protection
                  
                  "CMD", # drought
                  "AHM" # aridity
    ) %>% 
    dplyr::mutate(Prov = as.character(Prov)) %>% 
    left_join(pops_plot %>% 
                dplyr::select(Prov, cluster)) %>% 
    drop_na(cluster) %>% 
    pivot_longer(Latitude:AHM, names_to = "clim", values_to = "clim_value") %>% 
    ungroup() %>% 
    dplyr::mutate(cluster = factor(cluster),
           clim = factor(clim,
                         levels = c( "Latitude", "Elevation", 
                                    "MWMT", "MCMT",  
                                    "bFFP", "eFFP",
                                    "MSP_log", "PAS_log", 
                                    "CMD", "AHM"))) %>% 
    distinct() %>% 
  dplyr::mutate(clim = dplyr::recode(clim, 
                                     Latitude = "Latitude (°N)",
                                     Elevation = "Elevation (m)",
                                     MWMT = "MWMT (°C)",
                                     MCMT = "MCMT (°C)",
                                     MSP_log = "MSP (mm)",
                                     PAS_log = "PAS (mm)"))

pops_plot_clim_lab = pops_plot_clim %>% 
  dplyr::select(cluster, clim, clim_value) %>% 
  group_by(clim) %>% 
  dplyr::mutate(max_clim = max(clim_value),
         min_clim = min(clim_value),
         clim_value = max_clim + ((max_clim - min_clim) * .25)) %>% 
  distinct(cluster, clim_value)

library(tidytext)

cluster_clim_facet = pops_plot_clim %>% 
  ggplot(aes(x = clim_value, y = reorder_within(cluster, clim_value, clim), color = cluster)) +
  geom_point(size = 3, alpha = .3) +
    geom_point(data = pops_plot_clim_lab, shape = 16, size = 9) +
    scale_color_manual(values = clust_cols) +
    geom_text(data = pops_plot_clim_lab, aes(label = cluster),
              color = "black", fontface = "bold", size = 4.5) +
    scale_x_continuous(expand = expansion(mult = c(.07, .14))) +
    theme_bw(base_size = 12) +
    labs(x = expression(bold(~~~~~~~~~~~~~~~~~~~GEOGRAPHY
                             ~~~~~~~~~~~~~~~~~~~~~~~~~~
                               ~~~~~~~~~~~TEMPERATURE
                             ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               ~~~~GROWING~SEASON
                             ~~~~~~~~~~~~~~~~~~~~~~~~~~
                               ~~~~~~PRECIPITATION
                             ~~~~~~~~~~~~~~~~~~~~~~~~~~
                               ~~~~~~~~~~~~~~~ARIDITY))) +
    theme(strip.placement = "outside",
          #aspect.ratio = 2,
          strip.background = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          strip.text = element_text(size = 11.5),
          axis.title.x = element_text(hjust = 0, face = "bold")) +
  facet_grid2(. ~ clim,
             scales = "free",
             independent = "all",
             switch = "both") +
    facetted_pos_scales(
      x = list(
        clim == "Elevation (m)" ~ scale_x_continuous(expand = expansion(mult = c(.07, .16)),
                                                     breaks = seq(0, 3000, 1500)),
        clim == "MWMT (°C)" ~ scale_x_continuous(expand = expansion(mult = c(.07, .16)),
                                                 breaks = seq(12, 18, 3)),
        clim == "MCMT (°C)" ~ scale_x_continuous(expand = expansion(mult = c(.07, .15)),
                                                 breaks = seq(-30, 0, 10)),
        clim == "bFFP" ~ scale_x_continuous(expand = expansion(mult = c(.07, .16)),
                                            breaks = c(121, 152, 182),
                                            labels = c("May 1", "Jun 1", "Jul 1")),
        clim == "eFFP" ~ scale_x_continuous(expand = expansion(mult = c(.07, .16)),
                                            breaks = c(244, 274, 305),
                                            labels = c("Sep 1", "Oct 1", "Nov 1")),
        clim == "MSP (mm)" ~ scale_x_continuous(expand = expansion(mult = c(.07, .16)),
                                                labels = c(125, 250, 500),
                                                breaks = c(log(125), log(250), log(500))),
        clim == "PAS (mm)" ~ scale_x_continuous(expand = expansion(mult = c(.07, .16)),
                                                labels = c(125, 250, 500, 1000),
                                                breaks = c(log(125), log(250), log(500), log(1000))),
        clim == "CMD" ~ scale_x_continuous(expand = expansion(mult = c(.07, .16)),
                                           breaks = seq(100, 500, 200)),
        clim == "AHM" ~ scale_x_continuous(expand = expansion(mult = c(.07, .16)),
                                           breaks = seq(10, 40, 10))
      )
    )


# ggplot2::ggsave(plot = cluster_clim_facet,
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\cluster_clim_facet_9.jpeg",
#        device = jpeg,
#        width = 12.5,
#        height = 14,
#        units = 'in',
#        dpi = 300,
#        bg = 'white')

### ternary plot species

# new hybrid matches from Jon

hybrid_prop_new = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\SxGenecology_closestGenotypes_240208.csv") %>% 
  dplyr::mutate(Prov = prov,
                propGla = glauca,
                propSit = sitchensis,
                propEng = engelmannii)
  
pops_plot_species = pops_plot %>% 
  distinct(Prov, cluster) %>% 
  left_join(hybrid_prop_new %>% 
              dplyr::select(Prov, propEng, propGla, propSit) %>% 
              dplyr::mutate(Prov = as.character(Prov))) %>% 
  # left_join(pols_dat_skim %>% 
  #             dplyr::select(Prov, propEng, propGla, propSit, Class) %>% 
  #             dplyr::mutate(Prov = as.character(Prov))) %>% 
  distinct() %>% 
  dplyr::mutate(prop_total = propEng + propGla + propSit,
         propEng = propEng / prop_total, 
         propGla = propGla / prop_total, 
         propSit = propSit / prop_total)  %>% 
  left_join(sx_pops_clim %>% 
              dplyr::select(Prov, Latitude) %>% 
              dplyr::mutate(Prov = as.character(Prov)))%>% 
  group_by(cluster) %>% 
  dplyr::mutate(rank_prop = row_number(Latitude),
         rank_prop_mean = mean(rank_prop)) %>% 
  pivot_longer(propEng:propSit, names_to = "species", values_to = "prop") %>% 
  dplyr::mutate(species = factor(species, levels = c("propGla", "propEng", "propSit"))) %>% 
  distinct()

pops_lab_species = pops_plot_species %>% 
  distinct(rank_prop_mean, cluster, species)

library(ggnewscale)


(cluster_species = pops_plot_species  %>% 
    #drop_na(prop) %>% 
    ggplot(aes(x = factor(rank_prop), 
               y = prop, fill = species, group = cluster)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = cols_species,#c("#888888", "#555555", "#111111"),
                      labels = c("P. glauca", "P. engelmannii", "P. sitchensis"),
                      name = "Species") +
    theme_bw() +
    labs(fill = "Species") +
    geom_point(data = pops_lab_species,
               aes(x = rank_prop_mean, color = factor(cluster)),
               shape = 16, size = 9,
               y = 1.11,
               show.legend = FALSE) +
    scale_color_manual(values = clust_cols) +
    geom_text(aes(x = rank_prop_mean, label = cluster),
              color = "black", fontface = "bold", size = 4.5,
              y = 1.11, 
              show.legend = FALSE) +
    # new_scale_color() +
    # geom_point(aes(color = Class),
    #            shape = 15,
    #                size = 3,
    #            y = -.04) +
    # scale_color_manual(values = class_cols) +
    theme(strip.placement = "outside",
          #aspect.ratio = 20,
          strip.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          strip.text = element_blank(),
          legend.margin = margin(c(12, -8, 2, 2)),
          legend.spacing.y = unit(.2, 'cm'),
          legend.text = element_text(size = 12, face = "italic"),
          legend.title = element_blank(),
          legend.position = "left") + #c(.95, .7)) +
    coord_cartesian(expand = FALSE) +
    expand_limits(y = c(0, 1.2),
                  x = c(.08, 0)) +
    guides(fill = guide_legend(override.aes = list(shape = NA),
                               byrow = TRUE)) +
    facet_grid2(. ~ cluster,
                scales = "free_x",
                space = "free_x",
                switch = "both"))

legend_1 = get_legend(cluster_species)

# ggplot2::ggsave(plot = cluster_species,
#                 filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\cluster_species_new_matches.jpeg",
#                 device = jpeg,
#                 width = 14,
#                 height = 5,
#                 units = 'in',
#                 dpi = 300,
#                 bg = 'white')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rect = rectGrob(
  x = unit(.19, "in"),
  y = unit(1, "npc") - unit(.15, "in"),
  width = unit(.37, "in"),
  height = unit(.37, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))



font_s = 17.5

lab_a = textGrob(
  label = "(a)",
  x = unit(.22, "in"),
  y = unit(1, "npc") - unit(.23, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))


lab_b = textGrob(
  label = "(b)",
  x = unit(.22, "in"),
  y = unit(1, "npc") - unit(.23, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_c = textGrob(
  label = "(c)",
  x = unit(.22, "in"),
  y = unit(1, "npc") - unit(.23, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

map_fig = ggarrange((ggdraw(cluster_species) +
                       draw_grob(rect) +
                       draw_grob(lab_a)), #cluster_species,
                    (ggdraw(cluster_map_facet) +
                       draw_grob(rect) +
                       draw_grob(lab_b)),
                    (ggdraw(cluster_clim_facet) +
                       draw_grob(rect) +
                       draw_grob(lab_c)),
                    heights = c(.505, 1, .65),
                    ncol = 1)


ggplot2::ggsave(plot = map_fig,
                filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\map_fig.jpeg",
                device = jpeg,
                width = 14,
                height = 11.5,
                units = 'in',
                dpi = 300,
                bg = 'white')


################################################################################

library(scatterpie)


pca_plot = bind_cols(df_cluster, pca_spec$x[,1:4]) %>% 
  left_join(pops_plot %>% dplyr::select(Prov, cluster),
            by = "Prov") %>%
  distinct() %>% 
  dplyr::select(Prov, cluster, PC1, PC2, PC3) %>% 
  dplyr::mutate(PC1 = PC1 * -1,
                PC2 = PC2 * -1)


  
pca_plot_species = pca_plot %>% 
  left_join(hybrid_prop_new %>% mutate(Prov = as.character(Prov),
                                   prop_total = propEng + propGla + propSit,
                                   propEng = propEng / prop_total,
                                   propGla = propGla / prop_total,
                                   propSit = propSit / prop_total),
            by = "Prov")


loadings = pca_spec$rotation %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    dplyr::select(rowname:PC4) %>% 
    #pivot_longer(cols = PC1:PC3, names_to = "PC", values_to = "loading") %>% 
  dplyr::mutate(PC1_abs = abs(PC1),
                PC2_abs = abs(PC2),
                PC3_abs = abs(PC3),
                PC4_abs = abs(PC4)) %>% 
  # left_join(pca_plot %>% dplyr::select(PC1:PC3) %>% 
  #     pivot_longer(PC1:PC3, names_to = "PC", values_to = "value") %>% 
  #     dplyr::group_by(PC) %>% 
  #     dplyr::summarise(min = min(value),
  #                      max = max(value)) %>% 
  #     dplyr::mutate(range = max - min,
  #                   unit = range * .05)) %>% 
  # dplyr::group_by(PC) %>% 
  # dplyr::mutate(h_just = if_else(loading > 0, -.1, 1.1),
  #               sign_off = if_else(loading > 0, 1, 0)) %>% 
  # dplyr::group_by(PC, h_just) %>% 
  #dplyr::filter(loading_abs > .2) %>% 
  dplyr::mutate(name = str_replace(rowname, 
                                   "RE_upper",
                                   "RE['upper']"),
                name = str_replace(name,
                                   "NDRE3",
                                   "NDRE['740']"),
                # name = str_replace(name, 
                #                    "EWI9", "EWI['9']"),
                # name = str_replace(name, 
                #                    "Datt", "mDatt"),
                name = str_replace(name,
                                   "_early_summer",
                                   "^'early summer'"),
                name = str_replace(name,
                                      "_mid_summer",
                                      "^'mid summer'"),
                name = str_replace(name,
                                      "_late_summer",
                                   "^'late summer'"),
                name = str_replace(name,
                                      "_late_winter",
                                   "^'late winter'"),
                name = str_replace(name,
                                      "_decline",
                                      "^'decline'"),
                name = str_replace(name,
                                      "_green_up",
                                      "^'green-up'")
                # rowname = str_replace(rowname,
                #                    "_mid_summer",
                #                    "['mid summer']"),
                # rownname = str_replace(rowname,
                #                    "_late_summer",
                #                    "['late summer']"),
                # rownname = str_replace(rowname,
                #                    "_late_winter",
                #                    "['late winter']"),
                # rowname = str_replace(rowname,
                #                    "_decline",
                #                    "['decline']"),
                # rowname = str_replace(rowname,
                #                    "_greenup",
                #                    "['greenup']"),
                #name = if_else(PC == "PC2", rowname, name),
                ) %>% 
  mutate(index = str_split_i(rowname, "_", 1),
         season = str_split_i(rowname, "_", -1)) %>% 
  dplyr::group_by(index, season) %>%
  dplyr::mutate(max_PC1 = max(abs(PC1)),
                max_PC2 = max(abs(PC2)),
                max_PC3 = max(abs(PC3)),
           keep_PC1 = if_else(abs(PC1) == max_PC1, 1, 0),
           keep_PC2 = if_else(abs(PC2) == max_PC2, 1, 0),
           keep_PC3 = if_else(abs(PC3) == max_PC3, 1, 0)) %>% 
  dplyr::mutate(PC1_h = if_else(PC1 > 0 & name != "EWI['9']^'green-up'", .1, .8),
         PC2_v = if_else(PC2 > 0 & name != "BCC^'late winter'" & name != "ARI^'decline'", .1, 1),
         PC2_v = if_else(name == "RE['upper']^'late winter'",
                         0, PC2_v),
         PC2_h = if_else(PC2 > 0, .1, .8),
         PC2_h = if_else(name %in% c("mDatt^'late summer'", "NDRE['740']^'decline'"), .45, PC2_h),
         PC3_v = if_else(PC3 > 0 & name != "RE['upper']^'late summer'", 0, 1),
         name = paste0("bold(", name, ")")) %>% 
  dplyr::mutate(PC1_h = PC1_h - (PC1/2),
         PC2_v = PC2_v - PC2,
         PC3_v = PC3_v - PC3,
         PC2_h = PC2_h - (PC2/2))
  # distinct() %>% 
  # dplyr::group_by(PC) %>% 
  # dplyr::mutate(rank = min_rank(desc(loading_abs)))

pc1_range = (max(pca_plot$PC1) - min(pca_plot$PC1)) * .03
pc2_range = (max(pca_plot$PC2) - min(pca_plot$PC2)) * .03
pc3_range = (max(pca_plot$PC3) - min(pca_plot$PC3)) * .03

mult = 3
scale = .5

# (pc_plot_1 = ggplot() +
#     geom_point(data = pca_plot, aes(x = PC1, y = PC2, color = cluster), size = 4.5, alpha = .7) +
#     geom_segment(data=dplyr::filter(loadings, loading_abs > .2 & PC == "PC1"), 
#                  aes(x=0, 
#                      y= min(pca_plot$PC2) - (rank * pc2_range) - (pc2_range * mult) - (sign_off * (pc2_range/2)), 
#                      xend = loading * (max(pca_plot$PC1) - min(pca_plot$PC1)) * scale, 
#                      yend = min(pca_plot$PC2) - (rank * pc2_range) - (pc2_range * mult) - (sign_off * (pc2_range/2))), 
#                  color="red4") +
#     geom_text(data=dplyr::filter(loadings, loading_abs > .2 & PC == "PC1"), 
#                  aes(y= min(pca_plot$PC2) - (rank * pc2_range) - (pc2_range * mult) - (sign_off * (pc2_range/2)),  
#                      x = loading * (max(pca_plot$PC1) - min(pca_plot$PC1)) * scale,
#                      label = name,
#                      hjust = h_just),
#               size = 4.5, parse = TRUE) +
#     geom_segment(data=dplyr::filter(loadings, loading_abs > .2 & PC == "PC2"), 
#                  aes(y=0, 
#                      x= min(pca_plot$PC1) - (rank * pc1_range) - (pc1_range * mult) - (sign_off * (pc1_range/2)), 
#                      yend = loading * (max(pca_plot$PC2) - min(pca_plot$PC2)) * scale, 
#                      xend = min(pca_plot$PC1) - (rank * pc1_range) - (pc1_range * mult) - (sign_off * (pc1_range/2))), 
#                  color="red4") +
#     geom_text(data=dplyr::filter(loadings, loading_abs > .2 & PC == "PC2"), 
#               aes(x= min(pca_plot$PC1) - (rank * pc1_range) - (pc1_range * mult) - (sign_off * (pc1_range/2)),  
#                   y = loading * (max(pca_plot$PC2) - min(pca_plot$PC2)) * scale,
#                   label = name,
#                   hjust = h_just),
#               size = 4.5, parse = TRUE, angle = 90) +
#     theme_bw(20) +
#   labs(x = "PC1 — greenness (56.1%)", y = "PC2 — greenup (20.6%)") +
#     scale_color_manual(values = clust_cols, name = "Seedlot Class") +
#     scale_fill_manual(values = c("#666666", "#AAAAAA", "#333333"),
#                       labels = c("P. engelmannii", "P. glauca", "P. sitchensis")) +
#     scale_y_continuous(expand = expansion(mult = c(.05, .05))) +
#     scale_x_continuous(expand = expansion(mult = c(.05, .05))) +
#     coord_equal(clip = "off",
#                 ylim = c(min(pca_plot$PC2), max(pca_plot$PC2)),
#                 xlim = c(min(pca_plot$PC1), max(pca_plot$PC1))) +
#     theme(strip.text = element_blank(),
#           plot.background = element_rect(fill = NA, color = "black", size = 1),
#           aspect.ratio = 1,
#           axis.text = element_text(size = 10),
#           axis.title.y = element_text(margin = margin(t = 0, r = 90, b = 0, l = 0)),
#           axis.title.x = element_text(margin = margin(t = 90, r = 0, b = 0, l = 0)),
#           legend.position = "none"))

pca_plot_cluster = pca_plot %>% 
  dplyr::select(PC1:PC3, cluster) %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::summarise_all("mean") %>% 
  dplyr::mutate(PC2_y = if_else(cluster %in% c(7, 6), PC2 - .18, PC2),
                PC2_y = if_else(cluster %in% c(2), PC2 - .22, PC2_y),
                PC2_y = if_else(cluster %in% c(8, 5), PC2 + .04, PC2_y),
                PC1_x = if_else(cluster %in% c(1), PC1 + .15, PC1),
                PC1_x = if_else(cluster %in% c(5), PC1 - .11, PC1_x)) %>% 
  dplyr::mutate(PC3_y = if_else(cluster %in% c(7, 8), PC3 - .1, PC3),
                PC3_y = if_else(cluster %in% c(5), PC3 + .15, PC3_y),
                PC3_y = if_else(cluster %in% c(2), PC3 + .01, PC3_y),
                PC3_y = if_else(cluster %in% c(3), PC3 - .03, PC3_y),
                PC2_x = if_else(cluster %in% c(3, 10), PC2 - .05, PC2),
                PC2_x = if_else(cluster %in% c(2), PC2 + .03, PC2_x),
                PC2_x = if_else(cluster %in% c(1, 7, 8), PC2 + .12, PC2_x))

(pc_plot_1 = ggplot() +
    geom_point(data = pca_plot,
               aes(x = PC1, y = PC2, color = cluster, alpha = cluster), size = 5) +
    theme_bw(20) +
    labs(x = "PC1 — greenness (48.9%)", y = "PC2 — green-up (15.6%)") +
    scale_color_manual(values = clust_cols, name = "Seedlot Class") +
    scale_fill_manual(values = c("#666666", "#AAAAAA", "#333333"),
                      labels = c("P. engelmannii", "P. glauca", "P. sitchensis")) +
    scale_y_continuous(breaks = seq(-2, 1, 1),
                       labels = c("-2.0", "-1.0", "0.0", "1.0"),
                       expand = expansion(mult = c(.08, .05))) +
    scale_x_continuous(breaks = seq(-2, 1, 1),
                       labels = c("-2.0", "-1.0", "0.0", "1.0"),
                       expand = expansion(mult = c(.05, .05))) +
    #scale_alpha_manual(values = c(.8, .8, 0, 0, 0, 0, 0, 0)) +
    scale_alpha_manual(values = c(.8, .8, .8, .8, .8, .8, .8, .8)) +
    ggnewscale::new_scale_color() +
    geom_point(data = pca_plot_cluster, aes(x = PC1_x, y = PC2_y, color = cluster),
               shape = 16, size = 10.5) +
    scale_color_manual(values = clust_cols) +
    geom_text(data = pca_plot_cluster, aes(x = PC1_x, y = PC2_y, label = cluster),
              color = "black", fontface = "bold", size = 6) +
    annotate(geom = 'text',
             label = '   FOLIAR BIOMASS →\n   HEIGHT GROWTH →\n   PHOTOSYNTHETIC DOWNREGULATION →\n   ',
             x = -Inf, y = -Inf, hjust = 0, vjust = 0.1, fontface = "bold", color = "grey30", size = 4) +
    annotate(geom = 'text',
             label = '← NITROGEN STORAGE   \n← NEEDLE MATURATION   \n\n',
             x = Inf, y = Inf, hjust = 1, vjust = .4, angle = 90, fontface = "bold", color = "grey30", size = 4) +
    coord_equal(clip = "off",
                ylim = c(min(pca_plot$PC2), max(pca_plot$PC2)),
                xlim = c(min(pca_plot$PC1), max(pca_plot$PC1))) +
    theme(strip.text = element_blank(),
          panel.border = element_rect(linewidth = .5),
          aspect.ratio = 1,
          axis.text = element_text(size = 13),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
          legend.position = "none"))

# ggsave(plot = pc_plot_1,
#        filename = 'D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated\\PC1_2_1&2.jpeg',
#        device = jpeg,
#        width = 6,
#        height = 6,
#        units = 'in',
#        dpi = 300,
#        bg = 'white')

(pc_loadings_1 = loadings %>% 
    dplyr::filter((abs(PC1) > .195 & keep_PC1 == 1)
                  | (abs(PC2) > .21 & keep_PC2 == 1)) %>% 
    ggplot() +
    geom_segment(aes(x=0, 
                     y= 0,
                     xend = PC1,
                     yend = PC2),
                 color="red4") +
    geom_text(aes(x = PC1, 
                  y = PC2,
                  label = name,
                  hjust = PC1_h,
                  vjust = PC2_v),
              size = 4.5, parse = TRUE) +
    theme_bw(20) +
    labs(x = "PC1 — greenness (48.9%)", y = "PC2 — green-up (15.6%)") +
    scale_color_manual(values = clust_cols, name = "Seedlot Class") +
    scale_fill_manual(values = c("#666666", "#AAAAAA", "#333333"),
                      labels = c("P. engelmannii", "P. glauca", "P. sitchensis")) +
    scale_y_continuous(breaks = seq(-.4, .4, .2),
                       expand = expansion(mult = c(.15, .15))) +
    scale_x_continuous(breaks = seq(-.4, .4, .2),
                       expand = expansion(mult = c(.28, .28))) +
    coord_equal(clip = "off") +
    theme(strip.text = element_blank(),
          panel.border = element_rect(linewidth = .5),
          aspect.ratio = 1,
          axis.text = element_text(size = 13),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
          legend.position = "none"))

# S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
# IS_sqrt <- function(x){x^2*sign(x)}
# S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

(pc_plot_2 = ggplot() +
    geom_point(data = pca_plot,
               aes(x = PC2, y = PC3, color = cluster,
                   alpha = cluster), size = 5) +
    theme_bw(20) +
    labs(x = "PC2 — green-up (15.6%)", y = "PC3 — red edge stress (13.1%)") +
    scale_color_manual(values = clust_cols, 
                       name = "Seedlot Class") +
    scale_alpha_manual(values = c(.8, .8, .8, .8, .8, .8, .8, .8)) +
    # scale_fill_manual(values = c("#666666", "#AAAAAA", "#333333"),
    #                   labels = c("P. engelmannii", "P. glauca", "P. sitchensis")) +
    # scale_y_continuous(breaks = c((-2 + 2.05)^2, (-1 + 2.05)^2, (0 + 2.05)^2, (1 + 2.05)^2),
    #                    minor_breaks = c((-1.5 + 2.05)^2, (-.5 + 2.05)^2, (.5 + 2.05)^2, (1.5 + 2.05)^2),
    #                    labels = c("-2.0", "-1.0", "0.0", "1.0"),
    #                    expand = expansion(mult = c(.05, .09))) +
    scale_x_continuous(breaks = seq(-2, 1, 1),
                       labels = c("-2.0", "-1.0", "0.0", "1.0"),
                       expand = expansion(mult = c(.05, .05))) +
    ggnewscale::new_scale_color() +
    geom_point(data = pca_plot_cluster,
               aes(x = PC2_x, y = PC3_y, color = cluster),
               shape = 16, size = 10.5) +
    scale_color_manual(values = clust_cols) +
    geom_text(data = pca_plot_cluster,
              aes(x = PC2_x, y = PC3_y, label = cluster),
              color = "black", fontface = "bold", size = 6) +
    # annotate(geom = 'text', 
    #          label = '  ← NITROGEN STORAGE   \n  ← NEEDLE MATURATION   \n\n',
    #          x = -Inf, y = -Inf, hjust = 0, vjust = 0, fontface = "bold", color = "grey30", size = 4) +
    # annotate(geom = 'text', 
    #          label = 'DROUGHT AVOIDANCE →   \n← ANTHOCYANIN STRESS RESPONSE   \n\n',
    #          x = Inf, y = Inf, hjust = 1, vjust = 0, angle = 90, fontface = "bold", color = "grey30", size = 4) +
    # coord_equal(clip = "off",
    #             ylim = c((min(pca_plot$PC3) + 2.05)^2, (max(pca_plot$PC3) + 2.05)^2),
    #             xlim = c(min(pca_plot$PC2), max(pca_plot$PC2))) +
    theme(strip.text = element_blank(),
          panel.border = element_rect(linewidth = .5),
          aspect.ratio = 1,
          axis.text = element_text(size = 13),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
          legend.position = "none"))

# ggsave(plot = pc_plot_2,
#        filename = 'D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated\\PC2_3_7&8.jpeg',
#        device = jpeg,
#        width = 6,
#        height = 6,
#        units = 'in',
#        dpi = 300,
#        bg = 'white')

(pc_loadings_2 = loadings %>% 
    dplyr::filter((abs(PC2) > .21 & keep_PC2 == 1)
                  | (abs(PC3) > .21 & keep_PC3 == 1)) %>% 
    ggplot() +
    geom_segment(aes(x=0, 
                     y= 0,
                     xend = PC2,
                     yend = PC3),
                 color="red4") +
    geom_text(aes(x = PC2, 
                  y = PC3,
                  label = name,
                  hjust = PC2_h,
                  vjust = PC3_v),
              size = 4.5, parse = TRUE) +
    theme_bw(20) +
    labs(x = "PC2 — green-up (15.6%)", y = "PC3 — red edge stress (13.1%)") +
    scale_color_manual(values = clust_cols, name = "Seedlot Class") +
    scale_fill_manual(values = c("#666666", "#AAAAAA", "#333333"),
                      labels = c("P. engelmannii", "P. glauca", "P. sitchensis")) +
    scale_y_continuous(breaks = seq(-.4, .4, .2),
                       expand = expansion(mult = c(.15, .15))) +
    scale_x_continuous(breaks = seq(-.4, .4, .2),
                       expand = expansion(mult = c(.28, .28))) +
    coord_equal(clip = "off") +
    # annotate(geom = 'text', 
    #          label = '   ← GREATER NITROGEN STORAGE?\n   ← FASTER NEEDLE MATURATION?\n',
    #          x = -Inf, y = -Inf, hjust = 0, vjust = 0, fontface = "bold", color = "grey30") +
    theme(strip.text = element_blank(),
          panel.border = element_rect(linewidth = .5),
          aspect.ratio = 1,
          axis.text = element_text(size = 13),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
          legend.position = "none"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rect = rectGrob(
  x = unit(1.17, "in"),
  y = unit(1, "npc") - unit(.32, "in"),
  width = unit(.6, "in"),
  height = unit(.6, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

font_s = 30

lab_a = textGrob(
  label = "(a)",
  x = unit(1.21, "in"),
  y = unit(1, "npc") - unit(.41, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))


lab_b = textGrob(
  label = "(b)",
  x = unit(1.21, "in"),
  y = unit(1, "npc") - unit(.41, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_c = textGrob(
  label = "(c)",
  x = unit(1.21, "in"),
  y = unit(1, "npc") - unit(.41, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_d = textGrob(
  label = "(d)",
  x = unit(1.21, "in"),
  y = unit(1, "npc") - unit(.41, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))


pc_plots = ggarrange((ggdraw(pc_plot_1) + 
                       draw_grob(rect) +
                       draw_grob(lab_a)),
                     
                     (ggdraw(pc_loadings_1) + 
                        draw_grob(rect) +
                        draw_grob(lab_b)),
                     
                     (ggdraw(pc_plot_2) + 
                        draw_grob(rect) +
                        draw_grob(lab_c)),
                     
                     (ggdraw(pc_loadings_2) + 
                        draw_grob(rect) +
                        draw_grob(lab_d)),

                        ncol = 2, nrow = 2,
                     heights = c(1,1,1,1))

# pc_all = ggarrange(pc_plots, pc_loadings, nrow = 1,
#                    widths = c(1, 1.01))


ggsave(plot = pc_plots,
       filename = 'D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\figure_8.jpeg',
       device = jpeg,
       width = 15,
       height = 15,
       units = 'in',
       dpi = 300,
       bg = 'white')

################################################################################
library(PieGlyph)

pca_plot_species2 = pca_plot_species

(pc_species_1 = ggplot() +
   # geom_scatterpie(data = pca_plot_species, aes(x = PC1, y = PC2, r = PC1 / 30), 
   #                 color = NA,
   #                 cols = c("propGla", "propEng", "propSit")) +
    geom_pie_glyph(data = filter(pca_plot_species),# propEng != 1 & propGla != 1 & propSit != 1), 
                   aes(x = PC1, y = PC2),
                   colour = NA,#"black",
                   radius = .27,
                   linewidth = .01,
                   slices = c( "propGla", "propEng", "propSit")) +
    # geom_pie_glyph(data = filter(pca_plot_species, propEng == 1 | propGla == 1 | propSit == 1),
    #                aes(x = PC1, y = PC2),
    #                colour = NA,
    #            radius = .29,
    #            slices = c("propEng", "propGla", "propSit")) +
    #scale_radius_discrete(range=c(3,3), unit="cm") +
   theme_bw(20) +
    labs(x = "PC1 — greenness (48.9%)", y = "PC2 — green-up (15.6%)") +
   scale_color_manual(values = clust_cols, name = "Seedlot Class") +
    scale_fill_manual(values = cols_species,
                      labels = c("P. glauca", "P. engelmannii", "P. sitchensis"),
                      name = "Species") +
   scale_y_continuous(expand = expansion(mult = c(.05, .05)),
                      breaks = c(-1, 0, 1),
                      labels = c("-1.0", "0.0", "1.0")) +
   scale_x_continuous(expand = expansion(mult = c(.05, .05)),
                      breaks = c(-2, -1, 0, 1, 2),
                      labels = c("-2.0", "-1.0", "0.0", "1.0", "2.0")) +
   coord_equal(clip = "off",
               ylim = c(min(pca_plot$PC2), max(pca_plot$PC2)),
               xlim = c(min(pca_plot$PC1), max(pca_plot$PC1))) +
    theme(strip.text = element_blank(),
          panel.border = element_rect(linewidth = .5),
          axis.text = element_text(size = 13),
          aspect.ratio = 1,
          plot.margin = unit(c(0, 0, 0, 0), 
                             "cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 2)),
          axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 2, l = 0)),
          legend.position = "none"))


# pca_plot_species2 = pca_plot_species %>% 
#   dplyr::mutate(PC3 = (PC3 + 2.05)^2)

(pc_species_2 = ggplot() +
    # geom_scatterpie(data = pca_plot_species, aes(x = PC1, y = PC2, r = PC1 / 30), 
    #                 color = NA,
    #                 cols = c("propGla", "propEng", "propSit")) +
    geom_pie_glyph(data = pca_plot_species2,
                   aes(x = PC2, y = PC3),
                   colour = NA,
                   radius = .27,
                   slices = c( "propGla", "propEng", "propSit")) +
    #scale_radius_discrete(range=c(3,3), unit="cm") +
    theme_bw(20) +
    labs(x = "PC2 — green-up (15.6%)", y = "PC3 — red edge stress (13.1%)") +
    scale_color_manual(values = clust_cols, name = "Seedlot Class") +
    scale_fill_manual(values = cols_species,
                      labels = c("P. glauca", "P. engelmannii", "P. sitchensis"),
                      name = "Species") +
    # scale_y_continuous(breaks = c((-2 + 2.05)^2, (-1 + 2.05)^2, (0 + 2.05)^2, (1 + 2.05)^2),
    #                    minor_breaks = c((-1.5 + 2.05)^2, (-.5 + 2.05)^2, (.5 + 2.05)^2, (1.5 + 2.05)^2),
    #                    labels = c("-2.0", "-1.0", "0.0", "1.0"),
    #              expand = expansion(mult = c(.05, .09))) +
    scale_x_continuous(expand = expansion(mult = c(.05, .05)),
                       breaks = c(-1, 0, 1),
                       labels = c("-1.0", "0.0", "1.0")) +
    coord_equal(clip = "off",
                ylim = c(min(pca_plot_species2$PC3), max(pca_plot_species2$PC3)),
                xlim = c(min(pca_plot_species2$PC2), max(pca_plot_species2$PC2))) +
    theme(strip.text = element_blank(),
          panel.border = element_rect(linewidth = .5),
          axis.text = element_text(size = 13),
          aspect.ratio = 1,
          plot.margin = unit(c(0, 0, 0, 0), 
                             "cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 2)),
          axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 2, l = 0)),
          legend.position = "none"))


pops_plot_species_sf = pca_plot_species %>% 
  left_join(pops_plot_facet %>% 
              dplyr::select(Prov, geometry)) %>% 
  st_as_sf() %>% 
  st_transform(crs = 3005) %>% 
  bind_cols(st_coordinates(.))

saveRDS(pops_plot_species_sf, "D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_cluster_species_sf.rds")

# proj = "+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +gamma=12"
# 
# pops_sf_new = pca_plot_species %>% 
#   left_join(pops_plot_facet %>% 
#               dplyr::select(Prov, geometry)) %>% 
#   st_as_sf() %>% 
#   st_transform(crs = proj) %>% 
#   bind_cols(st_coordinates(.))
# 
# 
# #crs(pops_plot_species_sf, describe=TRUE, proj=TRUE)
# # rast_df2 = rast_df1 %>% 
# #   project(crs = proj)
# 
# library(elevatr)
# elevation_data1 = get_elev_raster(locations = countries, z = 4) %>% 
#   terra::crop((st_bbox(pops_plot_species_sf) + c(xmin = -7, ymin = -2, xmax = -13, ymax = 1))) %>% 
#   projectRaster(crs = proj)
# 
# rast_df1 = as.data.frame(elevation_data1, xy = TRUE)
# names(rast_df1)[3] <- "elev"
# 
# saveRDS(rast_df1, "D:\\Sync\\_Sites\\Sx_Genecology\\rast_df_rotate")

states_crop = states %>% 
  st_transform(crs = st_crs(3005)) %>% 
  dplyr::filter(!(name %in% c("Alaska", "Oregon", "Saskatchewan", "Nunavut"))) %>% 
  dplyr::mutate(#name = str_replace(name, " ", "\n"),
                name = str_replace(name, "Yukon", "Yukon                    \n\n"),
                name = str_replace(name, "Northwest Territories", "Northwest\nTerritories"),
                name = str_replace(name, "British Columbia", "\n       British\n       Columbia")) %>% 
  st_crop(st_bbox(pops_plot_species_sf %>% 
                    dplyr::filter(!(Prov %in% c(404:407)))) +
            c(+220000, -3000000, 40000, 50000))

states_crop_inset = states %>% 
  st_transform(crs = st_crs(3005)) %>% 
  dplyr::filter(name %in% c("Arizona", "New Mexico")) %>% 
  dplyr::mutate(name = str_replace(name, " ", "\n")) %>% 
  st_crop(st_bbox(pops_plot_species_sf %>% 
                    dplyr::filter(Prov %in% c(404:407))) +
            c(0, -20000, 0, -20000))

(map_species = ggplot() +
  geom_raster(data = rast_df1,  aes(x = x, y = y, fill = elev)) +
  scale_fill_gradient(high = "white", low = "grey20", guide = "none") +
  geom_sf(data = countries, fill = NA, color = "grey20") +
  geom_sf(data = states, fill = NA, color = "grey20") +
  geom_sf_text(data = states_crop, aes(label = name), 
               size = 5,
               angle = 0,
               color = "grey20") +
    geom_sf_text(data = states_crop_inset, aes(label = name), 
                 size = 5,
                 angle = 0,
                 color = "grey20") +
  geom_sf(data = ocean, fill = "white", color = "grey20") + 
    # annotate("text", x = 1900000, y = 1490000, label = "C    A    N    A    D    A",
    #          size = 5.5, color = "grey36", angle = 13) +
    # annotate("text", x = 2170000, y = -284500, label = "U  N  I  T  E  D    S  T  A  T  E  S\no f\nA  M  E  R  I  C  A",
    #          size = 4.1, color = "grey36", angle = 11) +
    annotate("text", x = 550000, y = 605000, label = "P   A   C   I   F   I   C\n\nO   C   E   A   N",
             size = 4, color = "grey36", angle = 0) +
  ggnewscale::new_scale_fill() +
    geom_pie_glyph(data = pops_plot_species_sf, aes(x = X, y = Y),
                   colour = NA,#"black",
                   radius = .3,
                   slices = c( "propGla", "propEng", "propSit")) +
  #geom_sf(data = pops_plot_facet, aes(color = factor(cluster)), shape = 16, size = 1.5) +
  #geom_sf(data = cluster_centroids, aes(color = factor(cluster)), shape = 19, size = 7) +
    scale_fill_manual(values = cols_species,
                      labels = c("P. glauca", "P. engelmannii", "P. sitchensis"),
                      name = "Species") +
  #geom_sf_text(data = cluster_centroids, aes(label = cluster), color = "black", fontface = "bold", size = 3) +
  coord_sf(xlim = c(375000, 2090000), ylim = c(175000,#-768333,
                                               2190000),
           crs = st_crs(3005)) +
  theme(panel.border = element_rect(fill = NA, color = "grey25", size = 1),
        legend.position = c(0.87, 0.91),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.key = element_blank(),
        legend.spacing.y = unit(.2, 'cm'),
        legend.text = element_text(size = 14, face = "italic"),
        legend.title = element_text(size = 18),
        legend.margin = margin(8, 8, 14, 8),
        plot.margin = unit(c(.25,.25,.25,.25), 'lines'),
        panel.background = element_rect(fill = "white", color = "black", linewidth = .5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank()))

#x
print((2090000 - 375000) * .4) # 686000

#y
print((2190000 - 175000) * .18) # 362700

map_inset = map_species +
  coord_sf(xlim = c(2250000, (2250000 + 686000)), 
           ylim = c(-898333, (-898333 + 362700)),
           crs = st_crs(3005)) +
    theme(legend.position = "none")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rect = rectGrob(
  x = unit(.9, "in"),
  y = unit(1, "npc") - unit(.17, "in"),
  width = unit(.5, "in"),
  height = unit(.5, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

rect2 = rectGrob(
  x = unit(.25, "in"),
  y = unit(1, "npc") - unit(.17, "in"),
  width = unit(.5, "in"),
  height = unit(.5, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

rect3 = rectGrob(
  x = unit(.01, "npc") + unit(.36, "in"),
  y = unit(.17, "npc") - unit(.18, "in"),
  width = unit(.5, "in"),
  height = unit(.5, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

font_s = 25

lab_a = textGrob(
  label = "(a)",
  x = unit(.93, "in"),
  y = unit(1, "npc") - unit(.25, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))


lab_b = textGrob(
  label = "(b)",
  x = unit(.93, "in"),
  y = unit(1, "npc") - unit(.25, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_c = textGrob(
  label = "(c)",
  x = unit(.3, "in"),
  y = unit(1, "npc") - unit(.25, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_d = textGrob(
  label = "(d)",
  x = unit(.01, "npc") + unit(.4, "in"),
  y = unit(.17, "npc") - unit(.27, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

map_all = ggdraw(map_species) +
    draw_grob(rect2) +
  #draw_grob(map_species) +
  draw_grob(as_grob(map_inset), hjust = 0, vjust = 1,
            y = .17, x = .01, width = .35, height = .15) +
    draw_grob(rect3) +
  draw_grob(lab_c) +
  draw_grob(lab_d)


leg_species = get_legend(map_species)

pc_species = ggarrange(ggdraw(pc_species_1) + 
                         draw_grob(rect) +
                         draw_grob(lab_a), 
                       ggdraw(pc_species_2) +
                         draw_grob(rect) +
                         draw_grob(lab_b),
                       vjust = 1,
                       nrow = 2,
                       widths = c(1,1))

# species_figure = ggarrange(pc_species, map_species,
#                            ncol = 2,
#                            widths = c(1,1.75))

species_figure = ggarrange(pc_species, map_all,
                           vjust = 1,
                           ncol = 2,
                           widths = c(1, 1.72))

ggsave(plot = species_figure,
       filename = 'D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\figure_9.jpeg',
       device = jpeg,
       width = 15,
       height = 11,
       units = 'in',
       dpi = 300,
       bg = 'white')


# Presentation

species_prez = ggarrange(pc_species_1, NULL, pc_species_2,
                           ncol = 3,
                           widths = c(1,.07, 1))

ggsave(plot = species_prez,
       filename = 'D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated\\PC_species.jpeg',
       device = jpeg,
       width = 16,
       height = 9,
       units = 'in',
       dpi = 300,
       bg = 'white')

################################################################################


############### traits unscaled


plot_band = 
  plot_change = mm_blups_change %>% 
  filter(index %in% R_list) %>% 
  dplyr::mutate(#Class = factor(Class, levels = c("Wildstand", "Orchard")),
    EST = timepoint_value) %>% 
  distinct(index, EST, timepoint, Prov) %>% 
  # dplyr::filter(Prov %in% c(404, 411, 338, 403) &
  #                 timepoint == "mid summer") %>% 
  left_join(pca_plot %>% dplyr::select(Prov, cluster))



plot_traits = readRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\BLUPS\\Updated\\df_cluster.rds") %>% 
  mutate(cat = if_else(timepoint %in% c("early summer", "mid summer", "late summer", "late winter"), "seasonal", timepoint)) %>% 
  group_by(index, timepoint) %>% 
  #dplyr::filter(index %in% keep_list) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-cat) %>% 
  unite("trait", c(index, timepoint)) %>% 
  # filter(trait %in% index_keep_list$trait
  #        #        #& !(trait %in% timepoint_drop_list$trait)
  # ) %>%
  pivot_wider(names_from = trait, values_from = ESTIMATE) %>% 
  drop_na() %>% 
  left_join(pops_plot_species_sf %>% dplyr::select(Prov, PC1:PC3, propGla, propEng, propSit, cluster, geometry)) %>% 
  st_as_sf()

# PCs
(pc_map = ggplot() +
    #geom_raster(data = rast_df1,  aes(x = x, y = y, fill = elev)) +
    #scale_fill_gradient(high = "white", low = "grey20", guide = "none") +
    geom_sf(data = countries, fill = "grey95", color = "grey20") +
    geom_sf(data = states, fill = NA) +
    geom_sf(data = ocean, fill = "white", color = "grey20") + 
    geom_sf(data = plot_traits %>% 
              dplyr::select(PC1:PC3) %>% 
              pivot_longer(PC1:PC3, names_to = "PC",
                           values_to = "PC_val"),
            aes(color = PC_val,
                                     # `NDVI_decline`
                                    ), shape = 16, size = 4) +
    #geom_sf(data = cluster_centroids, aes(color = factor(cluster)), shape = 19, size = 9) +
    # annotate("text", x = 1900000, y = 1490000, label = "C    A    N    A    D    A",
    #          size = 5.5, color = "grey36", angle = 13) +
    # annotate("text", x = 2170000, y = -354500, label = "U  N  I  T  E  D    S  T  A  T  E  S\no f\nA  M  E  R  I  C  A",
    #          size = 4.1, color = "grey36", angle = 11) +
    # annotate("text", x = 650000, y = -13000, label = "P   A   C   I   F   I   C\n\nO   C   E   A   N",
    #          size = 3, color = "grey36", angle = 0) +
    # annotate("text", vjust = 1,  x = 1966000, y = 494000, label = "   R  O  C  K  Y           M  O  U  N  T  A  I  N  S",
    #          size = 2.7, color = "grey36", angle = -45) +
    # annotate("text", vjust = 1,  x = 890000, y = 1100000, label = " C  O  A  S  T         M  O  U  N  T  A  I  N  S    ",
    #          size = 2.7, color = "grey36", angle = -61) +
  #scale_color_manual(values = clust_cols) +
  scale_color_viridis_c(option = "B") +
    #geom_sf_text(data = cluster_centroids, aes(label = cluster), color = "black", fontface = "bold", size = 4.5) +
    coord_sf(xlim = c(225000, 2800000), ylim = c(-768333, 2150000),
             crs = st_crs(3005)) +
    theme(panel.border = element_rect(fill = NA, color = "grey25", size = .5),
          #legend.position = "below",#c(0.94, 0.5),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 24),
          legend.margin = margin(8, 8, 14, 8),
          panel.background = element_rect(fill = "white", color = "black", linewidth = .5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          #strip.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_wrap(. ~ PC))

ggsave(plot = pc_map,
       filename = 'D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated\\SI\\PC_map.jpeg',
       device = jpeg,
       width = 15,
       height = 7,
       units = 'in',
       dpi = 300,
       bg = 'white')

################################################################################

rect = rectGrob(
  x = unit(.19, "in"),
  y = unit(1, "npc") - unit(.45, "in"),
  width = unit(.37, "in"),
  height = unit(.37, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

font_s = 17.5

lab_a = textGrob(
  label = "(a)",
  x = unit(.22, "in"),
  y = unit(1, "npc") - unit(.53, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))


lab_b = textGrob(
  label = "(b)",
  x = unit(.22, "in"),
  y = unit(1, "npc") - unit(.53, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_c = textGrob(
  label = "(c)",
  x = unit(.22, "in"),
  y = unit(1, "npc") - unit(.53, "in"),
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = font_s))

trait = plot_traits$`GCC_late winter`
index_title = expression(GCC^"late winter")
timepoint_ = "late winter"
bands = #c( "R668", 
           #"R717", "R740")
         # "R705", 
          #"R842")
  #c("R560", "R705")
c("R444", "R475", "R531", "R560", "R650", "R668")

# trait = plot_traits$`ARI_mid summer`
# index_title = expression(ARI^"mid summer")
# timepoint_ = "mid summer"
# bands = #c( "R668", "R705", "R717", "R740", "R842")
#   c("R560", "R705")
#   #c("R444", "R475", "R531", "R560", "R650", "R668", "R705")

# trait = plot_traits$`Datt_green-up`
# index_title = expression(mDatt^"green-up")
# timepoint_ = "green-up"
# bands = c( "R668", "R705", "R842")
#   #c("R560", "R705")
# #c("R444", "R475", "R531", "R560", "R650", "R668", "R705")

# trait = plot_traits$`GCC_green-up`
# index_title = expression(GCC^"green-up")
# timepoint_ = "green-up"
# bands = #c( "R668", "R705", "R717", "R740", "R842")
#   #c("R560", "R705")
#   c("R444", "R475", "R531", "R560", "R650", "R668")

# 
# (index_plot_species = plot_traits %>%
#     dplyr::mutate(cluster = reorder(cluster, `ARI_mid summer`)) %>%
#     ggplot() +
#     geom_pie_glyph(aes(x = cluster, y = `ARI_mid summer`),
#                    slices = c("propGla", "propEng", "propSit"),
#                    position = position_dodge(width = 2)) +
#     scale_fill_manual(values = c("#888888", "#555555", "#111111"),
#                       labels = c("P. glauca", "P. engelmannii", "P. sitchensis"),
#                       name = "Species") +
#     theme_bw(base_size = 22) +
#     labs(x = "Spectral band", y = index_title) +
#     #scale_y_discrete(limits = rev(levels(blups_r2_plot$clim))) +
#     #scale_fill_manual(values = clust_cols) +
#     #scale_y_continuous(expand = expansion(mult = c(.08, .18))) +
#     theme(#aspect.ratio = 1,
#       legend.position = "right",#c(.24, .79),
#       legend.text = element_text(size = 20),
#       axis.text.x = element_text(angle = 90, vjust = .5),
#       axis.title.x = element_blank(),
#       legend.background = element_rect(color = "black", linewidth = .5),
#       legend.title = element_blank()))



(index_map = ggplot() +
    #geom_raster(data = rast_df1,  aes(x = x, y = y, fill = elev)) +
    #scale_fill_gradient(high = "white", low = "grey20", guide = "none") +
    geom_sf(data = countries, fill = "grey95", color = "grey20") +
    geom_sf(data = states, fill = NA) +
    geom_sf(data = ocean, fill = "white", color = "grey20") + 
    geom_sf(data = plot_traits,
            aes(color = trait
                # `NDVI_decline`
            ), shape = 16, size = 3) +
  scale_color_gradient2(low = "steelblue", 
                        mid = "tan", 
                        high = "red3",
                        midpoint = max(trait) - ((max(trait) - min(trait))/2),
                       # breaks = seq(-100, 100, .5),
                        name = index_title) +
    coord_sf(xlim = c(225000, 2800000), ylim = c(-768333, 2150000),
             crs = st_crs(3005)) +
    theme(panel.border = element_rect(fill = NA, color = "grey25", size = .8),
          legend.position = c(0.14, 0.2),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.margin = margin(5, 5, 7, 5),
          panel.background = element_rect(fill = "white", color = "black", linewidth = .5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          #strip.text = element_blank(),
          axis.ticks = element_blank()))

(band_plot = plot_band %>%
    filter(index %in% bands) %>% 
    #filter(cluster %in% c(1, 3, 6)) %>% 
    #pivot_wider(names_from = index, values_from = EST) %>% 
    dplyr::filter(timepoint  == timepoint_) %>% 
    ggplot(aes(color = cluster,
              fill = cluster,
              group = cluster,
              x = index,
              y = EST
              #group = factor(Class, levels = c("Wildstand", "Orchard"))
   )) +
    stat_summary(geom = "line", linewidth = .6) +
    stat_summary(geom = "point", size = 2, alpha = .8) +
   theme_bw(base_size = 22) +
   labs(x = "Spectral band", y = "Reflectance") +
   #scale_y_discrete(limits = rev(levels(blups_r2_plot$clim))) +
   scale_fill_manual(values = clust_cols) +
   scale_color_manual(values = clust_cols) +
   #scale_y_continuous(expand = expansion(mult = c(.08, .18))) +
    theme(#aspect.ratio = 1.25,
          legend.position = "none",#c(.2, .7),
          legend.text = element_text(size = 10, hjust = 0),
          legend.margin = margin(c(-5, 4, 4, 4)),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14, margin = margin(r = 20, l = 0)),
          axis.text.x = element_text(vjust = .5, size = 12),
          axis.text.y = element_text(size = 12),
          #axis.title.x = element_blank(),
          legend.background = element_rect(color = "black", linewidth = .5),
          legend.title = element_blank()))


(index_plot = plot_traits %>%
    dplyr::mutate(cluster = reorder(cluster, trait)) %>% 
    ggplot(aes(x = cluster, y = trait, fill = cluster)) +
    geom_boxplot(aes(fill = cluster),
                 outlier.alpha = 0, width = .5, alpha = .7,
                 position = position_dodge(width = .75),
                 linewidth = .3, color = "black") +
    geom_point(size = 3, alpha = .8, shape = 21) +
    #scale_color_manual(values = clust_cols) +
    theme_bw(base_size = 22) +
    labs(x = "Cluster", y = index_title) +
    #scale_y_continuous(breaks = seq(.03, .05, .01)) +
    scale_fill_manual(values = clust_cols,
                      labels = c("1" = "1: Eng",
                                 "2" = "2: Eng",
                                 "3" = "3: Eng x Gla",
                                 "4" = "4: Eng x Gla",
                                 "5" = "5: Eng x Gla",
                                 "6" = "6: Sit",
                                 "7" = "7: Gla",
                                 "8" = "8: Gla")) +
    scale_y_continuous(expand = expansion(mult = c(.08, .22))) +
    theme(#aspect.ratio = 1.25,
          legend.position = #c(.75, .24),
            c(.24, .79),
          legend.text = element_text(size = 10, hjust = 0),
          legend.margin = margin(c(3, 4, 4, 4)),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14, margin = margin(r = 20, l = 0)),
          axis.text.x = element_text(vjust = .5, size = 12),
          axis.text.y = element_text(size = 12),
          #axis.title.x = element_blank(),
          legend.background = element_rect(color = "black", linewidth = .5),
          legend.title = element_blank()))

# index_viz = ggarrange(index_map, band_plot, index_plot, nrow = 1)

(index_viz = ggarrange((ggdraw(index_map) +
                       draw_grob(rect) +
                       draw_grob(lab_a)), #cluster_species,
                    (ggdraw(band_plot) +
                       draw_grob(rect) +
                       draw_grob(lab_b)),
                    (ggdraw(index_plot) +
                       draw_grob(rect) +
                       draw_grob(lab_c)),
                    widths = c(1, 1, 1),
                    ncol = 3))


ggsave(plot = index_viz,
       filename = 'D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\SI\\_index_viz_GCC_late_winter.jpeg',
       device = jpeg,
       width = 15,
       height = 6,
       units = 'in',
       dpi = 300,
       bg = 'white')


################################################################################
# mDatt vs GCC greenup

(regr_winter = plot_traits %>%
   dplyr::mutate(cluster = reorder(cluster, trait)) %>% 
   ggplot(aes(
              y = `CCI_late winter`,
              x = `PRI_late winter`,
              # y = `PRI_mid summer`,
              # x = `CCI_mid summer`,
              fill = cluster)) +
   geom_point(size = 3, alpha = .8, shape = 21) +
   #scale_color_manual(values = clust_cols) +
   theme_bw(base_size = 22) +
   #labs(x = "Spectral band", y = index_title) +
   #scale_y_discrete(limits = rev(levels(blups_r2_plot$clim))) +
   scale_fill_manual(values = clust_cols) +
   #scale_y_continuous(expand = expansion(mult = c(.08, .18))) +
   theme(#aspect.ratio = 1,
     legend.position = "right",#c(.24, .79),
     legend.text = element_text(size = 20),
     axis.text.x = element_text(angle = 90, vjust = .5),
     #axis.title.x = element_blank(),
     legend.background = element_rect(color = "black", linewidth = .5),
     legend.title = element_blank()))

ggsave(plot = greenup_regr,
       filename = 'D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated\\SI\\PRI_CCI_greenup_regression.jpeg',
       device = jpeg,
       width = 10,
       height = 10,
       units = 'in',
       dpi = 300,
       bg = 'white')

################################################################################

(pri_cci = plot_traits %>%
   dplyr::select(Prov, cluster, 
                 `PRI_late winter`,
                 `PRI_mid summer`,
                 `CCI_late winter`,
                 `CCI_mid summer`,
                 ) %>% 
   pivot_longer(`PRI_late winter`:`CCI_mid summer`, names_to = "index", values_to = "EST") %>% 
   separate(index, sep = "_", into = c("index", "timepoint")) %>% 
   pivot_wider(names_from = "index", values_from = "EST") %>% 
   ggplot(aes(#color = cluster,
              color = cluster,
              group = cluster,
              x = PRI,
              y = CCI
              #group = factor(Class, levels = c("Wildstand", "Orchard"))
   )) +
   geom_point(shape = 16, size = 4.5, alpha = .8) +
   geom_smooth(se = FALSE, method = "lm",
               aes(group = timepoint),
               show.legend = FALSE, color = "grey20", linetype = 3) +
   theme_bw(base_size = 22) +
   #labs(x = "Cluster", y = index_title) +
   #scale_y_continuous(breaks = seq(.03, .05, .01)) +
   scale_color_manual(values = clust_cols,
                     labels = c("1" = "1: Eng",
                                "2" = "2: Eng",
                                "3" = "3: Eng x Gla",
                                "4" = "4: Eng x Gla",
                                "5" = "5: Eng x Gla",
                                "6" = "6: Sit",
                                "7" = "7: Gla",
                                "8" = "8: Gla")) +
   #scale_y_continuous(expand = expansion(mult = c(.08, .18))) +
   theme(aspect.ratio = 1,
     legend.position = c(.38, .28),
     legend.text = element_text(size = 9, hjust = 0),
     legend.margin = margin(c(4, 4, 4, 4)),
     axis.title.x = element_text(size = 14),
     axis.title.y = element_text(size = 14, margin = margin(r = 20, l = 0)),
     axis.text.x = element_text(vjust = .5, size = 12),
     axis.text.y = element_text(size = 12),
     strip.background = element_rect(fill = "white"),
     legend.background = element_rect(color = "black", linewidth = .5),
     legend.title = element_blank()) +
   facet_wrap(. ~ timepoint,
              scales = "free")
   )

ggsave(plot = pri_cci,
       filename = 'D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\SI\\_PRI_CCI.jpeg',
       device = jpeg,
       width = 12,
       height = 6.5,
       units = 'in',
       dpi = 300,
       bg = 'white')



################################################################################


(plot_band %>% 
   filter(index %in% #c("R668", "R650", "R740", "R705", "R717", "R842")
           c("R444", "R475", "R531", "R560")
     #& cluster %in% c(4, 6, 7, 8)
          & timepoint %in% c("mid summer", "late winter", "late summer")) %>% 
    ggplot(aes(color = cluster,
               fill = cluster,
               group = cluster,
               x = index,
               y = EST
               #group = factor(Class, levels = c("Wildstand", "Orchard"))
    )) +
    stat_summary(geom = "line", linewidth = 1) +
    stat_summary(geom = "point", size = 5, alpha = .7) +
    theme_bw(base_size = 22) +
    labs(x = "Spectral band", y = "Reflectance") +
    #scale_y_discrete(limits = rev(levels(blups_r2_plot$clim))) +
    scale_fill_manual(values = clust_cols) +
    scale_color_manual(values = clust_cols) +
    #scale_y_continuous(expand = expansion(mult = c(.08, .18))) +
    #lims(y = c(0.02, .08)) +
    theme(aspect.ratio = 1,
          legend.position = "right",#c(.24, .79),
          legend.text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = .5),
          axis.title.x = element_blank(),
          legend.background = element_rect(color = "black", linewidth = .5),
          legend.title = element_blank()) +
    facet_wrap2(. ~ timepoint#, scales = "free_y"
                ))

################################################################################

pols_dat_skim %>% 
  dplyr::mutate(Prov = as.character(Prov),
                HT16 = as.numeric(HT16),
                DBH16 = as.numeric(DBH16),
                DBH18 = as.numeric(DBH18)) %>% 
  left_join(pca_plot %>% dplyr::select(Prov, cluster, PC1:PC3)) %>% 
  distinct() %>% 
  ggplot(aes(x = DBH18,
             y = PC3,
             color = cluster)) +
  geom_point(size = 3, alpha = .4, shape = 16) +
  theme_bw(base_size = 22) +
  #labs(x = "Spectral band", y = index_title) +
  #scale_y_discrete(limits = rev(levels(blups_r2_plot$clim))) +
  scale_color_manual(values = clust_cols) +
  #scale_y_continuous(expand = expansion(mult = c(.08, .18))) +
  theme(#aspect.ratio = 1,
    legend.position = "right",#c(.24, .79),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(angle = 90, vjust = .5),
    #axis.title.x = element_blank(),
    legend.background = element_rect(color = "black", linewidth = .5),
    legend.title = element_blank()) 


pols_dat_skim %>% 
  dplyr::mutate(Prov = as.character(Prov),
                HT16 = as.numeric(HT16),
                DBH16 = as.numeric(DBH16)) %>% 
  left_join(plot_band %>% dplyr::select(Prov, cluster)) %>% 
  distinct() %>% 
  pivot_longer(DBH16:HT16, names_to = "trait", values_to = "value") %>% 
  ggplot(aes(x = cluster,
             y = value,
             fill = cluster)) +
  theme_bw(base_size = 22) +
  geom_boxplot(aes(fill = cluster),
               outlier.alpha = 0, width = .5, alpha = .7,
               position = position_dodge(width = .75),
               linewidth = .3, color = "black") +
  geom_point(size = 3, alpha = .5, shape = 21) +
  scale_fill_manual(values = clust_cols) +
  #scale_y_continuous(expand = expansion(mult = c(.08, .18))) +
  theme(#aspect.ratio = 1,
    legend.position = "right",#c(.24, .79),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(angle = 90, vjust = .5),
    #axis.title.x = element_blank(),
    legend.background = element_rect(color = "black", linewidth = .5),
    legend.title = element_blank()) +
  facet_wrap(. ~ trait,
             scales = "free_y")



test = pols_dat_skim %>% 
  dplyr::mutate(Prov = as.character(Prov),
                HT16 = as.numeric(HT16),
                DBH16 = as.numeric(DBH16),
                DBH18 = as.numeric(DBH18)) %>% 
  group_by(Obs) %>% 
  dplyr::mutate(DBH16_18 = (DBH16 + DBH18) / 2) %>% 
  left_join(mm_dat_init %>% 
              filter(index %in% keep_list),
                     by = c("Prov", "Obs")) %>% 
  left_join(plot_traits %>%
              dplyr::select(Prov, cluster) %>%
              st_drop_geometry()) %>%
  distinct()

test %>% 
  filter(index == "PRI"
         & (Date == "2021-06-29"
         | Date == "2021-08-14")
         ) %>% 
  ggplot(aes(x = HT16,
             y = var_value,
             color = cluster)) +
  #labs(y = "CCI") +
  geom_point(size = 3, alpha = .3, shape = 16) +
  # stat_smooth(method = "lm", 
  #             aes(group = factor(Prov)),
  #             alpha = .5,
  #             se = FALSE) +
  stat_smooth(method = "lm", 
              #color = "black",
              #aes(group = "Prov"),
              alpha = .5,
              se = FALSE,
              linewidth = 1.4) +
  theme_bw(base_size = 22) +
  #labs(x = "Spectral band", y = index_title) +
  #scale_y_discrete(limits = rev(levels(blups_r2_plot$clim))) +
  scale_color_manual(values = clust_cols) +
  #scale_y_continuous(expand = expansion(mult = c(.08, .18))) +
  theme(aspect.ratio = 1,
    legend.position = "right",#c(.24, .79),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(angle = 90, vjust = .5),
    #axis.title.x = element_blank(),
    legend.background = element_rect(color = "black", linewidth = .5),
    legend.title = element_blank()) +
  facet_wrap(. ~  Date,
             nrow = 1)

################################################################################


# library(scatterpie)
# 
# (pc_species_1 = ggplot() +
#     geom_scatterpie(data = pca_plot_species, aes(x = PC1, y = PC2), 
#                     color = NA,
#                     cols = c("propGla", "propEng", "propSit"),
#                     pie_scale = .75) +
#     geom_segment(data=dplyr::filter(loadings, loading_abs > .15 & PC == "PC1"), 
#                  aes(x=0, 
#                      y= min(pca_plot$PC2) - (rank * pc2_range) - (pc2_range * mult) - (sign_off * (pc2_range/2)), 
#                      xend = loading * (max(pca_plot$PC1) - min(pca_plot$PC1)) * scale, 
#                      yend = min(pca_plot$PC2) - (rank * pc2_range) - (pc2_range * mult) - (sign_off * (pc2_range/2))), 
#                  color="red4") +
#     geom_text(data=dplyr::filter(loadings, loading_abs > .15 & PC == "PC1"), 
#               aes(y= min(pca_plot$PC2) - (rank * pc2_range) - (pc2_range * mult) - (sign_off * (pc2_range/2)),  
#                   x = loading * (max(pca_plot$PC1) - min(pca_plot$PC1)) * scale,
#                   label = name,
#                   hjust = h_just),
#               size = 4.5, parse = TRUE) +
#     geom_segment(data=dplyr::filter(loadings, loading_abs > .15 & PC == "PC2"), 
#                  aes(y=0, 
#                      x= min(pca_plot$PC1) - (rank * pc1_range) - (pc1_range * mult) - (sign_off * (pc1_range/2)), 
#                      yend = loading * (max(pca_plot$PC2) - min(pca_plot$PC2)) * scale, 
#                      xend = min(pca_plot$PC1) - (rank * pc1_range) - (pc1_range * mult) - (sign_off * (pc1_range/2))), 
#                  color="red4") +
#     geom_text(data=dplyr::filter(loadings, loading_abs > .15 & PC == "PC2"), 
#               aes(x= min(pca_plot$PC1) - (rank * pc1_range) - (pc1_range * mult) - (sign_off * (pc1_range/2)),  
#                   y = loading * (max(pca_plot$PC2) - min(pca_plot$PC2)) * scale,
#                   label = name,
#                   hjust = h_just),
#               size = 4.5, parse = TRUE, angle = 90) +
#     theme_bw(20) +
#     labs(x = "PC1 — greenness (56.1%)", y = "PC2 — green-up (20.6%)") +
#     scale_color_manual(values = clust_cols, name = "Seedlot Class") +
#     scale_fill_manual(values = c("propEng" = "#666666", "propGla" = "#AAAAAA", "propSit" = "#333333"),
#                       labels = c("propEng" = "P. engelmannii", "propGla" = "P. glauca", "propSit" = "P. sitchensis")
#     ) +
#     scale_y_continuous(expand = expansion(mult = c(.05, .05))) +
#     scale_x_continuous(expand = expansion(mult = c(.05, .05))) +
#     coord_equal(clip = "off",
#                 ylim = c(min(pca_plot$PC2), max(pca_plot$PC2)),
#                 xlim = c(min(pca_plot$PC1), max(pca_plot$PC1))) +
#     theme(strip.text = element_blank(),
#           plot.background = element_rect(fill = NA, color = "black", size = 1),
#           aspect.ratio = 1,
#           axis.text = element_text(size = 10),
#           axis.title.y = element_text(margin = margin(t = 0, r = 90, b = 0, l = 0)),
#           axis.title.x = element_text(margin = margin(t = 90, r = 0, b = 0, l = 0)),
#           legend.position = "none"))
# 
# 
# (pc_plot_2 = ggplot() +
#     geom_point(data = pca_plot, aes(x = PC2, y = PC3, color = cluster), size = 4.5, alpha = .7) +
#     geom_segment(data=dplyr::filter(loadings, loading_abs > .2 & PC == "PC2"), 
#                  aes(x=0, 
#                      y= min(pca_plot$PC3) - (rank * pc3_range) - (pc3_range * mult) - (sign_off * (pc3_range/2)), 
#                      xend = loading * (max(pca_plot$PC2) - min(pca_plot$PC2)) * scale, 
#                      yend = min(pca_plot$PC3) - (rank * pc3_range) - (pc3_range * mult) - (sign_off * (pc3_range/2))), 
#                  color="red4") +
#     geom_text(data=dplyr::filter(loadings, loading_abs > .2 & PC == "PC2"), 
#               aes(y= min(pca_plot$PC3) - (rank * pc3_range) - (pc3_range * mult) - (sign_off * (pc3_range/2)),  
#                   x = loading * (max(pca_plot$PC2) - min(pca_plot$PC2)) * scale,
#                   label = name,
#                   hjust = h_just),
#               size = 4.5, parse = TRUE) +
#     geom_segment(data=dplyr::filter(loadings, loading_abs > .2 & PC == "PC3"), 
#                  aes(y=0, 
#                      x= min(pca_plot$PC2) - (rank * pc2_range) - (pc2_range * mult) - (sign_off * (pc2_range/2)), 
#                      yend = loading * (max(pca_plot$PC3) - min(pca_plot$PC3)) * scale, 
#                      xend = min(pca_plot$PC2) - (rank * pc2_range) - (pc2_range * mult) - (sign_off * (pc2_range/2))), 
#                  color="red4") +
#     geom_text(data=dplyr::filter(loadings, loading_abs > .2 & PC == "PC3"), 
#               aes(x= min(pca_plot$PC2) - (rank * pc2_range) - (pc2_range * mult) - (sign_off * (pc2_range/2)),  
#                   y = loading * (max(pca_plot$PC3) - min(pca_plot$PC3)) * scale,
#                   label = name,
#                   hjust = h_just),
#               size = 4.5, parse = TRUE, angle = 90) +
#     theme_bw(20) +
#     labs(x = "PC2 — green-up (20.6%)", y = "PC3 — red edge shift (10.7%)") +
#     scale_color_manual(values = clust_cols, name = "Seedlot Class") +
#     scale_fill_manual(values = c("#666666", "#AAAAAA", "#333333"),
#                       labels = c("P. engelmannii", "P. glauca", "P. sitchensis")) +
#     scale_y_continuous(expand = expansion(mult = c(.05, .05))) +
#     scale_x_continuous(expand = expansion(mult = c(.05, .05))) +
#     coord_equal(clip = "off",
#                 ylim = c(min(pca_plot$PC3), max(pca_plot$PC3)),
#                 xlim = c(min(pca_plot$PC2), max(pca_plot$PC2))) +
#     theme(strip.text = element_blank(),
#           plot.background = element_rect(fill = NA, color = "black", size = 1),
#           aspect.ratio = 1,
#           axis.text = element_text(size = 10),
#           axis.title.y = element_text(margin = margin(t = 0, r = 90, b = 0, l = 0)),
#           axis.title.x = element_text(margin = margin(t = 90, r = 0, b = 0, l = 0)),
#           legend.position = "none"))
# 
# (pc_species_2 = ggplot() +
#     geom_scatterpie(data = pca_plot_species, aes(x = PC2, y = PC3), 
#                     color = NA,
#                     cols = c("propGla", "propEng", "propSit"),
#                     pie_scale = .75) +
#     geom_segment(data=dplyr::filter(loadings, loading_abs > .2 & PC == "PC2"), 
#                  aes(x=0, 
#                      y= min(pca_plot$PC3) - (rank * pc3_range) - (pc3_range * mult) - (sign_off * (pc3_range/2)), 
#                      xend = loading * (max(pca_plot$PC2) - min(pca_plot$PC2)) * scale, 
#                      yend = min(pca_plot$PC3) - (rank * pc3_range) - (pc3_range * mult) - (sign_off * (pc3_range/2))), 
#                  color="red4") +
#     geom_text(data=dplyr::filter(loadings, loading_abs > .2 & PC == "PC2"), 
#               aes(y= min(pca_plot$PC3) - (rank * pc3_range) - (pc3_range * mult) - (sign_off * (pc3_range/2)),  
#                   x = loading * (max(pca_plot$PC2) - min(pca_plot$PC2)) * scale,
#                   label = name,
#                   hjust = h_just),
#               size = 4.5, parse = TRUE) +
#     geom_segment(data=dplyr::filter(loadings, loading_abs > .2 & PC == "PC3"), 
#                  aes(y=0, 
#                      x= min(pca_plot$PC2) - (rank * pc2_range) - (pc2_range * mult) - (sign_off * (pc2_range/2)), 
#                      yend = loading * (max(pca_plot$PC3) - min(pca_plot$PC3)) * scale, 
#                      xend = min(pca_plot$PC2) - (rank * pc2_range) - (pc2_range * mult) - (sign_off * (pc2_range/2))), 
#                  color="red4") +
#     geom_text(data=dplyr::filter(loadings, loading_abs > .2 & PC == "PC3"), 
#               aes(x= min(pca_plot$PC2) - (rank * pc2_range) - (pc2_range * mult) - (sign_off * (pc2_range/2)),  
#                   y = loading * (max(pca_plot$PC3) - min(pca_plot$PC3)) * scale,
#                   label = name,
#                   hjust = h_just),
#               size = 4.5, parse = TRUE, angle = 90) +
#     theme_bw(20) +
#     labs(x = "PC2 — green-up (20.6%)", y = "PC3 — red edge shift (10.7%)") +
#     scale_color_manual(values = clust_cols, name = "Seedlot Class") +
#     scale_fill_manual(values = c("propEng" = "#666666", "propGla" = "#AAAAAA", "propSit" = "#333333"),
#                       labels = c("propEng" = "P. engelmannii", "propGla" = "P. glauca", "propSit" = "P. sitchensis")
#                       ) +
#     scale_y_continuous(expand = expansion(mult = c(.05, .05))) +
#     scale_x_continuous(expand = expansion(mult = c(.05, .05))) +
#     coord_equal(clip = "off",
#                 ylim = c(min(pca_plot$PC3), max(pca_plot$PC3)),
#                 xlim = c(min(pca_plot$PC2), max(pca_plot$PC2))) +
#     theme(strip.text = element_blank(),
#           plot.background = element_rect(fill = NA, color = "black", size = 1),
#           aspect.ratio = 1,
#           axis.text = element_text(size = 10),
#           legend.text = element_text(face = "italic", size = 12),
#           legend.title = element_blank(),
#           legend.position = c(.2, .88),
#           axis.title.y = element_text(margin = margin(t = 0, r = 90, b = 0, l = 0)),
#           axis.title.x = element_text(margin = margin(t = 90, r = 0, b = 0, l = 0))))
# 
# 
# pc_all = ggarrange(pc_plot_1, pc_species_1, pc_plot_2, pc_species_2,
#           ncol = 2, nrow = 2)
# 
# ggsave(plot = pc_all,
#        filename = 'D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated\\figure_8.jpeg',
#        device = jpeg,
#        width = 15,
#        height = 15,
#        units = 'in',
#        dpi = 300,
#        bg = 'white')
# 
# # 
# # (pc1_loadings = loadings %>% 
# #     dplyr::filter(PC == "PC1" & loading_abs > .2) %>% 
# #     dplyr::mutate(h_just = if_else(loading > 0, -.1, 1.1)) %>% 
# #   ggplot(aes(x = loading, y = reorder(rowname, loading_abs))) +
# #   geom_bar(stat = "identity") +
# #     geom_text(aes(label = rowname, hjust = h_just)) +
# #     theme_void() +
# #     scale_x_continuous(expand = expansion(mult = c(.1, .1))))
# # 
# # (pc_plot_2 = ggplot() +
# #     geom_point(data = pca_plot, aes(x = PC2, y = PC3, color = cluster), size = 4.5, alpha = .7) +
# #     # geom_text(data = sites2, 
# #     #           aes(x = PC1, y = PC2,, label = Site.x, hjust = hjust1, vjust = vjust1), 
# #     #           size = 6.5, color = "black", fontface = "bold",
# #     #           family = "Cambria") +
# #     theme_bw(20) +
# #     annotate(geom = "text", label = "G R E E N U P →",
# #              x = max(pca_plot$PC2),
# #              y = min(pca_plot$PC3),
# #              fontface = "bold",
# #              #family = "Cambria",
# #              size = 7,
# #              hjust = 1,
# #              vjust = 2.5,
# #              alpha = .5) +
# #     annotate(geom = "text", label = "R E D  E D G E →",
# #              x = min(pca_plot$PC2),
# #              y = max(pca_plot$PC3),
# #              fontface = "bold",
# #              size = 7,
# #              hjust = 1,
# #              vjust = -1,
# #              alpha = .5,
# #              angle = 90) +
# #     labs(x = "PC1 (42.4%)", y = "PC2 (28.2%)") +
# #     scale_color_manual(values = clust_cols, name = "Seedlot Class") +
# #     scale_fill_manual(values = c("#666666", "#AAAAAA", "#333333"),
# #                       labels = c("P. engelmannii", "P. glauca", "P. sitchensis")) +
# #     scale_y_continuous(expand = expansion(mult = c(.12, .07))) +
# #     scale_x_continuous(expand = expansion(mult = c(.12, .07))) +
# #     theme(strip.text = element_blank(),
# #           plot.background = element_rect(fill = NA, color = "black", size = 1),
# #           #aspect.ratio = 1,
# #           axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)),
# #           axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)),
# #           legend.position = "none"))
# # 
# # (pc_species_1 = ggplot() +
# #     #geom_point(data = pca_plot, aes(x = PC1, y = PC2, color = cluster), size = 4.5, alpha = .6) +
# #     geom_scatterpie(data = pca_plot_species, aes(x = PC2, y = PC3), 
# #                     color = NA,
# #                     cols = c("propGla", "propEng", "propSit"),
# #                     pie_scale = .75) +
# #     # geom_text(data = sites2, 
# #     #           aes(x = PC1, y = PC2,, label = Site.x, hjust = hjust1, vjust = vjust1), 
# #     #           size = 6.5, color = "black", fontface = "bold",
# #     #           family = "Cambria") +
# #     theme_bw(20) +
# #     # annotate(geom = "text", label = "C O L D →",
# #     #          x = max(pca_plot$PC1),
# #     #          y = min(pca_plot$PC2),
# #     #          fontface = "bold",
# #     #          family = "Cambria",
# #     #          size = 7,
# #     #          hjust = 1,
# #     #          vjust = 2.5,
# #     #          alpha = .5) +
# #     # annotate(geom = "text", label = "← W A R M",
# #     #          x = min(pca_plot$PC1),
# #   #          y = min(pca_plot$PC2),
# #   #          fontface = "bold",
# #   #          family = "Cambria",
# #   #          size = 7,
# #   #          hjust = .1,
# #   #          vjust = 2.5,
# #   #          alpha = .5) +
# #   # annotate(geom = "text", label = "M O I S T →",
# #   #          x = min(pca_plot$PC1),
# #   #          y = max(pca_p$PC2),
# #   #          fontface = "bold",
# #   #          family = "Cambria",
# #   #          size = 7,
# #   #          hjust = 1,
# #   #          vjust = -1,
# #   #          alpha = .5,
# #   #          angle = 90) +
# #   # annotate(geom = "text", label = "← A R I D",
# #   #          x = min(all_plot$PC1),
# #   #          y = min(all_plot$PC2),
# #   #          fontface = "bold",
# #   #          family = "Cambria",
# #   #          size = 7,
# #   #          hjust = .2,
# #   #          vjust = -1,
# #   #          alpha = .5,
# #   #          angle = 90) +
# #   #labs(x = "PC1 (42.4%)", y = "PC2 (28.2%)") +
# #     #scale_color_manual(values = clust_cols, name = "Seedlot Class") +
# #     scale_fill_manual(values = c("propEng" = "#666666", "propGla" = "#AAAAAA", "propSit" = "#333333"),
# #                       labels = c("propEng" = "P. engelmannii", "propGla" = "P. glauca", "propSit" = "P. sitchensis")) +
# #     scale_y_continuous(expand = expansion(mult = c(.12, .07))) +
# #     scale_x_continuous(expand = expansion(mult = c(.12, .07))) +
# #     theme(strip.text = element_blank(),
# #           plot.background = element_rect(fill = NA, color = "black", size = 1),
# #           aspect.ratio = 1,
# #           axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)),
# #           axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)),
# #           legend.position = "none"))
# # 
# # 
# # 
# # (pc_species_1 = ggplot() +
# #     #geom_point(data = pca_plot, aes(x = PC1, y = PC2, color = cluster), size = 4.5, alpha = .6) +
# #     geom_density(data = pca_plot_species, aes(x = PC2, y = PC3), 
# #                  color = NA,
# #                  cols = c("propGla", "propEng", "propSit"),
# #                  pie_scale = .8))
# # 
# # 
# # 
# # ################################################################################
# # 
# # # PC1: temperature
# # (pc_plot = pca_spec$rotation %>% 
# #     as.data.frame() %>% 
# #     rownames_to_column() %>% 
# #     dplyr::select(rowname:PC4) %>% 
# #     pivot_longer(cols = PC1:PC4, names_to = "PC", values_to = "loading") %>% 
# #     # dplyr::mutate(rowname = dplyr::recode(rowname, MAT = 'MAT (mean annual temp.)',
# #     #                                MWMT = 'MWMT (mean warmest-month temp.)',
# #     #                                MCMT = 'MCMT (mean coldest-month temp.)',
# #     #                                TD = 'TD (annual temp. difference)',
# #     #                                MAP = 'MAP (mean annual precipitation)',
# #     #                                AHM = "AHM (annual heat:moisture index)",
# #     #                                SHM = "SHM (summer heat:moisture index)",
# #     #                                DD_0 = 'DD_0 (degree-days < 0°C)',
# #     #                                DD5 = 'DD5 (degree-days > 5°C)',
# #     #                                DD18 = 'DD18 (degree-days < 18°C)',
# #     #                                DD_18 = 'DD_18 (degree-days > 18°C)',
# #   #                                NFFD = 'NFFD (number of frost-free days)',
# #   #                                FFP = 'FFP (length of frost-free period)',
# #   #                                bFFP = 'bFFP (beginning of frost-free period)',
# #   #                                eFFP = 'eFFP (end of frost-free period)',
# #   #                                PAS = 'PAS (precipitation as snow)',
# #   #                                EMT = 'EMT (extreme minimum temp.)',
# #   #                                EXT = 'EXT (extreme maximum temp.)',
# #   #                                Eref = "Eref (reference evapotranspiration)",
# #   #                                RH = 'RH (relative huimidity)',
# #   #                                CMI = 'CMI (climate moisture index)',
# #   #                                DD1040 = 'DD1040 (degree-days 10–40°C)',
# #   #                                MSP = 'MSP (mean summer precipitation)',
# #   #                                CMD = 'CMD (climatic moisture deficit)'),
# #   #        PC = dplyr::recode(PC, PC1 = "PC1 (cold–warm)",
# #   #                           PC2 = "PC2 (arid–moist)",
# #   #                           PC3 = "PC3 (lowland–montane)")) %>% 
# #   dplyr::filter(abs(loading) > .15) %>% 
# #   group_by(PC) %>% 
# #     dplyr::mutate(rowname = reorder_within(rowname, abs(loading), PC)) %>% 
# #     ggplot(aes(x = rowname, y = loading)) +
# #     geom_bar(stat = "identity") +
# #     labs(x = "",
# #          y = "ClimateNA variable loadings") +
# #     theme_bw(base_size = 15) +
# #     scale_x_reordered() +
# #     coord_flip() +
# #     theme(axis.text.x = element_text(family = "Cambria", size = 15),
# #           axis.text.y = element_text(family = "Cambria", hjust = 0, vjust = .5, size = 10),
# #           axis.title.x = element_text(family = "Cambria", size = 20, vjust = 1,
# #                                       margin = margin(t = 10, r = 0, b = 0, l = 0)),
# #           strip.background = element_rect(fill = "white"),
# #           strip.text = element_text(family = "Cambria", face = "bold", size = 14)) +
# #     facet_wrap(. ~ PC, scales = "free_y",
# #                nrow = 1))
# # 
# # ggsave(plot = pc_plot,
# #        filename = 'D:\\Sync\\Figures\\_FIGURES_GCB\\PC_loadings.tiff',
# #        device = tiff,
# #        width = 12,
# #        height = 15,
# #        units = 'in',
# #        dpi = 600,
# #        bg = 'white')
# # 
# # 
# # ################################################################################
# # 
# # (df_cluster_scaled = readRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\BLUPS\\df_cluster.rds") %>% 
# #     group_by(index, timepoint) %>% 
# #     dplyr::mutate(ESTIMATE = scale(ESTIMATE)) %>% 
# #     left_join(pops_plot_facet %>% 
# #                 dplyr::select(Prov, cluster)) %>% 
# #     distinct() %>% 
# #     pivot_wider(names_from = "timepoint", values_from = "ESTIMATE") %>% 
# #     ggplot(aes(x = greenup, y = decline, color = cluster)) +
# #     geom_point() +
# #     theme(aspect.ratio = 1) +
# #     facet_wrap(. ~ index,
# #                ncol = 3))
# # 
# # 
# # 
# # (df_cluster_clim = blups_summary %>% 
# #     left_join(pops_plot_facet %>% 
# #                 distinct(Prov, cluster, group_plot),
# #               by = "Prov") %>% 
# #     distinct(Prov, clim, cluster, group_plot, clim_value) %>% 
# #     filter(index %in% c("NIRvCCI", "MCARI", "NDRE3") &
# #              clim %in% c("PAS_log", "MSP_log") ) %>% 
# #     pivot_wider(names_from = clim, values_from = clim_value) %>% 
# #     ggplot(aes(x = PAS_log, y = MSP_log, color = cluster)) +
# #     geom_point() +
# #     theme(aspect.ratio = 1))
# # 
