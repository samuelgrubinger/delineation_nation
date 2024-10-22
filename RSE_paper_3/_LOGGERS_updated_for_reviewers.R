library(tidyverse)
library(zoo)
library(cowplot)
library(ggpubr)
library(myClim)
library(grid)
library(ggh4x)

timepoint_cols = c("late winter" = "darkslategray3",
                   "early summer" = "darkolivegreen4",
                   "mid summer" = "goldenrod",
                   "late summer" =  "goldenrod4")

###########################################
rect = rectGrob(
  y = unit(1, "npc") - unit(.19, "in"),
  x = unit(1, "npc") - unit(.19, "in"),
  width = unit(.4, "in"),
  height = unit(.4, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fill = "white", alpha = 1,
            color = "black"))

font_s = 20

lab_a = textGrob(
  label = "(a)",
  y = unit(1, "npc") - unit(.25, "in"),
  x = unit(1, "npc") - unit(.22, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = font_s))


lab_b = textGrob(
  label = "(b)",
  y = unit(1, "npc") - unit(.25, "in"),
  x = unit(1, "npc") - unit(.22, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_c = textGrob(
  label = "(c)",
  y = unit(1, "npc") - unit(.25, "in"),
  x = unit(1, "npc") - unit(.22, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_d = textGrob(
  label = "(d)",
  y = unit(1, "npc") - unit(.25, "in"),
  x = unit(1, "npc") - unit(.22, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = font_s))

lab_e = textGrob(
  label = "(e)",
  y = unit(1, "npc") - unit(.25, "in"),
  x = unit(1, "npc") - unit(.22, "in"),
  hjust = 1, vjust = 1,
  gp = gpar(fontsize = font_s))

################################################################################


date_display = scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 month",
                            date_labels = "%b %Y",
                            expand = c(0,0))
                            # limits = c(ymd("2020-01-20"), ymd("2023-04-10")))

date_null = theme(axis.text.x = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title.y = element_text(margin=margin(0,15,0,0),
                                              vjust = 0))

flight_dates = geom_vline(xintercept = c(ymd("2020-03-22"), ymd("2020-07-02"), ymd("2020-08-07"), ymd("2020-08-25"),
                                         ymd("2021-03-31"), ymd("2021-06-29"), ymd("2021-07-29"), ymd("2021-08-14"),
                                         ymd("2023-02-24")), 
           color = "grey20", linewidth = 1.1, alpha = .7)

met = read.csv("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\Met_station\\Skimikin_2019.csv") %>%
  mutate(Date = ymd(Date)) %>% 
  bind_rows(read.csv("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\Met_station\\Skimikin_2020.csv") %>% 
              mutate(Date = ymd(Date))) %>% 
  bind_rows(read.csv("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\Met_station\\Skimikin_2022.csv") %>% 
              mutate(Date = mdy(Date))) %>% 
  bind_rows(read.csv("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\Met_station\\Skimikin_2023.csv") %>% 
              mutate(Date = ymd(Date))) %>% 
  distinct()

temp = met %>%     
  group_by((Date)) %>% 
  mutate(Tmax = max(Temp, na.rm = TRUE),
         Tmin = min(Temp, na.rm = TRUE)) %>% 
  dplyr::select(Tmax, Tmin, Date) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(Tmax_smooth = zoo::rollmean(Tmax, k = 30, fill = NA),
         Tmin_smooth = zoo::rollmean(Tmin, k = 30, fill = NA)) %>% 
  filter(Date >= "2020-05-25" & Date <= "2023-04-01" & 
           (Date >= "2023-01-01" | Date <= "2021-12-31")) %>% 
  mutate(period = if_else(Date >= "2023-01-01", 2, 1))

temp2 = read.csv("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\Met_station\\EM40189_2023Sep13-1047.csv") %>% 
            mutate(Date = as_date(ymd_hm(Measurement.Time))) %>%     
  group_by((Date)) %>% 
  mutate(Tmax = max(Temp, na.rm = TRUE),
         Tmin = min(Temp, na.rm = TRUE)) %>% 
  dplyr::select(Tmax, Tmin, Date) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(Tmax_smooth = zoo::rollmean(Tmax, k = 30, fill = NA),
         Tmin_smooth = zoo::rollmean(Tmin, k = 30, fill = NA)) %>% 
  filter(Date >= "2020-05-25" & Date <= "2023-04-01" & 
           (Date >= "2023-01-01" | Date <= "2021-12-31")) %>% 
  mutate(period = if_else(Date >= "2023-01-01", 2, 1))

temp_plot = temp %>% 
  left_join(temp2, by = c("Date", "period")) %>% 
  mutate(Tmax = if_else(period == 1, Tmax.x, Tmax.y),
         Tmin = if_else(period == 1, Tmin.x, Tmin.y),
         Tmax_smooth = if_else(period == 1, Tmax_smooth.x, Tmax_smooth.y),
         Tmin_smooth = if_else(period == 1, Tmin_smooth.x, Tmin_smooth.y)) 


# (RH_met = met %>% 
#     filter(Date >= "2020-05-25" & Date <= "2023-04-01" & 
#              (Date >= "2023-01-01" | Date <= "2021-12-31")) %>% 
#     mutate(period = if_else(Date >= "2023-01-01", 2, 1)) %>% 
#     group_by(Date) %>% 
#     mutate(RH_avg = mean(RH, na.rm = TRUE)) %>% 
#     ungroup() %>% 
#     mutate(RH_new = zoo::rollmean(RH_avg, k = 40, fill = NA)) %>% 
#     # filter(Date > "2020-05-25" & 
#     #        (Date >= "2023-01-01" | Date <= "2021-12-31")) %>%
#     # drop_na() %>% 
#     ggplot(aes(x = ymd(Date), y = RH_new)) +
#     geom_line(color = "orchid4", alpha = .9, linewidth = 1.1, stat = "identity") +
#     scale_y_continuous(expand = expansion(mult = c(.02, .02))) +
#     expand_limits(y = c(0, 100)) +
#     theme_bw(base_size = 16) +
#     labs(y = "Relative humidty (%)") +
#     theme(axis.title.x=element_blank(), 
#           axis.text.x = element_text(angle = 90, vjust = 3.6, face = "bold"),
#           strip.background = element_blank(),
#           strip.text.x = element_blank()) +
#     date_display +
#     date_null +
#     flight_dates +
#     facet_grid(. ~ period,
#                scales = "free_x",
#                space = "free_x"))


temp_met_fig = temp_plot %>% 
  pivot_longer(Tmax:Tmin, names_to = "minmax", values_to = "temp") %>% 
    #mutate(RH_new = zoo::rollmean(RH_avg, k = 20, fill = NA)) %>% 
    #filter(Precip > 0) %>% 
    ggplot(aes(x = ymd(Date), y = temp, color = minmax, group = minmax)) +
    geom_hline(yintercept = 0, color = "grey20") +
    geom_line(size = 1.3, alpha = .4) +
    #geom_line(aes(y = Tmin), size = 1.3, color = "steelblue", alpha = .4) +
    geom_line(aes(y = Tmax_smooth), size = 1, color = "red3", alpha = .8) +
    geom_line(aes(y = Tmin_smooth), size = 1, color = "steelblue", alpha = .8) +
    scale_y_continuous(expand = expansion(mult = c(0.04,0.04)),
                       breaks = seq(-30, 40, 10)) +
  scale_color_manual(values = c(Tmax = "red3",
                                  Tmin = "steelblue"),
                     labels = c(Tmax = "Daily maximum",
                                Tmin = "Daily minimum")) +
    #expand_limits(y = c(0, 100)) +
    theme_bw(base_size = 16) +
    labs(y = "Temperature (°C)") +
    theme(axis.title.x=element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 3.6, face = "bold"),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    date_display +
    date_null +
    flight_dates +
    theme(axis.title.y = element_text(margin=margin(0,19,0,0)),
          legend.position = c(.3, .8),
          legend.title = element_blank(),
          legend.background = element_rect(color = "grey20", linewidth = .4),
          legend.margin = margin(2, 2, 7, 6)) +
    facet_grid(. ~ period,
               scales = "free_x",
               space = "free_x")

(temp_met = ggdraw(temp_met_fig) +
    draw_grob(rect) +
    draw_grob(lab_a))
    
soil_micro = mc_read_files(c("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TOMST\\data_94205651_2023_02_24_0.csv",
                        "D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TOMST\\data_94205652_2023_02_24_0.csv",
                        "D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TOMST\\data_94205654_2023_02_24_0.csv",
                        "D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TOMST\\data_94205655_2023_02_24_0.csv"),
                      dataformat_name = "TOMST", silent = T) %>% 
  mc_calc_vwc(soiltype = "sandy loam A") %>% 
  mc_reshape_wide() %>% 
  dplyr::select(datetime, contains("VWC_moisture")) %>% 
  mutate(Date = as_date(datetime)) %>% 
  group_by(Date) %>% 
  mutate(moisture = mean(
    c(`94205651_1_94205651_VWC_moisture`,
                           `94205652_1_94205652_VWC_moisture`,
                           `94205654_1_94205654_VWC_moisture`,
                           `94205655_1_94205655_VWC_moisture`), na.rm = TRUE)
         ) %>% 
  filter(Date >= "2020-05-25" & Date <= "2023-04-01" & 
           (Date >= "2023-01-01" | Date <= "2021-12-31")) %>% 
  mutate(period = if_else(Date >= "2023-01-01", 2, 1))


micro = read.delim("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TOMST\\data_94205651_2023_02_24_0.csv",
                   header = FALSE,
                   sep = ";") %>% 
  bind_rows(read.delim("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TOMST\\data_94205652_2023_02_24_0.csv",
                       header = FALSE,
                       sep = ";")) %>% 
  bind_rows(read.delim("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TOMST\\data_94205654_2023_02_24_0.csv",
                       header = FALSE,
                       sep = ";")) %>% 
  bind_rows(read.delim("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TOMST\\data_94205655_2023_02_24_0.csv",
                       header = FALSE,
                       sep = ";")) %>% 
  rename(reading = V1,
         date = V2,
         zone = V3,
         moisture = V7,
         shake = V8,
         error = V9) %>% 
  mutate(T1 = as.numeric(sub(",", ".", V4, fixed = TRUE)),
         T2 = as.numeric(sub(",", ".", V5, fixed = TRUE)),
         T3 = as.numeric(sub(",", ".", V6, fixed = TRUE))) %>% 
  separate(col = date, sep = " ", into = c("day", "time")) %>% 
  mutate(date = ymd(day)) %>% 
  filter(date >= as.Date("2020-07-03")) %>% 
  group_by(date) %>% 
  summarise(Tmax = max(T3, na.rm = TRUE),
            Tmean = mean(T3, na.rm = TRUE),
            Tmin = min(T3, na.rm = TRUE),
            moisture = mean(moisture, na.rm = TRUE)) %>% 
  filter(date >= "2020-05-25" & date <= "2023-04-01" & 
           (date >= "2023-01-01" | date <= "2021-12-31")) %>% 
  mutate(period = if_else(date >= "2023-01-01", 2, 1)) %>% 
  rename("Date" = "date") %>% 
  left_join(soil_micro, by = "Date")

min_soil = min(micro$moisture.y)
max_soil = max(micro$moisture.y)

date_labs = met %>% 
  filter(Date %in% c(ymd("2020-07-02"), ymd("2020-08-07"), ymd("2020-08-25"),
                       ymd("2021-03-31"), ymd("2021-06-29"), ymd("2021-07-29"), ymd("2021-08-14"),
                       ymd("2023-02-24"))) %>% 
  distinct(Date) %>% 
  mutate(period = if_else(Date >= "2023-01-01", 2, 1),
         Date_lab = as.character(Date))

soil_met_fig = met %>% 
    filter(Date >= "2020-05-25" & Date <= "2023-04-01" & 
             (Date >= "2023-01-01" | Date <= "2021-12-31")) %>% 
    full_join(dplyr::select(micro, Date, moisture.y),
              by = "Date") %>% 
    mutate(period = if_else(Date >= "2023-01-01", 2, 1)) %>% 
    group_by(Date) %>% 
  mutate(dummy = if_else(week(Date) < 10, min_soil, max_soil)) %>% 
  ggplot(aes(x = ymd(Date), y = dummy)) +
  geom_blank() +
  geom_line(aes(y = moisture.y), color = "chocolate", linewidth = 1.2, alpha = .7) +
    scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 month") +
    scale_y_continuous(expand = expansion(mult = c(.0, .03))) +
    expand_limits(y = 0) +
    theme_bw(base_size = 16) +
    labs(y = expression(Soil~moisture~(m^3/m^3))) +
    theme(axis.title.x=element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 3.6),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    date_display +
    date_null +
    flight_dates +
    theme(axis.title.y = element_text(margin=margin(0,14,0,0))) +
    facet_grid(. ~ period,
               scales = "free_x",
               space = "free_x")

(soil_met = ggdraw(soil_met_fig) +
    draw_grob(rect) +
    draw_grob(lab_b))



precip_met_fig = met %>% 
  mutate(period = if_else(Date >= "2023-01-01", 2, 1)) %>% 
  #group_by(Date) %>% 
  #filter(Rain_mm < 10) %>% 
  mutate(Week = round_date(Date, "week")) %>% 
  group_by(Week) %>% 
  mutate(Precip = sum(Rain_mm, na.rm = TRUE),
         period = median(period)) %>% 
  distinct(Date, Precip, period) %>% 
  filter(Date >= "2020-05-25" & Date <= "2023-04-01" &
           (Date >= "2023-01-01" | Date <= "2021-12-31")) %>%
  ggplot(aes(x = ymd(Date), y = Precip)) +
  geom_bar(fill = "steelblue4", width = 1, alpha = .7, stat = "identity") +
  expand_limits(y = c(0, 55)) +
  scale_y_continuous(expand = expansion(mult = c(.00, .00))) +
  geom_text(data = date_labs,
            aes(x = Date, y = 55, label = Date_lab, fontface = "bold"),
            angle = 90, vjust = -.6, hjust = 1.05, size = 3.9,
            color = "grey20") +
  theme_bw(base_size = 16) +
  labs(y = "Precipitation (mm)", x = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 3.4, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 0),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks = element_blank()) +
  date_display +
  #date_null +
  flight_dates +
  theme(axis.title.y = element_text(margin=margin(0,27,0,0))) +
  facet_grid(. ~ period,
             scales = "free_x",
             space = "free_x")

(precip_met = (ggdraw(precip_met_fig) +
                 draw_grob(rect) +
                 draw_grob(lab_c)))

 # (temp_met = ggplot(micro, aes(x = date, y = Tmax)) +
#   geom_line(size = 1, color = "red3", alpha = .5) +
#   #geom_point(shape = 16, size = 3, color = "red3") +
#   #geom_line(aes(y = Tmean), size = 1, color = "purple", alpha = .5) +
#   geom_line(aes(y = Tmin), size = 1, color = "steelblue", alpha = .5) +
#   geom_smooth(size = 1.2, color = "red3", method = "loess", span = .02, se = FALSE) +
#   #geom_smooth(aes(y = Tmean), size = 1.2, color = "purple", method = "loess", span = .02, se = FALSE) +
#   geom_smooth(aes(y = Tmin), size = 1.2, color = "steelblue", method = "loess", span = .02, se = FALSE) +
#   labs(x = "Date", y = "Daily temperature (C)") +
#     scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 month") +
#     theme_bw() +
#     labs(y = "Temperature (degrees C)") +
#     date_display +
#     theme(axis.text.x = element_text(angle = 90),
#           axis.title.x = element_blank()))

fl = read.csv("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\Chl_fl\\Skimikin_field.csv") %>% 
  mutate(fvfm_mean = rowMeans(.[9:15], na.rm = TRUE)) %>% 
  filter(is.finite(fvfm_mean)) %>% 
  group_by(Seedlot, Date) %>% 
  mutate(pop_mean = mean(fvfm_mean)) %>% 
  group_by(Date) %>% 
  mutate(all_mean = mean(fvfm_mean),
         period = 1) %>% 
  left_join(read.csv("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\Sx_MASTER_ver_21_ADDITIONAL_TRAITS.csv") %>% 
              filter(Site == "Skim") %>% 
              dplyr::select(Meas_seq, Rep, Prov, Tree) %>% 
              mutate(Tree = as.numeric(Tree)),
            by = c("Rep" = "Rep",
                   "Seedlot" = "Prov",
                   "Tree" = "Tree")) %>% 
  rename("Obs" = "Meas_seq") %>% 
  mutate(timepoint = case_when(substr(Date,7,7) %in% c("2", "3") ~ "late winter",
                                        Date %in% c("2020-07-02", "2021-06-28") ~ "early summer",
                                        Date %in% c("2020-08-07", "2021-07-29") ~ "mid summer",
                                        Date %in% c("2020-08-25", "2021-08-14") ~ "late summer"))


min_fl = min(fl$pop_mean)
max_fl = max(fl$pop_mean)

(fl_plot_fig = met %>% 
  filter(Date >= "2020-05-25" & #Date <= "2023-04-01" & 
           #(Date >= "2023-01-01" |
           Date <= "2021-09-30") %>% 
  mutate(period = if_else(Date >= "2023-01-01", 2, 1)) %>% 
  group_by(Date) %>% 
  #mutate(dummy = if_else(week(Date) < 10, min_refl, max_refl)) %>% 
  ggplot(aes(x = ymd(Date))) +
  geom_blank() +
  geom_boxplot(data = fl, aes(y = fvfm_mean, group = Date, fill = timepoint),
               alpha = .7) +
  geom_point(data = fl, 
             #fill = "steelblue4",
             aes(y = fvfm_mean, fill = timepoint),
             shape = 21, size = 3, alpha = .7) +
  scale_fill_manual(values = timepoint_cols) +
  scale_y_continuous(breaks = seq(.6, .85, .05),
                     expand = expansion(mult = c(.07, .07))) +
  #lims(y = c(-.2,-.02)) +
  labs(y = "Fv/Fm") +
  theme_bw(base_size = 16) +
  theme(axis.title.x=element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 3.5, face = "bold"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none") +
  date_display +
  #flight_dates +
  date_null +
  theme(axis.title.y = element_text(margin=margin(0,16,0,0))) +
  # facet_grid(. ~ period,
  #            scales = "free_x",
  #            space = "free_x") +
  rremove("x.text"))

(fl_plot = ggdraw(fl_plot_fig) +
    draw_grob(rect) +
    draw_grob(lab_a))

df_mm_blups = readRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_all_spectral.rds") %>% 
  dplyr::select(-`2023-02-24`) %>% 
  pivot_longer(`2020-07-02`:`2021-08-14`, names_to = "trait", values_to = "value") %>% 
  # #dplyr::select(-`2020-03-22`) %>% 
  # pivot_longer(`2020-07-02`:`2023-02-24`, names_to = "trait", values_to = "value") %>% 
  #dplyr::select(-(mean_summer:decline)) %>% 
  left_join(readRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\traits//pols_dat_skim.rds") %>% dplyr::select(Obs, Blk, Rep, Prov, Class)) %>% 
  filter(index == "PRI" & substr(trait,1,1) == "2") %>%  
  mutate(Date = trait,
         period = if_else(Date < "2023-01-01", 1, 2),
         TRAIT = value,
         Rep = as.numeric(Rep),
         Seedlot = as.numeric(Prov)) %>%
    dplyr::select(index, Date, Obs, Prov, TRAIT) %>% 
  filter(Obs %in% fl$Obs) %>% 
  distinct() %>% 
  mutate(period = if_else(Date <= "2021-12-31", 1, 2)) %>% 
  group_by(Prov, Date) %>% 
  mutate(pop_mean = mean(TRAIT)) %>% 
  group_by(Date) %>% 
  mutate(all_mean = mean(TRAIT),
         timepoint = case_when(substr(Date,7,7) %in% c("2", "3") ~ "late winter",
                                             Date %in% c("2020-07-02", "2021-06-29") ~ "early summer",
                                             Date %in% c("2020-08-07", "2021-07-29") ~ "mid summer",
                                             Date %in% c("2020-08-25", "2021-08-14") ~ "late summer"))
  

# dummy_blups = df_mm_blups %>% 
#   filter(Date > ymd("2021-08-13")) %>%
#   ungroup() %>% 
#   dplyr::select(Obs, Date, pop_mean, Prov) %>% 
#   pivot_wider(names_from = Date, values_from = pop_mean) %>% 
#   mutate(rate = (`2021-08-14` - `2023-02-24`) / 193,
#          `2021-12-31` = `2021-08-14` - (rate * 138),
#          `2023-01-01` = `2021-08-14` - (rate * 142)) %>% 
#   relocate(rate, .after = last_col()) %>% 
#   pivot_longer(`2021-08-14`:`2023-01-01`, names_to = "Date", values_to = "pop_mean") %>% 
#   mutate(period = if_else(Date <= "2021-12-31", 1, 2))

# dummy_blups_all = df_mm_blups %>% 
#   filter(Date > ymd("2021-08-13")) %>%
#   ungroup() %>% 
#   dplyr::select(Obs, Date, all_mean, Prov) %>% 
#   pivot_wider(names_from = Date, values_from = all_mean) %>% 
#   mutate(rate = (`2021-08-14` - `2023-02-24`) / 193,
#          `2021-12-31` = `2021-08-14` - (rate * 138),
#          `2023-01-01` = `2021-08-14` - (rate * 142)) %>% 
#   relocate(rate, .after = last_col()) %>% 
#   pivot_longer(`2021-08-14`:`2023-01-01`, names_to = "Date", values_to = "all_mean") %>% 
#   mutate(period = if_else(Date <= "2021-12-31", 1, 2))

min_refl = min(df_mm_blups$TRAIT)
max_refl = max(df_mm_blups$TRAIT)

(index_plot_fig = met %>% 
  filter(Date >= "2020-05-25" & #Date <= "2023-04-01" & 
           #(Date >= "2023-01-01" |
              Date <= "2021-09-30") %>% 
  mutate(period = if_else(Date >= "2023-01-01", 2, 1)) %>% 
  group_by(Date) %>% 
  mutate(dummy = if_else(week(Date) < 10, min_refl, max_refl)) %>% 
  ggplot(aes(x = ymd(Date))) +
  geom_blank() +
  geom_boxplot(data = df_mm_blups, aes(y = TRAIT, group = Date, fill = timepoint),
               alpha = .7) +
  geom_point(data = df_mm_blups, 
             #fill = "steelblue4",
             aes(y = TRAIT, group = Prov, fill = timepoint),
             shape = 21, size = 3, alpha = .8) +
    scale_fill_manual(values = timepoint_cols) +
  #lims(y = c(-.2,-.02)) +
  labs(y = "PRI") +
  theme_bw(base_size = 16) +
  theme(axis.title.x=element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 3.5, face = "bold"),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none") +
  date_display +
  #flight_dates +
  #date_null +
  scale_y_continuous(breaks = seq(-.2, 0, .05),
                       expand = expansion(mult = c(.07, .07))) +
  theme(axis.title.y = element_text(margin=margin(0,8,0,0))))
  # facet_grid(. ~ period,
  #            scales = "free_x",
  #            space = "free_x") +
  #rremove("x.text"))

index_plot = ggdraw(index_plot_fig) +
    draw_grob(rect) +
    draw_grob(lab_b)

(pri_fvfm = ggarrange(
  fl_plot,
  index_plot,
  ncol = 1, 
  heights = c(.83, 1),
  align = "h"))

p1 = ggarrange(#fl_plot,
               #index_plot,
               temp_met,
               soil_met,
               ncol = 1, align = "h")

(all_plot = ggarrange(p1,
                      precip_met,
                      heights = c(1, .6),
                      ncol = 1, align = "h"))

ggsave(plot = all_plot, 
       filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\figure_1.jpeg", 
       device = jpeg,
       width = 12,
       height = 12,
       units = "in",
       dpi = 300,
       bg = "white")


################################################################################


df_mm_blups = readRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Updated\\Skimikin_all_spectral_average_median.rds") %>% 
  #dplyr::select(-`2020-03-22`) %>% 
  pivot_longer(`2020-07-02`:`2023-02-24`, names_to = "trait", values_to = "value") %>% 
  #dplyr::select(-(mean_summer:decline)) %>% 
  left_join(readRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\traits//pols_dat_skim.rds") %>% dplyr::select(Obs, Blk, Rep, Prov, Class)) %>% 
  # filter(index %in% c(
  #   "NDVI", "mDatt",
  #   "ARI", "EWI9", 
  #   "GCC", "BCC",
  #   "PRI", "CCI",
  #   "RE_total", "NDRE3") & 
  filter(index %in% c(
    "PRI") & 
           substr(trait,1,1) == "2") %>%  
  mutate(Date = trait,
         period = if_else(Date < "2023-01-01", 1, 2),
         TRAIT = value,
         Rep = as.numeric(Rep),
         Seedlot = as.numeric(Prov)) %>%
  dplyr::select(index, Date, Obs, Prov, TRAIT) %>% 
  filter(Obs %in% fl$Obs) %>% 
  distinct() %>% 
  mutate(period = if_else(Date <= "2021-12-31", 1, 2)) %>% 
  group_by(Prov, Date, index) %>% 
  mutate(pop_mean = mean(TRAIT)) %>% 
  group_by(Date) %>% 
  mutate(all_mean = mean(TRAIT))


pri_dat = fl %>% 
  rename("pop_mean_fv_fm" = "pop_mean",
         "all_mean_fv_fm" = "all_mean") %>% 
  dplyr::mutate(Date = case_match(Date, "2021-06-28" ~ "2021-06-29", .default = Date)) %>% 
  left_join(df_mm_blups) %>% 
  dplyr::mutate(timepoint = case_when(substr(Date,7,7) %in% c("2", "3") ~ "late winter",
                                      Date %in% c("2020-07-02", "2021-06-29") ~ "early summer",
                                      Date %in% c("2020-08-07", "2021-07-29") ~ "mid summer",
                                      Date %in% c("2020-08-25", "2021-08-14") ~ "late summer"),
                Season = if_else(substr(Date,7,7) %in% c("2", "3"), "Winter", "Summer"))

lm_pri = lm(fvfm_mean ~ TRAIT, data = pri_dat)

lm_pri_date = pri_dat %>% 
  group_by(Date) %>% 
  nest() %>% 
  mutate(model = map(data,  ~ lm(fvfm_mean ~ TRAIT, data = .)),
         summary = map(model, broom::glance),
         summary = map(summary, dplyr::select, adj.r.squared)) %>% 
  unnest(summary) %>% 
  mutate(r.squared = round(adj.r.squared, 2))

r2_pri = summary(lm_pri)$adj.r.squared %>% 
  round(2) %>% 
  str_remove("^0+")

(fv_fm_pri = pri_dat %>% 
  ggplot(aes(x = TRAIT, y = fvfm_mean, color = timepoint)) +
  geom_point(size = 4, alpha = .8) +
    geom_line(stat = "smooth", method = "lm", se = FALSE, alpha = .7, linewidth = 1) +
    geom_line(stat = "smooth", method = "lm", se = FALSE, color = "black",
                linewidth = 1, alpha = .7) +
    labs(x = "PRI", y = "Fv/Fm") +
    scale_color_manual(values = timepoint_cols) +
    scale_y_continuous(breaks = seq(.6, .85, .05),
                       expand = expansion(mult = c(.04, .17))) +
    geom_text(aes(x = -.164, y = .8705),
              label = expression(R^2),
              size = 8,
              alpha = 1,
              parse = TRUE,
              color = "black") +
    geom_text(aes(x = -.15, y = .87), 
              label = paste0(" = ", r2_pri),
              size = 8,
              alpha = 1,
              parse = FALSE,
              color = "black") +
  theme_bw(base_size = 18) +
    theme(aspect.ratio = 1.1,
          legend.position = c(.8, .14),
          legend.text = element_text(size = 20),
          legend.background = element_rect(color = "black", linewidth = .5),
          legend.title = element_blank())
  # +
  #   facet_wrap2(. ~ index,
  #               nrow = 2,
  #               strip.position = "left",
  #               scales = "free_x")
  )

fv_fm_lab = ggdraw(fv_fm_pri) +
  draw_grob(rect) +
  draw_grob(lab_c)

(all_pri_fvfm = ggarrange(pri_fvfm,
                      fv_fm_lab,
                      widths = c(1, .89),
                      ncol = 2, align = "h"))

ggsave(plot = all_pri_fvfm, 
       filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\Updated_median\\SI\\timepoints_PRI_fvfm.jpeg", 
       device = jpeg,
       width = 16,
       height = 8,
       units = "in",
       dpi = 300,
       bg = "white")

# (fv_fm_pri = pri_dat %>% 
#     left_join(lm_pri_date, by = "Date") %>% 
#     ggplot(aes(x = TRAIT, y = fvfm_mean, color = timepoint)) +
#     geom_point(size = 5, alpha = .7) +
#     geom_line(stat = "smooth", method = "lm", se = FALSE, alpha = .7, linewidth = 1) +
#     geom_line(stat = "smooth", method = "lm", se = FALSE, color = "black",
#               linewidth = 1, alpha = .7) +
#     labs(x = "PRI", y = "Fv/Fm") +
#     scale_color_manual(values = timepoint_cols) +
#     # geom_text(aes(x = -.164, y = .8505),
#     #           label = expression(R^2),
#     #           size = 8,
#     #           alpha = 1,
#     #           parse = TRUE,
#     #           color = "black") +
#     geom_text(aes(x = -.15, y = .85, label = r.squared),
#               size = 8,
#               alpha = 1,
#               parse = FALSE,
#               color = "black") +
#     theme_bw(base_size = 22) +
#     theme(aspect.ratio = 1,
#           legend.position = c(.8, .14),
#           legend.text = element_text(size = 20),
#           legend.background = element_rect(color = "black", linewidth = .5),
#           legend.title = element_blank())
#   +
#     facet_wrap2(. ~ Date,
#                 nrow = 2,
#                 strip.position = "left",
#                 scales = "free_y")
# )
# 
# 
# cluster_species_sf = readRDS("D:\\Sync\\_Sites\\Skimikin_spectral\\CSV\\TRAITS\\Skimikin_cluster_species_sf.rds")
# 
# sx_pops_clim = read.csv("D:\\Sync\\_Sites\\Sx_Genecology\\_Climate\\output\\Sx_CC_Seedlots_Normal_1961_1990_Y_S.csv") %>% 
#   dplyr::rename("Prov" = "Seedlot") %>% 
#   dplyr::select(-X, -Longitude, -MAR) %>% 
#   dplyr::mutate(MSP_log = log(MSP),
#                 PAS_log = log(PAS))
# 
# clust_cols = c("1" = "#7E1A0A", "2" = "#cd6f37", "3" = "#CA9A28", 
#                "4" = "#7e9633", "5" = "#4d8d57", "6" = "#14748f",
#                "7" = "#9070b0", "8" = "#b06a70", "9" = "#8A2665")
# 
# (fv_fm_pri = pri_dat %>% 
#     dplyr::mutate(Prov = as.character(Prov)) %>% 
#     left_join(cluster_species_sf %>% 
#                 dplyr::select(Prov, cluster, propGla:propEng, elev_m)) %>% 
#     left_join(sx_pops_clim %>% 
#                 dplyr::mutate(Prov = as.character(Prov))) %>% 
#     drop_na(cluster) %>% 
#     ggplot(aes(x = factor(timepoint,
#                           levels = c("late winter", "mid summer", "late summer")),
#                y = TRAIT, color = cluster, group = Prov)) +
#     stat_summary(geom = "line", linewidth = 1) +
#     stat_summary(geom = "point", size = 5, alpha = .7) +
#     #labs(x = "PRI", y = "Fv/Fm") +
#     scale_color_manual(values = clust_cols) +
#     # geom_text(aes(x = -.164, y = .8505),
#     #           label = expression(R^2),
#     #           size = 8,
#     #           alpha = 1,
#     #           parse = TRUE,
#     #           color = "black") +
#     # geom_text(aes(x = -.15, y = .85), 
#     #           label = paste0(" = ", r2_pri),
#     #           size = 10,
#     #           alpha = 1,
#     #           parse = FALSE,
#     #           color = "black") +
#     theme_bw(base_size = 22) +
#     theme(aspect.ratio = 1,
#           legend.position = c(.8, .25),
#           legend.text = element_text(size = 20),
#           legend.background = element_rect(color = "black", linewidth = .5),
#           legend.title = element_blank())
#   # +
#   #   facet_wrap2(. ~ timepoint,
#   #               scales = "free")
# )
# 
# ggsave(plot = fv_fm_pri, 
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\SI\\timepoints_PRI_Fv_Fm.jpeg", 
#        device = jpeg,
#        width = 10,
#        height = 10,
#        units = "in",
#        dpi = 300,
#        bg = "white")
# 
# (fv_fm_pri = pri_dat %>% 
#     dplyr::mutate(Prov = as.character(Prov)) %>% 
#     left_join(cluster_species_sf %>% 
#                 dplyr::select(Prov, cluster, propGla:propEng, elev_m)) %>% 
#     left_join(sx_pops_clim %>% 
#                 dplyr::mutate(Prov = as.character(Prov))) %>% 
#     left_join(pri_dat %>% 
#                 dplyr::select(Name, Height) %>% 
#                 group_by(Name) %>% 
#                 mutate(Height2 = mean(Height, na.rm = TRUE)) %>% 
#                 dplyr::select(-Height) %>% 
#                 ungroup()) %>% 
#     drop_na(cluster) %>% 
#     ggplot(aes(x = TRAIT,
#                y = fvfm_mean,
#                color = cluster,
#                size = 3)) +
#     geom_point() +
#     geom_smooth(method = "lm", color = "black") +
#     #labs(x = "PRI", y = "Fv/Fm") +
#     scale_color_manual(values = clust_cols) +
#     # geom_text(aes(x = -.164, y = .8505),
#     #           label = expression(R^2),
#     #           size = 8,
#     #           alpha = 1,
#     #           parse = TRUE,
#     #           color = "black") +
#     # geom_text(aes(x = -.15, y = .85), 
#     #           label = paste0(" = ", r2_pri),
#     #           size = 10,
#     #           alpha = 1,
#     #           parse = FALSE,
#   #           color = "black") +
#   theme_bw(base_size = 22) +
#     theme(aspect.ratio = 1,
#           legend.position = c(.8, .25),
#           legend.text = element_text(size = 20),
#           legend.background = element_rect(color = "black", linewidth = .5),
#           legend.title = element_blank()) +
#     facet_wrap2(. ~ factor(timepoint,
#                            levels = c("late winter", "mid summer", "late summer")),
#                 scales = "free")
# )
# 
# ggsave(plot = fv_fm_pri, 
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\SI\\timepoints_PRI_Fv_Fm.jpeg", 
#        device = jpeg,
#        width = 10,
#        height = 10,
#        units = "in",
#        dpi = 300,
#        bg = "white")
# 
# ################################################################################
# 
# (fv_fm_mwmt = pri_dat %>% 
#     dplyr::mutate(Prov = as.character(Prov)) %>% 
#     left_join(cluster_species_sf %>% 
#                 dplyr::select(Prov, cluster, propGla:propEng, elev_m)) %>% 
#     left_join(sx_pops_clim %>% 
#                 dplyr::mutate(Prov = as.character(Prov))) %>% 
#     drop_na(cluster) %>% 
#     ggplot(aes(x = MWMT, y = fvfm_mean, color = cluster)) +
#     geom_point(size = 5, alpha = .7) +
#     #geom_line(stat = "smooth", method = "lm", se = FALSE, alpha = .7, linewidth = 1) +
#     geom_line(stat = "smooth", method = "lm", se = FALSE, color = "black",
#               linewidth = 1, alpha = .7) +
#     labs(x = "MWMT", y = "Fv/Fm") +
#     scale_color_manual(values = clust_cols) +
#     # geom_text(aes(x = -.164, y = .8505),
#     #           label = expression(R^2),
#     #           size = 8,
#     #           alpha = 1,
#     #           parse = TRUE,
#     #           color = "black") +
#     # geom_text(aes(x = -.15, y = .85), 
#     #           label = paste0(" = ", r2_pri),
#     #           size = 10,
#     #           alpha = 1,
#     #           parse = FALSE,
#   #           color = "black") +
#   theme_bw(base_size = 22) +
#     theme(aspect.ratio = 1,
#           legend.position = c(.8, .14),
#           legend.text = element_text(size = 20),
#           legend.background = element_rect(color = "black", linewidth = .5),
#           legend.title = element_blank()) 
#   +
#     facet_wrap(. ~ timepoint,
#                 scales = "free")
# )
# 
# 
# (fv_fm_elev = pri_dat %>% 
#     dplyr::mutate(Prov = as.character(Prov)) %>% 
#     left_join(cluster_species_sf %>% 
#                 dplyr::select(Prov, cluster, propGla:propEng, elev_m)) %>% 
#     left_join(sx_pops_clim %>% 
#                 dplyr::mutate(Prov = as.character(Prov))) %>% 
#     drop_na(cluster) %>% 
#     ggplot(aes(x = PAS, y = fvfm_mean, color = cluster)) +
#     geom_point(size = 5, alpha = .7) +
#     #geom_line(stat = "smooth", method = "lm", se = FALSE, alpha = .7, linewidth = 1) +
#     geom_line(stat = "smooth", method = "lm", se = FALSE, color = "black",
#               linewidth = 1, alpha = .7) +
#     labs(x = "Elevation", y = "Fv/Fm") +
#     scale_color_manual(values = clust_cols) +
#     # geom_text(aes(x = -.164, y = .8505),
#     #           label = expression(R^2),
#     #           size = 8,
#     #           alpha = 1,
#     #           parse = TRUE,
#     #           color = "black") +
#     # geom_text(aes(x = -.15, y = .85), 
#     #           label = paste0(" = ", r2_pri),
#     #           size = 10,
#     #           alpha = 1,
#     #           parse = FALSE,
#   #           color = "black") +
#   theme_bw(base_size = 22) +
#     theme(aspect.ratio = 1,
#           legend.position = c(.8, .14),
#           legend.text = element_text(size = 20),
#           legend.background = element_rect(color = "black", linewidth = .5),
#           legend.title = element_blank()) 
#   +
#     facet_wrap(. ~ timepoint,
#                scales = "free")
# )
# 
# #(fv_fm_clim = ggarrange(fv_fm_mwmt, fv_fm_elev, ncol = 1))
# 
# ggsave(plot = fv_fm_clim, 
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\SI\\figure_fvfm_clim.jpeg", 
#        device = jpeg,
#        width = 10,
#        height = 10,
#        units = "in",
#        dpi = 300,
#        bg = "white")
# 
# (fv_fm_pri = pri_dat %>% 
#     dplyr::mutate(Prov = as.character(Prov)) %>% 
#     left_join(cluster_species_sf %>% 
#                 dplyr::select(Prov, cluster, propGla:propEng)) %>% 
#     drop_na(cluster) %>% 
#     group_by(timepoint) %>% 
#     dplyr::mutate(cluster = reorder(cluster, fvfm_mean)) %>% 
#     ggplot(aes(x = cluster, y = fvfm_mean, fill = cluster)) +
#     geom_boxplot(aes(fill = cluster),
#                  outlier.alpha = 0, width = .5, alpha = .7,
#                  position = position_dodge(width = .75),
#                  linewidth = .3, color = "black") +
#     geom_point(size = 3, alpha = .8, shape = 21) +
#     #geom_line(stat = "smooth", method = "lm", se = FALSE, alpha = .7, linewidth = 1) +
#     geom_line(stat = "smooth", method = "lm", se = FALSE, color = "black",
#               linewidth = 1, alpha = .7) +
#     labs(x = "PRI", y = "Fv/Fm") +
#     scale_fill_manual(values = clust_cols) +
#     # geom_text(aes(x = -.164, y = .8505),
#     #           label = expression(R^2),
#     #           size = 8,
#     #           alpha = 1,
#     #           parse = TRUE,
#     #           color = "black") +
#     # geom_text(aes(x = -.15, y = .85), 
#     #           label = paste0(" = ", r2_pri),
#     #           size = 10,
#     #           alpha = 1,
#     #           parse = FALSE,
#   #           color = "black") +
#   theme_bw(base_size = 22) +
#     theme(aspect.ratio = 1,
#           legend.position = c(.8, .14),
#           legend.text = element_text(size = 20),
#           legend.background = element_rect(color = "black", linewidth = .5),
#           legend.title = element_blank()) +
#     facet_wrap(. ~ timepoint,
#                scales = "free_x")
# )
# 
# ggsave(plot = fv_fm_pri, 
#        filename = "D:\\Sync\\Figures\\_FIGURES_RSE\\submit\\SI\\figure_pri_fvfm.jpeg", 
#        device = jpeg,
#        width = 10,
#        height = 10,
#        units = "in",
#        dpi = 300,
#        bg = "white")
