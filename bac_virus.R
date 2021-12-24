    
    library(lubridate);library(scales);library(ggpubr);library(tidyverse);
    library(grid);library(lemon)
    
    Auto <- read_csv("Data/Picoauto.csv") %>% 
      drop_na(Pro) %>% 
      select(Station, Id, Id1, Pro)
    bac_virus <- read_csv("Data/Bac_virus.csv")
    Sampletime <- read_csv("Data/Sampletime.csv") %>% 
      mutate(Time = ymd_hms(paste(Date0, Time0))) %>%
      mutate(Gap = case_when(
        Station == "SEATS" ~ 0,
        Station == "M4"    ~ 4 * 3600 *24,
        Station == "K11"   ~ 10 * 3600 *24)) %>%
      mutate(Time = Time - Gap) %>% 
      select(Id, Time)
    
    
    bac_virus %>% count(Id, Id1) %>% filter(n > 1)
    Auto %>% anti_join(bac_virus, by = c("Id","Id1"))
    
    Sampletime %>% count(Id) %>% filter(n > 1)
    Auto %>% anti_join(Sampletime, by = "Id")
    
    dat <- Auto %>% 
      left_join(bac_virus,  by = c("Id","Id1")) %>% 
      left_join(Sampletime, by = "Id") %>% 
      arrange(Time) %>% 
      mutate(hBac = (Bac - Pro)*9.1/1000,
             Virus = Virus   * 0.2/1000) %>% 
      select(Station,Id,Time,hBac,Virus) %>% 
      gather(Para, Value, hBac, Virus) %>% 
      mutate(Station = factor(Station,levels = c("SEATS","M4","K11")))
    lab.p   <- c(hBac = "Bacteria",Virus = "Virus")

    # use station SEATS for night-time: 20190620-21
    Night1 <- ymd_hms(c("2019-06-19 18:53:33", "2019-06-20 18:53:47",
                        "2019-06-21 18:54:01", "2019-06-22 18:54:14"));
    Night2 <- ymd_hms(c("2019-06-20 05:41:11", "2019-06-21 05:41:24",
                        "2019-06-22 05:41:37", "2019-06-23 05:41:51"));
    
    Icolor <- c("#ee7e32","#1fb050","#6f3996") 
    lwd_pt <- 1/(.pt*72.27/96) 
    theme_set(theme_bw(base_size = 10, base_line_size = .4*lwd_pt, base_rect_size = .4*lwd_pt) +
                theme(plot.title = element_text(face = "bold", hjust = 0),
                      panel.border = element_rect(fill = NA, colour = "black", size = rel(2)),
                      axis.text  = element_text(colour = "black"),
                      axis.ticks = element_line(colour = "black"),
                      axis.title = element_text(size = rel(0.9)),
                      strip.text = element_text(face="bold"),
                      strip.background = element_rect( size = rel(1.5)),
                      plot.background  = element_rect(colour = NA),
                      
                      panel.spacing.y=unit(-.2, "lines"),
                      panel.grid =element_blank(),
                      legend.position = "bottom",
                      legend.title = element_blank(),
                      legend.background = element_blank(),
                      legend.text = element_text(size = rel(0.7)),
                      legend.key.width  = unit(.6, "cm"),
                      legend.key.height = unit(0.24, "cm"),
                      legend.spacing.x  = unit(0.25, 'cm')))
    

    ggplot(dat, aes(Time, Value, color = Station)) +
      annotate("rect", xmin = Night1, xmax = Night2, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") +
      geom_line(stat = "summary", fun = mean, size=.8*lwd_pt) +
      stat_summary(fun.data = 'mean_sd', geom = "errorbar", width=2400, size=.5*lwd_pt,color="black") +
      geom_point( stat = "summary", fun = mean,  size =1.5*lwd_pt) +
      scale_x_datetime(breaks = seq(ymd_hms("2019-06-20 06:00:00"), ymd_hms("2019-06-22 18:00:00"),
                                    "6 hours"), labels = date_format("%H")) +
      xlab("Time of day (h)") + ylab(NULL) + ggtitle("Biomass (µg C L-1)") +
      facet_rep_grid(Para~.,scales = "free_y",labeller=labeller(Para = as_labeller(lab.p))) +
      coord_cartesian(xlim = c(ymd_hms("2019-06-20 06:50:00"), ymd_hms("2019-06-22 17:30:00"))) +
      scale_color_manual(values = Icolor)  #pdf 6.8 * 3.6

    
   
      