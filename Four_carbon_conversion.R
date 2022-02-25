  # Four approaches were used for carbon/biomass conversion.

    library(lubridate); library(tidyverse); library(ggpubr); library(scales); library(lemon);
    
    #' FSC to size, then biovolume and carbon per cell (equation parameters refer to Calibration_FSC.R)
    #' 
    #' @author Changlin Li 
    #' 
    #' @param FSC normalized to standard beads (0.97-µm)
    #' 
    #' @return Dia cell size (µm, equivalent spherical diameter, ESD)
    #' @return Carbon (μg C per cell)
    #' 
    FSC_Vol <- function(FSC) {
      Dia <- 10 ^ (log10(FSC) * 0.250378257 + 0.135092677)
      Vol <- 4/3*pi*((Dia)/2)^3
      return(Vol)
    }
    #' Four approaches were used for carbon (per cell) conversion.
    #'  
    #' @note our study, with a conversion factor of 280 fg C µm^-3
    #' @references Elemental composition of single cells of various strains of marine Prochlorococcus and Synechococcus using X-ray microanalysis.
    #' @references Seasonal variability of picoplankton in the Northern South China Sea at the SEATS station.
        Carbon_c280   <- function(Vol) { Vol * 280 / 1000 }  
        
    #' @note $log carbon (pg cell^-1) = 0.94 * log volume (μm^3) - 0.6$ (Origin: Eppley et al. 1970)
    #' @references Assessing the dynamics and ecology of marine picophytoplankton: the importance of the eukaryotic component.
        Carbon_log    <- function(Vol) { 10^(log10(Vol) *.94 -.6) }
        
    #' @note $C (pg) = 0.433 * (Biovolume)^0.863$ (Origin: Verity et al. 1992)
    #' @references Assessing the dynamics and ecology of marine picophytoplankton: the importance of the eukaryotic component.    
        Carbon_e0.433 <- function(Vol) { Vol^.863 *.433 }
    
    #' @note $C (pg) = 0.261 * (Biovolume)^0.860$
    #' @references SeaFlow data v1, high-resolution abundance, size and biomass of small phytoplankton in the North Pacific.
    #' @references Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton. 
        Carbon_e0.261 <- function(Vol) { Vol^.860 *.261 }
    
    #rm(list = ls())
    Pico <- read_csv("Data/Picoauto.csv")
    Sampletime <- read_csv("Data/Sampletime.csv") %>% 
      mutate(Time = ymd_hms(paste(Date0, Time0))) %>%
      mutate(Gap = case_when(
        Station == "SEATS" ~ 0,
        Station == "M4"    ~ 4 * 3600 *24,
        Station == "K11"   ~ 10 * 3600 *24)) %>%
      mutate(Time = Time - Gap) %>% 
      select(Id, Time)
    
    # join the data
    Sampletime %>% count(Id) %>% filter(n > 1)
    Pico %>% anti_join(Sampletime, by = "Id")


    Pico <- Pico %>% left_join(Sampletime, by = "Id") %>%
      arrange(Time) %>%
      mutate(Syn_c280    = Carbon_c280  (FSC_Vol(SynFSC)) * Syn,
             Pro_c280    = Carbon_c280  (FSC_Vol(ProFSC)) * Pro,
             Peuk_c280   = Carbon_c280  (FSC_Vol(PeukFSC))* Peuk,
             Syn_log     = Carbon_log   (FSC_Vol(SynFSC)) * Syn,
             Pro_log     = Carbon_log   (FSC_Vol(ProFSC)) * Pro,
             Peuk_log    = Carbon_log   (FSC_Vol(PeukFSC))* Peuk,
             Syn_e0.433  = Carbon_e0.433(FSC_Vol(SynFSC)) * Syn,
             Pro_e0.433  = Carbon_e0.433(FSC_Vol(ProFSC)) * Pro,
             Peuk_e0.433 = Carbon_e0.433(FSC_Vol(PeukFSC))* Peuk,
             Syn_e0.261  = Carbon_e0.261(FSC_Vol(SynFSC)) * Syn,
             Pro_e0.261  = Carbon_e0.261(FSC_Vol(ProFSC)) * Pro,
             Peuk_e0.261 = Carbon_e0.261(FSC_Vol(PeukFSC))* Peuk) %>% 
      select(Station,Id, Id1, Time, contains("_")) %>% 
      pivot_longer(cols      = Syn_c280:Peuk_e0.261,
                   names_to  = c("pico", "equation"),
                   names_sep = "_",
                   values_to = "biomass") %>% 
      drop_na(biomass) %>% 
      group_by(Station, pico, equation, Id) %>% 
      summarise(Time = mean(Time), u = mean(biomass), SD = sd(biomass)) %>% 
      mutate(pico = factor(pico, levels = c("Pro","Syn","Peuk")),
             Station = factor(Station, levels = c("SEATS","M4","K11")),
             equation = factor(equation, levels = c("c280", "log", "e0.433", "e0.261")))
    
    library(suncalc) # Get Sunlight times # use station SEATS for night-time
    # The sampling area is located in time zone 8
    Sun <- getSunlightTimes(as.Date(c("2019-06-19", "2019-06-20", "2019-06-21", "2019-06-22", "2019-06-23")),
                            keep = c("sunrise", "sunset"), lat = 18, lon = 116)
    Night1 <- Sun$sunset[-5] + 8 * 3600
    Night2 <- Sun$sunrise[-1]+ 8 * 3600
    
    lab.b   <- c(Pro = "Prochlorococcus",Syn = "Synechococcus",Peuk = "Picoeukaryotes") 
    Icolor <- c("#ee7e32","#1471b8","#6f3996","#1fb050") 
    lwd_pt <- 1/(.pt*72.27/96)
    theme_set(theme_bw(base_size = 10, base_line_size = .4*lwd_pt, base_rect_size = .4*lwd_pt) +
                theme(plot.title = element_text(face = "bold", hjust = 0.5),
                      panel.border = element_rect(fill = NA, colour = "black", size = rel(2)), 
                      axis.text  = element_text(colour = "black"),
                      axis.ticks = element_line(colour = "black"),
                      strip.text = element_text(face="bold"),
                      strip.background = element_rect(fill = "transparent", size = rel(1.5)),
                      plot.background  = element_rect(colour = NA),
                      
                      panel.spacing.y=unit(-.2, "lines"),
                      panel.grid =element_blank(),
                      legend.position = "bottom",
                      legend.background = element_rect(colour = "black", linetype = 2),
                      legend.key.width  = unit(.6, "cm"),
                      legend.key.height = unit(0.24, "cm"),
                      legend.spacing.x  = unit(0.25, 'cm')))
    
    library(ggh4x)
    ggplot(Pico, aes(Time, u, colour = equation)) +
      annotate("rect", xmin = Night1, xmax = Night2, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") +
      geom_line(size=.3) +
      geom_errorbar(aes(ymin = u - SD, ymax = u + SD), width=2400, size = .2, color="black") +
      geom_point(size =.7) +
      scale_x_datetime(breaks = seq(ymd_hms("2019-06-20 06:00:00"), ymd_hms("2019-06-22 18:00:00"),
                                    "6 hours"), labels = date_format("%H")) +
      xlab(NULL) +  ylab(NULL) +
      labs(title = "Biomass (µg C L-1) was estimated using four different approaches") +
      facet_grid2(vars(pico), vars(Station), 
                  scales = "free", independent = "y",
                  axes = "all", remove_labels = "x",
                  labeller=labeller(pico = as_labeller(lab.b))) +
      coord_cartesian(xlim = c(ymd_hms("2019-06-20 06:50:00"), ymd_hms("2019-06-22 17:30:00"))) +
      scale_color_manual(values = Icolor) +
      theme(panel.spacing.y=unit(.2, "lines"))
    
    SEATS <- Pico %>% filter(Station == "SEATS")
    M4    <- Pico %>% filter(Station == "M4")
    K11   <- Pico %>% filter(Station == "K11")
    
    c280   <- expression("fgC" ~ cell^-1 == 280 %*% Volume)
    `log`  <- expression(log (pgC ~ cell^-1) == .94 %*% log(Volume) - .6)
    e0.433 <- expression("pgC" ~ cell^-1 == .433 %*% Volume^.863)
    e0.261 <- expression("pgC" ~ cell^-1 == .261 %*% Volume^.860)

    # stat_summary(fun.data = 'mean_sd', geom = "errorbar", color="black") did not work;
    p1 <- ggplot(SEATS, aes(Time, u, color = equation)) + 
      annotate("rect", xmin = Night1, xmax = Night2, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") + 
      geom_line(size = .8*lwd_pt) +
      geom_errorbar(aes(ymin = u - SD, ymax = u + SD), width=2400, size = .5*lwd_pt, color="black") +
      geom_point(size = 1.5*lwd_pt) +
      scale_x_datetime(breaks = seq(ymd_hms("2019-06-20 06:00:00"), ymd_hms("2019-06-22 18:00:00"), 
                                    "6 hours"), labels = date_format("%H")) +
      xlab("Time of day (local time)") + ylab("Biomass (µg C L-1)") +
      facet_rep_grid(pico~Station,scales = "free_y")+
      coord_cartesian(xlim = c(ymd_hms("2019-06-20 06:30:00"), ymd_hms("2019-06-22 07:15:00"))) +
      scale_color_manual(values = Icolor) +
      theme(strip.background.y = element_blank(),
            strip.text.y = element_blank()) +
      guides(color = "none")
    
    p2 <- ggplot(M4, aes(Time, u, color = equation)) +
      annotate("rect", xmin = Night1, xmax = Night2, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") +
      geom_line(size = .8*lwd_pt) +
      geom_errorbar(aes(ymin = u - SD, ymax = u + SD), width=2400, size = .5*lwd_pt, color="black") +
      geom_point(size = 1.5*lwd_pt) +
      scale_x_datetime(breaks = seq(ymd_hms("2019-06-20 06:00:00"), ymd_hms("2019-06-22 18:00:00"),
                                    "6 hours"), labels = date_format("%H")) +
      xlab(NULL) +  ylab(NULL) +
      facet_rep_grid(pico~Station,scales = "free_y")+
      coord_cartesian(xlim = c(ymd_hms("2019-06-20 18:00:00"), ymd_hms("2019-06-22 18:30:00"))) +
      scale_color_manual(name="Equations:", 
                         values = Icolor,
                         breaks=c("c280","log","e0.433","e0.261"), 
                         labels=c(c280,log,e0.433,e0.261)) +
      theme(strip.background.y = element_blank(),
            strip.text.y = element_blank()) +
      guides(color=guide_legend(nrow=2,byrow=TRUE,label.hjust = 0))

    p3 <- ggplot(K11, aes(Time, u, color = equation)) + 
      annotate("rect", xmin = Night1, xmax = Night2, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") +
      geom_line(size = .8*lwd_pt) +
      geom_errorbar(aes(ymin = u - SD, ymax = u + SD), width=2400, size = .5*lwd_pt, color="black") +
      geom_point(size = 1.5*lwd_pt) +
      scale_x_datetime(breaks = seq(ymd_hms("2019-06-20 06:00:00"), ymd_hms("2019-06-22 18:00:00"), 
                                    "6 hours"), labels = date_format("%H")) +
      xlab(NULL) +  ylab(NULL) +
      facet_rep_grid(pico~Station,scales = "free_y",labeller=labeller(pico = as_labeller(lab.b)))+
      coord_cartesian(xlim = c(ymd_hms("2019-06-20 17:15:00"), ymd_hms("2019-06-21 20:40:00"))) +
      scale_color_manual(values = Icolor) +
      guides(color = "none")
    
    # use same scale
    Ran1 <- ggplot_build(p1)$layout$panel_params[[1]]$x.range; ran1 <- Ran1[2] - Ran1[1];#ran1
    Ran2 <- ggplot_build(p2)$layout$panel_params[[1]]$x.range; ran2 <- Ran2[2] - Ran2[1];#ran2
    Ran3 <- ggplot_build(p3)$layout$panel_params[[1]]$x.range; ran3 <- Ran3[2] - Ran3[1];#ran3
    
    library(patchwork)
    p1 + p2 + p3 + plot_layout(widths = c(1,ran2/ran1,ran3/ran1))
    #pdf 7.1 * 5
    
    