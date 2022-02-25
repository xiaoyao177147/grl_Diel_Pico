    
    library(tidyverse); library(lubridate); library(ggpubr); library(scales);
    #' size and biomass calculation (equation parameters refer to Calibration_FSC.R)
    #' 
    #' @author Changlin Li 
    #' 
    #' @param FSC normalized to standard beads (0.97-µm)
    #' @param Abu abudance (10^3 cells mL-1)
    #'
    #' @return Dia cell size (µm, equivalent spherical diameter, ESD)
    #' @return Biomass  (μg C L-1)
    #' 
    FSC_Dia <- function(FSC){
      Dia <- 10 ^ (log10(FSC) * 0.250378257 + 0.135092677)
      return(Dia)
    }
    
    FSC_Abu_Biomass <- function(FSC, Abu) {
      Dia <- 10 ^ (log10(FSC) * 0.250378257 + 0.135092677)
      Vol <- 4/3*pi*((Dia)/2)^3
      Carbon  <- Vol * 280
      Biomass <- Carbon * Abu / 1000
      return(Biomass)
    }
    
    # the average peak time was approximately 04:30 for Prochlorococcus in the SCS study
    df <- read_csv("SeaFlow/Data/result.csv") %>% 
      filter(cruise %in% c("SCS2019_M4","SCS2019_SEATS","SCS2019_K11") & para == "Abu_Pro") %>% 
      mutate(peak_time = phi/2/pi * 24) 
    # Assume that cell division does not occur during the photoperiod,
    # cell division occurs only during the several hours of darkness
    (peak_time <- mean(df$peak_time)) # ~ 4.5
    
    
    # the median values of simulated data were the same as that of the observed data
    Sampletime <- read_csv("Data/Sampletime.csv") %>% 
      mutate(Time = ymd_hms(paste(Date0, Time0))) %>%
      select(Id, Time)

    Pico <- read_csv("Data/Picoauto.csv") %>% 
      left_join(Sampletime, by = "Id") %>% 
      drop_na(Pro) %>%
      mutate(Pro_Dia     = FSC_Dia(ProFSC),
             Pro_Biomass = FSC_Abu_Biomass(ProFSC, Pro))
    median(Pico$Pro); median(Pico$Pro_Biomass); median(Pico$Pro_Dia)
    
    # net primary production rate (p) was the mean biomass-based net growth rate during the daytime 
    # calculated by three methods (Figure 1c); the tibble fig3 are from Field_SCS.R
    Rates <- fig3 %>%
      filter(Pico == "Pro" & Para == "bio") %>% 
      group_by(Id1) %>% 
      summarise(u = mean(rate))
    Rates$u[1] # 0.0415103

# Start the model ---------------------------------------------------------

    #' @param cn cell numbers (10^3 cells mL-1)
    #' @param b  biomass
    #' @param cs cell size
    #' @param cc carbon per cell
    #' @param g  grazing rate per hour
    #' @param p  net primary production rate per hour
    #' @param r  respiration rate in dark per hour
    #' @param d  rate of cell division per hour
    
    p <- 0.0415103
    r <- p/4
    g <- 0.0155383 # Grazing occurs at a constant rate throughout the period of 24 hours.
    d <- 0.03554575
    # The simulation was conducted on 5-minute time steps for 24 h.
    p <- c(rep(0, 6*12), rep(p, 12*12), rep(0, 6*12)); # Net primary production occurs at a constant rate during the day.
    r <- c(rep(r, 6*12), rep(0, 12*12), rep(r, 6*12)); # Respiration occurs at a constant rate throughout the night.
    d <- c(rep(d, 4.5*12), rep(0, 13.5*12), rep(d, 6*12)); # cells divide according to this schedule, peak_time = 4.5
    
    b0 <- vector(mode = "numeric", length = 0);
    cn <- b0
    cc <- b0
    b  <- b0
    # the input data for starting the model at 00:00
    cn[1] <- 167.2394
    b[1]  <- 2.580935
    cc[1] <- b[1] / cn[1]
    
    for (k in 2:289){
      #k=2
      cn[k] = cn[k-1] * (1 + d[k-1]/12 - g/12)
      b[k]  = b[k-1]  * (1 + p[k-1]/12 - r[k-1]/12 - g/12)
      cc[k] = b[k]/cn[k]
    }
    
    cs = (cc * 1000 / 280 * 3 / 4 / pi) ^ (1/3) * 2 # cell size
    
    cn[289]/cn[1];
    cc[289]/cc[1];
    b[289]/b[1];
    median(cn);median(b);median(cs)
    
    #max(cn);min(cn);max(cn)-min(cn);
    #max(cs);min(cs);
    #max(b) ;min(b);
    
    # draw plots
    dat1 <- tibble(t = seq(0,24,1/12),
                   k = t * 12 + 1,
                   b  = b,
                   cc = cc,
                   cn = cn,
                   cs = cs) 
    
    #' ggplot2 - adding another y-axis on right side of a plot 
    #' 
    #' @source https://stackoverflow.com/questions/36754891/ggplot2-adding-secondary-y-axis-on-top-of-a-plot
    #' 
    #' @param p1/p2 ggplot plot
    #'
    #' @return An ggplot grobs with two or three (p1 use sec.axis) y-axis. 
    #' The two/three layers are completely proportionally overlaid, which is exactly the same as if they were drawn separately.
    #' 
    plus_y_axis <- function(p1, p2) {
      # Get the ggplot grobs
      g1 <- ggplotGrob(p1)
      g2 <- ggplotGrob(p2)
      
      # Get the location of the plot panel in g1.
      # These are used later when transformed elements of g2 are put back into g1
      pp <- c(subset(g1$layout, name == "panel", se = t:r))
      
      # ggplot contains many labels that are themselves complex grob; 
      # usually a text grob surrounded by margins.
      # When moving the grobs from, say, the left to the right of a plot,
      # make sure the margins and the justifications are swapped around.
      # The function below does the swapping.
      # Taken from the cowplot package:
      # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
      hinvert_title_grob <- function(grob){
        # Swap the widths
        widths <- grob$widths
        grob$widths[1] <- widths[3]
        grob$widths[3] <- widths[1]
        grob$vp[[1]]$layout$widths[1] <- widths[3]
        grob$vp[[1]]$layout$widths[3] <- widths[1]
        
        # Fix the justification
        grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
        grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
        grob$children[[1]]$x <- unit(1, "npc") - grob$children[[1]]$x
        grob
      }
      
      # Get the y axis title from g2 - "Elevation (ft)" 
      index <- which(g2$layout$name == "ylab-l") # Which grob contains the y axis title?
      ylab <- g2$grobs[[index]]                # Extract that grob
      ylab <- hinvert_title_grob(ylab)         # Swap margins and fix justifications
      
      # Put the transformed label on the right side of g1
      g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
      g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = "off", name = "ylab-r")
      
      # Get the y axis from g2 (axis line, tick marks, and tick mark labels)
      index <- which(g2$layout$name == "axis-l")  # Which grob
      yaxis <- g2$grobs[[index]]                  # Extract the grob
      
      # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
      # The relevant grobs are contained in axis$children:
      #   axis$children[[1]] contains the axis line;
      #   axis$children[[2]] contains the tick marks and tick mark labels.
      
      # First, move the axis line to the left
      yaxis$children[[1]]$x <- unit.c(unit(0, "npc"), unit(0, "npc"))
      
      # Second, swap tick marks and tick mark labels
      ticks <- yaxis$children[[2]]
      ticks$widths <- rev(ticks$widths)
      ticks$grobs <- rev(ticks$grobs)
      
      # Third, move the tick marks
      ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, "npc") + unit(3, "pt")
      
      # Fourth, swap margins and fix justifications for the tick mark labels
      ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
      
      # Fifth, put ticks back into yaxis
      yaxis$children[[2]] <- ticks
      
      # Put the transformed yaxis on the right side of g1
      g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
      g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = "off", name = "axis-r")
      
      # Draw it
      # grid.newpage(); grid.draw(g1)
      return(g1)
    }
    
    library(grid);  library(gtable)
    lwd_pt <- 1/(.pt*72.27/96)
    theme_set(theme_bw(base_size = 10, base_line_size = .7*lwd_pt, base_rect_size = .7*lwd_pt))

    Theme1 <- theme(axis.text.y.left   = element_text(color = "#1fb050"),
                    axis.ticks.y.left  = element_line(color = "#1fb050"),
                    axis.line.y.left   = element_line(color = "#1fb050"),
                    axis.title.y.left  = element_text(color = "#1fb050", size = rel(0.9)),

                    axis.text.y.right  = element_text(color = "#ee7e32"),
                    axis.ticks.y.right = element_line(color = "#ee7e32"),
                    axis.line.y.right  = element_line(color = "#ee7e32"),
                    axis.title.y.right = element_text(color = "#ee7e32", angle = 90, hjust = 0.5, size = rel(0.9)),
                    axis.title.x       = element_text(size  = rel(0.9)),
                    panel.grid         = element_blank(),
                    panel.border       = element_blank(),
                    plot.title         = element_text(size = rel(1.1)))
    Theme2 <- theme(axis.text.y   = element_text(color = "#6f3996"),
                    axis.ticks.y  = element_line(color = "#6f3996"),
                    axis.line.y   = element_line(color = "#6f3996"),
                    axis.title.y  = element_text(color = "#6f3996", size = rel(0.9)))
    

    # Set up an equal scale combination plot - dual or three y-axis
    p1 <- ggplot(dat1, aes(x = t, y = cn)) + geom_point() + xlab(NULL) + guides(colour = "none");#p1
    p2 <- ggplot(dat1, aes(x = t, y = cs)) + geom_point() + xlab(NULL) + guides(colour = "none");#p2
    p3 <- ggplot(dat1, aes(x = t, y = b))  + geom_point() + xlab(NULL) + guides(colour = "none");#p3
    range1 <- ggplot_build(p1)$layout$panel_params[[1]]$y.range
    range2 <- ggplot_build(p2)$layout$panel_params[[1]]$y.range
    range3 <- ggplot_build(p3)$layout$panel_params[[1]]$y.range
    k2 <- (range2[2] - range2[1]) / (range1[2] - range1[1])
    b2 <- (range1[2] * range2[1] - range1[1] * range2[2]) / (range1[2] - range1[1])
    k3 <- (range3[2] - range3[1]) / (range1[2] - range1[1])
    b3 <- (range1[2] * range3[1] - range1[1] * range3[2]) / (range1[2] - range1[1])
    
    
    plot_1a <- ggplot(dat1, aes(x = t)) +
      annotate("rect", xmin = -Inf, xmax = 6, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") +
      annotate("rect", xmin = 18, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") +
      geom_line(aes(y = cn, color = "y1"), size = lwd_pt) +
      #geom_point(aes(y = cn, color = "y1"), size = .7) +
      scale_y_continuous(sec.axis = sec_axis(trans = ~ .* k2 + b2, name = "Cell size (mean ESD, µm)\n",
                                             breaks = c(.45,.47,.49,.51)),
                         breaks = breaks_pretty(4)) +
      geom_line( aes(y = (cs - b2) / k2, color = "y2"), size = lwd_pt) +
      #geom_point(aes(y = (cs - b2) / k2, color = "y2"), size = .7) +
      geom_line( aes(y = (b  - b3) / k3, color = "y3"), size = lwd_pt) +
      #geom_point(aes(y = (b  - b3) / k3, color = "y3"), size = .7) +
      annotate(geom='segment',y=-Inf,yend=-Inf,x=as.POSIXct(-Inf),xend=as.POSIXct(Inf), size = .7*lwd_pt) +
      annotate(geom='segment',y= Inf,yend= Inf,x=as.POSIXct(-Inf),xend=as.POSIXct(Inf), size = .7*lwd_pt) +
      scale_x_continuous(breaks = seq(0,24,3), expand = expansion(mult = .03)) +
      scale_color_manual(values = c("y1" = "#1fb050", "y2" = "#ee7e32", "y3" = "#6f3996")) +
      guides(colour = "none") +
      labs(x = NULL, y = "Abundance (103 cells mL-1)", title = "Constant grazing") +
      Theme1; #plot_1a
    
    plot_1b <- ggplot(dat1, aes(x = t, y = b)) +
      geom_point() +
      scale_y_continuous(breaks = c(2.3,2.5,2.7,2.9)) +
      xlab(NULL) + ylab("Biomass (µg C L-1)") +
      Theme2; # plot_1b
    g1 <- plus_y_axis(plot_1a, plot_1b)
    #grid.newpage(); grid.draw(g1)
    
# Grazing rate = 0 --------------------------------------------------------
    g <- 0
    
    cn2 <- b0
    cc2 <- b0
    b2  <- b0
    # the input data for starting the model at 00:00
    cn2[1] <- 167.2394
    b2[1]  <- 2.580935
    cc2[1] <- b2[1] / cn2[1]
    
    for (k in 2:289){
      #k=2
      cn2[k] = cn2[k-1] * (1 + d[k-1]/12 - g/12)
      b2[k]  = b2[k-1]  * (1 + p[k-1]/12 - r[k-1]/12 - g/12)
      cc2[k] = b2[k]/cn2[k]
    }
    
    cs2 = (cc2 * 1000 / 280 * 3 / 4 / pi) ^ (1/3) * 2 # cell size
    
    dat2 <- tibble(t = seq(0,24,1/12),
                   k = t * 12 + 1,
                   b  = b2,
                   cc = cc2,
                   cn = cn2,
                   cs = cs2) 
    
    # Set up an equal scale combination plot - dual or three y-axis
    p1 <- ggplot(dat2, aes(x = t, y = cn)) + geom_point() + xlab(NULL) + guides(colour = "none");#p1
    p2 <- ggplot(dat2, aes(x = t, y = cs)) + geom_point() + xlab(NULL) + guides(colour = "none");#p2
    p3 <- ggplot(dat2, aes(x = t, y = b))  + geom_point() + xlab(NULL) + guides(colour = "none");#p3
    range1 <- ggplot_build(p1)$layout$panel_params[[1]]$y.range
    range2 <- ggplot_build(p2)$layout$panel_params[[1]]$y.range
    range3 <- ggplot_build(p3)$layout$panel_params[[1]]$y.range
    k2 <- (range2[2] - range2[1]) / (range1[2] - range1[1])
    b2 <- (range1[2] * range2[1] - range1[1] * range2[2]) / (range1[2] - range1[1])
    k3 <- (range3[2] - range3[1]) / (range1[2] - range1[1])
    b3 <- (range1[2] * range3[1] - range1[1] * range3[2]) / (range1[2] - range1[1])
    
    plot_2a <- ggplot(dat2, aes(x = t)) +
      annotate("rect", xmin = -Inf, xmax = 6, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") +
      annotate("rect", xmin = 18, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") +
      geom_line(aes(y = cn, color = "y1"), size = lwd_pt, linetype = 2) +
      #geom_point(aes(y = cn, color = "y1"), pch = 15, size = .7) +
      scale_y_continuous(sec.axis = sec_axis(trans = ~ .* k2 + b2, name = "Cell size (mean ESD, µm)\n",
                                             breaks = c(.45,.47,.49,.51))) +
      geom_line( aes(y = (cs - b2) / k2, color = "y2"), size = lwd_pt, linetype = 2) +
      #geom_point(aes(y = (cs - b2) / k2, color = "y2"), pch = 15, size = .7) +
      geom_line( aes(y = (b  - b3) / k3, color = "y3"), size = lwd_pt, linetype = 2) +
      #geom_point(aes(y = (b  - b3) / k3, color = "y3"), pch = 15, size = .7) +
      annotate(geom='segment',y=-Inf,yend=-Inf,x=as.POSIXct(-Inf),xend=as.POSIXct(Inf), size = .7*lwd_pt) +
      annotate(geom='segment',y= Inf,yend= Inf,x=as.POSIXct(-Inf),xend=as.POSIXct(Inf), size = .7*lwd_pt) +
      scale_x_continuous(breaks = seq(0,24,3), expand = expansion(mult = .03)) +
      scale_color_manual(values = c("y1" = "#1fb050", "y2" = "#ee7e32", "y3" = "#6f3996")) +
      guides(colour = "none") +
      labs(x = "Time of day (normalized)", y = "Abundance (103 cells mL-1)", title = "Grazing = 0") +
      Theme1; #plot_2a
    
    plot_2b <- ggplot(dat2, aes(x = t, y = b)) +
      geom_point() +
      xlab("Time of day (normalized)") + ylab("Biomass (µg C L-1)") +
      Theme2; #plot_2b
    g2 <- plus_y_axis(plot_2a, plot_2b)
    #grid.newpage(); grid.draw(g2)
    
    library(cowplot)
    nn <- length(list.files("D:/Demo/1/", pattern=".pdf"))
    pdf(file = paste0("D:/Demo/1/", "plot_", nn, ".pdf"),  width = 3.74, height = 4.8)
    plot_grid(g1, g2, labels = c("a", "b"), align = 'h', ncol = 1)
    dev.off()
    
    