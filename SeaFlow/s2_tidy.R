# ### The second step of data analysis of SeaFlow data v1.3 ####
# 
# Filtering and smoothing of raw data
# 
# 1. Data cleaning: (i) Cruises without at least one complete set of 24-h continuous observations;
#                   (ii) Discard Prochlorococcus abundance mean value < 80 * 10^3 cells/mL data;
#                   (iii) The abundance of Synechococcus or picoeukaryotes changes more than 5 times in 12 hours;
#                   (iv) Some obvious outliers in 5 cruises.
# 
# 2. Mesor-based normalization 
# 
# 3. Outlier identification and replacement 
# 
# 4. Applying rolling median functions based on time
# 
# 5. Reducing resolution of the data: To one value per 0.5h.

    
# setting ---------------------------------------------------
    
    library(tidyverse); library(lubridate); library(scales); 
    library(grid);  library(gtable);  library(egg);

    #' Time-Based Rolling Mesor (computed using cosinor method)
    #' 
    #' align = "center"
    #' 
    #' rolling window: at least 24 hours. The entire rolling window is divided into 48 time periods. If all values are NA 
    #' in a whole time period, the scrolling window will be expanded until 48 time periods containing data are included.
    #' 
    #' @author Changlin Li 
    #' 
    #' @param t time
    #' @param h time (half an hour)
    #' @param y value
    #' 
    #'
    #' @return An object of the same class as y with the rolling mesor
    #' 
    norm_mesor = function(t, h, y) {
      #t = df$t; h = df$h; y = df$value;
      mesor = numeric(length(y))
      
      for( i in seq_along(y) ) {
        #if(is.na(y[i])) {mesor[i] <- NA;len[i] <- NA; next}
        for(n in 24:240){
          #i = 1; n = 24
          window = h >= (h[i] - n/2) & h <= (h[i] + n/2); #h[window]
          ifelse(length(unique(h[window])) < 48, next, break)
        }
        lm1 <- lm(y[window] ~ cos(pi/12*t[window]) + sin(pi/12*t[window]))
        mesor[i] <- coef(lm1)[1] # summary(lm1)
      }
      return(mesor)
    }
    

    #' Time-Based Rolling Median
    #' 
    #' align = "center"
    #' 
    #' @author Changlin Li 
    #' 
    #' @param t time
    #' @param y value
    #' @param width time width (single) of the rolling window: minute
    #'
    #' @return An object of the same class as y with the rolling median
    #' 
    rollmedian = function(t, y, width) {
      out = numeric(length(t))
      for( i in seq_along(t) ) {
        window = t >= (t[i] - width * 60) & t <= (t[i] + width * 60); #window
        out[i] = median( y[window], na.rm = TRUE)
      }     
      return(out)
    }
    
    
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

    
    setwd("./SeaFlow")
    Sun  <- read_csv('Data/Solar.csv') 
    dat0 <- read_csv('Data/SeaFlow_data.csv') %>% 
      select(-contains("Qc_"), -contains("croco"), -depth, -lat, -lon) %>% 
      rename(Abu_Pro = abundance_prochloro, Abu_Syn = abundance_synecho, Abu_Euk = abundance_picoeuk,
             Dia_Pro = diam_prochloro,      Dia_Syn = diam_synecho,      Dia_Euk = diam_picoeuk,
             Bio_Pro = biomass_prochloro,   Bio_Syn = biomass_synecho,   Bio_Euk = biomass_picoeuk)

    # left_join two tibble
    Sun %>% count(cruise, time) %>% filter(n > 1)
    dat0 %>% anti_join(Sun, by = c("cruise", "time"))
# Data-cleaning #(i) ------------------------------------------------
    
    #Data-cleaning: filter 2 cruises without Prochlorococcus data
    dat1 <- dat0 %>% 
      left_join(Sun, by = c("cruise", "time")) %>% 
      filter(!cruise %in% c("TN248", "TN280")) %>% 
      arrange(cruise, time) 
   
    # draw plots
    Sys.setlocale(category = 'LC_ALL', locale = 'English')
    theme_set(theme_bw(base_size = 10))

    # color refer to red-green-black.pdf (www.nature.com/doifinder/10.1038/nature07832)
    Theme1 <- theme(axis.text.y.left   = element_text(color = "#ee1f24"),
                    axis.ticks.y.left  = element_line(color = "#ee1f24"),
                    axis.line.y.left   = element_line(color = "#ee1f24", size = rel(2)), # the left axis will be cover half when gtable_frame
                    axis.title.y.left  = element_text(color = "#ee1f24"),

                    axis.text.y.right  = element_text(color = "#00a94a"),
                    axis.ticks.y.right = element_line(color = "#00a94a"),
                    axis.line.y.right  = element_line(color = "#00a94a"),
                    axis.title.y.right = element_text(color = "#00a94a"),
                    panel.grid         = element_blank(),
                    panel.border       = element_blank(),
                    plot.title         = element_text(size = rel(.9), margin = margin(-3,0,1,0)))
    Theme2 <- theme(axis.text.y   = element_text(color = "black"),
                    axis.ticks.y  = element_line(color = "black"),
                    axis.line.y   = element_line(color = "black"),
                    axis.title.y  = element_text(color = "black", size = rel(.15)))


    Plot <- function(i, Para) {
      # i = 2; Para = "Dia"
      df <- dat_p %>% filter(cruise == Icruise[i] & para == Para)
      xlim1 <- min(df$time_n)
      xlim2 <- max(df$time_n)
      Night1 <- seq.POSIXt(date(xlim1) - hours(6), date(xlim2) + hours(18), by = "24 hours")
      Night2 <- seq.POSIXt(date(xlim1) + hours(6), date(xlim2) + hours(30), by = "24 hours")

      # Set up an equal scale combination plot - dual or three y-axis
      p1 <- ggplot(df, aes(x = time_n, y = Pro)) + geom_point() + xlab(NULL) + guides(colour = "none");#p1
      p2 <- ggplot(df, aes(x = time_n, y = Syn)) + geom_point() + xlab(NULL) + guides(colour = "none");#p2
      p3 <- ggplot(df, aes(x = time_n, y = Euk)) + geom_point() + xlab(NULL) + guides(colour = "none");#p3
      range1 <- ggplot_build(p1)$layout$panel_params[[1]]$y.range
      range2 <- ggplot_build(p2)$layout$panel_params[[1]]$y.range
      range3 <- ggplot_build(p3)$layout$panel_params[[1]]$y.range
      k2 <- (range2[2] - range2[1]) / (range1[2] - range1[1])
      b2 <- (range1[2] * range2[1] - range1[1] * range2[2]) / (range1[2] - range1[1])
      k3 <- (range3[2] - range3[1]) / (range1[2] - range1[1])
      b3 <- (range1[2] * range3[1] - range1[1] * range3[2]) / (range1[2] - range1[1])


      plot_1a <- ggplot(df, aes(x = time_n)) +
        annotate("rect", xmin = Night1, xmax = Night2, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") +
        geom_point(aes(y = Pro, color = "y1"), size= .5, alpha = .8) +
        scale_y_continuous(sec.axis = sec_axis(trans = ~ .* k2 + b2)) +  #, name = "Syn"
        geom_point(aes(y = (Syn - b2) / k2, color = "y2"), size= .5, alpha = .8) +
        geom_point(aes(y = (Euk - b3) / k3, color = "y3"), size= .5, alpha = .8) +
        annotate(geom='segment',y=-Inf,yend=-Inf,x=as.POSIXct(-Inf),xend=as.POSIXct(Inf), size = rel(.7)) +
        annotate(geom='segment',y= Inf,yend= Inf,x=as.POSIXct(-Inf),xend=as.POSIXct(Inf), size = rel(.7)) +
        scale_x_datetime(labels = date_format("%b %d"), expand = expansion(mult = .03)) +
        scale_color_manual(values = c("y1" = "#ee1f24", "y2" = "#00a94a", "y3" = "black")) +
        guides(colour = "none") +
        coord_cartesian(xlim = c(xlim1, xlim2)) +
        labs(x = NULL, y = NULL) +
        Theme1

      plot_1b <- ggplot(df, aes(x = time_n, y = Euk)) +
        geom_point() +
        xlab(NULL) + ylab("") +
        coord_cartesian(xlim = c(xlim1, xlim2)) +
        Theme2; #plot_1b
      g1 <- plus_y_axis(plot_1a, plot_1b)
      # grid.arrange(g1)
      g1 <- gtable_frame(g1)
    }

    
    dat_p <- dat1 %>%
      pivot_longer(cols = Abu_Pro:Bio_Euk,
                   names_to = c("para", ".value"),
                   names_pattern = "(.*)_(.*)")
    # Icruise <- unique(dat_p$cruise)
    # P.abu <- lapply(1:length(Icruise), Plot, Para = "Abu")
    # P.dia <- lapply(1:length(Icruise), Plot, Para = "Dia")
    # P.bio <- lapply(1:length(Icruise), Plot, Para = "Bio")
    # 
    # for (i in seq_along(Icruise)){
    #   #i = 3
    #   pdf(file = paste0("./Fig/s2-0_original/", Icruise[i], ".pdf"), width = 16, height = 8)
    #   grid.newpage()
    #   combined <- gtable_rbind(P.abu[[i]], P.dia[[i]], P.bio[[i]])
    #   grid.draw(combined)
    #   dev.off()
    # }


# Data-cleaning #(i) & #(ii) ---------------------------------------------------
    
    # Data-cleaning #2: discard 7 cruises
    # Coastal (*Prochlorococcus* abundance mean value < 80 data): CN12ID, CN13ID, (TN248,TN280 -- see code before)
    # cruises without a complete 24-h continuous observation:     KOK1512, KM1510,KM1603,KOK1803, KM1915
    
    dat1 <- dat1 %>% 
      filter(!cruise %in% c("CN12ID","CN13ID","KM1510","KM1603","KOK1512","KM1915","KOK1803"))
   
    # Data-cleaning #3: Discard *Prochlorococcus* abundance mean value < 80 data segment
    Icruise <- unique(dat1$cruise)
    i=1
    df <- dat1 %>% filter(cruise == Icruise[i]) %>%
      select(cruise, time, time_n, Abu_Pro, Abu_Syn, Abu_Euk) %>%
      arrange(time); min(df$Abu_Pro);i =i+1
    
    
    dat1 <- dat1 %>% 
      filter(!( (cruise == "CN11ID" &(time < ymd_hms("2011-10-02 18:50:00")  | time > ymd_hms("2011-10-07 18:59:00"))) |
                (cruise == "KM1502" & time < ymd_hms("2015-03-22 03:22:00")) |
                (cruise == "KM1712" & time > ymd_hms("2017-08-16 16:35:00")) |
                (cruise == "KM1713" & time < ymd_hms("2017-09-19 03:01:00")) |
                (cruise == "KM1906" & time < ymd_hms("2019-04-26 21:32:00")) |
                (cruise == "KOK1606"& time > ymd_hms("2016-04-26 08:50:00")  & time < ymd_hms("2016-04-28 15:00:00"))|
                (cruise == "MGL1704"& time > ymd_hms("2017-05-31 16:09:00")  & time < ymd_hms("2017-06-08 06:12:00"))|
                (cruise == "TN271"  & time < ymd_hms("2011-10-26 22:02:00")) |
                (cruise == "TN292"  & time < ymd_hms("2013-03-08 22:33:00")) |
                (cruise == "Tokyo_3"& time < ymd_hms("2011-09-20 09:35:00")) |  
                (cruise == "Tokyo_3"& time > ymd_hms("2011-09-25 12:02:00")  & time < ymd_hms("2011-09-30 07:32:00")) ))
    
    
    # Pro value of cruise KM1709 is a mixture of two groups, one is ~200 and the other is ~150; the latter to NA.
    KM1709 <- dat1 %>% 
      filter(cruise == "KM1709") %>% 
      mutate(pr_abu = predict(loess(Abu_Pro ~ as.numeric(time), span = 0.1))) 
    #ggplot(KM1709, aes(time, Abu_Pro)) + geom_point() + geom_line(aes(y = pr_abu))
    KM1709 <- KM1709 %>% 
      mutate(Abu_Pro = ifelse(Abu_Pro > pr_abu, Abu_Pro, NA),
             Bio_Pro = ifelse(Abu_Pro > pr_abu, Bio_Pro, NA)) %>% 
      select(-pr_abu)
    
    
    dat1 <- dat1 %>% 
      filter(cruise != "KM1709") %>% 
      bind_rows(KM1709) %>% 
      arrange(cruise, time) 

    # # check plot again
    # dat_p <- dat1 %>%
    #   pivot_longer(cols = Abu_Pro:Bio_Euk,
    #                names_to = c("para", ".value"),
    #                names_pattern = "(.*)_(.*)")
    # Icruise <- unique(dat_p$cruise)
    # 
    # P.abu <- lapply(1:length(Icruise), Plot, Para = "Abu")
    # P.dia <- lapply(1:length(Icruise), Plot, Para = "Dia")
    # P.bio <- lapply(1:length(Icruise), Plot, Para = "Bio")
    # 
    # for (i in seq_along(Icruise)){
    #   pdf(file = paste0("./Fig/s2-1_clean/", Icruise[i], ".pdf"), width = 16, height = 8)
    #   grid.newpage()
    #   combined <- gtable_rbind(P.abu[[i]], P.dia[[i]], P.bio[[i]])
    #   grid.draw(combined)
    #   dev.off()
    # }

# Data-cleaning #(iii) % #(iv) ---------------------------------------------------    
    
    dat1_check <- dat1 %>%
      select(cruise, time_n,  contains("Abu")) %>% 
      group_by(cruise) %>% 
      mutate(t = round_date(time_n, "12 hours")) %>% 
      ungroup() %>% 
      group_by(cruise, t) %>% 
      summarise(r_syn = max(Abu_Syn, na.rm = TRUE) / min(Abu_Syn, na.rm = TRUE),
                r_euk = max(Abu_Euk, na.rm = TRUE) / min(Abu_Euk, na.rm = TRUE))

    i = 1;
    dat1_check1 <- dat1_check %>% filter(cruise == Icruise [i]); i = i + 1
    
    # (iii) Dramatic fluctuations in the abundance of Synechococcus or picoeukaryotes.
    # We speculated that the dramatic fluctuation might arise from instrument malfunction or the passage into a distinctly different water mass.
    # We defined the condition of dramatic fluctuation as that the abundance of Synechococcus or picoeukaryotes changed more than 5 times within 12 hours.
    # However, most of the excluded time periods were longer than the duration with five times change, and they were determined by comprehensively 
    # considering the time series variations of the three parameters (abundance, cell size, and biomass) of three picophytoplankton populations.
    
    dat2 <- dat1 %>% 
      filter(!( (cruise == "CN11ID"    &(time < ymd_hms("2011-10-03 19:20:00") | time > ymd_hms("2011-10-06 07:33:00")))|
                (cruise == "FK180310-1"& time > ymd_hms("2018-03-23 03:57:00"))|
                (cruise == "KM1502"    & time < ymd_hms("2015-03-24 08:39:00"))|
                (cruise == "KM1508"    &(time < ymd_hms("2015-05-23 02:42:00") | time > ymd_hms("2015-05-26 09:36:00")))|
                (cruise == "KM1512"    & time < ymd_hms("2015-07-19 04:27:00"))|
                (cruise == "KM1513"    & time < ymd_hms("2015-07-25 15:15:00"))|
                (cruise == "KM1518"    & time > ymd_hms("2015-11-14 19:00:00") & time < ymd_hms("2015-11-14 21:00:00")) |
                (cruise == "KM1601"    & time < ymd_hms("2016-01-12 04:32:00") & Abu_Syn > 1.3) |
                (cruise == "KM1602"    & time > ymd_hms("2016-02-09 02:50:00") & time < ymd_hms("2016-02-09 03:10:00")) |
                (cruise == "KM1708"    & time > ymd_hms("2017-06-22 16:00:00"))| 
                (cruise == "KM1709"    & time < ymd_hms("2017-06-26 23:21:00"))| 
                (cruise == "KM1712"    & time > ymd_hms("2017-08-15 23:32:00"))|
                (cruise == "KM1712"    & time > ymd_hms("2017-08-12 20:00:00") & time < ymd_hms("2017-08-13 03:17:00")) |
                (cruise == "KM1713"    & time < ymd_hms("2017-09-19 22:16:00"))|
                (cruise == "KM1717"    & time < ymd_hms("2017-11-08 07:46:00"))|
                (cruise == "KM1821"    & time > ymd_hms("2018-11-19 06:50:00"))|
                (cruise == "KM1823"    & time > ymd_hms("2018-12-12 08:03:00") & time < ymd_hms("2018-12-12 13:00:00")) |
                (cruise == "KM1901"    & time > ymd_hms("2019-01-18 15:03:00"))|
                (cruise == "KM1903"    & time > ymd_hms("2019-02-19 02:32:00") & time < ymd_hms("2019-02-19 04:48:00")) |
                (cruise == "KM1906"    & time < ymd_hms("2019-04-27 02:16:00"))|
                (cruise == "KM1909"    & time < ymd_hms("2019-06-11 07:03:00"))|
                (cruise == "KM1912"    &(time < ymd_hms("2019-07-01 06:30:00") | time > ymd_hms("2019-07-03 16:00:00")))|
                (cruise == "KM1917"    &(time < ymd_hms("2019-09-04 04:23:00") | time > ymd_hms("2019-09-07 05:25:00")))|
                (cruise == "KN210-04"  & time > ymd_hms("2013-03-30 04:07:00") & time < ymd_hms("2013-03-31 00:43:00")) |
                (cruise == "KN210-04"  & time > ymd_hms("2013-04-04 00:00:00") & time < ymd_hms("2013-04-05 02:54:00")) |
                (cruise == "KN210-04"  & time > ymd_hms("2013-04-17 06:11:00") & time < ymd_hms("2013-04-18 03:31:00")) |
                (cruise == "KN210-04"  & time > ymd_hms("2013-04-22 11:23:00") & time < ymd_hms("2013-04-23 00:50:00")) |
                (cruise == "KN210-04"  & time > ymd_hms("2013-04-24 23:43:00"))|
                (cruise == "KOK1515"   & time > ymd_hms("2015-10-14 08:40:00"))|
                (cruise == "KOK1604"   &(time < ymd_hms("2016-04-13 23:31:00") | time > ymd_hms("2016-04-16 17:01:00")))|
                (cruise == "KOK1606"   & time < ymd_hms("2016-04-20 00:50:00"))|
                (cruise == "KOK1606"   & time > ymd_hms("2016-04-21 19:20:00") & time < ymd_hms("2016-04-22 06:15:00")) |
                (cruise == "KOK1606"   & time > ymd_hms("2016-04-24 00:52:00") & time < ymd_hms("2016-05-01 10:24:00")) |
                (cruise == "KOK1607"   & time < ymd_hms("2016-05-12 18:46:00"))|
                (cruise == "KOK1608"   & time < ymd_hms("2016-07-11 04:11:00"))|
                (cruise == "KOK1806"   &(time < ymd_hms("2018-07-15 13:30:00") | time > ymd_hms("2018-07-17 04:15:00")))|
                (cruise == "KOK1807"   & time < ymd_hms("2018-07-24 03:07:00"))|
                (cruise == "SR1917"    &(time < ymd_hms("2019-11-07 07:40:00") | time > ymd_hms("2019-11-26 02:34:00")))|
                (cruise == "TN271"     &(time < ymd_hms("2011-10-28 21:20:00") | time > ymd_hms("2011-11-02 21:30:00")))|
                (cruise == "TN292"     & time < ymd_hms("2013-03-09 21:16:00")))) %>% 
      filter(!(cruise %in% c("MGL1704", "Tokyo_3"))) # After data cleaning (ii) and (iii), data of the two cruises is very short
    
    
    # removed some obvious outliers.
    dat2 <- dat2 %>% 
      filter(!( (cruise == "FK180310-2"& time > ymd_hms("2018-04-04 11:37:00") & time < ymd_hms("2018-04-04 11:44:00")) |
                (cruise == "KM1709"    & time > ymd_hms("2017-07-08 22:49:00") & time < ymd_hms("2017-07-08 23:20:00")) |
                (cruise == "KM1802"    & time > ymd_hms("2018-01-17 01:43:00") & time < ymd_hms("2018-01-17 03:05:00")) |
                (cruise == "KM1903"    & time > ymd_hms("2019-02-20 23:00:00") & time < ymd_hms("2019-02-21 00:10:00"))))

    
    # some Peuk value of cruise KM1508 is very irregular, to NA.
    dat2 <- dat2 %>% 
      mutate(Dia_Euk = ifelse(cruise == "KM1508" & (time < ymd_hms("2015-05-23 09:01:00") | time > ymd_hms("2015-05-26 01:33:00")),
                                        NA, Dia_Euk),
             Bio_Euk = ifelse(cruise == "KM1508" & (time < ymd_hms("2015-05-23 09:01:00") | time > ymd_hms("2015-05-26 01:33:00")),
                                        NA, Bio_Euk))
    dat2[dat2$cruise == "KM1508" & is.na(dat2$Dia_Euk),]    # n = 271 

    # some Peuk value of cruise KM1709 is very irregular, to NA.
    dat2 <- dat2 %>% 
      mutate(Abu_Euk = ifelse(cruise == "KM1709" & time > ymd_hms("2017-06-28 11:00:00") & time < ymd_hms("2017-06-28 19:10:00"),
                              NA, Abu_Euk))
    dat2[dat2$cruise == "KM1709" & is.na(dat2$Abu_Euk),]  # n = 161 
    

    # # check plot again
    # dat_p <- dat2 %>%
    #   pivot_longer(cols = Abu_Pro:Bio_Euk,
    #                names_to = c("para", ".value"),
    #                names_pattern = "(.*)_(.*)")
    # Icruise <- unique(dat_p$cruise)
    # 
    # P.abu <- lapply(1:length(Icruise), Plot, Para = "Abu")
    # P.dia <- lapply(1:length(Icruise), Plot, Para = "Dia")
    # P.bio <- lapply(1:length(Icruise), Plot, Para = "Bio")
    # 
    # for (i in seq_along(Icruise)){
    #   pdf(file = paste0("./Fig/s2-2_clean/", Icruise[i], ".pdf"), width = 16, height = 8)
    #   grid.newpage()
    #   combined <- gtable_rbind(P.abu[[i]], P.dia[[i]], P.bio[[i]])
    #   grid.draw(combined)
    #   dev.off()
    # }

# Mesor-based normalization -----------------------------------------------
    
    # 14 min
    dat3 <- dat2 %>% 
      select(-time, -lat, -lon) %>% 
      gather(para, value, Abu_Pro:Bio_Euk) %>% 
      drop_na(value) %>% 
      #filter(cruise == "") %>% 
      group_by(cruise, para) %>% 
      mutate(h0 = round_date(time_n, "30 minutes"),
             t  = time_length(interval(date(min(time_n)), time_n), "hour"),
             h  = time_length(interval(date(min(time_n)), h0), "minute") / 60,
             mesor = norm_mesor(t, h, value))  %>% 
      ungroup() %>% 
      mutate(final = value/mesor)
    
    # # draw plots
    # Plot1 <- function(df) {
    #   #df <- dat_p[[3]][[1]]
    #   
    #   # Set up an equal scale combination plot - dual y-axis
    #   p1 <- ggplot(df, aes(x = time_n, y = value)) + geom_point() + xlab(NULL) + guides(colour = "none");#p1
    #   p2 <- ggplot(df, aes(x = time_n, y = final)) + geom_point() + xlab(NULL) + guides(colour = "none");#p2
    #   range1 <- ggplot_build(p1)$layout$panel_params[[1]]$y.range
    #   range2 <- ggplot_build(p2)$layout$panel_params[[1]]$y.range
    #   k2 <- (range2[2] - range2[1]) / (range1[2] - range1[1])
    #   b2 <- (range1[2] * range2[1] - range1[1] * range2[2]) / (range1[2] - range1[1])
    #   y1 <- (range1[1] + (1-range2[1])/(range2[2]-1)*range1[2])*(range2[2]-1)/(range2[2]-range2[1]) # y = 1
    #   
    #   ggplot(df, aes(x = time_n)) +
    #     geom_point(aes(y = value), size= 1, color = "gray") +
    #     geom_line(aes(y = mesor), color = "red", lwd = 1.2) +
    #     geom_hline(yintercept = y1, color = "green", lwd = 1.2, linetype=2) +
    #     scale_y_continuous(sec.axis = sec_axis(trans = ~ .* k2 + b2)) +
    #     geom_point(aes(y = (final - b2) / k2), size= 1, color = "blue", alpha = .5) + #007e7e  #b997fb
    #     scale_x_datetime(labels = date_format("%b %d"), expand = expansion(mult = .03)) +
    #     guides(colour = "none") +
    #     labs(x = NULL, y = NULL)
    # }
    # 
    # 
    # # import the SCS2019 data   
    # dat_SCS <- read_csv("Data/dat_SCS.csv") %>% 
    #   mutate(final = value/mesor) %>% 
    #   select(-t) 
    # 
    # dat_p <- dat3 %>%
    #   select(-h0, -t, -h) %>% 
    #   bind_rows(dat_SCS) %>%
    #   group_by(cruise, para) %>%
    #   group_nest() %>%
    #   mutate(gplot = map(data, Plot1))
    # Icruise <- unique(dat_p$cruise)
    # 
    # library(cowplot);
    # for (i in seq_along(Icruise)){
    #   df <- dat_p %>% filter(cruise == Icruise[i])
    #   
    #   plot_grid(df$gplot[[1]], df$gplot[[2]], df$gplot[[3]],
    #             df$gplot[[4]], df$gplot[[5]], df$gplot[[6]],
    #             df$gplot[[7]], df$gplot[[8]], df$gplot[[9]],
    #             labels = c("Abu_Euk","Abu_Pro","Abu_Syn",
    #                        "Bio_Euk","Bio_Pro","Bio_Syn",
    #                        "Dia_Euk","Dia_Pro","Dia_Syn"), ncol = 1, align = 'v')
    #   ggsave(filename = paste0("./Fig/s2-3_norm/", Icruise[i], ".pdf"),  width = 19, height = 25)
    # }
    

# Identify and replace outliers in a time series --------------------------
    
    #tsoutliers: Cruises with large time intervals should be divided into segments and finally merged.
    dat3 <- dat3 %>% 
      mutate(cruise = case_when(
        cruise == "KOK1606"  & time_n < ymd_hms("2016-04-21 12:00:00") ~ "KOK1606_a",
        cruise == "KOK1606"  & time_n > ymd_hms("2016-04-21 12:00:00") &
                               time_n < ymd_hms("2016-04-26 12:00:00") ~ "KOK1606_b",
        cruise == "KOK1606"  & time_n > ymd_hms("2016-04-26 12:00:00") ~ "KOK1606_c",
        cruise == "KN210-04" & time_n < ymd_hms("2013-03-30 12:00:00") ~ "KN210-04_a",
        cruise == "KN210-04" & time_n > ymd_hms("2013-03-30 12:00:00") &
                               time_n < ymd_hms("2013-04-04 00:00:00") ~ "KN210-04_b",
        cruise == "KN210-04" & time_n > ymd_hms("2013-04-04 00:00:00") &
                               time_n < ymd_hms("2013-04-17 12:00:00") ~ "KN210-04_c",
        cruise == "KN210-04" & time_n > ymd_hms("2013-04-17 12:00:00") &
                               time_n < ymd_hms("2013-04-22 18:00:00") ~ "KN210-04_d",
        cruise == "KN210-04" & time_n > ymd_hms("2013-04-22 18:00:00") ~ "KN210-04_e",
        TRUE                                                           ~ cruise
      ))  %>%
      mutate(time1 = round_date(time_n, "3 minutes"))

    
    ## main loop over the data
    Icruise <- unique(dat3$cruise)
    Ipara   <- unique(dat3$para)

    for (i in seq_along(Icruise)){
      for(j in seq_along(Ipara)){
        #i=25; j = 4
        df0 <- dat3 %>% filter(cruise == Icruise[i] & para == Ipara[j])
        df1 <- tibble(time1 = with(df0, seq.POSIXt(min(time1), max(time1), by = "3 mins"))) %>%
          left_join(df0, by = "time1") %>%
          mutate(index = row_number())
        
        t <- zoo::zoo(df1$final, order.by = seq(min(df1$time1), max(df1$time1), by ="3 mins"))
        o <- forecast::tsoutliers(t) %>% as_tibble(); o
        
        if(dim(o)[1] > 0){
          df2 <- df1 %>% left_join(o, by = "index") %>% drop_na(final, replacements) 
          
          # replace outliers in dat3
          sel = dat3$cruise == Icruise[i] & dat3$para == Ipara[j] & dat3$time_n %in% df2$time_n
          dat3$final[sel] <- df2$replacements
          
          ggplot() +
            geom_point(data = df0, aes(time_n, final), color='gray') +
            geom_point(data = df2, aes(time_n, replacements), col = "red", pch = 17, size = 2.5) +
            geom_point(data = df2, aes(time_n, final), col = "blue", pch = 6) +
            xlab(NULL) + ylab(NULL) +
            labs(title = paste0(Icruise[i], "_", Ipara[j], ": outliers = ", dim(df2)[1]))
          # ggsave(filename = paste0("./Fig/s2-4_tsoutlier/", Icruise[i],"_", Ipara[j], ".pdf"), width = 14, height = 8)
        }
      }
    }
    
    dat3 <- dat3 %>%
      mutate(cruise = case_when(
        cruise %in% c("KOK1606_a", "KOK1606_b", "KOK1606_c")   ~ "KOK1606",
        cruise %in% c("KN210-04_a", "KN210-04_b","KN210-04_c","KN210-04_d","KN210-04_e") ~ "KN210-04",
        TRUE                                      ~ cruise
      )) %>%
      select(-h0, -t, -h, -mesor, -time1) 
   
# rolling smooth & Data resolution reduction -----------------------------------------------------
    
    # applying rolling median functions based on time (50 min rolling window) -- 3 min
    dat4 <- dat3 %>% 
      arrange(cruise, para, time_n) %>% 
      group_by(cruise, para) %>% 
      mutate(Final = rollmedian(time_n, final, 25)) %>%
      ungroup()
    
    dat4 %>% ungroup() %>%  count(time_n) %>% filter(n > 9)
    # A tibble: 0 x 2
    
    # Data resolution reduction: to one value per 0.5h
    dat_seaflow <- dat4 %>% 
      mutate(time_n = round_date(time_n, "30 minutes")) %>% 
      group_by(cruise, para, time_n) %>%
      summarise_at(vars(Final), median) %>% 
      ungroup() %>% 
      spread(para, Final)

    # # draw plots
    # dat_p <- dat_seaflow %>%
    #   pivot_longer(cols = Abu_Euk:Dia_Syn,
    #                names_to = c("para", ".value"),
    #                names_pattern = "(.*)_(.*)")
    # Icruise <- unique(dat_p$cruise)
    # 
    # P.abu <- lapply(1:length(Icruise), Plot, Para = "Abu")
    # P.dia <- lapply(1:length(Icruise), Plot, Para = "Dia")
    # P.bio <- lapply(1:length(Icruise), Plot, Para = "Bio")
    # 
    # for (i in seq_along(Icruise)){
    #   pdf(file = paste0("./Fig/s2-5_final/", Icruise[i], ".pdf"), width = 16, height = 8)
    #   grid.newpage()
    #   combined <- gtable_rbind(P.abu[[i]], P.dia[[i]], P.bio[[i]])
    #   grid.draw(combined)
    #   dev.off()
    # }
    write_csv(dat_seaflow,"Data/dat_seaflow.csv")  
  
  