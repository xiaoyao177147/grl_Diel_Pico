#### The third step of data analysis of SeaFlow data v1.3 ####
#
# Summarize the results and draw plots.
#
# 1. Plot: average half-hourly values
#
# 2. Cosinor analysis: 
#            (1) set the r2 value threshold at 0.36; 
#            (2) The diel pattern of cell size of cruise *KOK1604*/*KOK1515* is different with others, the raw data is wrong?
#                Email to Francois Ribalet (The SeaFlow dataset Maintainer), 
#                "Problem 1: It is possible that the time zone for this cruise (KOK1604) was not GMT like we assumed, 
#                we will investigate." The examination results have not yet been returned. 
#         
# 2. Summarize and plot. 
#
##
    
    library(tidyverse); library(lubridate); library(scales);
    
    setwd("./SeaFlow")
    dat3 <- read_csv("Data/dat_seaflow.csv")
    # import the SCS2019 data
    dat_SCS <- read_csv("Data/dat_SCS.csv") %>%
      mutate(final = value/mesor)
    
    
    dat3_1 <- dat3 %>% 
      select(cruise, time_n) %>% 
      group_nest(cruise) %>% 
      mutate(tday = map(data, ~seq.Date(date(min(.$time_n))- 1, date(max(.$time_n))+ 1, by = "days"))) %>% 
      select(-data) %>% 
      unnest(tday) %>% 
      mutate(time_n = map(tday, function(t) t + minutes(30) * 0:47)) %>% # Fill in the time gap: one value per half hour
      unnest(time_n) %>% 
      arrange(cruise, time_n) %>% 
      mutate(time_n = as_datetime(time_n))
    
    # left_join two tibble
    dat3 %>% count(cruise, time_n) %>% filter(n > 1)
    dat3_1 %>% anti_join(dat3, by = c("cruise", "time_n"))
    
    dat3 <- dat3_1 %>% 
      left_join(dat3, by = c("cruise", "time_n")) %>% 
      group_by(cruise) %>% 
      mutate(t = time_length(interval(min(time_n), time_n), "hour")) %>% 
      ungroup() %>% 
      gather(para, value, Abu_Euk:Dia_Syn) %>% 
      arrange(cruise, para, time_n)
    

# Draw mean + sd for half-hourly measurements -------------------------------
    # transform the SCS2019 data
    dat_SCS1 <- dat_SCS %>%
      mutate(time_n = round_date(time_n, "30 minutes")) %>%
      group_by(cruise) %>%
      mutate(t  = time_length(interval(date(min(time_n)), time_n), "hour")) %>%
      ungroup() %>%
      select(-value, -mesor) %>% 
      rename(value = final)
    
    dat4 <- dat3 %>% 
      bind_rows(dat_SCS1) %>% 
      mutate(t1 = t %% 24) %>% 
      drop_na(value) %>% 
      group_by(t1, para) %>% 
      summarise(n  = length(value),
                u  = mean(value),
                SD = sd(value)) %>% 
      ungroup()
    
    
    dat4_1 <- dat4 %>% filter(t1 == 0) %>% mutate(t1 = 24)
    dat4 <- dat4 %>% 
      bind_rows(dat4_1) %>% 
      separate(para, c("Para", "Pico")) %>% 
      mutate(Pico = factor(Pico, levels = c("Pro","Syn","Euk")),
             Para = factor(Para, levels = c("Abu","Dia","Bio")))
    
    
    library(ggh4x)
    lab.p   <- c(Abu = "Abundance",Dia = "Cell size (mean ESD)",Bio = "Biomass")
    lab.b1 <- c(Pro = "Prochlorococcus",Syn = "Synechococcus",Euk = "Picoeukaryotes")
    lwd_pt <- 1/(.pt*72.27/96)
    theme_set(theme_bw(base_size = 10, base_line_size = .5*lwd_pt, base_rect_size = .5*lwd_pt))
    
    p1 <- ggplot(dat4, aes(t1, u)) +
      annotate("rect", xmin = -Inf, xmax = 6, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") +
      annotate("rect", xmin = 18, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") +
      geom_errorbar(aes(ymin = u-SD, ymax = u+SD), width = 0, size = .5*lwd_pt, alpha = .5) +
      geom_point(aes(color = Pico), size = 1.5*lwd_pt) +
      geom_line(aes(color = Pico), size = .8*lwd_pt) +
      facet_grid2(vars(Pico), vars(Para),
                  scales = "free_y", independent = "y", switch = "y",
                  axes = "all", remove_labels = "x",
                  labeller=labeller(Para = as_labeller(lab.p),
                                    Pico = as_labeller(lab.b1))) +
      scale_x_continuous(breaks = c(0,6,12,18,24), expand = expansion(mult = .03)) +
      scale_y_continuous(breaks = breaks_pretty(3)) +
      ylab(NULL) + xlab('Time of day ("standard day")') +
      scale_color_manual(values = c("#1fb050","#1471b8","#ee7e32")) +
      guides(color = "none") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            panel.border = element_rect(fill = NA, colour = "black", size = rel(2)),
            axis.text  = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            axis.title = element_text(size = rel(0.9)),
            strip.text = element_text(size = rel(.9)),
            strip.text.x = element_text(face = "bold"),
            strip.background = element_rect(fill = "transparent", size = rel(1.5)),
            strip.background.y = element_rect(colour = NA),
            strip.placement.y = "outside",
            plot.background  = element_rect(colour = NA),
            panel.spacing.x=unit(.5, "lines"),
            panel.spacing.y=unit(.45, "lines"),
            panel.grid =element_blank(),
            plot.margin = unit(c(5,6,0,6), "pt"));p1
    # ggsave(filename = paste0("D:/Demo/mesor.pdf"),  width = 7.41, height = 4.8)
    
# analysis the data at "daily" level ----------------------------------------------

    Icruise <- unique(dat3$cruise)
    Ipara   <- unique(dat3$para)
    
    result  <- tibble(cruise= character(), tday  = POSIXct(), para  = character(),
                      na.d1 = numeric(),   na.d2 = numeric(), na.d3 = numeric(),
                      p     = numeric(),   a     = numeric(), b     = numeric(), 
                      Mesor = numeric(),   phi   = numeric(), r2    = numeric())
    
    for (i in seq_along(Icruise)){
      #i=26
      df0 <- dat3 %>% filter(cruise == Icruise[i])
      Iday <- unique(df0$tday)
      
      for(j in 2:(length(Iday)- 1)){
        #j=26
        df1 <- df0 %>% filter(tday %within% interval(Iday[j] - days(1), Iday[j] + days(1)))
        
        #pdf(file = paste0("./Fig/s3-1_cosinor/", Icruise[i],"_", Iday[j], ".pdf"), width = 34, height = 16.5, pointsize = 24)
        par(mfrow = c(3,3), mar = c(4,4,2,4))
        
        for(p in seq_along(Ipara)){
          #p= 7
          df <- df1 %>% filter(para == Ipara[p])
          na_num <- df %>% group_by(tday) %>% summarise(na = sum(is.na(value))) # missing value
          
          if(!all(is.na(df$value))) {
            Ylim <- range(df$value, na.rm = TRUE)
            Xlim <- range(df$t,     na.rm = TRUE)
            # lm fit
            lm1 <- with(df, lm(value ~ cos(pi/12*t) + sin(pi/12*t)))
            #summary(lm1)
            ft    <- summary(lm1)$fstatistic
            p_val <- unname(1 - pf(ft[1], ft[2], ft[3]))
            
            coef1 <- coef(lm1)
            a <- coef1[2]
            b <- coef1[3]
            phase <- atan2(b, a) %% (2*pi) # refer to package HarmonicRegression
            #amplitude <- sqrt(a^2 + b^2)
            r2 <- summary(lm1)$r.squared
            Mesor <- coef1[1]
            
            df <- df %>% mutate(y1 = a * cos(pi/12 * t) + b * sin(pi/12 * t) + Mesor) 
            
            plot(value ~ t, df
                 ,"p", cex=.6, pch = 5
                 ,xaxt = "n", xlab = "", ylab = "")
            lines(y1 ~ t, df, col = "blue", lty = 2, lwd = 4)
            
            axis(1, seq(min(df$t), max(df$t) + 0.5, 4), seq(0, 72, 4))
            title(main = paste0(Icruise[i], "_", Iday[j], "_", Ipara[p]), adj = 0)
            text(x = (Xlim[1] + Xlim[2])/2, y = (Ylim[1] * 2 + Ylim[2] * 8) / 10
                 ,paste0("r2 = "   ,  round(r2, 2) 
                         ,", p = "  , format(signif(p_val, 3), scientific = TRUE),"\n"
                         ,"pt = ", round(phase /2/pi * 24, 1)), col = "green", cex = 2)
            
            result <- result %>%
              add_row(cruise= Icruise[i],   tday  = Iday[j],      para  = Ipara[p],
                      na.d1 = na_num$na[1], na.d2 = na_num$na[2], na.d3 = na_num$na[3],
                      p     = p_val,        a     = a,            b     = b,    
                      Mesor = Mesor,        phi   = phase,        r2    = r2)
          } else {
            t <- 1:72
            plot(t, xaxt = "n", xlab = "", ylab = "")
            title(main = paste0(Icruise[i], "_", Iday[j], "_", Ipara[p]), adj = 0)
            text(x = 28, y = 40
                 ,paste0("NA"), col = "red", cex = 2)
            
            result <- result %>%
              add_row(cruise= Icruise[i],   tday     = Iday[j],      para  = Ipara[p],
                      na.d1 = na_num$na[1], na.d2    = na_num$na[2], na.d3 = na_num$na[3])
          }
        }
        #dev.off()
      }
    }
    

# analysis the SCS data at daily level similar the SeaFlow dataset --------
     
    seats_1 <- dat_SCS %>% 
      filter(time_n < ymd_hms("2019-06-22 01:00:00") & cruise == "SCS2019_SEATS") %>%
      mutate(tday = "2019-06-20")
    seats_2 <- dat_SCS %>% filter(cruise == "SCS2019_SEATS") %>%
      mutate(tday = "2019-06-21")
    m4_1    <- dat_SCS %>% 
      filter(cruise == "SCS2019_M4") %>% 
      mutate(tday = "2019-06-25")
    m4_2    <- dat_SCS %>% 
      filter(time_n > ymd_hms("2019-06-24 23:30:00") & cruise == "SCS2019_M4") %>% 
      mutate(tday = "2019-06-26")
    k11     <- dat_SCS %>% filter(cruise == "SCS2019_K11") %>%
      mutate(tday = "2019-07-01")
    
    # merge file
    dat_SCS2 <- seats_1 %>% 
      bind_rows(seats_2) %>% 
      bind_rows(m4_1) %>% 
      bind_rows(m4_2) %>% 
      bind_rows(k11) %>% 
      select(-mesor, -value) %>% 
      rename(value = final) %>% 
      mutate(tday = as_date(tday))
    
    
    Iday2 <- unique(dat_SCS2$tday)
    
    for(i in seq_along(Iday2)){
      #i=1
      df0 <- dat_SCS2 %>% filter(tday == Iday2[i])
      
      #pdf(file = paste0("./Fig/s3-1_cosinor/", df0$cruise[1],"_", Iday2[i], ".pdf"), width = 34, height = 16.5, pointsize = 24)
      par(mfrow = c(3,3), mar = c(4,4,2,4))
      
      for(p in seq_along(Ipara)){
        #p= 7
        df <- df0 %>% filter(para == Ipara[p])
        
        Ylim <- range(df$value, na.rm = TRUE)
        Xlim <- range(df$t,     na.rm = TRUE)
        # lm fit
        lm1 <- with(df, lm(value ~ cos(pi/12*t) + sin(pi/12*t)))
        #summary(lm1)
        ft    <- summary(lm1)$fstatistic
        p_val <- unname(1 - pf(ft[1], ft[2], ft[3]))
        
        coef1 <- coef(lm1)
        a <- coef1[2]
        b <- coef1[3]
        phase <- atan2(b, a) %% (2*pi) # refer to package HarmonicRegression
        #amplitude <- sqrt(a^2 + b^2)
        r2 <- summary(lm1)$r.squared
        Mesor <- coef1[1]
        
        df <- df %>% mutate(y1 = a * cos(pi/12 * t) + b * sin(pi/12 * t) + Mesor) 
        
        plot(value ~ t, df
             ,"p", cex=.6, pch = 5
             ,xlab = "", ylab = "")
        lines(y1 ~ t, df, col = "blue", lty = 2, lwd = 4)
        
        title(main = paste0(df0$cruise[1], "_", Iday2[i], "_", Ipara[p]), adj = 0)
        text(x = (Xlim[1] + Xlim[2])/2, y = (Ylim[1] * 2 + Ylim[2] * 8) / 10
             ,paste0("r2 = "   ,  round(r2, 2) 
                     ,", p = "  , format(signif(p_val, 3), scientific = TRUE),"\n"
                     ,"pt = ", round(phase /2/pi * 24, 1)), col = "green", cex = 2)
        
        result <- result %>%
          add_row(cruise= df0$cruise[1], tday = Iday2[i], para = Ipara[p],
                  na.d2 = 0, 
                  p     = p_val,         a    = a,        b    = b,
                  Mesor = Mesor,         phi  = phase,    r2   = r2)
      }
      #dev.off()
    }
    #write_csv(result,"Data/result.csv")

# summarise the resule and draw plots ---------------------------------------------------
    
    #library(lubridate); library(scales); library(tidyverse); 
    #setwd("./SeaFlow")
    #result <- read_csv("Data/result.csv")
    
    result1 <- result %>% 
      drop_na(p) %>% 
      filter(na.d2 <= 20) %>% 
      filter(cruise != "KOK1604") %>% 
      filter(cruise != "KOK1515") %>% 
      mutate(peak_time = phi/2/pi * 24) %>% 
      mutate(diel = if_else(round(r2, 2) >= .36, 1, 0))
    
    # Calculate the cruises proportion with significant 24h periodicity
    out1 <- result1 %>% 
      group_by(para) %>% 
      summarise(n       =length(diel),
                n.diel  = sum(diel),
                n.ratio = round(n.diel / n * 100, 1))

    message("draw plots") 
    library(circular);library(bpDir);# Attaching package: ‘circular’ The following objects are masked from ‘package:stats’: sd, var
    # Firstly, draw by R package circular and bpDir;
    # Extract main statistics and outliers from CircularBoxplot for ggplot2, subsequently check them by plots and outliers (number); 
    # Refer the method: https://stackoverflow.com/questions/30078797/how-to-mimic-geom-boxplot-with-outliers-using-geom-boxplotstat-identity
    
    cir_df <- result1 %>% 
      filter(diel == 1) %>% 
      select(para, peak_time)
    cir_df <- split(cir_df, cir_df$para)
    # calculate a data.frame with the main statistics and a separated data.frame with the outliers
    # whether the number of outliers is equal to the CircularBoxplot function
    cir_bp <- tibble(para          = character(),
                     lower.whisker = numeric(),
                     lower.hinge   = numeric(),
                     median        = numeric(),
                     upper.hinge   = numeric(),
                     upper.whisker = numeric(),
                     Outliers      = logical())
    
    for(d in cir_df){
      #d  <- cir_df$Abu_Pro
      pt <- d$peak_time
      A  <- circular(pt, units = "hours", template = "clock24", modulo = "2pi")
      
      # output images -- draw by R package circular and bpDir
      #pdf(file = paste0("./Fig/s3-2_circular/", d$para[1], ".pdf"), width = 15, height = 9, pointsize = 16)
      par(mfrow = c(1,2), mar = c(4,4,2,4))
      plot(A, col = "blue", main = d$para[1])
      mtext(side=1, line=2, at=0.07, adj=0, cex=1.7, "draw by R package circular and bpDir")
      bp.cir <- CircularBoxplot (A, template = NULL, H = TRUE, stack = TRUE,
                                 place = "outside", constant = 1.5)
      #dev.off()
      

      quan.cir <- quantile.circular(A);
      stats    <-  c(quan.cir[[5]],quan.cir[[4]],quan.cir[[3]],quan.cir[[2]],quan.cir[[1]])

      iqr      <- diff(stats[c(2, 4)])  
      outliers <- pt < (stats[2] - 1.5 * iqr) | pt > (stats[4] + 1.5 * iqr)
       sort(pt);quan.cir;stats
      
      if (any(outliers)) stats[c(1, 5)] <- range(c(stats[2:4], pt[!outliers]), na.rm=TRUE)
      
      outlier_values = pt[outliers]; #outlier_values
      if (length(outlier_values) == 0) outlier_values <- NA_real_
      
      cir_bp <- cir_bp %>% 
        add_row(para          = d$para[1],
                lower.whisker = stats[1],
                lower.hinge   = stats[2],
                median        = stats[3],
                upper.hinge   = stats[4],
                upper.whisker = stats[5],
                Outliers      = length(outlier_values) == length(bp.cir$farout))
      # for Abu_Pro, should deal with it specially
      if(d$para[1] == "Abu_Pro") { 
        lw <- stats[2] - 1.5 * iqr + 24
        cir_bp$lower.whisker[cir_bp$para == "Abu_Pro"] <- min(pt[pt>lw])
      }
    }
    # for Bio_euk, the lower.whisker derived from bpDir is different from quan.cir, because:
    # bp.cir$statistics$CounterClockwiseHinge, all are integers, so I think the IQR from 
    #  (bp.cir$statistics$ClockwiseHinge - bp.cir$statistics$CounterClockwiseHinge) / 360 * 24 is wrong;
    # in the Package source -- CircularBoxplot: # defining the whiskers
    # d <- (rad(round(deg(range(circular(IQR, modulo = "2pi")))))) # I think radian --> degree --> round --> radian!!!
    # email to the maintainer Davide Buttarazzi, he agree
    # "I'm afraid it is indeed an unfortunate rounding effect which causes this. By transforming hours to radians/degrees 
    # back and forth there is a loss of precision. I could not avoid using this workaround of multiple conversions because 
    # in some special cases R gave me erratic behaviour. It's likely that this is fixed in the meantime compared to 
    # last time I worked on the package."
    
    options(digits.secs = 2)
    out3 <- cir_bp %>% 
      select(para, lower.hinge, median, upper.hinge) %>% 
      mutate(lower.hinge1 = strftime(Sys.Date() + seconds(lower.hinge) *3600, format = "%H:%M:%OS", tz = "UTC"),
             median1      = strftime(Sys.Date() + seconds(median)      *3600, format = "%H:%M:%OS", tz = "UTC"),
             upper.hinge1 = strftime(Sys.Date() + seconds(upper.hinge) *3600, format = "%H:%M:%OS", tz = "UTC"))
    
    # Create a new virtual tibble to determine the relative positions of the three groups
    X <- tibble(pico = c("Euk","Syn","Pro"),
                x    = c(34,26,18))
    result2 <- result1 %>% 
      filter(diel == 1) %>% 
      select(para, peak_time) %>% 
      separate(para,c("para","pico")) %>% 
      left_join(X, by = "pico") %>% 
      mutate(para  = factor(para, levels = c("Abu", "Dia", "Bio")),
             pico  = factor(pico, levels = c("Pro", "Syn", "Euk"))) 
    cir_bp <- cir_bp %>% 
      separate(para,c("para","pico")) %>% 
      left_join(X, by = "pico") %>% 
      mutate(para  = factor(para, levels = c("Abu", "Dia", "Bio")),
             pico  = factor(pico, levels = c("Pro", "Syn", "Euk"))) 
    
    
    # For Abu_Pro, the peak_time crosses 24 o'clock and needs special treatment.
    pro_abu <- cir_bp %>% filter(para == "Abu" & pico == "Pro")
    cir_bp$lower.whisker[cir_bp$para == "Abu" & cir_bp$pico == "Pro"] <- NA
    #cir_bp$lower.whisker[cir_bp$para == "Abu" & cir_bp$pico == "Pro"] <- pro_abu$lower.whisker
    
    library(ggbeeswarm)
    lwd_pt <- 1/(.pt*72.27/96) 
    # theme_set(theme_bw(base_size = 10, base_line_size = .4*lwd_pt, base_rect_size = .4*lwd_pt))
    
    ann_text <- tibble(para = c("Abu", "Dia", "Bio")) %>% 
      mutate(para  = factor(para, levels = c("Abu", "Dia", "Bio")))
    
   p2 <- ggplot() + 
      annotate("segment", x=39, y=seq(0,24,2), xend=37.8, yend=seq(0,24,2),   col="black", size = .5*lwd_pt)+
      annotate("segment", x=39, y=seq(1,24,2), xend=38.2, yend=seq(1,24,2),   col="black", size = .5*lwd_pt)+
      geom_vline(xintercept  = c(14,39), colour = "black", size = .5*lwd_pt) +
      geom_vline(xintercept  = c(22,30), colour = "black", size = .5*lwd_pt, linetype = "dashed") +
      geom_quasirandom(data = result2, aes(x, peak_time, color = pico),
                       size = .6, shape = 20, alpha = .5, varwidth=TRUE) +
      geom_boxplot(data = cir_bp,
                   aes(x, ymin  = lower.whisker, lower = lower.hinge, middle = median,
                       upper = upper.hinge,   ymax  = upper.whisker, color = pico),
                   stat = "identity" , lwd = .3, fill = "transparent", width = 6) +
      geom_segment(data = cir_bp, aes(x = x-1.5, xend = x+1.5, y = lower.whisker, yend = lower.whisker, color=pico),
                   size = .3) +
      geom_segment(data = cir_bp, aes(x = x-1.5, xend = x+1.5, y = upper.whisker, yend = upper.whisker, color=pico),
                   size = .3) +
      geom_segment(data = pro_abu,aes(x = x-1.5, xend = x+1.5, y = lower.whisker, yend = lower.whisker, color=pico),
                   size = .3) +
      geom_segment(data = pro_abu, aes(x = x, xend = x, y = lower.whisker, yend = 24, color = pico), size = .3) +
      geom_segment(data = pro_abu, aes(x = x, xend = x, y = 0, yend = lower.hinge, color = pico), size = .3) +
      geom_text(data = ann_text,x=0, y=0, label = c("Abundance", "Cell size \n(mean ESD)", "Biomass"),size = 9/.pt) +
      geom_text(data = ann_text,x=41.8,y=0, label = "24/0",size = 8/.pt) +
      geom_text(data = ann_text,x=41.8,y=2, label = "2",   size = 8/.pt) +
      geom_text(data = ann_text,x=41.8,y=4, label = "4",   size = 8/.pt) +
      geom_text(data = ann_text,x=41.8,y=6, label = "6",   size = 8/.pt) +
      geom_text(data = ann_text,x=41.8,y=8, label = "8",   size = 8/.pt) +
      geom_text(data = ann_text,x=41.8,y=10,label = "10",  size = 8/.pt) +
      geom_text(data = ann_text,x=41.8,y=12,label = "12",  size = 8/.pt) +
      geom_text(data = ann_text,x=41.8,y=14,label = "14",  size = 8/.pt) +
      geom_text(data = ann_text,x=41.8,y=16,label = "16",  size = 8/.pt) +
      geom_text(data = ann_text,x=41.8,y=18,label = "18",  size = 8/.pt) +
      geom_text(data = ann_text,x=41.8,y=20,label = "20",  size = 8/.pt) +
      geom_text(data = ann_text,x=41.8,y=22,label = "22",  size = 8/.pt) +
      facet_grid(.~ para) +
      xlab(NULL) + ylab(NULL) +
      coord_polar(theta = "y") +
      scale_x_continuous(limits = c(0, 39.5)) +
      scale_y_continuous(limits = c(0, 24), breaks = seq(0,24,4)) +
      scale_color_manual(values = c("#1fb050","black","#ee7e32"),
                         breaks = c("Pro","Syn","Euk"),
                         labels = c("Prochlorococcus","Synechococcus","Picoeukaryotes")) +
      theme(plot.background  = element_rect(colour = NA),
            axis.text = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(colour = "black"),
            strip.background = element_blank(),
            strip.text = element_blank(),
            panel.spacing.x = unit(-0.8, "lines"),
            panel.border  = element_blank(),
            panel.grid = element_blank(),
            legend.position = c(.67, .12),
            legend.title = element_blank(),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.key.width  = unit(.33, "cm"),
            legend.key.height = unit(.33, "cm"),
            plot.margin = unit(c(0,-6,0,-6), "pt")) #;p2 # pdf 7.48 * 4.3

    
    library(cowplot)
    plot_grid(p1, p2, labels = c("a", "b"), align = 'V', ncol = 1, rel_widths = c(1, 2))

    nn <- length(list.files("D:/Demo/1/", pattern=".pdf"))
    ggsave(filename = paste0("D:/Demo/1/", "plot_", nn, ".pdf"),  width = 7, height = 10)
    # At 100% magnification, the figure below first deletes 1+3 outer borders and shift++up 11 times.
    