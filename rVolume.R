  
    library(lubridate); library(scales);library(tidyverse); library(ggpubr);
    
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
    
    Vol <- function(Dia) 4/3*pi*((Dia)/2)^3 # biovolume
    
    
    Pico0 <- read_csv("Data/Picoauto.csv")
    Sampletime <- read_csv("Data/Sampletime.csv") %>% 
      mutate(Time = ymd_hms(paste(Date0, Time0))) %>%
      select(Id, Time)
    
    Sampletime %>% count(Id) %>% filter(n > 1)
    Pico0 %>% anti_join(Sampletime, by = "Id")
    
    Pico <- Pico0 %>% 
      left_join(Sampletime, by = "Id") %>%
      arrange(Time) %>%
      group_by(Station) %>%
      mutate(t = time_length(interval(date(min(Time)), Time), "hour")) %>% 
      ungroup() %>%
      mutate(Syn_Biomass  = FSC_Abu_Biomass(SynFSC, Syn),
             Pro_Biomass  = FSC_Abu_Biomass(ProFSC, Pro),
             Peuk_Biomass = FSC_Abu_Biomass(PeukFSC, Peuk),
             Syn_Dia      = FSC_Dia(SynFSC),
             Pro_Dia      = FSC_Dia(ProFSC),
             Peuk_Dia     = FSC_Dia(PeukFSC)) %>% 
      select(Station, Id, Id1, Time, t, everything()) %>% 
      select(-contains("FSC")) %>% 
      rename(Syn_Abu = Syn, Pro_Abu = Pro, Peuk_Abu = Peuk) %>%
      gather(Pico, Value, Syn_Abu:Peuk_Dia) %>%
      drop_na(Value)
    
    
    
    Istation <- unique(Pico$Station)
    Ipico    <- unique(Pico$Pico)
    result  <- tibble(station= character(), pico      = character(),
                      p      = numeric(),   amplitude = numeric(),
                      Mesor  = numeric(),   phi       = numeric(), r2 = numeric())
    
    for(i in seq_along(Istation)){
      #i=1
      df0 <- Pico %>% filter(Station == Istation[i])
      par(mfrow = c(3,3), mar = c(4,4,2,4))
      
      for(j in seq_along(Ipico)){
        #j= 7
        df <- df0 %>% filter(Pico == Ipico[j])
        
        Ylim <- range(df$Value, na.rm = TRUE)
        Xlim <- range(df$t,     na.rm = TRUE)
        # lm fit
        lm1 <- with(df, lm(Value ~ cos(pi/12*t) + sin(pi/12*t)))
        #summary(lm1)
        ft    <- summary(lm1)$fstatistic
        p_val <- unname(1 - pf(ft[1], ft[2], ft[3]))
        
        coef1 <- coef(lm1)
        a <- coef1[2]
        b <- coef1[3]
        phase <- atan2(b, a) %% (2*pi) # refer to package HarmonicRegression
        amplitude <- sqrt(a^2 + b^2)
        r2 <- summary(lm1)$r.squared
        Mesor <- coef1[1]
        
        df <- df %>% mutate(y1 = a * cos(pi/12 * t) + b * sin(pi/12 * t) + Mesor) 
        
        plot(Value ~ t, df
             ,"p", cex=.6, pch = 5
             ,xlab = "", ylab = "")
        lines(y1 ~ t, df, col = "blue", lty = 2, lwd = 4)
        title(main = paste0(Istation[i], "_", Ipico[j]), adj = 0)
        text(x = (Xlim[1] + Xlim[2])/2, y = (Ylim[1] * 2 + Ylim[2] * 8) / 10
             ,paste0("r2 = "   ,  round(r2, 2) 
                     ,", p = "  , format(signif(p_val, 3), scientific = TRUE),"\n"
                     ,"pt = ", round(phase /2/pi * 24, 1)), col = "green", cex = 2)
        
        result <- result %>%
          add_row(station= Istation[i], pico      = Ipico[j],
                  p      = p_val,       amplitude = amplitude,        
                  Mesor  = Mesor,       phi       = phase, r2 = r2)
      }
    }
    
    
    
    field <- result %>% 
      filter(pico != "Peuk_Biomass") %>%
      separate(pico, c("pico", "para")) %>% 
      mutate(peak      = Mesor + amplitude,
             trough    = Mesor - amplitude,
             rIncrease = 2* amplitude / trough,
             rVol      = ifelse(para == "Dia", Vol(peak) / Vol(trough) - 1, NA))

    # summarise average daily percent increase
    out2 <- field %>% 
      group_by(pico, para) %>% 
      summarise(m.rIncrease  = round(mean(rIncrease) * 100, 1),
                sd.rIncrease = round(sd(rIncrease) * 100, 1))
    
    # # biovolume daily percent increase
    # out3 <- field %>% 
    #   drop_na(rVol) %>% 
    #   group_by(pico) %>% 
    #   summarise(m.rVol  = round(mean(rVol) * 100, 1),
    #             sd.rVol = round(sd(rVol) * 100, 1))
  

    # the culture diel increase
    Lab  <- read_csv("Data/FSC_size.csv") %>% 
      select(Pico,Culture,Id,Size) %>% 
      pivot_wider(names_from  = Id,
                  values_from = Size) %>% 
      mutate(rSize = `3`/`1` - 1,
             rVol = Vol(`3`)/Vol(`1`) - 1) %>% 
      select(Pico, rVol)%>% 
      rename(pico = Pico) %>% 
      mutate(Id = "culture")

    
    dat <- field %>% 
      drop_na(rVol) %>% 
      select(pico,rVol) %>% 
      mutate(Id = "field") %>% 
      bind_rows(Lab) %>% 
      mutate(pico = factor(pico, levels = c("Pro", "Syn", "Peuk")))


    Icolor <- c("#f5e855","#6e5691")
    lwd_pt <- 1/(.pt*72.27/96)
    theme_set(theme_bw(base_size = 10, base_line_size = .4*lwd_pt, base_rect_size = .4*lwd_pt) +
                theme(plot.title = element_text(face = "bold", hjust = 0.5),
                      panel.border = element_rect(fill = NA, colour = "black", size = rel(2)),
                      axis.text  = element_text(colour = "black"),
                      axis.ticks = element_line(colour = "black"),
                      plot.background  = element_rect(colour = NA),
                      panel.grid =element_blank(),
                      legend.title = element_blank()))
    
    p <- ggplot(dat, aes(x = pico, y = rVol * 100, fill = Id)) +
      geom_bar(stat = "summary", fun = mean,  position = position_dodge(width = .82),width = .75) +
      geom_point(position = position_dodge(.75), size = 1.2, shape = 18) +
      stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
                   width = 0.25,position = position_dodge( .75), size = .4) +
      geom_hline(aes(yintercept=0),color = "grey", size = .3) +
      scale_fill_manual( values = Icolor) +
      scale_x_discrete(breaks = c("Pro", "Syn", "Peuk"),
                       labels = c("Prochlorococcus","Synechococcus","Picoeukaryotes")) +
      xlab(NULL) + ylab("Daily percent increase in biovolume") +
      ylim(-10.45,74) +
      theme(legend.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) 

    # Levene test
    library(car); 
    dat %>% split(.$pico) %>% map(~leveneTest(rVol ~ Id, data = .x, center = mean))
    # all p value > 0.05, so (var.equal = TRUE)
    
    p +　stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE), size=2.4)
    ggsave(filename = "Figure/rVol.pdf",  width = 4.3, height = 3.3)
    
    
    