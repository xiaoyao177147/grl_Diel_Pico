#### Same data preprocessing of SeaFlow data v1.3 was applied to SCS2019 cruise ####
#
# 1. s1_Ntime: Calculate sunrise and sunset times; Time normalization;
#
# 2. s2_tidy: Mesor-based normalization;
# 
##

    library(tidyverse); library(lubridate); 
    
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
    
    #' Calculate the corresponding Sunrise and Sunset Times for each time point  
    #' 
    #' @author Changlin Li
    #' 
    #' @param Lat/Lon longitude and latitude
    #' @param t time
    #'
    #' @return variable sunrise, sunset times of each data
    #' 
    library(suncalc)
    rise_set <- function(t, Lat, Lon){
      #t = Sun$time[1]; Lat = Sun$lat[1]; Lon = Sun$lon[1]
      ts <- getSunlightTimes(date = seq.Date(date(t) - 1, date(t) + 1, by = 1), Lat, Lon, 
                             keep = c("sunrise", "sunset"))  # yesterday; today; tomorrow
      
      # https://stackoverflow.com/questions/50495809/r-if-else-in-function
      sunrise <- sapply(t, function(el) if(el > ts$sunrise[1] & el <= ts$sunset[1]){ 
        as.character(ts$sunrise[1])
        } else if (el > ts$sunset[1] & el <= ts$sunset[2]) {
          as.character(ts$sunrise[2])
          } else if (el > ts$sunset[2] & el <= ts$sunset[3]) {
            as.character(ts$sunrise[3])
            } else { NA})
      
      sunset <- sapply(t, function(el) if(el > ts$sunrise[1] & el <= ts$sunrise[2]){ 
        as.character(ts$sunset[1])
        } else if (el > ts$sunrise[2] & el <= ts$sunrise[3]) {
          as.character(ts$sunset[2])
          } else if (el > ts$sunrise[3] & el <= ts$sunset[3]) {
            as.character(ts$sunset[3])
            } else { NA})
      data.frame(sunrise, sunset)
    }
    
    #' Time-Based Rolling Mesor (computed using cosinor method)
    #' 
    #' rolling window: at least 51 points (Ignore missing values) (SCS ~1.5h per data)
    #' 
    #' @author Changlin Li 
    #' 
    #' @param t time
    #' @param y value  
    #'
    #' @return An object of the same class as y with the rolling mesor
    #' 
    norm_mesor = function(t, Num, y) {
      #t = df$t; y = df$value; Num = df$Num
      mesor = numeric(length(y))
      
      for( i in seq_along(y) ) {
        for(n in 8:24){
          #i = 27; n = 8
          window = Num >= (Num[i] - n) & Num <= (Num[i] + n); y[window]
          ifelse(length(na.omit(y[window])) < 51, next, break)
        }
        lm1 <- lm(y[window] ~ cos(pi/12*t[window]) + sin(pi/12*t[window])) # lm fit
        mesor[i] <- coef(lm1)[1] # summary(lm1)
      }
      return(mesor)
    }
    
    
    
    
    Pico <- read_csv("Data/Picoauto.csv") %>% select(-Station)
    # The sampling area is located in time zone 8
    Sampletime <- read_csv("Data/Sampletime.csv") %>% 
      mutate(time = ymd_hms(paste(Date0, Time0)),
             Time = time - 8 * 3600) %>% 
      select(-Date0, -Time0, -Sample)
    
    
    # Get Sunlight times
    Sun <- tibble(Station = c("SEATS", "M4", "K11"),
                  lat     = c(18,  19.75,  21.83),
                  lon     = c(116, 115.12, 119.1)) %>% 
      left_join(Sampletime, by = "Station") %>% 
      mutate(Sun  = pmap(list(t = Time, Lat = lat, Lon = lon), rise_set)) %>% 
      unnest(Sun) %>% 
      select(-time)

    
    # normalization
    Sun <- Sun %>%
      mutate(sunrise= ymd_hms(sunrise),
             sunset = ymd_hms(sunset)) %>%
      mutate(Stime  = time_length(interval(sunrise, Time),  "second"), 
             Sdur   = abs(time_length(interval(sunrise, sunset),"second")),
             time_n = date(sunrise) + seconds(6 + Stime/Sdur*12) * 3600 + days(1)) %>%
      select(Station, Id, time_n)
    
    
    # left_join two tibble
    Sun %>% count(Id) %>% filter(n > 1)
    Pico %>% anti_join(Sun, by = "Id")
    
    Pico <- Pico %>% 
      left_join(Sun, by = "Id") %>%
      group_by(Station) %>% 
      mutate(t = time_length(interval(date(min(time_n)), time_n), "hour")) %>% 
      ungroup() %>% 
      mutate(cruise  = paste0("SCS2019_", Station),
             Num     = str_sub(Id, start = 2),
             Num     = as.numeric(Num),
             Dia_Syn = FSC_Dia(SynFSC),
             Dia_Pro = FSC_Dia(ProFSC),
             Dia_Euk = FSC_Dia(PeukFSC),
             Bio_Syn = FSC_Abu_Biomass(SynFSC, Syn),
             Bio_Pro = FSC_Abu_Biomass(ProFSC, Pro),
             Bio_Euk = FSC_Abu_Biomass(PeukFSC, Peuk)) %>%
      select(cruise, Id1, time_n, t, Num, Syn:Peuk, Dia_Syn:Bio_Euk) %>% 
      rename(Abu_Syn = Syn, Abu_Pro = Pro, Abu_Euk = Peuk) %>% 
      gather(para, value, Abu_Syn:Bio_Euk) %>% 
      drop_na(value) %>% 
      arrange(cruise, para, time_n)
    
    # Mesor-based normalization  
    Pico <- Pico %>% 
      group_by(cruise, para) %>% 
      mutate(mesor = norm_mesor(t, Num, value)) %>% 
      select(-Num)
    
    write_csv(Pico, "./SeaFlow/Data/dat_SCS.csv")

    