#### The first step of data analysis of SeaFlow data v1.3 ####
#
# SeaFlow data v1: High-resolution abundance, size and biomass of small phytoplankton in the North Pacific.
#
# Ref: Ribalet, F. et al. SeaFlow data v1: High-resolution abundance, size and biomass of 
#      small phytoplankton in the North Pacific. Zenodo, https://doi.org/10.5281/zenodo.2678021 (2019).
#
# Version: SeaFlow_allstats_v1.3_2020-08-21
#
# 1. Some values of latitude and longitude are wrong, fix them manually;
#    
# 2. Calculate sunrise, sunset times; 
#     
# 3. Time normalization: adjusted to a “standard day” with sunrise at 06:00 and sunset at 18:00
#    
##

    library(tidyverse); library(lubridate); library(data.table) #need shift function
    setwd("./SeaFlow")
    dat0 <- read_csv("Data/SeaFlow_data.csv")
    
    # check parameters of dataset (time \ lat \ lon); summary(dat1$depth) #all depths equal 5m
    dat1 <- dat0 %>% 
      arrange(cruise,time) %>% 
      select(-(abundance_prochloro:biomass_croco),-depth) %>% 
      group_by(cruise) %>% 
      mutate(diff_time  = as.numeric(difftime(time,shift(time,fill = NA), units = "mins")),
             diff_lat   = abs(lat - shift(lat,fill = NA)),
             diff_lon   = abs(lon - shift(lon,fill = NA)),
             ship_v_lat = diff_lat / diff_time * 3600,
             ship_v_lon = diff_lon / diff_time * 3600)
    
    diff_time0  <- dat1 %>% split(.$cruise) %>% map(~summary(.$diff_time))
    diff_lat0   <- dat1 %>% split(.$cruise) %>% map(~summary(.$diff_lat))
    diff_lon0   <- dat1 %>% split(.$cruise) %>% map(~summary(.$diff_lon))
    ship_v_lat0 <- dat1 %>% split(.$cruise) %>% map(~summary(.$ship_v_lat))
    ship_v_lon0 <- dat1 %>% split(.$cruise) %>% map(~summary(.$ship_v_lon))
    i = 1;
    dat1_check <- dat1 %>% filter(cruise == unique(dat1$cruise)[i]); diff_time0[i];diff_lat0[i];diff_lon0[i];ship_v_lat0[i];ship_v_lon0[i]; i=i+1
    

    # Some values of latitude and longitude are wrong, fix them manually; 
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-22 11:28:23")] <- -.0009
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-03-30 02:45:33")] <- -37.948
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-25 05:25:18")] <- 2.831
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-03-28 10:19:13")] <- -37.654
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-23 10:43:08")] <- 0.651
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-03-31 15:07:42")] <- -34.549
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-02 23:04:13")] <- -28.233
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-06 05:03:44")] <- -22.4956
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-06 15:53:33")] <- -22.499
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-07 21:56:37")] <- -22.476
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-10 16:38:35")] <- -15.38
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-22 12:32:09")] <- -.0012
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-22 12:35:09")] <- -.0014
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-22 12:38:09")] <- -.0005
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-22 12:41:09")] <- -.0000
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-22 12:44:09")] <- -.0002
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-22 12:47:10")] <- -.0000
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-22 12:50:10")] <- -.0002
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-22 12:53:10")] <- -.0013
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-27 07:48:56")] <- 5.9205
    dat0$lat [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-06 12:25:32")] <- -22.4983
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-03-31 16:41:28")] <- -42.435
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-08 20:37:02")] <- -32.463
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-08 23:07:14")] <- -32.1
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-09 03:57:38")] <- -31.318
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-15 12:59:04")] <- -25.996
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-18 04:36:42")] <- -28.5001
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-21 10:56:49")] <- -30.929
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-22 23:31:24")] <- -33.297
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-25 19:21:46")] <- -39.0001
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-25 22:07:39")] <- -39.0001
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-26 11:21:30")] <- -40.37
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-26 19:55:24")] <- -41.249
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-27 18:29:31")] <- -41.335
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-27 23:54:38")] <- -41.777
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-04-28 08:23:28")] <- -43.474
    dat0$lon [dat0$cruise == "KN210-04" & dat0$time == ymd_hms("2013-05-02 14:39:34")] <- -52.325
    
    dat0$lat [dat0$cruise == "KOK1606" & dat0$time == ymd_hms("2016-05-02 07:00:59")] <- 27.5108
    dat0$lat [dat0$cruise == "KOK1606" & dat0$time == ymd_hms("2016-04-22 07:00:58")] <- 26.9938
    dat0$lat [dat0$cruise == "KOK1606" & dat0$time == ymd_hms("2016-04-28 07:00:29")] <- 36.5481
    dat0$lon [dat0$cruise == "KOK1606" & dat0$time == ymd_hms("2016-04-28 07:00:29")] <- -158
    
    dat0$lon [dat0$cruise == "SR1917" & dat0$time == ymd_hms("2019-11-16 00:45:05")]  <- -179.9962
    
    dat0$lat [dat0$cruise == "Tokyo_3" & dat0$time == ymd_hms("2011-09-26 22:17:56")] <- 43.7498
    dat0$lon [dat0$cruise == "Tokyo_3" & dat0$time == ymd_hms("2011-09-26 22:17:56")] <- 177.8073
    dat0$lon [dat0$cruise == "Tokyo_3" & dat0$time == ymd_hms("2011-09-26 22:20:56")] <- 177.8303
    dat0$lat [dat0$cruise == "Tokyo_3" & dat0$time == ymd_hms("2011-09-26 10:35:05")] <- 42.9722
    dat0$lon [dat0$cruise == "Tokyo_3" & dat0$time == ymd_hms("2011-09-26 10:35:05")] <- 172.5734
    dat0$lat [dat0$cruise == "Tokyo_3" & dat0$time == ymd_hms("2011-09-26 12:02:11")] <- 43.085
    dat0$lat [dat0$cruise == "Tokyo_3" & dat0$time == ymd_hms("2011-09-26 12:05:11")] <- 43.085
    dat0$lon [dat0$cruise == "Tokyo_3" & dat0$time == ymd_hms("2011-09-26 12:02:11")] <- 173.229
    dat0$lon [dat0$cruise == "Tokyo_3" & dat0$time == ymd_hms("2011-09-26 12:05:11")] <- 173.229
    dat0$lon [dat0$cruise == "Tokyo_3" & dat0$time == ymd_hms("2011-09-26 12:26:13")] <- 173.362
    
    
    # Some latitude and longitude values of two cruises (KOK1512 & TN248) are NA, calculated using adjacent values.
    # The two cruises will be deleted in the subsequent data-cleaning step.
    # dd <- dat0 %>% group_by(cruise) %>% summarise(n1=sum(is.na(lat)),n2=sum(is.na(lon)),n3=sum(is.na(time)))
    dat0$lat [dat0$cruise == "KOK1512" & dat0$time < ymd_hms("2015-09-27 03:12:00")]  <- 21.3162
    dat0$lon [dat0$cruise == "KOK1512" & dat0$time < ymd_hms("2015-09-27 03:12:00")]  <- -157.8862
    dat0$lat [dat0$cruise == "KOK1512" & dat0$time == ymd_hms("2015-09-27 07:00:59")] <- 21.3436
    dat0$lon [dat0$cruise == "KOK1512" & dat0$time == ymd_hms("2015-09-27 07:00:59")] <- -158.2735
    dat0$lat [dat0$cruise == "TN248"   & dat0$time < ymd_hms("2010-05-04 23:13:00")]  <- 49.9424
    dat0$lon [dat0$cruise == "TN248"   & dat0$time < ymd_hms("2010-05-04 23:13:00")]  <- -136.3091
    
    
    #' Calculate the nearest Sunrise and Sunset Times for each time point  
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
    
    # Get Sunlight times, 20 min
    Sun <- dat0 %>%
      select(time,lat,lon,cruise) %>% 
      #slice_sample(n = 100) %>% 
      mutate(Sun  = pmap(list(t = time, Lat = lat, Lon = lon), rise_set)) %>% 
      unnest(Sun)
    
    # normalized time
    Sun1 <- Sun %>%
      mutate(sunrise= ymd_hms(sunrise),
             sunset = ymd_hms(sunset),
             id     = if_else(sunrise < sunset, 1, -1)) %>% # 1: day; -1: night
      mutate(Stime  = time_length(interval(sunrise, time),   "second"), 
             Sdur   = time_length(interval(sunrise, sunset), "second"),
             time_n = date(sunrise) + seconds(6 + Stime/Sdur*12* id) * 3600) %>% 
      arrange(cruise, time)  #%>%
    #   group_by(cruise) %>% 
    #   mutate(diff_time   = as.numeric(difftime(time,   shift(time,   fill = NA), units = "mins")),
    #          diff_time_n = as.numeric(difftime(time_n, shift(time_n, fill = NA), units = "mins")),
    #          diff        = abs(diff_time/diff_time_n - 1))
    # boxplot(Sun1_n$diff)
    
    
    # # check
    # check <- Sun1 %>%
    #   select(cruise, sunrise, sunset, id, time,time_n,Sdur) %>%
    #   mutate(Sdur = round(Sdur / 3600, 2),
    #          lag  = time_length(interval(time, time_n),"hour"),
    #          lag  = round(lag, 2)) %>%
    #   arrange(cruise, time)
    # Icruise <- unique(check$cruise)
    # library(patchwork); library(scales)
    # for (i in 1:length(Icruise)){
    #   #i=25
    #   df   <- check %>% filter(cruise == Icruise[i])
    #   lim1 <- date(min(df$time_n)) + hours(6)
    #   lim2 <- date(max(df$time_n)) + hours(18)
    #   t618 <- seq.POSIXt(lim1, lim2, by = "12 hours")
    # 
    #   p1 <- ggplot(df, aes(time, id *.5)) +
    #     geom_point(position = position_jitter(w = 0, h = 0.35), size = .5) +
    #     scale_x_datetime(date_breaks = "1 day",labels = date_format("%H")) +
    #     xlab(NULL) + ylab("Original") +
    #     theme(axis.title = element_text(size = rel(1.5), face = "bold"))
    #   p2 <- ggplot(df, aes(time_n, id *.5)) +
    #     geom_point(position = position_jitter(w = 0, h = 0.35), color = "blue", size = .5) +
    #     scale_x_datetime(date_breaks = "1 day",labels = date_format("%H")) +
    #     annotate("segment", x = t618, xend = t618, y = -Inf, yend = Inf, size = .1) +
    #     xlab(NULL) + ylab("Normalized") +
    #     theme(axis.title = element_text(size = rel(1.5), face = "bold"));p2
    #   p1/p2
    #   ggsave(filename = paste(Icruise[i], ".pdf", sep = ""), path = "./Fig/s1_nLST", width = 18, height = 8)
    # }
    
    # hours of daytime
    Day_length <- Sun1 %>% 
      filter(id == 1) %>% 
      group_by(cruise) %>% 
      summarise(Sdur = mean(Sdur) / 3600)
    summary(Day_length$Sdur); 
    
    
    Sun1 <- Sun1 %>%
      select(cruise, lat, lon, time, time_n)
    write_csv(Sun1, "Data/Solar.csv")
    
    