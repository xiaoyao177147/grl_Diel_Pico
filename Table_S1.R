    rm(list = ls())
    library(lubridate);library(tidyverse); 
    
    Auto <- read_csv("Data/Picoauto.csv") #%>% drop_na(Pro) 
    bac_virus <- read_csv("Data/Bac_virus.csv")
    
    bac_virus %>% count(Id, Id1) %>% filter(n > 1)
    Auto %>% anti_join(bac_virus, by = c("Id","Id1"))
    
    dat <- Auto %>% 
      left_join(bac_virus, by = c("Id","Id1")) %>%
      group_by(Station,Id) %>% 
      summarise_at(vars(Syn:Virus), mean, na.rm = TRUE) %>% 
      mutate(hBac = Bac - Pro) %>%
      select(Station:Peuk,hBac,Virus) 
    
    specify_decimal <- function(x, k = 2) trimws(format(round(x, k), nsmall=k)) #set two decimal
   
    Table1 <- dat %>% 
      gather(Pico, Abu, Syn:Virus) %>% 
      group_by(Station,Pico) %>% 
      summarise(Mean = specify_decimal(mean(Abu)),
                SD = specify_decimal(sd(Abu)))
    
    Chla <- read_csv("Data/Chla.csv") %>%  
      group_by(Station) %>% 
      summarise(Mean = specify_decimal(mean(Chla)),
                SD = specify_decimal(sd(Chla)))

    
