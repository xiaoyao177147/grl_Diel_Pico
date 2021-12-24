
    library(tidyverse); library(basicTrendline);
    library(investr) #get confidence interval
    dat  <- read_csv("Data/FSC_size.csv")

    
    y <- dat$Size;  x <- dat$FSC
    fIT2P <- nls(Size ~ SSpower2P(FSC, a, b), data = dat)
    #trendline(x = dat$FSC, y = dat$Size,  model="power2P", ePos.x = "topleft", summary=TRUE, eDigit=5)
    Pred  <- data.frame(FSC = seq(min(dat$FSC),max(dat$FSC),length.out = 300)) 
    Perd2 <- predFit(fIT2P, newdata=Pred, se.fit=T, interval="confidence")$fit  %>% cbind(Pred)
    
    trendline(x = log10(dat$FSC), y = log10(dat$Size), model="line2P", ePos.x = "topleft", summary=TRUE, eDigit=5)
    
    lm1 <- lm(log10(dat$Size)~log10(dat$FSC))
    options(digits=9) 
    (Pare <- summary(lm1))
    (a <- Pare$coefficients[1,1]); (b <- Pare$coefficients[2,1])
    
    
    lwd_pt <- 1/(.pt*72.27/96)
    options(scipen=9)
    theme_set(theme_bw(base_size = 10, base_line_size = .4*lwd_pt, base_rect_size = .4*lwd_pt) +
                theme(plot.title = element_text(face = "bold", hjust = 0.5),
                      panel.border = element_rect(fill = NA, colour = "black", size = rel(2)),
                      axis.text  = element_text(colour = "black"),
                      axis.ticks = element_line(colour = "black"),
                      plot.background  = element_rect(colour = NA),
                      panel.grid = element_blank()))
    
    ggplot() +
      geom_ribbon(data = Perd2, aes(FSC,fit,  ymin=lwr, ymax=upr), alpha=0.6, fill="grey") +
      geom_line(  data = Perd2, aes(FSC,fit), size=2*lwd_pt) +
      geom_point( data = dat, aes(FSC, Size), size=4*lwd_pt, shape = 21, fill = "#b4d893") +
      ylab("ESD (µm)") + xlab("FSC") +
      annotation_logticks(size = .4*lwd_pt)+
      scale_x_log10(limits = c(10^-1.4,10^0.23), breaks = c(0.05,0.1,0.5,1))+
      scale_y_log10(limits = c(10^-.935,10^.38), breaks = c(0.1,0.5,1,2)) +
      annotate("text", x=0.45,  y=.35,parse=TRUE, label="log~italic(ESD)== .2504 ~ log ~italic(FSC) + .1351",size=3) +
      annotate("text", x=0.26,y=0.25,size=3,
               parse=TRUE, label=as.character(expression(italic(R)^2 == 0.93*","~italic(p)<~.0001)))  #pdf 3.42 * 3.42

    
    