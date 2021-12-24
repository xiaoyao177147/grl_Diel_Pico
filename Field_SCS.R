       
     library(plyr); library(tidyverse); library(ggpubr);
     library(grid);  library(lubridate); library(scales);library(cowplot); library(ggh4x)

# Functions ---------------------------------------------------------------
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
    
    #' equation of regression 
    lm_eqn <- function(lm1){
      formula<-sprintf("italic(y) == %.4f %+.4f * italic(x)",coef(lm1)[1], coef(lm1)[2]) 
      eq <- substitute(~italic(R)^2~"="~r2~","~italic(p)~"="~P,
                       list( r2 = format(summary(lm1)$r.squared, digits = 1, nsmall = 4),
                             P = format(summary(lm1)$coefficients[2,4], digits = 1, nsmall = 4)))
      data.frame(formula = formula, eq = as.character(as.expression(eq)), stringsAsFactors = FALSE)
    }
    
    
    #' Dilution plots of net growth rate (9h−1) versus the seawater fraction in parallel experiments were performed.
    #' 
    #' @author Changlin Li 
    #' 
    plot1 <- function(df){
      G <- lm(Rate ~ Dilu, data = df[df$Treat == 'um0.1', ])
      V <- lm(Rate ~ Dilu, data = df[df$Treat == 'kDa30', ])
      lm1 <- lm(Rate ~ Dilu * Treat, data = df) #anova(lm1)
      # library(HH); df <- df %>% mutate(Treat = factor(Treat)); ancova(Rate ~ Dilu * Treat, data = df)
      p <- ggplot(df, aes(Dilu, Rate, col = Treat)) +
        geom_point(size = 2) + stat_smooth(method = lm, se = FALSE) +
        labs(title = paste(df$ExpID[1], df$Pico[1],df$Incubation[1], sep = "."), face = "bold") + xlab(NULL) + ylab(NULL) +
        theme_gray(base_size = 18) +
        geom_text(aes(x = 0.5, y = max(Rate, na.rm = T) + 0.05, label = lm_eqn(G)$formula),
                  parse = TRUE,  size = 6, col = "#00BFC4") +
        geom_text(aes(x = 0.6, y = max(Rate, na.rm = T), label = lm_eqn(G)$eq),
                  parse = TRUE,  size = 6, col = "#00BFC4") +
        geom_text(aes(x = 0.5, y = min(Rate, na.rm = T) - 0.05, label = lm_eqn(V)$formula),
                  parse = TRUE,  size = 6, col = "#F8766D") +
        geom_text(aes(x = 0.6, y = min(Rate, na.rm = T) - 0.1,  label = lm_eqn(V)$eq),
                  parse = TRUE,  size = 6, col = "#F8766D") +
        geom_text(aes(x = 0.76, y = min(Rate[Id == "A"], na.rm = TRUE) - 0.05,
                      label = sprintf("v == %.4f", (-coef(V)[2] + coef(G)[2]))), parse = TRUE, size = 6, col = 'black') +
        geom_text(aes(x = 0.9, y = min(Rate[Id == "A"], na.rm = TRUE) - 0.05,
                      label = sprintf("italic(p) == %.4f", anova(lm1)[2, 5])), parse = TRUE, size = 6, col = 'black') 
      
      if (df$Pico[1] == 'Syn'){
        p <- p +  scale_x_continuous(labels = NULL) +
          scale_colour_hue(breaks = c("kDa30", "um0.1"), labels=c("30 kDa ", "0.1 um")) +
          theme(legend.position  = c(0.7, 1.04),
                legend.direction = "horizontal",
                legend.title     = element_blank())
      } else {
          p <- p +  scale_x_continuous(labels = NULL) + guides(color = "none")}
      return(p)
      }
    
    #' The rates (d-1) of gross growth (µgross), total mortality (TM), grazing (G) and viral lysis (V) 
    #' of different phytoplankton groups inferred from modified dilution experiments.
    #'
    output <- function(df) {
      G <- lm(Rate ~ Dilu, data = df[df$Treat == 'um0.1', ])
      V <- lm(Rate ~ Dilu, data = df[df$Treat == 'kDa30', ])
      lm1 <- lm(Rate ~ Dilu * Treat, data = df)  #anova(lm1)
      op <- data.frame(u1 =  coef(G)[1],     p.u1 = summary(G)$coefficients[1, 4],
                       m1 = -coef(G)[2],     p.m1 = summary(G)$coefficients[2, 4],
                       u2 =  coef(V)[1],     p.u2 = summary(V)$coefficients[1, 4],
                       m2 = -coef(V)[2],     p.m2 = summary(V)$coefficients[2, 4],
                       v  = -coef(V)[2] + coef(G)[2],      p.v = anova(lm1)[2, 5])
    }
    
    
    #'Make some modifications to ggplot2::geom_errorbar
    #' @export
    #' @rdname geom_linerange
    geom_uperrorbar <- function(mapping = NULL, data = NULL,
                                stat = "identity", position = "identity",
                                ...,
                                na.rm = FALSE,
                                orientation = NA,
                                show.legend = NA,
                                inherit.aes = TRUE) {
      layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomUperrorbar,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
          na.rm = na.rm,
          orientation = orientation,
          ...
        )
      )
    }
    
    #' @rdname ggplot2-ggproto
    #' @format NULL
    #' @usage NULL
    #' @export
    GeomUperrorbar <- ggproto("GeomUperrorbar", Geom,
                              default_aes = aes(colour = "black", size = 0.5, linetype = 1, width = 0.5,
                                                alpha = NA),
                              
                              draw_key = draw_key_path,
                              
                              required_aes = c("x|y", "ymin|xmin", "ymax|xmax"),
                              
                              setup_params = function(data, params) {
                                GeomLinerange$setup_params(data, params)
                              },
                              
                              extra_params = c("na.rm", "orientation"),
                              
                              setup_data = function(data, params) {
                                data$flipped_aes <- params$flipped_aes
                                data <- flip_data(data, params$flipped_aes)
                                data$width <- data$width %||%
                                  params$width %||% (resolution(data$x, FALSE) * 0.9)
                                data <- transform(data,
                                                  xmin = x - width / 2, xmax = x + width / 2, width = NULL
                                )
                                flip_data(data, params$flipped_aes)
                              },
                              
                              draw_panel = function(data, panel_params, coord, width = NULL, flipped_aes = FALSE) {
                                data <- flip_data(data, flipped_aes)
                                #x <- as.vector(rbind(data$xmin, data$xmax, NA, data$x,    data$x,    NA, data$xmin, data$xmax))
                                #y <- as.vector(rbind(data$ymax, data$ymax, NA, data$ymax, data$ymin, NA, data$ymin, data$ymin))
                                sel <- data$y < 0 
                                data$ymax[sel] <- data$ymin[sel]
                                x <- as.vector(rbind(data$xmin, data$xmax, NA, data$x,    data$x))
                                y <- as.vector(rbind(data$ymax, data$ymax, NA, data$ymax, data$y))
                                data <- new_data_frame(list(
                                  x = x,
                                  y = y,
                                  colour = rep(data$colour, each = 5),
                                  alpha = rep(data$alpha, each = 5),
                                  size = rep(data$size, each = 5),
                                  linetype = rep(data$linetype, each = 5),
                                  group = rep(1:(nrow(data)), each = 5),
                                  row.names = 1:(nrow(data) * 5)
                                ))
                                data <- flip_data(data, flipped_aes)
                                GeomPath$draw_panel(data, panel_params, coord)
                              }
    )
    # needed for GeomUperrorbar
    new_data_frame <- function(x = list(), n = NULL) {
      if (length(x) != 0 && is.null(names(x))) {
        abort("Elements must be named")
      }
      lengths <- vapply(x, length, integer(1))
      if (is.null(n)) {
        n <- if (length(x) == 0 || min(lengths) == 0) 0 else max(lengths)
      }
      for (i in seq_along(x)) {
        if (lengths[i] == n) next
        if (lengths[i] != 1) {
          abort("Elements must equal the number of rows or 1")
        }
        x[[i]] <- rep(x[[i]], n)
      }
      
      class(x) <- "data.frame"
      
      attr(x, "row.names") <- .set_row_names(n)
      x
    }
    
    
    #' Build a function to find better breaks--adjust
    #' 
    #' @author Changlin Li 
    #' 
    #' @param x the range
    #'
    #' @return  b adjusted breaks
    #' 
    Adj_break <- function(x){
      #x <- c(0, 12.5)
      # scale_continuous expand defaults are 0.05
      #x <- c(x[1]-(x[2]-x[1])*.05, x[2]+(x[2]-x[1])*.05)
      d <- tibble(i = 3:9) %>% 
        mutate(b  = map(i, function(i) get_breaks(i)(x)),
               x1 = x[1],
               x2 = x[2]) %>% 
        unnest_longer(b) %>% 
        filter(b > x1 & b < x2) %>% 
        group_by(i) %>% 
        mutate(n    = length(b),
               Rank = rank(b)) %>% 
        filter(n > 1 & n < 8) %>% 
        filter(!(n %in% c(5, 7) & Rank %in% c(2, 4, 6))) %>%
        mutate(n = ifelse(n == 6 & Rank %in% c(2, 4, 6), 9, n)) %>% 
        group_by(i, n) %>% 
        mutate(Num = length(b),
               sd1 = sd(c(min(b) - x1, x2 - max(b)))) %>% 
        ungroup() %>%
        select(-Rank, -x1, -x2) %>% 
        group_by(Num) %>%
        filter(sd1 == min(sd1)) %>% 
        mutate(n.digits = count_decimals(b), 
               n.digits = max(n.digits)) %>% 
        ungroup() %>% 
        nest(b = b) %>%
        distinct(b, .keep_all = TRUE) %>%
        filter(n.digits == min(n.digits)) %>%
        filter(Num == max(Num))
      d$b[[1]]$b
    }
    
    #'calculate the number of digits after zero
    #'@source https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r
    count_decimals = function(x) {
      #length zero input
      if (length(x) == 0) return(numeric())
      
      #count decimals
      x_nchr = x %>% abs() %>% as.character() %>% nchar() %>% as.numeric()
      x_int = floor(x) %>% abs() %>% nchar()
      x_nchr = x_nchr - 1 - x_int
      x_nchr[x_nchr < 0] = 0
      x_nchr
    }
    
    # A function designed to set the legend spacing increase 
    #(https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2/26971729)
    draw_key_polygon3 <- function(data, params, size) {
      lwd <- min(data$size, min(size) / 4)
      
      grid::rectGrob(
        width = grid::unit(0.6, "npc"),
        height = grid::unit(0.6, "npc"),
        gp = grid::gpar(
          col = data$colour,
          fill = alpha(data$fill, data$alpha),
          lty = data$linetype,
          lwd = lwd * .pt,
          linejoin = "mitre"
        ))
    }
    # register new key drawing function, 
    # the effect is global & persistent throughout the R session
    GeomBar$draw_key = draw_key_polygon3

# 1. Dilution --------------------------------------------------------------
    
    Dilu0 <- read_csv("Data/Dilution.csv") %>% 
    # As large variation and NA values between parallel FSC data, the method of averaging is adopted  
      group_by(ExpID, Id) %>% 
      mutate_at(c("SynFSC", "ProFSC", "PeukFSC"), mean, na.rm = TRUE) %>% 
      ungroup() %>% 
      mutate(SynBiomass = FSC_Abu_Biomass(SynFSC, Syn),
             ProBiomass = FSC_Abu_Biomass(ProFSC, Pro),
             PeukBiomass= FSC_Abu_Biomass(PeukFSC, Peuk)) %>% 
      select(-SynFSC,-ProFSC,-PeukFSC, -Id1)
    
    T0 <- Dilu0 %>% filter(Treat == "T0") %>%
      rename(T0Syn = Syn, T0Pro = Pro, T0Peuk = Peuk,
             T0SynBio = SynBiomass, T0ProBio = ProBiomass, T0PeukBio = PeukBiomass) %>% 
      select(ExpID, contains("T0"))
    
    Dilu <- Dilu0 %>% 
      filter(Treat != "T0") %>%
      left_join(T0, by = "ExpID") %>%
      mutate(RateSyn_abu  = log(Syn  / Dilu / T0Syn),
             RatePro_abu  = log(Pro  / Dilu / T0Pro),
             RatePeuk_abu = log(Peuk / Dilu / T0Peuk),
             RateSyn_bio  = log(SynBiomass  / Dilu / T0SynBio),
             RatePro_bio  = log(ProBiomass  / Dilu / T0ProBio),
             RatePeuk_bio = log(PeukBiomass / Dilu / T0PeukBio)) %>% 
      select(-(Syn:T0PeukBio),-Station)
    
    #Add 4 columns of Id == "A" dilution data--30kda
    Dilu_A <- filter(Dilu, Id == "A") %>% mutate(Treat = "kDa30")
    Dilu <- Dilu %>% 
      mutate(Treat = if_else(Id == "A", "um0.1", Treat)) %>% 
      bind_rows(Dilu_A) %>% 
      gather(Pico, Rate, RateSyn_abu:RatePeuk_bio) %>% 
      separate(Pico, c("Pico", "Para")) %>% 
      mutate(Pico = str_sub(Pico, 5))
    
    # df <- Dilu %>%
    #   filter(Incubation == "Day-1" & ExpID == "S1", Pico == "Pro", Para == "bio")
    # p1 <- dlply(Dilu, .(Incubation, ExpID, Pico, Para), plot1)
    # 
    # png(file = paste("Figure/Abu_seats-1.png"), width = 1600, height = 950)
    # ggarrange(p1$`Day-1.S1.Peuk.abu`,p1$`Day-1.S1.Syn.abu`, p1$`Day-1.S1.Pro.abu`,
    #           p1$`Night-1.S2.Peuk.abu`,p1$`Night-1.S2.Syn.abu`, p1$`Night-1.S2.Pro.abu`,
    #           ncol = 3, nrow = 2, widths = c(1,1), heights = c(1,1))
    # dev.off()
    # png(file = paste("Figure/Abu_seats-2.png"), width = 1600, height = 950)
    # ggarrange(p1$`Day-2.S3.Peuk.abu`,p1$`Day-2.S3.Syn.abu`, p1$`Day-2.S3.Pro.abu`,
    #           p1$`Night-2.S4.Peuk.abu`,p1$`Night-2.S4.Syn.abu`, p1$`Night-2.S4.Pro.abu`,
    #           ncol = 3, nrow = 2, widths = c(1,1), heights = c(1,1))
    # dev.off()
    # png(file = paste("Figure/Abu_m4-1.png"), width = 1600, height = 950)
    # ggarrange(p1$`Night-1.M1.Peuk.abu`,p1$`Night-1.M1.Syn.abu`, p1$`Night-1.M1.Pro.abu`,
    #           p1$`Day-1.M2.Peuk.abu`,p1$`Day-1.M2.Syn.abu`, p1$`Day-1.M2.Pro.abu`,
    #           ncol = 3, nrow = 2, widths = c(1,1), heights = c(1,1))
    # dev.off()
    # png(file = paste("Figure/Abu_m4-2.png"), width = 1600, height = 950)
    # ggarrange(p1$`Night-2.M3.Peuk.abu`,p1$`Night-2.M3.Syn.abu`, p1$`Night-2.M3.Pro.abu`,
    #           p1$`Day-2.M4.Peuk.abu`,p1$`Day-2.M4.Syn.abu`, p1$`Day-2.M4.Pro.abu`,
    #           ncol = 3, nrow = 2, widths = c(1,1), heights = c(1,1))
    # dev.off()
    # 
    # png(file = paste("Figure/Abu_k11.png"), width = 1600, height = 950)
    # ggarrange(p1$`Night-1.K1.Peuk.abu`,p1$`Night-1.K1.Syn.abu`, p1$`Night-1.K1.Pro.abu`,
    #           p1$`Day-1.K2.Peuk.abu`,p1$`Day-1.K2.Syn.abu`, p1$`Day-1.K2.Pro.abu`,
    #           ncol = 3, nrow = 2, widths = c(1,1), heights = c(1,1))
    # dev.off()
    # 
    # png(file = paste("Figure/Bio_seats-1.png"), width = 1600, height = 950)
    # ggarrange(p1$`Day-1.S1.Peuk.bio`,p1$`Day-1.S1.Syn.bio`, p1$`Day-1.S1.Pro.bio`,
    #           p1$`Night-1.S2.Peuk.bio`,p1$`Night-1.S2.Syn.bio`, p1$`Night-1.S2.Pro.bio`,
    #           ncol = 3, nrow = 2, widths = c(1,1), heights = c(1,1))
    # dev.off()
    # png(file = paste("Figure/Bio_seats-2.png"), width = 1600, height = 950)
    # ggarrange(p1$`Day-2.S3.Peuk.bio`,p1$`Day-2.S3.Syn.bio`, p1$`Day-2.S3.Pro.bio`,
    #           p1$`Night-2.S4.Peuk.bio`,p1$`Night-2.S4.Syn.bio`, p1$`Night-2.S4.Pro.bio`,
    #           ncol = 3, nrow = 2, widths = c(1,1), heights = c(1,1))
    # dev.off()
    # png(file = paste("Figure/Bio_m4-1.png"), width = 1600, height = 950)
    # ggarrange(p1$`Night-1.M1.Peuk.bio`,p1$`Night-1.M1.Syn.bio`, p1$`Night-1.M1.Pro.bio`,
    #           p1$`Day-1.M2.Peuk.bio`,p1$`Day-1.M2.Syn.bio`, p1$`Day-1.M2.Pro.bio`,
    #           ncol = 3, nrow = 2, widths = c(1,1), heights = c(1,1))
    # dev.off()
    # png(file = paste("Figure/Bio_m4-2.png"), width = 1600, height = 950)
    # ggarrange(p1$`Night-2.M3.Peuk.bio`,p1$`Night-2.M3.Syn.bio`, p1$`Night-2.M3.Pro.bio`,
    #           p1$`Day-2.M4.Peuk.bio`,p1$`Day-2.M4.Syn.bio`, p1$`Day-2.M4.Pro.bio`,
    #           ncol = 3, nrow = 2, widths = c(1,1), heights = c(1,1))
    # dev.off()
    # png(file = paste("Figure/Bio_k11.png"), width = 1600, height = 950)
    # ggarrange(p1$`Night-1.K1.Peuk.bio`,p1$`Night-1.K1.Syn.bio`, p1$`Night-1.K1.Pro.bio`,
    #           p1$`Day-1.K2.Peuk.bio`,p1$`Day-1.K2.Syn.bio`, p1$`Day-1.K2.Pro.bio`,
    #           ncol = 3, nrow = 2, widths = c(1,1), heights = c(1,1))
    # dev.off()

    Result <- ddply(Dilu, .(Incubation, ExpID, Pico, Para), output) %>% as_tibble()
    

## Table_S3 ----------------------------------------------------------------
    
    table_s3 <- Result %>% 
      select(Incubation, ExpID, Pico, Para, u2, p.u2, m1, p.m1, v, p.v) %>% 
      rename(Period = Incubation) %>% 
      filter(Pico != "Peuk") %>% 
      mutate(u2 = sprintf("%0.3f", u2/9),
             m1 = sprintf("%0.3f", m1/9),
             v  = sprintf("%0.3f", v/9),
             p.u2 = case_when(
               p.u2 < .01               ~ "p < 0.01",
               p.u2 < .05 & p.u2 >= .01 ~ "p < 0.05",
               TRUE                     ~ paste0("p = ", sprintf("%0.2f", p.u2))),
             p.m1 = case_when(
               p.m1 < .01               ~ "p < 0.01",
               p.m1 < .05 & p.m1 >= .01 ~ "p < 0.05",
               TRUE                     ~ paste0("p = ", sprintf("%0.2f", p.m1))),
             p.v = case_when(
               p.v < .01               ~ "p < 0.01",
               p.v < .05 & p.v >= .01  ~ "p < 0.05",
               TRUE                     ~ paste0("p = ", sprintf("%0.2f", p.v)))) 
    
    Table_s3 <- table_s3 %>% 
      mutate(u2 = paste0(u2, " (", p.u2, ")"),
             m1 = paste0(m1, " (", p.m1, ")"),
             v = paste0(v, " (", p.v, ")")) %>% 
      select(-p.u2, -p.m1, -p.v) %>% 
      pivot_wider(names_from = Para, 
                  values_from = c(u2, m1, v)) %>% 
      mutate(Letter = str_sub(ExpID,1,1),
             Station = case_when(
               Letter == "S" ~ "SEATS",
               Letter == "M" ~ "M4",
               Letter == "K" ~ "K11")) %>% 
      #rename(Picophytoplankter = Pico) %>% 
      select(Pico, Station, Period, contains("abu"), contains("bio")) %>% 
      arrange(Pico, Station, Period)
    

## Figure_part2 ------------------------------------------------------------

    fig2 <- Result %>% 
      filter(Pico != "Peuk") %>% 
      select(Incubation:Para, m1, u2, v) %>% 
      gather(Rate, Value, m1:v) %>% 
      separate(Incubation, c("Id","Id1"))  %>% 
      mutate(Rate = factor(Rate, levels = c("u2", "m1","v")))
    
    
    lwd_pt <- 1/(.pt*72.27/96)
    lab.b  <- c(abu = "Abundance",bio = "Biomass")
    lab.b1 <- c(Pro = "Prochlorococcus",Syn = "Synechococcus",Peuk = "Picoeukaryotes")
    Rate_label <- c('Grazer-mediated \nmortality (h-1)', 'Intrinsic growth \nrate (h-1)', 'Virus-mediated \nmortality (h-1)')
    names(Rate_label) <- c("m1", "u2", "v")
    
    p2 <- ggplot(fig2, aes(Pico, Value/9, fill = Id)) +
      geom_hline(aes(yintercept = 0),color = "grey", size = .5*lwd_pt) +
      geom_bar(stat = "summary", fun = mean, width = 0.6, size = .3*lwd_pt,
               position = position_dodge( .7)) +
      stat_summary(fun.data = 'mean_sd', geom = "uperrorbar", colour = "black", width = .2, size = .3*lwd_pt,
                   position = position_dodge( .7)) +
      stat_compare_means(aes(label =..p.signif..), 
                         hide.ns = TRUE, 
                         label.y.npc = .84,
                         symnum.args = list(cutpoints = c(0,0.05,1), symbol = c("*","ns")),
                         size = 4) +
      xlab(NULL) + ylab(NULL)  +
      facet_grid2(vars(Rate), vars(Para), 
                  scales = "free_y", switch = "y", #independent = "y", 
                  axes = "all", remove_labels = "all",
                  labeller=labeller(Rate = Rate_label,
                                    Para = as_labeller(lab.b))) +
      scale_x_discrete(labels = c(Pro = "PRO", Syn = "SYN")) +
      scale_fill_manual(values = c("#b8529e","#cacac9")) +
      theme_classic(base_size = 10, base_line_size = .5*lwd_pt, base_rect_size = .5*lwd_pt) +
      theme(plot.title = element_text(face = "bold", hjust = .5),
            axis.text  = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            axis.title = element_text(size = rel(.9)),
            plot.background  = element_rect(colour = NA),
            panel.grid   = element_blank(),
            legend.text  = element_text(size = rel(.8)),
            legend.title = element_blank(),
            strip.placement.y = "outside",
            strip.background = element_rect(colour = NA),
            strip.text = element_text(size = rel(.9)),
            legend.background = element_rect(fill = "white", color = "black"),
            legend.key.width  = unit(.4, "cm"),
            legend.key.height = unit(.4, "cm")) #; p2
    

# Time-series -------------------------------------------------------------
    
    Pico0 <- read_csv("Data/Picoauto.csv")
    Sampletime <- read_csv("Data/Sampletime.csv") %>% 
      mutate(Time = ymd_hms(paste(Date0, Time0))) %>%
      mutate(Gap = case_when(
        Station == "SEATS" ~ 0,
        Station == "M4"    ~ 4 * 3600 *24, # For the drawing, the time is pushed forward by 4 days
        Station == "K11"   ~ 10 * 3600 *24)) %>%
      mutate(Time = Time - Gap) %>% 
      select(Id, Time)
    
    Sampletime %>% count(Id) %>% filter(n > 1)
    Pico0 %>% anti_join(Sampletime, by = "Id")
    
    Pico <- Pico0 %>% 
      left_join(Sampletime, by = "Id") %>%
      arrange(Time) %>%
      mutate(Syn_Biomass  = FSC_Abu_Biomass(SynFSC, Syn),
             Pro_Biomass  = FSC_Abu_Biomass(ProFSC, Pro),
             Peuk_Biomass = FSC_Abu_Biomass(PeukFSC, Peuk),
             Syn_Dia      = FSC_Dia(SynFSC),
             Pro_Dia      = FSC_Dia(ProFSC),
             Peuk_Dia     = FSC_Dia(PeukFSC)) %>% 
      select(Station,Id, Id1, Time, contains("Syn"), contains("Pro"), contains("Peuk")) %>% 
      rename(Syn_Abu = Syn, Pro_Abu = Pro, Peuk_Abu = Peuk,
             Syn_FSC = SynFSC, Pro_FSC = ProFSC, Peuk_FSC = PeukFSC) %>%
      gather(Pico, Value, Syn_Abu:Peuk_Dia) %>%
      separate(Pico, c("Pico", "Para")) %>% 
      mutate(Pico = factor(Pico, levels = c("Pro","Syn","Peuk")),
             Para = factor(Para, levels = c("Abu","Dia","Biomass")),
             Station = factor(Station, levels = c("SEATS","M4","K11"))) %>% 
      drop_na(Value) %>% 
      filter(Para != "FSC")

## Figure_part1 ------------------------------------------------------------
    
    library(suncalc) # Get Sunlight times # use station SEATS for night-time
    # The sampling area is located in time zone 8
    Sun <- getSunlightTimes(as.Date(c("2019-06-19", "2019-06-20", "2019-06-21", "2019-06-22", "2019-06-23")),
                            keep = c("sunrise", "sunset"), lat = 18, lon = 116)
    Night1 <- Sun$sunset[-5] + 8 * 3600
    Night2 <- Sun$sunrise[-1]+ 8 * 3600
    lab.p   <- c(Abu = "Abundance (103 cells m-1)",Dia = "Cell size (mean ESD, µm)",Biomass = "Biomass (µg C L-1)")
    
    p1 <- ggplot(Pico, aes(Time, Value, color = Station)) + 
      annotate("rect", xmin = Night1, xmax = Night2, ymin = -Inf, ymax = Inf, alpha = .4, fill = "gray") + 
      geom_line(stat = "summary", fun = mean, size = .8*lwd_pt) +
      stat_summary(fun.data = 'mean_sd', geom = "errorbar", width=2400, size = .5*lwd_pt,color="black") +
      geom_point(stat = "summary", fun = mean, size = 1.5*lwd_pt) +
      scale_x_datetime(breaks = seq(ymd_hms("2019-06-20 06:00:00"), ymd_hms("2019-06-22 18:00:00"), 
                                    "6 hours"), labels = date_format("%H")) +
      scale_y_continuous(breaks = Adj_break) +
      ylab(NULL) + xlab("Time of day (local time)") +
      facet_grid2(vars(Pico), vars(Para), 
                  scales = "free_y", independent = "y", switch = "y",
                  axes = "all", remove_labels = "x",
                  labeller=labeller(Para = as_labeller(lab.p),
                                    Pico = as_labeller(lab.b1))) +
      coord_cartesian(xlim = c(ymd_hms("2019-06-20 06:50:00"), ymd_hms("2019-06-22 17:30:00"))) +
      scale_color_manual(values = c("#ee7e32","#1fb050","#6f3996")) +
      theme_bw(base_size = 10, base_line_size = .5*lwd_pt, base_rect_size = .5*lwd_pt) +
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
            legend.title = element_blank(),
            legend.background = element_rect(fill = "white", color = "black"), # a little big
            legend.text = element_text(size = rel(0.8)),
            #legend.key.width  = unit(.6, "cm"),
            legend.key.height = unit(0.27, "cm")) #;p1

## Figure_part3 ------------------------------------------------------------
    
    Sampletime <- read_csv("Data/Sampletime.csv") %>% 
      mutate(Time = ymd_hms(paste(Date0, Time0))) %>%
      select(Id, Time)
    
    Pico1 <- Pico0 %>% 
      left_join(Sampletime, by = "Id") %>%
      arrange(Time) %>%
      mutate(Syn_bio  = FSC_Abu_Biomass(SynFSC, Syn),
             Pro_bio  = FSC_Abu_Biomass(ProFSC, Pro),
             Peuk_bio = FSC_Abu_Biomass(PeukFSC, Peuk)) %>% 
      select(Station,Id, Id1, Time, Syn:Peuk, contains("bio")) %>% 
      rename(Syn_abu = Syn, Pro_abu = Pro, Peuk_abu = Peuk) %>%
      gather(Pico, Value, Syn_abu:Peuk_bio) %>%
      separate(Pico, c("Pico", "Para")) %>% 
      drop_na(Value)  %>% 
      filter(Id %in% c("S1", "S9", "S17", "S25", "S33", "M1", "M9", "M17", "M25", "M33", "K1", "K9", "K17"),
             Pico != "Peuk") %>% 
      group_by(Station, Id, Pico, Para) %>% 
      summarise_at(c("Time", "Value"), mean, na.rm = TRUE)
    
    Pico1 <- Pico1 %>%
      ungroup() %>% 
      mutate(Id1 = case_when(
        Id %in% c("S1")              ~ "dawn1",
        Id %in% c("S17", "M9", "K9") ~ "dawn2",
        Id %in% c("S33","M25")       ~ "dawn3",
        Id %in% c("S9", "M1", "K1")  ~ "dusk1",
        Id %in% c("S25","M17","K17") ~ "dusk2",
        Id %in% c("M33")             ~ "dusk3")) %>% 
      select(-Id) %>% 
      pivot_wider(names_from  = Id1,
                  values_from = c(Value,Time)) %>% 
      mutate(
        r1 = log(Value_dusk1 / Value_dawn1) / as.numeric(difftime(Time_dusk1, Time_dawn1, units = "hours")),
        r2 = log(Value_dawn2 / Value_dusk1) / as.numeric(difftime(Time_dawn2, Time_dusk1, units = "hours")),
        r3 = log(Value_dusk2 / Value_dawn2) / as.numeric(difftime(Time_dusk2, Time_dawn2, units = "hours")),
        r4 = log(Value_dawn3 / Value_dusk2) / as.numeric(difftime(Time_dawn3, Time_dusk2, units = "hours")),
        r5 = log(Value_dusk3 / Value_dawn3) / as.numeric(difftime(Time_dusk3, Time_dawn3, units = "hours"))) %>% 
      select(Station, Pico, Para, r1:r5) %>% 
      gather(Rate, net_s, r1:r5) %>% 
      drop_na(net_s) %>% 
      mutate(Id1 = case_when(
        Rate == "r2" ~ "Night-1",
        Rate == "r4" ~ "Night-2",
        Rate == "r1" ~ "Day-1",
        Rate == "r5" ~ "Day-2",
        Rate == "r3" & Station == "SEATS" ~ "Day-2",
        Rate == "r3" & Station != "SEATS" ~ "Day-1")) %>% 
      select(-Rate)
    

    # dilution result    
    dilution2 <- Result %>% 
      select(Incubation, ExpID, Pico, Para, m1, u2, v) %>% 
      rename(Id1 = Incubation) %>% 
      filter(Pico != "Peuk") %>% 
      mutate(m1_ = if_else(m1 < 0, 0, m1),
             v_ =  if_else(v  < 0, 0, v)) %>% 
      mutate(net_d  = (u2 - m1 -v) / 9,
             net_d2 = (u2 - m1_ -v_) / 9,
             Letter = str_sub(ExpID,1,1),
             Station = case_when(
               Letter == "S" ~ "SEATS",
               Letter == "M" ~ "M4",
               Letter == "K" ~ "K11"))
    
    Pico1 %>% count(Station, Pico, Para, Id1) %>% filter(n > 1)
    dilution2 %>% anti_join(Pico1, by = c("Station", "Pico", "Para", "Id1"))
    
    fig3 <- dilution2 %>% 
      left_join(Pico1, by = c("Station", "Pico", "Para", "Id1")) %>% 
      select(Station, Pico, Para, Id1, net_d, net_d2, net_s) %>% 
      gather(net, rate, net_d:net_s)  %>%
      mutate(Id1 = str_split_fixed(Id1, pattern = "-", n = 2)[, 1]) # %>%
    # group_by(Id1, Pico, Para, net) %>% 
    # summarise(rate = mean(rate))


    p3 <- ggplot(fig3, aes(Id1, rate, fill = net)) +
      geom_hline(aes(yintercept=0), color = "grey", size = .5*lwd_pt) +
      geom_bar(stat = "summary",fun = mean, width = 0.8, size = .3*lwd_pt,
               position = position_dodge( .9)) +
      stat_summary(fun.data = 'mean_sd', geom = "uperrorbar", colour = "black", width = .3, size = .3*lwd_pt,
                   position = position_dodge( .9)) +
      xlab("Net growth rates (h-1)") + ylab(NULL) +
      facet_grid2(vars(Pico), vars(Para), 
                  scales = "free_y", switch = "y", #independent = "y", 
                  axes = "all", remove_labels = "all",
                  labeller=labeller(Para = as_labeller(lab.b),
                                    Pico = as_labeller(lab.b1))) +
      scale_fill_manual(values = c("#5fc8ee","#fab015","#8ac540"), 
                        breaks = c("net_d", "net_d2", "net_s"),
                        labels = c("Incubation", "Incubation-corrected", "Time series")) +
      theme_classic(base_size = 10, base_line_size = .5*lwd_pt, base_rect_size = .5*lwd_pt) +
      theme(plot.title = element_text(face = "bold", hjust = .5),
            axis.text  = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            axis.title = element_text(size = rel(.9)),
            plot.background  = element_rect(colour = NA),
            panel.grid   = element_blank(),
            legend.text  = element_text(size = rel(.8)),
            legend.title = element_blank(),
            strip.placement.y = "outside",
            strip.background = element_rect(colour = NA),
            strip.text = element_text(size = rel(.9)),
            legend.background = element_rect(fill = "white", color = "black"),
            legend.key.width  = unit(.4, "cm"),
            legend.key.height = unit(.4, "cm")) #;p3

# Layouts in ggplot -------------------------------------------------------
    
    lg  <- plot_grid(get_legend(p2), get_legend(p3))
    
    # remove the some y_axis
    g2 <- ggplotGrob(p2 + theme(legend.position="none"))
    g2[["grobs"]][[20]] <- zeroGrob() # "axis-l-1-2"
    g2[["grobs"]][[21]] <- zeroGrob() # "axis-l-2-2"
    g2[["grobs"]][[22]] <- zeroGrob() # "axis-l-3-2"
    
    # remove the some y_axis
    g3 <- ggplotGrob(p3 + theme(legend.position="none"))
    g3[["grobs"]][[14]] <- zeroGrob() # "axis-l-1-2"
    g3[["grobs"]][[15]] <- zeroGrob() # "axis-l-2-2"

    p30 <- plot_grid(g3, lg, ncol = 1, rel_heights = c(5.6, 2))#; p30
    p23 <- plot_grid(g2, p30, labels = c("b", "c"), align = 'h', nrow = 1, rel_widths = c(7, 6))
    
    plot_grid(p1, p23, labels = c("a", ""), ncol = 1, rel_heights = c(7, 6)) +
      draw_plot_label(x = .65, y = .12, 
                      label =  "PRO: ProchlorococcusSYN: Synechococcus",
                      hjust = 0, size = 4)
    
    nn <- length(list.files("D:/Demo/1/", pattern=".pdf"))
    ggsave(filename = paste0("D:/Demo/1/", "plot_", nn, ".pdf"),  width = 7.4, height = 8.05)
    
    