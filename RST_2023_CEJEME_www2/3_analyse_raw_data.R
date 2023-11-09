if(!require("RColorBrewer")) {install.packages("RColorBrewer"); library(RColorBrewer)}
if(!require("ggsci")) {install.packages("ggsci"); library(ggsci)}
if(!require("viridis")) {install.packages("viridis"); library(viridis)}
if(!require("wesanderson")) {install.packages("wesanderson"); library(wesanderson)}

theme.novpadding <-
  list(axis.line = 
         list(col =  'transparent'),
       layout.heights =
         list(top.padding = 0.2,
              main.key.padding = 0.2,
              key.axis.padding = 1,
              axis.xlab.padding = 0.5,
              xlab.key.padding = 1,
              key.sub.padding = 1,
              bottom.padding = 1),
       layout.widths =
         list(left.padding = 0.5,
              key.ylab.padding = 2,
              ylab.axis.padding = 0,
              axis.key.padding = 2,
              right.padding =0))
par(mfrow=c(5,2))
par()

# PL UE 2011, 2013, 2015, 2017, 2019, 2021
lista_PL_UE <-list()
lista_PL_UE[[1]] <- c("Miasto Kraków", "Miasto Szczecin", "Miasto Wrocław", "Miasto Warszawa", "Miasto Poznań", "Miasto Łódź") #6
lista_PL_UE[[2]] <- c("Kaliski", "Poznański", "Trójmieski", "Gliwicki", "Katowicki", "Leszczyński", "Lubelski", "Tyski", 
                      "Bielski", "Rybnicki", "Warszawski Zachodni") #11
lista_PL_UE[[3]] <- c("Bydgosko-toruński", "Nowotarski", "Opolski", "Siedlecki", "Skierniewicki", "Wrocławski", 
                      "Legnicko-głogowski", "Krakowski", "Żyrardowski") #9
lista_PL_UE[[4]] <- c("Łomżyński", "Puławski", "Sandomiersko-jędrzejowski", "Oświęcimski", "Bytomski", 
                      "Koniński", "Łódzki", "Piotrkowski", "Rzeszowski") #9 
lista_PL_UE[[5]] <- c("Sieradzki", "Tarnobrzeski","Tarnowski", "Warszawski Wschodni", "Częstochowski",
                      "Białostocki", "Pilski", "Gdański", "Zielonogórski", "Sosnowiecki") #10
lista_PL_UE[[6]] <- c("Gorzowski", "Bialski", "Kielecki", "Ostrołęcki", 
                      "Nowosądecki", "Suwalski", "Chełmsko-zamojski", "Ciechanowski", "Przemyski", "Krośnieński") #10
lista_PL_UE[[7]] <- c("Nyski", "Olsztyński", "Słupski", "Chojnicki", "Elbląski", "Jeleniogórski", "Koszaliński", "Starogardzki", 
                      "Grudziącki", "Płocki", "Świecki") #11
lista_PL_UE[[8]] <- c("Szczeciński", "Wałbrzyski","Ełcki", "Radomski", "Włocławski", "Szczeciniecko-pyrzycki", "Inowrocławski") #7

# PL GDP 2000, 2004, 2008, 2012, 2016, 2018
lista_PL_GDP <-list()
lista_PL_GDP[[1]] <- c("Miasto Kraków", "Miasto Szczecin", "Miasto Wrocław", "Miasto Warszawa", "Miasto Poznań", "Miasto Łódź") #6
lista_PL_GDP[[2]] <- c("Chojnicki", "Suwalski", "Bialski", "Ełcki", "Łomżyński", "Nowotarski", "Nyski",
                       "Ostrołęcki", "Puławski", "Siedlecki") #10
lista_PL_GDP[[3]] <- c("Słupski", "Skierniewicki", "Włocławski", "Ciechanowski", "Grudziądzki", "Inowrocławski", "Przemyski", "Świecki", 
                       "Krośnieński", "Żyrardowski", "Gorzowski") # 11
lista_PL_GDP[[4]] <- c("Bytomski", "Elbląski", "Jeleniogórski", "Koszaliński", "Łódzki", "Nowosądecki", "Pilski", "Sieradzki", "Starogardzki", 
                       "Szczeciniecko-pyrzycki", "Szczeciński", "Tarnowski") #12
lista_PL_GDP[[5]] <- c("Radomski", "Chełmsko-zamojski", "Sandomiersko-jędrzejowski","Warszawski Wschodni", 
                       "Białostocki", "Gdański", "Olsztyński", "Opolski", "Oświęcimski", "Zielonogórski") #10
lista_PL_GDP[[6]] <- c("Gliwicki", "Rybnicki", "Koniński", "Leszczyński", "Piotrkowski", "Rzeszowski", "Tarnobrzeski", 
                       "Tyski", "Wałbrzyski", "Częstochowski", "Krakowski", "Płocki") #12
lista_PL_GDP[[7]] <- c("Legnicko-głogowski", "Lubelski", "Bydgosko-toruński", "Kaliski", "Kielecki", "Poznański", "Wrocławski", "Sosnowiecki", 
                       "Warszawski Zachodni", "Bielski") #10
lista_PL_GDP[[8]] <- c("Katowicki", "Trójmiejski") # 2

####### drawing timelines with grouped/selected regions
draw_timelines_selected<-function(dataset, lista, rows, columns, text, variable, 
                                  div=1, width = 8.27, height = 11.69){
  draw <- list()
  f <- "Times"
  for (i in 1:length(lista)){
    data <- dataset
    pal <- brewer.pal(20, "Greys") #wes_palette("Zissou1", 100, type="continuous")
    a <- mean(data$MEAN/div)
    d <- data[which(data$Name %in% lista[[i]]),] %>%
      dplyr::filter(Value!=0)
    
    draw[[i]] <- ggplot(data=d,
                        aes(x=Date,
                            y=Value/div, 
                            group=Name, 
                            colour = Name)) +
      geom_line(size=0.5) +
      theme_ipsum() +
      scale_colour_grey()+#scale_color_viridis(discrete=TRUE, option="D") +
      xlab('date') +
      ylab(paste("value", variable)) +
      scale_y_continuous(name=" ") +
      labs(color = "Subregions:", face="bold") +
      theme(plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"),
            legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
            axis.title.x = element_text(size = 8, family = f),
            axis.title.y = element_text(size = 8, family = f),
            axis.text.x = element_text(size = 8, family = f),
            axis.text.y = element_text(size = 8, family = f),
            legend.text = element_text(size = 6, family = f),
            legend.title = element_text(size = 7, family = f),
            legend.position = "right",
            legend.key.size = unit(0.5, "cm"),
            legend.key.width = unit(0.5,"cm"))
    plot(draw[[i]])
    
  }
  png(file = paste0(text, ".png"), 
      width = width, 
      height = height, 
      units ="in", 
      res = 300)
  p <- marrangeGrob(draw, nrow=rows, ncol=columns, top = NULL)
  print(p)
  dev.off()
}
  
####### drawing maps in matrix
draw_raw_matrixMaps <- function(data, map, var, nclr, div, kolorki, i, key, st, years){
  data$Month <- as.character(month(data$Date))
  data$Year <- as.character(year(data$Date))
  ranges<-data[which(data$Year %in% years),] %>%
    dplyr::select(ID, Name, Date,Month,Year,Value) %>%
    dplyr::filter(Month=="1")
  pal <- brewer.pal(nclr, kolorki) # we select 7 colors from the palette
  breaks_qt <- classIntervals(ranges$Value/div, n=nclr, style=st)
  br <- breaks_qt$brks
  
  d<-data %>%
    dplyr::select(ID, Name, Date,Month,Year,Value) %>%
    dplyr::filter(Month=="1") %>%
    dplyr::filter(Year==i)
  sp <- merge(x = map, y = d, by.x = "ID", by.y = "ID")
  sp@data$bracket <- cut(sp@data$Value/div, breaks_qt$brks)
  spplot(sp, "bracket", lwd=0.1, col.regions=pal,
         colorkey=key, 
         main =list(label=i, cex=0.8, fontfamily=f),
         par.settings = theme.novpadding)
}
