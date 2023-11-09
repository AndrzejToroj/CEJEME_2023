############################### LIBRARIES ######################################
if(!require("RColorBrewer")) {install.packages("RColorBrewer"); library(RColorBrewer)}
if(!require("classInt")) {install.packages("classInt"); library(classInt)}
if(!require("grid")) {install.packages("grid"); library(grid)}
if(!require("ggplot2")) {install.packages("ggplot2"); library(ggplot2)}
if(!require("lattice")) {install.packages("lattice"); library(lattice)}
if(!require("gridExtra")) {install.packages("gridExtra"); library(gridExtra)}
if(!require("invgamma")) {install.packages("invgamma"); library(invgamma)}
if(!require("ggsci")) {install.packages("ggsci"); library(ggsci)}
if(!require("viridisLite")) {install.packages("viridisLite"); library(viridisLite)}
if(!require("wesanderson")) {install.packages("wesanderson"); library(wesanderson)}
if(!require("rstudioapi")) {install.packages("rstudioapi"); library(rstudioapi)}
if(!require("coda")) {install.packages("coda"); library(coda)}

############################### PARAMETERS #####################################
print("Setting parameters for simulation")
# the path with main folder, where all functions and result folders are located
main_path <- "C:/Andrzej/OneDrive - SGH/moje_papery/MG_MP/replication_pack_CEJEME"
simulation_path <- paste0(main_path, "/post_simul")
variable <- "GDP" # UE for unemployment rate, GDP for Gross Domestic Product
country <- "PL" # PL for Poland, USA for the United States
eval(parse(text = paste0("simulation_file <- paste0(simulation_path, '/posteriorB_", country, "_", variable, ".RData')")))
load(simulation_file) 
chainA$simul <- chainA$simul[1:16000,]
chainB$simul <- chainB$simul[1:16000,]
main_path <- "C:/Andrzej/OneDrive - SGH/moje_papery/MG_MP/replication_pack_CEJEME"

analysis_path <- paste0(main_path, "/output/analysis_", format(Sys.time(), "%b_%d"))
# the path with saved data (matrix W and datasets)
basic_file <- paste0(main_path, "/basic_data.RData")
# the path to the file with two chains
# the path to the file with classification of regions
file <- paste0(main_path, "/Regions_classification.csv")
# identifiers of dataset


# list of 8 regions which represents macroregions
list <- c("Miasto Kraków","Poznański", "Opolski", "Chojnicki", "Żyrardowski", 
          "Kielecki", "Białostocki", "Miasto Warszawa") 

############################### SETTINGS #######################################
setwd(main_path)
source("0_MC_MG_HF_conv_functions.R") # functions
source("3_analyse_raw_data.R") 
table_reg <- read_delim(file, delim = ";")
load(basic_file)
# settings <- list()
# if(variable == "UE") {
#   settings <- settings_list$PL_UE
# } else {
#   settings <- settings_list$PL_GDP
# }
hyperpar0 = list(alpha_prior=settings$alpha_prior,
                 v_prior=settings$v_prior,
                 delta_prior=settings$delta_prior,
                 m_prior=settings$m_prior,
                 M_prior=settings$M_prior)
attach(hyperpar0)
#ef_impulse <- ifelse(variable =="UE", -1, 1)
if(variable =="UE") {
  dataset <- PL_UE
} else {
  dataset <- PL_GDP
}
if(country =="PL") {
  country_map <- PL_map
} else {
  country_map <- USA_map
} # PL_map for Poland, USA_map for the United States
N <- ifelse(country =="PL", n_regions, n_states)

posterior_a <- chainA$simul
posterior_b <- chainB$simul
posterior <- rbind(posterior_a, posterior_b)
theta_0 <- colMeans(posterior)
theta0 <- relist(theta_0, N)

print("Settings and parameters are loaded")

########################### PARAMETERS FOR PLOTS ################################
n <- N
main_colour <- "grey"
main_colour2<- "black"
cex <- 1
n_col <- 4
n_row <- 7
m <- n_col*n_row
#colnames(Y) <- unique(dataset$Name)
names <- colnames(Y)
names[58] <- "Krosnieński"
names[35] <- "Legnicko-Głogowski"
names[4] <- "Łomżynski"
names[21] <- "Bydgosko-Torunski"
names[32] <- "Miasto Łódz"

pages <- ceiling(N/m)
labels_hamilton <- c("2012","2013","2014","2015","2016","2017","2018","2019","2020","2021")
labels_hamilton_short <- c("2012","2014","2016","2018","2020")
dates <- unique(dataset$Period)
br <- dates[seq(1,110,24)]
f <- "Sans"
#parameters for function which draw spatial effect, it influences the number and colors of breaks
neg_color <- "Greys"
middle_point <- 6
max_point <- 9
print("Parameters for plots are set")

regions_to_plot <- c()
for (i in list){
  regions_to_plot <- append(regions_to_plot, match(i, names))
}
print("List of 8 regions is created")

########################### ANALYSIS OF SIMULATION ###############################

# check if sub directory exists 
if (file.exists(analysis_path)){
  # specifying the working directory
  setwd(file.path(analysis_path))
} else {
  # create a new sub directory inside
  # the main path
  dir.create(file.path(analysis_path))
  # specifying the working directory
  setwd(file.path(analysis_path))
}


########################### ANALYSIS OF RAW DATA ###############################
## TIME SERIES
draw_timelines_selected(PL_UE, lista_PL_UE, 4, 2,"Figure_4","'stopa bezrobocia'", 
                        div=1, width = 8.27, height = 11.69)
draw_timelines_selected(PL_GDP, lista_PL_GDP, 4, 2,"Figure_1","'PKB'", 
                        div=1000, width = 8.27, height = 11.69)

## MAPS
#### UNEMPLOYMENT RATE IN POLAND
require(gridExtra)
temp <- unique(year(PL_UE$Date))
#png(file = paste0("mapM_BEZR_PL.png"), width = 8.27, height = 8, units ="in",res=300)
temp <- c(2011, 2013, 2015, 2017, 2019, 2021)
png(file = paste0("Figure_3.png"), width = 8.27, height = 8, units ="in",res=300)
plots = lapply(temp, function(.x) draw_raw_matrixMaps(PL_UE, PL_map,"'stopa bezrobocia'",8,1,"Greys",.x,FALSE,"equal",temp))
do.call(grid.arrange, c(plots, ncol=3))
dev.off()
png(file = paste0("Figure_3_leg.png"), width = 4, height =3, units ="in",res=300)
draw_raw_matrixMaps(PL_UE, PL_map,"'stopa bezrobocia'",8,1,"Greys",2018,list(space='left',height = 1,width =3,cex=30),"equal",temp)
dev.off()

#### GDP IN POLAND
temp <- unique(year(PL_GDP$Date))
#png(file = paste0("Figure_3.png"), width = 8.27, height = 8, units ="in",res=300)
temp <- c(2000, 2004, 2008, 2012, 2016, 2018)
png(file = paste0("Figure_2.png"), width = 8.27, height = 8, units ="in",res=300)
plots = lapply(temp, function(.x) draw_raw_matrixMaps(PL_GDP,PL_map,"'PKB'",8,1000,"Greys",.x,FALSE,"quantile",temp))
do.call(grid.arrange, c(plots, ncol=3))
dev.off()
# legend
png(file = paste0("Figure_2_leg.png"), width = 4, height =3, units ="in",res=300)
draw_raw_matrixMaps(PL_GDP,PL_map,"'PKB'",8,1000,"Greys",2018,list(space='left',height = 1,width =3),"quantile",temp)
dev.off()


#----------------------------- CONVERGENCE TEST -------------------------------#
print("Convergence test - Gelman, Geweke")

posterior_chA <- posterior_a$rho
posterior_chB <- posterior_b$rho
mh.draws1<-mcmc(posterior_chA)
mh.draws2<-mcmc(posterior_chB)
combined.chains <- mcmc.list(mh.draws1, mh.draws2)
gelman.diag(combined.chains)
png(file=paste0("Figure_5_6_", country, variable, ".png"), 
    width=10, 
    height=5, 
    units="in",
    res=300)
gelman.plot(combined.chains, 
            bin.width=10, 
            max.bins=50,
            confidence=0.95, 
            transform=FALSE, 
            auto.layout=TRUE,
            col=c(main_colour, main_colour2),  
            xlab='', 
            ylab='')
dev.off()
print("Geweke plot is saved")

png(file = paste0("Figure_7a_", country, variable, ".png"), 
    width = 6, 
    height = 5, 
    units ="in",
    res=300)
traceplot(mh.draws1, 
     type="s", 
     xpd=NA, 
     ylab="Rho", 
     xlab="Sample", 
     las=1,  
     main=paste0("Chain 1", " - ", variable))
dev.off()

png(file = paste0("Figure_7b_", country, variable, ".png"), 
    width = 6, 
    height = 5, 
    units ="in",
    res=300)
traceplot(mh.draws2, 
          type="s", 
          xpd=NA, 
          ylab="Rho", 
          xlab="Sample", 
          las=1,  
          main=paste0("Chain 2", " - ", variable))
dev.off()

png(file=paste0("Figure_8_", country, variable, ".png"), 
    width=7, 
    height=5, 
    units="in",
    res=300)
# Calculate the density for each histogram
density1 <- density(posterior_chA)
density2 <- density(posterior_chB)
# Plot only the density line
plot(density1, main = paste0("Density Plot of Two Chains for ", variable), 
     xlab = "Value of Rho", ylab = "Density", 
     col = "grey", lty = 1)
lines(density2, col = "black", lty = 3)
# Add a legend
legend("topright", legend = c("Chain 1", "Chain 2"), fill = c("gray", "black"), border = NA)
dev.off()

png(file=paste0("2ch_autocorr_plot_rho_", country, variable, ".png"), 
    width=7, 
    height=5, 
    units="in",
    res=300)
autocorr.plot(combined.chains)
dev.off()

#------------------------------------------------------------------------------#
sigma_domain <- seq(from = 0, to = max(posterior[,(2*n+2):(3*n+1)]), by = 0.01)
sigma_prior <- dinvgamma(sigma_domain, shape = v_prior/2, scale = delta_prior/2)
m1_domain <- seq(from = min(posterior[,2:(n+1)]), to = max(posterior[,2:(n+1)]), by = 0.01)
m0_domain <- seq(from = min(posterior[,(n+2):(2*n+1)]), to = max(posterior[,(n+2):(2*n+1)]), by = 0.01)
m_domain <- seq(from = min(c(m0_domain, m1_domain)), to = max(c(m1_domain, m0_domain)), by = 0.01)
m1_prior <- dnorm(m_domain, mean=m_prior[2], sd=M_prior[2,2]^0.5)
m0_prior <- dnorm(m_domain, mean=m_prior[1], sd=M_prior[1,1]^0.5)
p_domain <- seq(from=0, to=1, by=0.01)
p11_prior <- dbeta(p_domain, alpha_prior[2,2], alpha_prior[2,1])
p00_prior <- dbeta(p_domain, alpha_prior[1,1], alpha_prior[1,2])

lowerbound_rho <- 1/min(eigen(W)$values)
lowerbound_rho2 <- 0.8
rho_domain <- seq(from=lowerbound_rho2, to=1, by=0.01)
rho_prior <- rep(1/(1-lowerbound_rho2), length(rho_domain))

v_m1<-posterior[,2:(n+1)]
v_m0<-posterior[,(n+2):(2*n+1)]
v_omega<-posterior[,(2*n+2):(3*n+1)]
v_p0<-posterior[,(3*n+2):(4*n+1)]
v_p1<-posterior[,(4*n+2):(5*n+1)]
v_rho<-posterior[,1]

#---------------------------------TABLES---------------------------------------#
rho.sum <- matrix(NA, nrow=4, ncol=4)
colnames(rho.sum) <- c("post mean", "post SD", "HPDI 95 L", "HPDI 95 U")
post <- posterior
post.sum <- matrix(NA, nrow=4, ncol=ncol(posterior))

v_m1<-posterior[,2:(n+1)]
v_m0<-posterior[,(n+2):(2*n+1)]
v_omega<-posterior[,(2*n+2):(3*n+1)]
v_p0<-posterior[,(3*n+2):(4*n+1)]
v_p1<-posterior[,(4*n+2):(5*n+1)]
post.sum[1,] <- colMeans(post)
post.sum[2,] <- vapply(post, 2, FUN = sd)
post.sum[3,] <- vapply(post, 2, FUN = quantile, probs = c(0.025))
post.sum[4,] <- vapply(post, 2, FUN = quantile, probs = c(0.975))
post2 <- cbind(t(post.sum[,(n+2):(2*n+1)]), 
               t(post.sum[, 2:(n+1)]),
               t(post.sum[,(3*n+2):(4*n+1)]),
               t(post.sum[,(4*n+2):(5*n+1)]),
               t(post.sum[,(2*n+2):(3*n+1)]))
rownames(post2) <- names
colnames(post2) <- paste0(rep(c("m0", "m1", "p00", "p11", "sigma"), each = 4), " ", rep(c("post mean", "post SD", "HPDI 95 L", "HPDI 95 U"), 5))
post2 <- round(post2,3)
write.table(post2, 
            file = paste0(country, variable, "_results.csv"), 
            sep = ";", 
            dec = ",")
rho.sum <- post.sum[,1]
print("All results are saved in .csv format")

theta_posterior_means <- list(rho = post.sum[1,1],
                              mu_1 = as.vector(post2[,5]),
                              mu_0 = as.vector(post2[,1]),
                              omega_d = as.vector(post2[,17]),
                              p_00 = as.vector(post2[,9]),
                              p_11 = as.vector(post2[,13]))
p_Hamilton <- Hamilton_filter(Y, theta_posterior_means, W)
p_Hamilton <- p_Hamilton$p_1

#-------------------------------M ILLUSTRATION----------------------------------#
for (i in 1:pages){
    page<-i
    png(file=paste0("m1m0_2chains_", country, variable,"_", page,".png"), 
        width=8.27, 
        height=11.69, 
        units="in",
        res=300)
    par(mfrow=c(n_row, n_col), 
        family="serif",
        mar=c(3, 2, 2, 0) + 0.1,
        mgp=c(1.5,0.2,0))
    for (pp in 1:m) {
      pp<-pp+(page-1)*m
      if (pp<=N){    
        hist(v_m1[,pp], freq=FALSE, main=colnames(Y)[pp], border=rgb(1, 1, 1, 0, maxColorValue=1),
             xlab=NULL, ylab=NULL, nclass=20, col="white",#col = rgb(0, 0, 0, 0.5, maxColorValue = 1),
             xlim=c(min(m_domain), max(m_domain)), #ylim = c(0,1.5), 
             cex.main=cex, cex.axis=cex/1.2, tck=-0.02)
        hist(v_m0[,pp], freq = FALSE, main=colnames(Y)[pp], border=rgb(1, 1, 1, 0, maxColorValue = 1),
             xlab = NULL, ylab = NULL, nclass=20, col=rgb(1, 0, 0, 0.5, maxColorValue = 1),
             xlim = c(min(m_domain), max(m_domain)), #ylim = c(0,1.5), 
             cex.main=cex, cex.axis=cex/1.2, add=TRUE)
        lines(x=m_domain, y=m1_prior, lwd=2, col="darkgrey")
        lines(x=m_domain, y=m0_prior, lwd=2, col=main_colour2)
        legend(x="topleft", legend=c("m1 prior", "m1 posterior", "m0 prior", "m0 posterior"), 
               fill = c("darkgrey", "white", main_colour2, rgb(1, 0, 0, 0.5, maxColorValue=1)), 
               bty="n", cex=cex/1.4)}}
    dev.off()
}
print("Plots of m_00 and m_11 distribution are saved")

#-------------------------------P ILLUSTRATION----------------------------------#
for (i in 1:pages){
  page<-i
  png(file = paste0("p1p0_2chains_", country,variable, "_",page, ".png"), 
      width = 8.27, 
      height = 11.69, 
      units ="in", 
      res=300)
  par(mfrow = c(n_row, n_col),
      family="serif",
      mar=c(3, 1, 2, 1)+ 0.1,
      mgp=c(1.5,0.2,0))
  for (pp in 1:m) {
    pp<-pp+(page-1)*m
    if (pp<=N){
      hist(v_p1[,pp], freq=FALSE, main=names[pp], col="white",
           xlab=NULL, ylab=NULL, nclass=10, #col = rgb(0, 0, 0, 0.5, maxColorValue = 1),
           xlim=c(min(p_domain), max(p_domain)), #ylim = c(0, 8), 
           cex.main=cex, cex.axis = cex/1.2,tck=-0.02)
      hist(v_p0[,pp], freq=FALSE, main=names[pp], border=main_colour2,
           xlab=NULL, ylab=NULL, nclass=10, col=rgb(1, 0, 0, 0.5, maxColorValue=1),
           xlim=c(min(p_domain), max(p_domain)), #ylim = c(0, 8), 
           add=TRUE, cex.main = cex, cex.axis = cex/1.2)
      lines(x=p_domain, y=p11_prior, lwd = 2, col = "darkgrey")
      lines(x=p_domain, y=p00_prior, lwd = 2, col = main_colour2)
      legend(x="topleft", legend = c("p11 prior", "p11 posterior", "p00 prior", "p00 posterior"), 
             fill = c("darkgrey", "white", main_colour2, rgb(1, 0, 0, 0.5, maxColorValue = 1)), 
             bty = "n", cex = cex/1.4)
    }}
  dev.off()
}
print("Plots of p_00 and p_11 distribution are saved")

#-------------------------------RHO ILLUSTRATION----------------------------------#
png(file = paste0("Figure_15_", country, "_", variable, ".png"), 
    width=5, 
    height=5, 
    units="in", 
    res=300)

hist(v_rho, 
     freq=FALSE, 
     main="", 
     xlab=NULL, 
     ylab=NULL, 
     nclass=20, 
     col=rgb(0, 0, 0, 0.5, maxColorValue=1),
     xlim = c(lowerbound_rho2, 1), 
     #ylim = c(0, 15), 
     cex.main = cex, 
     cex.axis = cex/1.7)
par(new = TRUE)
lines(x = rho_domain, y = rho_prior, lwd = 2, col = "grey")
legend(x = "topleft", legend = c("rho prior", "rho posterior"), 
       fill = c("grey", rgb(0, 0, 0, 0.5, maxColorValue = 1)), 
       bty = "n", cex = cex/2)
dev.off()
print("Histogram of rho is saved")

#---------------------------------HAMILTON-------------------------------------# 
for (i in 1:pages){
  page<-i
  if (variable == "GDP") etykiety <- lata[2:19]
  eval(parse(text = paste0("if (variable == 'UE') etykiety <- year(unique(", country, "_UE_ch$Period))")))
  png(file = paste0("Hamilton_2chains_", country,variable,"_",page,".png"), 
      width = 8.27, 
      height = 11.69, 
      units = "in",
      res = 300)
  par(mfrow = c(n_row, n_col), 
      family="serif",
      mar=c(2, 2, 2, 0.5) + 0.1,
      mgp=c(1,0,0))
  eval(parse(text = paste0("dates <- unique(PL_", variable, "_ch$Period)")))
  for (pp in 1:m) {
    pp <- pp+(page-1)*m
    if (pp <= N){ 
      plot(x=1:length(dates), 
           y=p_Hamilton[,pp],
           type="l", 
           lwd=2,
           xlab="",
           ylab="probability of expansion",
           main=names[pp],
           col="black",
           cex.axis=cex/2,
           cex.main=cex,
           cex.lab=cex/2,
           xaxt="none",
           yaxt="none")
      axis(2, cex.axis=cex/1.5, tck=-0.015)
      axis(1, seq(1,nrow(p_Hamilton),1), cex.axis=cex/1.5, srt=45,tck=-0.015, #col.axis="red",
           labels=etykiety)
    }
  }
  dev.off()
}
print("Plots of the probability of expansion for every region are saved")

#-------------------------------SPATIAL EFFECT---------------------------------# 
middle_point <- 3
max_point <- 8
e <- c()
impulse <- abs(theta_posterior_means$mu_1- theta_posterior_means$mu_0)
for (pp in 1:N) {
  impulse2 <- as.matrix(rep(0,N))
  impulse2[pp] <- impulse[pp]
  effect <- solve(diag(N) - theta_posterior_means$rho * W) %*% impulse2
  e <- cbind(e, effect)}
effect_mean <- mean(e)
vec_e <- c(e)
breaks_qt <- classIntervals(vec_e, 2, n=max_point, style="jenks")
r <- breaks_qt$brks 
n_col <- 4
n_row <- 5
m <- n_col*n_row
pages <- ceiling(N/m)
pal <- c()

png(file = paste0("legend_effect_", country, "_", variable, ".png"), 
    width=8.27, 
    height=11.69, 
    units="in",
    res=300)
draw_impulse2(country_map, N, names, theta_posterior_means, W, e, r, TRUE, 1)
dev.off()

for (page in 1:pages){
  if (m+(page-1)*(m-1)>N){
    dif <- m - (N-(m+(page-1)*(m-1))+1)
    temp <- seq(1+(page-1)*(m-1),N)
    png(file = paste0("effect_", country, variable, "_", page, ".png"), 
        width=8.27, height=11.69, units="in",res=300)
    plots = lapply(temp, function(.x) draw_impulse2(country_map, 73, names, 
                                                    theta_posterior_means, 
                                                    W, e, r, FALSE, .x))
    p<-marrangeGrob(plots, nrow=n_row, ncol=n_col, top = NULL)
    print(p)
    dev.off()
  }else{
    temp <- seq(1+(page-1)*(m-1), m+(page-1)*(m-1)-1)
    png(file=paste0("effect_", country, variable, "_", page, ".png"), 
        width=8.27, height=11.69, units ="in",res=300)
    plots = lapply(temp, function(.x) draw_impulse2(country_map, 73, names, 
                                                    theta_posterior_means, 
                                                    W, e, r, FALSE, .x))
    do.call(grid.arrange, plots)
    dev.off()
  }
}
print("Plots of spatial effects are saved")

############################### 8 REGIONS #######################################
#-------------------------------M ILLUSTRATION----------------------------------#  

png(file=paste0("Figure_9_10_", country, variable, ".png"), 
    width=9, 
    height=5, 
    units="in", 
    res=300)
par(mfrow=c(2, 4), 
    family="serif", 
    mar=c(3, 2, 2, 0)+ 0.1, 
    mgp=c(1.5,0.2,0))
up_limit <- ifelse(variable == "UE", 5, 1)
for (pp in regions_to_plot) {
    hist(v_m1[,pp], freq=FALSE, main=colnames(Y)[pp], border="black",
         xlab=NULL, ylab=NULL, nclass=20, col=rgb(0,0,0,0.5),#col = rgb(0, 0, 0, 0.5, maxColorValue = 1),
         xlim=c(min(m_domain), max(m_domain)), ylim=c(0,up_limit), 
         cex.main=cex, cex.axis = cex/1.2,tck=-0.02)
    hist(v_m0[,pp], freq=FALSE, main = colnames(Y)[pp], border=rgb(0.5,0.5,0.5,0.5),
         xlab=NULL, ylab=NULL, nclass = 20, col = rgb(0.5,0.5,0.5,0.5),
         xlim=c(min(m_domain), max(m_domain)), ylim=c(0,up_limit), 
         cex.main = cex, cex.axis=cex/1.2, add=TRUE)
    lines(x=m_domain, y=m1_prior, lwd=2, col=main_colour2)#, lty = "dashed")
    lines(x=m_domain, y=m0_prior, lwd=2, col=rgb(0.5,0.5,0.5,0.5))
    legend(x="topleft", 
           legend = c("m1 prior (line)", "m1 posterior (bars)", "m0 prior (line)", "m0 posterior (bars)"), 
           fill = c("black", rgb(0,0,0,0.5), rgb(0.5,0.5,0.5,0.5), rgb(0.5,0.5,0.5,0.5)), 
           #density = c(20, NA, NA, NA),
           bty = "n", 
           cex=cex/1.4)
  }
dev.off()

#-------------------------------P ILLUSTRATION----------------------------------#
png(file=paste0("Figure_11_12_", country, variable, ".png"), 
    width=9, 
    height=4, 
    units="in", 
    res=300)
par(mfrow=c(2, 4), 
    family="serif",
    mar=c(3, 1, 2, 1)+ 0.1, 
    mgp=c(1.5,0.2,0))
for (pp in regions_to_plot) {
    hist(v_p1[,pp], freq=FALSE, main=names[pp], col="white",
         xlab=NULL, ylab=NULL, nclass=10, border="black",
         xlim=c(min(p_domain), max(p_domain)), ylim = c(0, 18), 
         cex.main=cex, cex.axis=cex/1.2, tck=-0.02)
    hist(v_p0[,pp], freq=FALSE, main=names[pp], border=main_colour2,
         xlab=NULL, ylab=NULL, nclass=10, col=rgb(0.5,0.5,0.5,0.5),
         xlim=c(min(p_domain), max(p_domain)), ylim=c(0, 18), 
         add=TRUE, cex.main=cex, cex.axis=cex/1.2)
    if (variable == "GDP") {
      lines(x=p_domain, y=p11_prior, lwd=2, col="black")
      legend(x="topleft", 
             legend = c("p11 posterior", "p00 posterior", "p11 & p00 prior"), 
             fill = c("white", "darkgrey", "black"), 
             bty="n", 
             cex=cex/1.4)
    }
    if (variable == "UE") {
      lines(x=p_domain, y=p11_prior, lwd=2, col=main_colour2, lty = "dashed")
      lines(x=p_domain, y=p00_prior, lwd=2, col=main_colour2)
      legend(x="topleft", 
             legend = c("p11 posterior", "p00 posterior", "p11 prior", "p00 prior"), 
             fill = c("white", "darkgrey", "black", "black"), 
             density = c(NA, NA, 40, NA),
             bty="n", 
             cex=cex/1.4)
    }
  }
dev.off()

#-------------------------------SPATIAL EFFECT---------------------------------# 
png(file=paste0("Figure_16_17_", country, "_", variable, ".png"), 
    width=9, 
    height=5, 
    units="in", 
    res=300)
plots = lapply(regions_to_plot, function(.x) draw_impulse2(country_map, 73, names, 
                                                           theta_posterior_means, 
                                                           W, e, r, FALSE, .x))
marrangeGrob(plots, 
             nrow=2, 
             ncol=4, top = NULL)
dev.off()


regions_to_plot <- list()
for (macro in unique(table_reg$WOJEWODZTWO)){
  temp_table_reg <- table_reg %>% filter(WOJEWODZTWO == macro)
  print(temp_table_reg)
  list_reg <- unique(temp_table_reg$Name)
  for (i in list_reg){
    regions_to_plot[[macro]] <- append(regions_to_plot[[macro]], match(i, names))
  }
}

#---------------------------------HAMILTON-------------------------------------#
par(mfrow=c(4, 4), 
    family="serif", 
    mar=c(2, 2, 2, 0.5) + 0.1, 
    mgp=c(1, 0, 0))

dH <- data.frame(p_Hamilton)
names(dH) <- names
i <- 0
list_macro <- list()
for (macro in unique(table_reg$WOJEWODZTWO)) {
  i<-i+1
  num_reg <- regions_to_plot[[macro]]
  temp_dH <- dH %>% 
    select(num_reg) %>% 
    mutate(year=dates) %>%
    pivot_longer(!year, names_to="region", values_to="value")
  print(temp_dH)
  temp_dH$macro <- macro
  list_macro[[i]] <- temp_dH
}
final_table <- data.frame(do.call(rbind, list_macro))
if (variable == "GDP") {
  final_table$year <- as.numeric(year(final_table$year))
} else {
  final_table$year <- paste0(year(final_table$year), "-", month(final_table$year))
}
do.call(rbind, 
        list_macro)


p <- ggplot(final_table, aes(year, value, group=region, colour=factor(region))) + 
  geom_line(size=0.5) +
  theme_classic() +
  scale_colour_grey() +
  xlab('date') +
  ylab("probability") +
  facet_wrap(vars(macro)) +
  scale_x_discrete(breaks=br, 
                   labels=labels_hamilton_short) +
  theme(legend.position="none")
fname <- paste0("Figure_13_14_", country, variable, ".png")
ggsave(fname, 
       plot=p, 
       device = "png",
       width=28, 
       height=22, 
       units="cm", 
       dpi=300)

#---------------------------OTHER CONVERGENCE TESTS----------------------------#
main_colour <- "darkgray"
main_colour2<- "black"

mh.draws <- mcmc(posterior)
tested <- seq(from=1,to=16000,by=1)#1:16690
combined_chains <- mcmc.list(mcmc(posterior_a[tested,]), mcmc(posterior_b[tested,]))
gelman.diag(combined_chains)
a <- geweke.diag(mh.draws)

conv_test <- summary(mh.draws)

gelman.diag(combined_chains, 
            confidence = 0.95, 
            transform = FALSE, 
            autoburnin = TRUE,
            multivariate = TRUE)
a <- geweke.diag(mh.draws)

n <- 0
for (i in a$z){
  n<-n+1
  ifelse(abs(i) > 1.96, print(a$z[n]), "FALSE")
}

table(abs(a$z) > 1.96)
df<-as.data.frame(cbind(a$z))
png(file = paste0("Z_plot_2CH", country, "_", variable, ".png"), 
    width = 10, 
    height = 5, 
    units ="in",
    res=300)
barplot(df[,1], 
        main="Z-score for given parameters", 
        xlab="parameter", 
        border =FALSE, 
        col=main_colour)
abline(h=c(-1.96, 1.96), col=main_colour2)
abline(h=0, col='black')
dev.off()

#---------------------------DECOMPOSITION----------------------------#

if(variable == "UE") {
  decomp_s_Hamilton <- p_Hamilton>0.5
  decomp_local_cycle <- decomp_s_Hamilton * kronecker(matrix(1, nrow = nrow(decomp_s_Hamilton), ncol = 1), t(as.matrix(theta_posterior_means$mu_1))) +
                        (1-decomp_s_Hamilton) * kronecker(matrix(1, nrow = nrow(decomp_s_Hamilton), ncol = 1), t(as.matrix(theta_posterior_means$mu_0)))
  decomp_systematic <- t(solve(diag(73) - theta_posterior_means$rho) %*% t(decomp_local_cycle))
  decomp_external_cycle <- decomp_systematic - decomp_local_cycle
  decomp_epsilon <- Y - decomp_systematic
  
  fit_l <- data.frame(names = PL_map@data$NAME_LATN, local_share = rep(NA, 73), external_share = rep(NA, 73), residual_share = rep(NA, 73))
  for (ii in 1:ncol(decomp_local_cycle)) {
    temp <- lm(Y[,ii] ~ decomp_systematic[,ii])
    fit_l$residual_share[ii] <- 1 - summary(temp)$r.squared
    temp2 <- lm(decomp_systematic[,ii] ~ decomp_local_cycle[,ii])
    fit_l$local_share[ii] <- summary(temp2)$r.squared * (1 - fit_l$residual_share[ii])
    fit_l$external_share[ii] <- 1 - fit_l$residual_share[ii] - fit_l$local_share[ii]
  }
  #View(fit_l)
  fit_l[,2:4] <- round(fit_l[,2:4], 3)
  PL_map <- merge(x = PL_map, y = fit_l, by.x = "NAME_LATN", by.y = "names")
  pal <- colorRampPalette(c("white", "black"), bias = 1)
  png("Figure_18.png")
  spplot(PL_map, zcol = "local_share", colorkey = TRUE, col.regions = pal(100), cuts = 99,
         par.settings = list(axis.line = list(col =  'transparent')),
         main = "")
  dev.off()
  png("Figure_19.png")
  spplot(PL_map, zcol = "external_share", colorkey = TRUE, col.regions = pal(100), cuts = 99,
         par.settings = list(axis.line = list(col =  'transparent')),
         main = "")
  dev.off()
}
