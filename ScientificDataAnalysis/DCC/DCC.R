library(geoChronR) #devtools::install_github("nickmckay/geoChronR")
library(lipdR) #devtools::install_github("nickmckay/lipd-utilities",subdir = "R")
library(purrr)
library(magrittr)
library(ggplot2)
library(compositeR)#devtools::install_github("nickmckay/compositeR")
library(foreach)
library(doParallel)
#load database
D <- readLipd("../lipdFilesWithEnsembles/")

#extract timeseries
TS <- extractTs(D)


sg <- pullTsVariable(TS,variable = "interpretation1_seasonalityGeneral")
ic <- pullTsVariable(TS,"paleoData_inCompilation")
u <- pullTsVariable(TS,"paleoData_units")


#filter by compilation and seasonality
tu <- which(tolower(ic) == "temp12kensemble" & (tolower(sg) == "annual" | tolower(sg) == "summeronly" | tolower(sg) == "winteronly") & tolower(u) == "degc")

fTS <- TS[tu]

ls <- map_dbl(fTS,function(x) sum(!is.na(x$paleoData_values) & !is.na(x$age)))
ls2 <- map_dbl(fTS,function(x) length(x$paleoData_values))

fTS <- fTS[which(ls > 10 & ls2 >10)]


#bin the TS
binvec <-  seq(-50, to = 12050, by = 100)
binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

#setup ensemble
nens <- 500


#composite lat bins
latbins <- seq(-90,90,by = 30)
lat <- geoChronR::pullTsVariable(fTS,"geo_latitude")

scaled <- comps <- counts <- c()
ensOut <- vector(mode = "list",length = nens)
library(foreach)
library(doParallel)
registerDoParallel(4)

allLatBins <- vector(mode = "list",length = 4)
allLatBins[[1]] <- seq(-90,90,by = 30)
allLatBins[[2]] <- c(-90,-30,0,30,90)
allLatBins[[3]] <- c(-90,0,90)
allLatBins[[4]] <- c(-90,90)

#load in scaling data
targets <- list.files(".",pattern = "PAGES2k",full.names = TRUE)
targetsShort <- list.files(".",pattern = "PAGES2k",full.names = FALSE)
sw <- 100 #pages 2k scaling window
targ <- purrr::map(targets,read.csv)

#for(alb in 1){
alb <- 1
sw <- 100
latbins <- allLatBins[[alb]]


#for(alb in 1){
alb <- 1
latbins <- allLatBins[[alb]]

ensOut[[alb]] <- foreach(i = 1:nens) %dopar% {
  if(i%%10 == 0){
    write.csv(i,file = "~/Desktop/nscprog.csv")
  }

  scaled <- c()
  for(lb in 1:(length(latbins)-1)){
    fi <- which(lat > latbins[lb] & lat <= latbins[lb+1])
    tc <- compositeEnsembles(fTS[fi],binvec,spread = TRUE,duration = 3000, searchRange = c(0,7000),gaussianizeInput = FALSE,normalizeVariance = FALSE)
    #tc <- compositeEnsembles(fTS[fi],binvec,spread = spread,...)
    comps <- cbind(comps,tc$composite)
    counts <- cbind(counts,tc$count)

  #  thisTarget <- which(grepl(targets,pattern = paste0(latbins[lb],"to",latbins[lb+1])))
    thisTarget <- which(stringr::str_starts(string = targetsShort, paste0(latbins[lb] ,"to", latbins[lb+1],"-scaleWindow",sw,"-PAGES2k.csv")))

    if(length(thisTarget) != 1){
      stop("target matching problem")
    }

    thisScaled <- scaleComposite(composite = tc$composite,binvec = binvec,scaleYears = 1950-targ[[thisTarget]][,1],scaleData = targ[[thisTarget]][,-1],scaleWindow = 1950-c(0,2000),scaleVariance =FALSE)


    scaled <- cbind(scaled,thisScaled)
  }

  #weight by areas
  zonalWeights <- sin(latbins[-1]*pi/180)-sin(latbins[-length(latbins)]*pi/180)
  zonalWeights <- zonalWeights/sum(zonalWeights)


  zonalNames <- stringr::str_c(latbins[-1]," to ",latbins[-length(latbins)])
  scaledDf <- as.data.frame(scaled)
  names(scaledDf) <- zonalNames
  scaledDf$year <- binAges
  scaledDf$GlobalMean <- rowSums(t(t(scaled)*zonalWeights))
  #scaledDf$counts <- scaled[,ncol(scaled)]



  return(scaledDf)
  # ensOut[[i]] <- scaledDf
  # print(i)
}

save(list = c("ensOut"),file = "12kensOutCNS500.RData")


#plotting!
library(magrittr)
library(purrr)
library(tidyverse)
allLatMeans <- cbind(binAges, map_dfc(ensOut[[alb]],extract2,"GlobalMean"))

#write out data
settings <- paste0(nens,"-",length(latbins)-1,"bands-",sw,"yr2kwindow")
write_csv(path = paste0("globalMean",settings,".csv"),x = allLatMeans,col_names = FALSE)



library(ggplot2)
globMean <- plotTimeseriesEnsRibbons(X = binAges,Y = as.matrix(allLatMeans),x.bin = seq(-1,12000,by = 10))+
  scale_x_reverse(name = "Year (BP)",breaks = seq(0,12000,2000),oob = scales::squish)+
  scale_y_continuous(name = "Temperature (deg C) (wrt 1000-2000 AD)",oob = scales::squish)+
  theme_bw()+
  ggtitle("Global Mean Temperature (Composite Plus Scale)")
globMean

ggsave(filename = paste0("DCC-GlobalMean12k-2k-",alb,".pdf"),globMean )

#plot bands:


colorsHi <- RColorBrewer::brewer.pal(6,"Dark2")

plot12k <- ggplot()

for(lb in 1:(length(latbins)-1)){

  out <- as.matrix(map_dfc(ensOut[[alb]],extract2,lb))
  write_csv(path = paste0("latBand",latbins[lb],"to",latbins[lb+1],".csv"),x = as.data.frame(cbind(binAges,out)),col_names = FALSE)

  plot12k <- plotTimeseriesEnsRibbons(plot12k,X = binAges, Y = out,alp = .5,colorHigh = colorsHi[lb],lineColor = colorsHi[lb],lineWidth = 1,x.bin = seq(-1,12000,by = 10))+
    geom_text(aes(x = 6000), y = (lb * .25) - 3 ,label = paste(latbins[lb],"to",latbins[lb+1]),color = colorsHi[lb])

}
plot12k <- plot12k +
  scale_x_reverse(name = "Year (BP)",breaks = seq(0,12000,2000))+
  scale_y_continuous(name = "Temperature (deg C) (wrt 1000-2000 AD)",oob = scales::squish)+
  ggtitle("Zonal Mean Temperature (Composite Plus Scale)")+
  theme_bw()

ggsave(filename = paste0("DCC-LatBands12k-2k-",alb,".pdf"),plot12k )
plot12k




#plot 2k reconstructions
targets <- list.files(".",pattern = "PAGES",full.names = TRUE)
targetsShort <- list.files(".",pattern = "PAGES",full.names = FALSE)

targ <- purrr::map(targets,read.csv)

plot2k <- ggplot()


for(lb in 1:(length(latbins)-1)){
plotlb <- ggplot()
  thisTarget <- which(stringr::str_starts(string = targetsShort, paste0(latbins[lb],"to",latbins[lb+1])))

  if(length(thisTarget) != 1){
    stop("target matching problem")
  }

  out <- as.matrix(targ[[thisTarget]])
  out2 <- as.matrix(map_dfc(ensOut[[alb]],extract2,lb))
  ba2 <- 1950-binAges
  out2 <- out2[which(ba2 > 0), ]
  ba2 <- ba2[which(ba2 > 0) ]


  #plot this band
  plotlb <- plotTimeseriesEnsRibbons(plotlb,X = out[,1], Y = scale(out[,-1],scale = FALSE),alp = .8,colorHigh = colorsHi[lb],lineColor = colorsHi[lb],lineWidth = 1,x.bin = seq(0,2000,by=2)) %>%
      plotTimeseriesEnsLines(X = ba2, Y = scale(out2,scale = FALSE))+
    scale_x_continuous(name = "Year (AD)",breaks = seq(0,2000,500))+
    scale_y_continuous(name = "Temperature (deg C) (wrt 1-2000 AD)",oob = scales::squish)+
    ggtitle( paste(latbins[lb],"to",latbins[lb+1]))+
    theme_bw()

  ggsave(filename = paste0("12k2kcompLat_",latbins[lb],"to",latbins[lb+1],"-",alb,".pdf"),plot = plotlb)

  #plot all of them
  plot2k <- plotTimeseriesEnsRibbons(plot2k,X = out[,1], Y = out[,-1],alp = .5,colorHigh = colorsHi[lb],lineColor = colorsHi[lb],lineWidth = 1,x.bin = seq(0,2000,by=2))+
    geom_text(aes(x = 1500), y = (lb * .35), label = paste(latbins[lb],"to",latbins[lb+1]),color = colorsHi[lb])

  out <- as.matrix(map_dfc(ensOut[[alb]],extract2,lb))



}


