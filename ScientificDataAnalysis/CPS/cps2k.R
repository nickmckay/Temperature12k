
library(geoChronR) #devtools::install_github("nickmckay/geoChronR")
library(lipdR) #devtools::install_github("nickmckay/lipd-utilities",subdir = "R")
library(purrr)
library(magrittr)
library(ggplot2)
library(compositeR)#devtools::install_github("nickmckay/compositeR")
library(foreach)
library(doParallel)
#load database
#load database
D <- readLipd("PAGES2kTemp_v2_1_0/")

#extract timeseries
TS <- extractTs(D)

#filter by compilation
fTS <- filterTs(TS,"paleoData_useInGlobalTemperatureAnalysis ==  TRUE")

#bin the TS
binvec <-  seq(0, to = 2000, by = 1)
binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

dur <- map_dbl(fTS,function(x){abs(diff(range(x$year)))})


#setup latbins
allLatBins <- vector(mode = "list",length = 4)
allLatBins[[1]] <- seq(-90,90,by = 30)
allLatBins[[2]] <- c(-90,-30,0,30,90)
allLatBins[[3]] <- c(-90,0,90)
allLatBins[[4]] <- c(-90,90)
library(foreach)
library(doParallel)
registerDoParallel(4)
for(alb in 1){#:length(allLatBins)){

#composite lat bins
latbins <- allLatBins[[alb]]
lat <- geoChronR::pullTsVariable(fTS,"geo_latitude")

#load in scaling data
targets <- list.files(".",pattern = "ERA",full.names = TRUE)
targetsShort <- list.files(".",pattern = "ERA",full.names = FALSE)
targ <- purrr::map(targets,read.csv)

#scaling window
sw <- 100

#set up ensembles?
nens <- 100

scaled <- comps <- counts <- c()
for(lb in 1:(length(latbins)-1)){
  scaleEns <- c()
  fi <- which(lat > latbins[lb] & lat <= latbins[lb+1])

  dur <- map_dbl(fTS[fi],function(x){abs(diff(range(x$year)))})
sw <- floor(min(dur))

scaledList <- foreach(n = 1:nens) %dopar% {
  tc <- compositeEnsembles(fTS[fi],binvec,ageVar = "year",spread = TRUE,duration = sw, searchRange = c(0,2000),binFun = simpleBinTs)
  comps <- cbind(comps,tc$composite)
  counts <- cbind(counts,tc$count)

  thisTarget <- which(stringr::str_starts(string = targetsShort, paste0(latbins[lb],"to",latbins[lb+1])))
  if(length(thisTarget) != 1){
    stop("target matching problem")
  }



  thisScaled <- scaleComposite(composite = tc$composite,binvec = binvec,scaleYears = targ[[thisTarget]][,1],scaleData = targ[[thisTarget]][,-1])

  return(thisScaled)

}

  scaleEns <-  as.matrix(purrr::map_dfc(scaledList,extract))

  out <- cbind(binAges,scaleEns)
  write.csv(x = out,file = paste0(".",paste0(latbins[lb] ,"to", latbins[lb+1],"-scaleWindow",sw,"-PAGES2k.csv")), row.names = FALSE,col.names = FALSE)
}



#plot 2k reconstructions
sw <- 200
targets <- list.files(".",pattern = "PAGES",full.names = TRUE)
targ <- purrr::map(targets,read.csv)

colorsHi <- RColorBrewer::brewer.pal(6,"Set3")
plot2k <- ggplot()

for(lb in 1:(length(latbins)-1)){
  fi <- which(lat > latbins[lb] & lat <= latbins[lb+1])
  dur <- map_dbl(fTS[fi],function(x){abs(diff(range(x$year)))})
  print(floor(min(dur)))
  thisTarget <- which(grepl(targets,pattern = paste0(latbins[lb] ,"to", latbins[lb+1],"-scaleWindow",sw,"-PAGES2k.csv")))
  if(length(thisTarget) != 1){
    stop("target matching problem")
  }

  out <- as.matrix(targ[[thisTarget]])

  plot2k <- plotTimeseriesEnsRibbons(plot2k,X = out[,1], Y = out[,-1],alp = .5,colorHigh = colorsHi[lb],lineColor = colorsHi[lb],lineWidth = 1)+
    geom_text(aes(x = 1500), y = (lb * .5) - 0 ,label = paste(latbins[lb],"to",latbins[lb+1]),color = colorsHi[lb])

}
plot2k





scaleEns <- as.matrix(out[,-1])
plotTimeseriesEnsLines(X = binAges,Y = scaleEns,maxPlotN = 100)





#weight by areas
zonalWeights <- sin(latbins[-1]*pi/180)-sin(latbins[-length(latbins)]*pi/180)
zonalWeights <- zonalWeights/sum(zonalWeights)


zonalNames <- stringr::str_c(latbins[-1]," to ",latbins[-length(latbins)])
scaledDf <- as.data.frame(scaled)
names(scaledDf) <- zonalNames
scaledDf$year <- binAges
GlobalMean <- rowSums(t(t(scaled)*zonalWeights))

tidyScale <- tidyr::pivot_longer(scaledDf,cols = -year,names_to = "Latitude Band")
tidyScale$`Latitude Band` <- factor(tidyScale$`Latitude Band`,levels = rev(zonalNames))

library(ggplot2)
ggplot()+geom_line(data = tidyScale,aes(x = year, y = value, colour = `Latitude Band`))+
  geom_line(aes(x = binAges,y = GlobalMean),color = "black",size = 1.5)+
  scale_x_continuous(name = "Year BP")+
  scale_y_continuous(name = "Temperature (wrt 1850-2000) (deg C)")+
  theme_bw()


}



