library(geoChronR) #devtools::install_github("nickmckay/geoChronR")
library(lipdR) #devtools::install_github("nickmckay/lipd-utilities",subdir = "R")
library(purrr)
library(magrittr)
library(ggplot2)
library(compositeR)#devtools::install_github("nickmckay/compositeR")
library(foreach)
library(doParallel)
#load database
setwd("~/GitHub/Temperature12k/ScientificDataAnalysis/DCC/")
td <- getwd()
#load("~/Dropbox/Temp12kSerialization/Temp12k/expandedMasterDatabase/newEnsemblesIn.RData")
D <- readLipd("../lipdFilesWithEnsembles/")
setwd(td)

#extract timeseries
TS <- extractTs(D)


sg <- pullTsVariable(TS,variable = "interpretation1_seasonalityGeneral")
ic <- pullTsVariable(TS,"paleoData_inCompilation")
u <- pullTsVariable(TS,"paleoData_units")


#filter by compilation and seasonality
tu <- which(tolower(ic) == "temp12kensemble" & (tolower(sg) == "annual" | tolower(sg) == "summeronly" | tolower(sg) == "winteronly") & tolower(u) == "degc")

fTS <- TS[tu]

# ls <- map_dbl(fTS,function(x) sum(!is.na(x$paleoData_values) & !is.na(x$age)))
# ls2 <- map_dbl(fTS,function(x) length(x$paleoData_values))
#
# fTS <- fTS[which(ls > 10 & ls2 >10)]


#bin the TS
binvec <-  seq(-50, to = 12050, by = 100)
binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

#setup ensemble
nens <- 500


#composite lat bins
latbins <- seq(-90,90,by = 30)
lat <- lipdR::pullTsVariable(fTS,"geo_latitude")

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
    tc <- compositeEnsembles(fTS[fi],binvec,spread = TRUE,duration = 3000, searchRange = c(0,7000),gaussianizeInput = FALSE,ageVar = "ageEnsemble",normalizeVariance = FALSE)
    #tc <- compositeEnsembles(fTS[fi],binvec,spread = spread,...)
    comps <- cbind(comps,tc$composite)
    counts <- cbind(counts,tc$count)

  # #  thisTarget <- which(grepl(targets,pattern = paste0(latbins[lb],"to",latbins[lb+1])))
  #   thisTarget <- which(stringr::str_starts(string = targetsShort, paste0(latbins[lb] ,"to", latbins[lb+1],"-scaleWindow",sw,"-PAGES2k.csv")))
  #
  #   if(length(thisTarget) != 1){
  #     stop("target matching problem")
  #   }
  #
  #   thisScaled <- scaleComposite(composite = tc$composite,binvec = binvec,scaleYears = 1950-targ[[thisTarget]][,1],scaleData = targ[[thisTarget]][,-1],scaleWindow = 1950-c(0,2000),scaleVariance =FALSE)


    scaled <- cbind(scaled,tc$composite)
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

allLatMeans <- as.matrix(cbind(binAges, map_dfc(ensOut[[alb]],extract2,"GlobalMean")))
plotTimeseriesEnsRibbons(X = allLatMeans[,1],Y = allLatMeans[,-1])
#write out data
settings <- paste0(nens,"-",length(latbins)-1,"bands")
readr::write_csv(path = paste0("globalMean",settings,".csv"),x = as.data.frame(allLatMeans),col_names = FALSE)

for(lb in 1:(length(latbins)-1)){
  lbn <- paste0(latbins[lb],"to",latbins[lb+1])
  out <- cbind(binAges,as.matrix(map_dfc(ensOut[[alb]],extract2,lb)))
  write.csv(file = paste0(lbn,settings,".csv"),x = out)

}


