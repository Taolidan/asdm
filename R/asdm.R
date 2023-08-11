#' @title Approximate Species Distribution Model
#' @description
#' Developing Approximate Species Distribution Model.
#' @param xylist Thinned coordinates of the target species, dataframe with two columns (X, Y)
#' @param map1 Predicted distribution map of species A (Raster* classes, one layer)
#' @param map2 Predicted distribution map of species B (Raster* classes, one layer)
#' @param map3 Predicted distribution map of combination data (Raster* classes, one layer)
#' @param valuescale The theoretical maximum of inputted predicted distribution map (1 or 1000)
#' @param overlap The overlap value matrix/dataframe(3*3) of target species, A and B (follow this order of species)
#' @param k1 Initial value of k1 (a float or "random" which would create a float between 0.7 and 1.3)
#' @param k2 Initial value of k2 (a float or "random" which would create a float between 0.7 and 1.3)
#' @param datasplit Vector with 3 numbers in range of 0-1, which referred to: the proportion of train, validation and test data
#' @param valiMaxTech Maximized criterion by validation, only support "entropy" and "maxtss"
#' @param k_rep Replicate number of k calculation
#' @param tss_rep Replicate number of sampling when calculating maxtss value
#' @param ACC_rep Replicate number of sampling when calculating ACC value
#' @param REC_rep Replicate number of sampling when calculating REC value
#' @param F1_rep Replicate number of sampling when calculating F1 value
#' @param dJ Gradient
#' @param kE.p Initial kE
#' @param step Step
#' @param threshold Exit threshold
#' @param k1.train Whether k1 is trainable, TRUE or FALSE
#' @param k2.train Whether k2 is trainable, TRUE or FALSE
#' @importFrom stats na.omit runif
#' @export
asdm <- function(xylist, map1, map2, map3, valuescale=1000, overlap, k1=1.0, k2=1.0, datasplit=c(1,0,0), k_rep=10, tss_rep=10, ACC_rep=10, REC_rep=10, F1_rep=10, valiMaxTech="entropy", dJ=1, kE.p=999, step=0.05, threshold=0.001, k1.train=T, k2.train=T){
  print("###########ASDM Starting##########")
  if (valuescale!=1000){
    map1 <- raster::calc(map1, fun=function(x){x*1000/valuescale})
    map2 <- raster::calc(map2, fun=function(x){x*1000/valuescale})
    map3 <- raster::calc(map3, fun=function(x){x*1000/valuescale})
  }
  w1 <- 1-(overlap[1,1]+overlap[2,1]-overlap[3,1])
  w2 <- 1-(overlap[1,1]+overlap[3,1]-overlap[2,1])
  w3 <- 1-(overlap[2,1]+overlap[3,1]-overlap[1,1])

  cellnum <- sum(which(as.matrix(map1)>=0))
  map.one <- raster::overlay(map1, fun=function(x){return(x*999.99/x)})
  map.zero <- raster::overlay(map1, fun=function(x){return(x/x-0.9999)})
  if (sum(datasplit) != 1){
    stop("The sum of datasplit should be 1.")
  }
  if (datasplit[3] == c(0)){
    k_rep <- 1
  }
  k1.list <- c()
  k2.list <- c()
  E.train.list <- c()
  E.vali.list <- c()
  E.test.list <- c()
  T.train.list <- c()
  T.vali.list <- c()
  T.test.list <- c()
  T.threshold.train.list <- c()
  print("###########Calculating K values##########")
  k1.initial <- k1
  k2.initial <- k2
  index <- sample(c(1:length(xylist[,1])), datasplit[3]*length(xylist[,1]))
  xylist.test <- xylist[index,]
  index <- setdiff(c(1:length(xylist[,1])), index)
  for (k_rep_i in c(1:k_rep)){
    print(paste("Replicate: ",k_rep_i,sep=""))
    index2 <- sample(c(1:length(index)), datasplit[1]/(datasplit[1]+datasplit[2])*length(index))
    xylist.train <- xylist[index[index2],]
    index2 <- setdiff(c(1:length(index)), index2)
    xylist.validation <- xylist[index[index2],]
    k1 <- k1.initial
    k2 <- k2.initial
    if (k1.initial=="random"){k1 <- runif(1, min=0.7,max=1.3)}
    if (k2.initial=="random"){k2 <- runif(1, min=0.7,max=1.3)}
    i <- 0
    dJ <- dJ
    kE.p <- kE.p
    step <- step
    threshold <- threshold
    k1.train <- k1.train
    k2.train <- k2.train
    repeat{
      i <- i+1
      map.overlay <- raster::overlay(map3, map1, map2, fun=function(x,y,z){return(x*k1/w1 - y*k2*w2/w1 - z*k2*w3/w1)})
      map.overlay <- raster::overlay(map.overlay, map.zero, fun="max")
      map.overlay <- raster::overlay(map.overlay, map.one, fun="min")
      xydata <- raster::extract(map.overlay, xylist.train)
      entropy <- unlist(lapply(xydata*0.001, function(x){return(-x*log(x,2)-(1-x)*log(1-x, 2))}))
      entropy <- na.omit(entropy)
      kE <- sum(entropy)/length(entropy)

      map.overlay <- raster::overlay(map3, map1, map2, fun=function(x,y,z){return(x*(k1+dJ*step)/w1 - y*k2*w2/w1 - z*k2*w3/w1)})
      map.overlay <- raster::overlay(map.overlay, map.zero, fun="max")
      map.overlay <- raster::overlay(map.overlay, map.one, fun="min")
      xydata <- raster::extract(map.overlay, xylist.train)
      entropy <- unlist(lapply(xydata*0.001, function(x){return(-x*log(x,2)-(1-x)*log(1-x, 2))}))
      entropy <- na.omit(entropy)
      dE1 <- sum(entropy)/length(entropy)

      map.overlay <- raster::overlay(map3, map1, map2, fun=function(x,y,z){return(x*k1/w1 - y*(k2+dJ*step)*w2/w1 - z*(k2+dJ*step)*w3/w1)})
      map.overlay <- raster::overlay(map.overlay, map.zero, fun="max")
      map.overlay <- raster::overlay(map.overlay, map.one, fun="min")
      xydata <- raster::extract(map.overlay, xylist.train)
      entropy <- unlist(lapply(xydata*0.001, function(x){return(-x*log(x,2)-(1-x)*log(1-x, 2))}))
      entropy <- na.omit(entropy)
      dE2 <- sum(entropy)/length(entropy)

      k1_out <- k1
      k2_out <- k2
      kE_out <- kE
      print(paste("k1 = ",k1_out,", k2 = ",k2_out,", entropy = ",kE_out, sep=""))
      if (i > 30 | abs(kE.p-kE) < threshold){
        break
      }
      if (k1.train){
        dJ1 <- (dE1-kE)/step
        k1 <- k1 + step*dJ1
      }
      if (k2.train){
        dJ2 <- (dE2-kE)/step
        k2 <- k2 + step*dJ2
      }
      kE.p <- kE
    }
    k1.list <- c(k1.list, k1_out)
    k2.list <- c(k2.list, k2_out)

    map.made <- raster::overlay(map3, map1, map2, fun=function(x,y,z){return(x*k1/w1 - y*k2*w2/w1 - z*k2*w3/w1)})
    map.made <- raster::overlay(map.made, map.zero, fun="max")
    map.made <- raster::overlay(map.made, map.one, fun="min")
    E.train.list <- c(E.train.list, calcu.map.entropy(map.made, xylist.train))
    T.train.out <- getRepCalcuMAXTSS(xylist.train, map.made, tss_rep)
    T.train.list <- c(T.train.list, T.train.out[["maxtss.mean"]])
    T.threshold.train.list <- c(T.threshold.train.list, T.train.out[["threshold.mean"]])

    output <- list()
    output[["map.made"]] <- NULL
    output[["k1"]] <- NULL
    output[["k2"]] <- NULL
    output[["entropy.train"]] <- NULL
    output[["entropy.train.list"]] <- NULL
    output[["maxtss.train"]] <- NULL
    output[["maxtss.train.list"]] <- NULL
    output[["maxtss.threshold.train.list"]] <- NULL
    output[["entropy.validation"]] <- NULL
    output[["entropy.validation.list"]] <- NULL
    output[["maxtss.validation"]] <- NULL
    output[["maxtss.validation.list"]] <- NULL
    output[["maxtss.threshold.validation.list"]] <- NULL
    output[["entropy.test"]] <- NULL
    output[["entropy.test.list"]] <- NULL
    output[["maxtss.test"]] <- NULL
    output[["maxtss.test.list"]] <- NULL
    output[["maxtss.threshold.test.list"]] <- NULL

    if (datasplit[3]==0){
      if (datasplit[2]==0){
        output[["map.made"]] <- map.made
        output[["k1"]] <- k1
        output[["k2"]] <- k2
        output[["entropy.train"]] <- E.train.list[1]
        output[["entropy.train.list"]] <- E.train.list
        output[["maxtss.train"]] <- T.train.list[1]
        output[["maxtss.train.list"]] <- T.train.list
        output[["maxtss.threshold"]] <- T.threshold.train.list[1]
        output[["calcu.map.meva"]] <- calcu.map.meva(xylist.train, map.made, output[["maxtss.threshold"]])
        output[["calcu.map.ACC"]] <- calcu.map.ACC(xylist.train, map.made, rep=10, output[["maxtss.threshold"]])
        output[["calcu.map.REC"]] <- calcu.map.REC(xylist.train, map.made, rep=10, output[["maxtss.threshold"]])
        output[["calcu.map.F1"]] <- calcu.map.F1(xylist.train, map.made, rep=10, output[["maxtss.threshold"]])
        print("###########Done##########")
        return(output)
      } else{
        E.vali.list <- c(E.vali.list, calcu.map.entropy(map.made, xylist.validation))
        T.vali.out <- getRepCalcuMAXTSS(xylist.validation, map.made, tss_rep)
        T.vali.list <- c(T.vali.list, T.vali.out[["maxtss.mean"]])
        T.threshold.vali.list <- c(T.vali.list, T.vali.out[["threshold.mean"]])
        output[["map.made"]] <- map.made
        output[["k1"]] <- k1
        output[["k2"]] <- k2
        output[["entropy.train"]] <- E.train.list[1]
        output[["entropy.train.list"]] <- E.train.list
        output[["maxtss.train"]] <- T.train.list[1]
        output[["maxtss.train.list"]] <- T.train.list
        output[["maxtss.threshold"]] <- T.threshold.train.list[1]
        output[["entropy.test"]] <- E.vali.list
        output[["maxtss.test"]] <- T.vali.list
        output[["calcu.map.meva"]] <- calcu.map.meva(xylist.validation, map.made, output[["maxtss.threshold"]])
        output[["calcu.map.ACC"]] <- calcu.map.ACC(xylist.validation, map.made, rep=ACC_rep, output[["maxtss.threshold"]])
        output[["calcu.map.REC"]] <- calcu.map.REC(xylist.validation, map.made, rep=REC_rep, output[["maxtss.threshold"]])
        output[["calcu.map.F1"]] <- calcu.map.F1(xylist.validation, map.made, rep=F1_rep, output[["maxtss.threshold"]])
        print("###########Done##########")
        return(output)
      }
    } else{
      E.vali.list <- c(E.vali.list, calcu.map.entropy(map.made, xylist.validation))
      T.vali.out <- getRepCalcuMAXTSS(xylist.validation, map.made, tss_rep)
      T.vali.list <- c(T.vali.list, T.vali.out[["maxtss.mean"]])
    }
  }
  E.test.list <- c(E.test.list, calcu.map.entropy(map.made, xylist.test))
  T.test.out <- getRepCalcuMAXTSS(xylist.test, map.made, tss_rep)
  T.test.list <- c(T.test.list, T.test.out[["maxtss.mean"]])
  if (valiMaxTech == "entropy"){bestk.index <- which.max(E.vali.list)} else {
    if (valiMaxTech == "maxtss"){bestk.index <- which.max(T.vali.list)} else {stop("Unsupported valiMaxTech !")}
  output[["k1"]] <- k1.list[bestk.index]
  output[["k2"]] <- k2.list[bestk.index]
  map.made <- raster::overlay(map3, map1, map2, fun=function(x,y,z){return(x*output[["k1"]]/w1 - y*output[["k2"]]*w2/w1 - z*output[["k2"]]*w3/w1)})
  map.made <- raster::overlay(map.made, map.zero, fun="max")
  map.made <- raster::overlay(map.made, map.one, fun="min")
  output[["map.made"]] <- map.made
  output[["entropy.train"]] <- E.train.list[bestk.index]
  output[["entropy.train.list"]] <- E.train.list
  output[["maxtss.train"]] <- T.train.list[bestk.index]
  output[["maxtss.train.list"]] <- T.train.list
  output[["maxtss.threshold"]] <- mean(T.threshold.train.list)
  output[["entropy.validation"]] <- E.vali.list[bestk.index]
  output[["entropy.validation.list"]] <- E.vali.list
  output[["maxtss.validation"]] <- T.vali.list[bestk.index]
  output[["maxtss.validation.list"]] <- T.vali.list
  output[["entropy.test"]] <- E.test.list
  output[["maxtss.test"]] <- T.test.list
  output[["calcu.map.meva"]] <- calcu.map.meva(xylist.test, map.made, output[["maxtss.threshold"]])
  output[["calcu.map.ACC"]] <- calcu.map.ACC(xylist.test, map.made, rep=10, output[["maxtss.threshold"]])
  output[["calcu.map.REC"]] <- calcu.map.REC(xylist.test, map.made, rep=10, output[["maxtss.threshold"]])
  output[["calcu.map.F1"]] <- calcu.map.F1(xylist.test, map.made, rep=10, output[["maxtss.threshold"]])
  print("###########Done##########")
  return(output)
}}



#' @title Calculate MaxTSS
#' @description
#' Calculate MaxTSS of model.
#' @param xylist Coordinates (dataframe), with two columns (X and Y)
#' @param map A Raster* class, the output of ASDM
#' @param rep Replications
#' @importFrom stats na.omit runif
#' @export
repCalcuMAXTSS<-function(xylist, map, rep=1){
  rep <- rep
  reslist <- list()
  for (i in c(1:rep)){
    samplesize <- nrow(xylist)
    sp.map <- raster::rasterize(cbind(xylist[,c(1,2)]), map, field=1)
    sp.map <- raster::overlay(sp.map, fun=function(x){return(!is.na(x))})
    map <- raster::overlay(sp.map, map, fun=function(x, y){return(y-x*y)})
    raster::NAvalue(map) <- 0
    # obtain psuedo-absent points randomly, sample size = 3*presenceSize
    Pred <- raster::sampleRandom(map, size=3*samplesize, na.rm=TRUE, ext=NULL,
                                 cells=FALSE, rowcol=FALSE, xy=FALSE, sp=FALSE, asRaster=FALSE)
    Sp.occ <- c(rep(0, 3*samplesize),rep(1, samplesize))
    Pred <- c(Pred, xylist[,3])
    Pred <- Pred/1000
    res <- ecospat::ecospat.max.tss(Pred, Sp.occ)
    reslist[[i]] <- list(max.TSS = res$max.TSS, max.threshold = res$max.threshold)
  }
  return(reslist)
}

#' @title Calculate R2
#' @description
#' Calculate R2 of model.
#' @param map1 Map for evaluation (ASDM output)
#' @param map2 Map for reference ('true' SDM output)
#' @importFrom stats na.omit runif
#' @export
calcu.map.R2 <- function(map1, map2){
  map.diff <- raster::overlay(map1, map2, fun=function(x,y){return((x-y)^2)})
  map.mean <- raster::cellStats(map1, "mean")
  map.diff2 <- raster::overlay(map1, fun=function(x){return((x-map.mean)^2)})
  R2 <- 1-raster::cellStats(map.diff, "sum")/raster::cellStats(map.diff2, "sum")
  print(R2)
}

#' @title Calculate entropy
#' @description
#' Calculate entropy of model.
#' @param xylist Coordinates of the target species (dataframe), with two columns (X and Y)
#' @param map The predicted distribution map of the target species (Raster* class)
#' @importFrom stats na.omit runif
#' @examples
#' # example code
#'
#' @export
calcu.map.entropy <- function(map, xylist){
  xydata <- raster::extract(map, xylist)
  entropy <- unlist(lapply(xydata*0.001, function(x){return(-x*log(x,2)-(1-x)*log(1-x, 2))}))
  entropy <- na.omit(entropy)
  entropy <- sum(entropy)/length(entropy)
  return(entropy)
}

#' @title Calculate MaxTSS (easy way)
#' @description
#' Calculate MaxTSS of model by a high-level approach.
#' @param xylist Coordinates of the target species, dataframe with 2 columns(X and Y)
#' @param map The predicted distribution map of the target species (Raster* class)
#' @param tss_rep Replication time for calculating MaxTSS
#' @importFrom stats na.omit runif
#' @examples
#' # example code
#'
#' @export
getRepCalcuMAXTSS <- function(xylist, map, tss_rep){
  if (length(xylist) == 0){return(NULL)}
  pred <- cbind(xylist, raster::extract(map, xylist))
  colnames(pred) <- c("x","y","asdm.pred")
  maxtss <- repCalcuMAXTSS(pred[,c(1,2,3)], map, rep=tss_rep)
  maxtss <- as.data.frame(maxtss)
  maxtss.mean <- apply(maxtss[,grep("max.TSS", colnames(maxtss))], 1, mean) # 0.93125
  threshold.mean <- apply(maxtss[,grep("max.threshold", colnames(maxtss))], 1, mean) # 0.352
  # threshold <- maxtss[,grep("max.threshold", colnames(maxtss))]
  # maxtss <- maxtss[,grep("max.TSS", colnames(maxtss))]
  # output <- list(pred=pred, maxtss=maxtss, threshold=threshold, maxtss.mean=maxtss.mean, threshold.mean=threshold.mean)
  output <- list(maxtss.mean=maxtss.mean, threshold.mean=threshold.mean)
  return(output)
}

#' @title Calculate RMSE
#' @description
#' Calculate RMSE of model.
#' @param map1 Map for evaluation (ASDM output)
#' @param map2 Map for reference ('true' SDM output)
#' @param cellnum Total number of grids.
#' @importFrom stats na.omit runif
#' @examples
#' # example code
#'
#' @export
calcu.map.RMSE <- function(map1, map2, cellnum){
  overlay <- raster::overlay(map1, map2, fun=function(x,y){return((x-y)^2)})
  RMSE <- (raster::cellStats(overlay, "sum")/cellnum)^0.5
  return(RMSE)
}

#' @title Calculate MAE
#' @description
#' Calculate MAE of model.
#' @param map1 Map for evaluation (ASDM output)
#' @param map2 Map for reference ('true' SDM output)
#' @param cellnum Total number of grids.
#' @importFrom stats na.omit runif
#' @examples
#' # example code
#'
#' @export
calcu.map.MAE <- function(map1, map2, cellnum){
  overlay <- raster::overlay(map1, map2, fun=function(x,y){return(abs(x-y))})
  MAE <- raster::cellStats(overlay, "sum")/cellnum
  return(MAE)
}


#' @title Calculate Model Evaluation for a Given Threshold Value
#' @description
#' Calculate Model Evaluation for a Given Threshold Value
#' @param xylist Coordinates (dataframe), with two columns (X and Y)
#' @param map A Raster* class, the output of ASDM
#' @param threshold Threshold of suitability, a float between 0 and 1
#' @importFrom stats na.omit runif
#' @examples
#' # example code
#'
#' @export
calcu.map.meva <- function(xylist, map, threshold){
  xylist <- cbind(xylist, raster::extract(map, xylist))
  samplesize <- nrow(xylist)
  sp.map <- raster::rasterize(cbind(xylist[,c(1,2)]), map, field=1)
  sp.map <- raster::overlay(sp.map, fun=function(x){return(!is.na(x))})
  map <- raster::overlay(sp.map, map, fun=function(x, y){return(y-x*y)})
  raster::NAvalue(map) <- 0
  # obtain psuedo-absent points randomly, sample size = 3*presenceSize
  Pred <- raster::sampleRandom(map, size=3*samplesize, na.rm=TRUE, ext=NULL,
                               cells=FALSE, rowcol=FALSE, xy=FALSE, sp=FALSE, asRaster=FALSE)
  Sp.occ <- c(rep(0, 3*samplesize),rep(1, samplesize))
  Pred <- c(Pred, xylist[,3])
  Pred <- Pred/1000
  res <- ecospat::ecospat.meva.table(Pred, Sp.occ, threshold)
  return(res)
}


#' @title Calculate Precision
#' @description
#' Calculate Precision of model.
#' @param xylist Coordinates (dataframe), with two columns (X and Y)
#' @param map A Raster* class, the output of ASDM
#' @param rep Replications
#' @param threshold Threshold of suitability
#' @importFrom stats na.omit runif
#' @examples
#' # example code
#'
#' @export
calcu.map.ACC <- function(xylist, map, rep=1, threshold){
  xylist <- cbind(xylist, raster::extract(map, xylist))
  samplesize <- nrow(xylist)
  reslist <- list()
  for (i in c(1:rep)){
    sp.map <- raster::rasterize(cbind(xylist[,c(1,2)]), map, field=1)
    sp.map <- raster::overlay(sp.map, fun=function(x){return(!is.na(x))})
    map <- raster::overlay(sp.map, map, fun=function(x, y){return(y-x*y)})
    raster::NAvalue(map) <- 0
    # obtain psuedo-absent points randomly, sample size = 3*presenceSize
    Pred <- raster::sampleRandom(map, size=3*samplesize, na.rm=TRUE, ext=NULL,
                               cells=FALSE, rowcol=FALSE, xy=FALSE, sp=FALSE, asRaster=FALSE)
    FP <- sum(Pred/1000>=threshold, na.rm=T)
    TP <- sum(xylist[,3]/1000>=threshold, na.rm=T)
    reslist[[i]] <- TP/(TP+FP)
  }

  return(unlist(reslist))
}

#' @title Calculate Recall Rate
#' @description
#' Calculate recall rate of model.
#' @param xylist Coordinates (dataframe), with two columns (X and Y)
#' @param map A Raster* class, the output of ASDM
#' @param rep Replications
#' @param threshold Threshold of suitability
#' @importFrom stats na.omit runif
#' @examples
#' # example code
#'
#' @export
calcu.map.REC <- function(xylist, map, rep=1, threshold){
  xylist <- cbind(xylist, raster::extract(map, xylist))
  samplesize <- nrow(xylist)
  reslist <- list()
  for (i in c(1:rep)){
    sp.map <- raster::rasterize(cbind(xylist[,c(1,2)]), map, field=1)
    sp.map <- raster::overlay(sp.map, fun=function(x){return(!is.na(x))})
    map <- raster::overlay(sp.map, map, fun=function(x, y){return(y-x*y)})
    raster::NAvalue(map) <- 0
    # obtain psuedo-absent points randomly, sample size = 3*presenceSize
    Pred <- raster::sampleRandom(map, size=3*samplesize, na.rm=TRUE, ext=NULL,
                                 cells=FALSE, rowcol=FALSE, xy=FALSE, sp=FALSE, asRaster=FALSE)
    FN <- sum(xylist[,3]/1000<threshold, na.rm=T)
    TP <- sum(xylist[,3]/1000>=threshold, na.rm=T)
    reslist[[i]] <- TP/(TP+FN)
  }

  return(unlist(reslist))
}

#' @title Calculate F1 Score
#' @description
#' Calculate F1 score of model.
#' @param xylist Coordinates (dataframe), with two columns (X and Y)
#' @param map A Raster* class, the output of ASDM
#' @param rep Replications
#' @param threshold Threshold of suitability
#' @importFrom stats na.omit runif
#' @examples
#' # example code
#'
#' @export
calcu.map.F1 <- function(xylist, map, rep=1, threshold){
  xylist <- cbind(xylist, raster::extract(map, xylist))
  samplesize <- nrow(xylist)
  reslist <- list()
  for (i in c(1:rep)){
    sp.map <- raster::rasterize(cbind(xylist[,c(1,2)]), map, field=1)
    sp.map <- raster::overlay(sp.map, fun=function(x){return(!is.na(x))})
    map <- raster::overlay(sp.map, map, fun=function(x, y){return(y-x*y)})
    raster::NAvalue(map) <- 0
    # obtain psuedo-absent points randomly, sample size = 3*presenceSize
    Pred <- raster::sampleRandom(map, size=3*samplesize, na.rm=TRUE, ext=NULL,
                                 cells=FALSE, rowcol=FALSE, xy=FALSE, sp=FALSE, asRaster=FALSE)
    FP <- sum(Pred/1000>=threshold, na.rm=T)
    FN <- sum(xylist[,3]/1000<threshold, na.rm=T)
    TP <- sum(xylist[,3]/1000>=threshold, na.rm=T)
    ACC <- TP/(TP+FP)
    REC <- TP/(TP+FN)
    reslist[[i]] <- 2*ACC*REC/(ACC+REC)
  }
  return(unlist(reslist))
}
