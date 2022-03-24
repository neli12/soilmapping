####################Cubist in R#################################################
##Load required packages####
library(caret)
library(raster)

# goof function from ithir package ----
goof <- function(observed,predicted, plot.it = FALSE, type="DSM"){
  # Coefficient of determination
  rLM <- lm(predicted ~ observed)
  R2 <- as.matrix(summary(rLM)$adj.r.squared)
  
  # Standard error of prediction ^2
  SEP2 <- mean((observed - predicted)^2)
  
  # Standard error of prediction
  SEP <- sqrt(SEP2)
  
  #Bias
  bias <- mean(predicted) - mean(observed)
  
  # residual  variance
  SEP2c <- sum(((predicted - bias - observed)^2) / length(observed))
  SEPc <- sqrt(SEP2c)
  
  # ratio of performance to deviation
  RPD <- sd(observed) / SEP
  
  # Ratio of performance to interquartile distance
  IQ <- c(quantile(observed))[3] - c(quantile(observed))[2]
  RPIQ <- IQ / SEP
  
  # Concordance
  mx <- mean(observed)
  my <- mean(predicted)
  s2x <- var(observed)
  s2y <- var(predicted)
  sxy <- mean((observed-mx) * (predicted-my))
  ccc <- 2 * sxy / (s2x + s2y + (mx - my)^2)
  
  if (plot.it==TRUE){eqscplot(observed, predicted)
    abline(a = 0, b = 1, col = "brown4")}
  
  if (type == "DSM"){ gf <- data.frame(R2=R2, concordance=ccc, MSE=SEP2, RMSE=SEP, bias=bias, row.names=NULL)}
  else if (type == "spec"){ gf <- data.frame(R2=R2, concordance=ccc, MSE=SEP2, RMSE=SEP, bias=bias, 
                                             MSEc=SEP2c,RMSEc=SEPc, RPD=RPD, RPIQ=RPIQ, row.names=NULL)}
  else {stop("ERROR: Revise the type of output you require. Select from either DSM or spec")} 
  
  return(gf)
}

###Set working directory###
setwd("C:/Users/FREY/Google Drive/GEOCIS/IC/Merylin/dados_2nd_revision")

###Load raster covariates, in our case, satellite image###

##SYSI Combined##
SYSI.Combined <- stack("SYSIs/SYSI_Combined.tif") ##before the / is the folder name
SYSI.Combined[SYSI.Combined == 0] <- NA  ##set 0 values as NA
SYSI.Combined[SYSI.Combined < 0] <- NA  ##set 0 values as NA
names(SYSI.Combined) <- c("SYSI.Combined1", "SYSI.Combined2","SYSI.Combined3", "SYSI.Combined4", "SYSI.Combined5", "SYSI.Combined6")


##SYSI L8 OLI##
SYSI.L8OLI <- stack("SYSIs/SYSI_L8-OLI.tif") ##before the / is the folder name
SYSI.L8OLI[SYSI.L8OLI == 0] <- NA  ##set 0 values as NA
SYSI.L8OLI[SYSI.L8OLI < 0] <- NA  ##set 0 values as NA
names(SYSI.L8OLI) <- c("SYSI.L8OLI1", "SYSI.L8OLI2","SYSI.L8OLI3", "SYSI.L8OLI4", "SYSI.L8OLI5", "SYSI.L8OLI6") ##renaming covariates

##SYSI S2-MSI##
SYSI.S2MSI <- stack("SYSIs/SYSI_S2-MSI.tif") ##before the / is the folder name
SYSI.S2MSI[SYSI.S2MSI == 0] <- NA  ##set 0 values as NA
SYSI.S2MSI[SYSI.S2MSI < 0] <- NA  ##set 0 values as NA
names(SYSI.S2MSI) <- c("SYSI.S2MSI1", "SYSI.S2MSI2","SYSI.S2MSI3", "SYSI.S2MSI4", "SYSI.S2MSI5", "SYSI.S2MSI6", "SYSI.S2MSI7", "SYSI.S2MSI8", "SYSI.S2MSI9") ##renaming covariates


###Read dataset###
dat1 <- read.csv("tables/dados.csv", sep=";")
summary(dat1)   #summary of the dataset

###Extract covariates values###
coordinates(dat1) <- ~ x+y  #set your data as spatial object
plot(dat1)                  #plot your data

SYSI.Combined.df <- as.data.frame(extract(SYSI.Combined, dat1))
SYSI.L8OLI.df <- as.data.frame(extract(SYSI.L8OLI, dat1))
SYSI.S2MSI.df <- as.data.frame(extract(SYSI.S2MSI, dat1))


##Join covariate values in a data frame###
SYSI.covs <- cbind(SYSI.Combined.df, SYSI.L8OLI.df, SYSI.L8OLI.df)


##Select 30% of the data
val_rows <- sample(nrow(dat1), 685)


###Split dataset and covariates into training and validation sets###
dat1_train <- dat1[-val_rows, ]   #dat1 corresponds to the data loaded in R (all of it)
covariates_train <- dat1[-val_rows, ]  #The SYSI.cov data contains six columns corresponding to the SYSI vis-NIR-SWIR bands

dat1_test <- dat1[val_rows, ]
covariates_test <- SYSI.covs[val_rows, ]


###Fit a cubist model to your data  ---  Example for clay content and SYSI Combined --- Repit for the remaining soil properties and covariates###
##The hyperparameters committees and neighbors were used by default. No hyperparameters were tunned##
clay.fit <- train(x = covariates_train[,1:6], y = dat1_train$Clay, method = "cubist", trControl = trainControl(method = "cv"))
summary(clay.fit)   ##to see the models and rules created
clay.fit            ##to see the performance of the models

###Predictions validation###
predicted.clay.test <- predict(clay.fit, covariates_test[,1:6])  ##validation


##Performance of validation set
goof(dat1_test$clay, predicted.clay.test, type="spec")


###Predict to the entire study area using the raster covariates###
clay.raster <- raster::predict(SYSI.Combined, clay.fit)

crs(clay.raster) <- "+init=epsg:4326"   ##Set coordinates reference system to WGS 1984(degrees)
plot(clay.raster)


##Save the predicted map##
writeRaster(clay.raster, "raster/clay_content.tif", overwrite=TRUE)



