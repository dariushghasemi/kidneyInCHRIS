#=========================================#
#          Data Transformation
#=========================================#


#-----------------------------------------------------#
#---------------- My standardization -----------------
#-----------------------------------------------------#

#First way of controlling Instruments' effect for UACR

#Mixed regression
#install.packages("lmer")
as.factor(chris$Sex)

#Serum Creatinine & eGFR controlling Instruments 
#------------------------------------------------------#
std <- function(df, y, x) {
  df <- chris
  #df <- df #
  #attach(df)
  x <- as.name(deparse(substitute(x))) #quoting and enquoting 
  y <- as.name(deparse(substitute(y)))
  #m <- list()
  formula = cat(paste(y, sep = "~", x)) #or print(varname, quote=FALSE)
  #data = paste(df)
  m <- lm(formula, df)
  #m <- lm(SerumCreatinine ~ SerumCreatinine.Ins, data = chris)
  myresidm <- data.frame(resm = m$residuals)
  #var <- character()
  #var <- append(var, paste(names(chris)[which(colnames(chris)==y)], "Res", sep = "."))#"chris$",
  var <- paste(names(df)[which(colnames(df) == y)], "Res", sep = ".")
  var <- paste(var)
  df[,var] <- NA
  rn <- as.numeric(rownames(myresidm))
  #df[rn, var] <- myresidm$resm
  standres <- data.frame(rn, myresidm$resm, check.rows = TRUE)#paste(var) = 
  colnames(standres) <- c("id", paste(var))
  #assign(paste(df), myresidm$resm, envir=.GlobalEnv)
  #return(head(standres))
  return(m)#return(df[,var])
}
#------------------------------------------------------#
std(chris, "SerumCreatinine", "SerumCreatinine.Unit")
std(chris, SerumCreatinine, SerumCreatinine.Ins)

#var<-array()
var <- character()
append(var, "res")

#get an index of a variable
names(chris)[3]
which(colnames(df)=="B")
grep("B", colnames(df))
#------------------------------------------------------#
standardization <- 
  function(df, trait, instrument){
    #trait <- df$SerumCreatinine
    #week <- df$week
    #instrument <- df$SerumCreatinine.Ins #noquote(paste0("trait",".Ins"))
    library(lme4)
    mSCr.Std <- lmer(trait ~ instrument + (1|week), data = df, na.action = "na.exclude")
    #noquote(paste0("df","$","trait",".Std")) <- mean(df$trait, na.rm = TRUE) + resid(mSCr.Std)
    trait.Std <- mean(trait, na.rm = TRUE) + resid(mSCr.Std)
    #df2 <- cbind(df, SerumCreatinine.Std)
    return(name())
  }
#------------------------------------------------------#
#Standardizing for change of the instrument and week of participation
chris$SerumCreatinine.Std   <- standardization(chris, chris$SerumCreatinine,   chris$SerumCreatinine.Ins)
#noquote(paste("df","$","trait"))



#-----------------------------------------------------#
#---------------- Min-Max Normalization --------------
#-----------------------------------------------------#

# Min-Max standardization generally and for instrument

minMaxNorm <- function(x){#data, trait, instrument
  #TSH0 = data[instrument == 0, "TSH"]
  #TSH1 = data[instrument == 1, "TSH"]
  #min0 = min(TSH0, na.rm = TRUE)
  #min1 = min(TSH1, na.rm = TRUE)
  #max0 = max(TSH0, na.rm = TRUE)
  #max1 = max(TSH1, na.rm = TRUE)
  #Std  = ifelse(instrument == 0, (trait - min0)/ (max0 - min0), (trait - min1)/ (max1 - min1))
  Std = (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  return(Std)
}

summary(minMaxNorm(chris, chris$TSH, chris$TSH.Ins))

ifelse(chris$TSH.Ins == 0, chris$TSH - min(chris$TSH) / (max0 - min0), trait - min1 / (max1 - min1))
ifelse(chris$TSH.Ins == 0, mean(chris[chris$TSH.Ins == 0, "TSH"], na.rm = TRUE), mean(chris[chris$TSH.Ins == 1, "TSH"], na.rm = TRUE))


#-----------------------------------------------------#
#--------------- Quantile Normalization --------------
#-----------------------------------------------------#

#Quantile normalization

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("preprocessCore")
library(preprocessCore)

normalize.quantiles <- function(x,copy=TRUE){
  
  rows <- dim(x)[1]
  cols <- dim(x)[2]
  
  if (!is.matrix(x)){
    stop("Matrix expected in normalize.quantiles")
  }
  
  if (is.integer(x)){
    x <- matrix(as.double(x),rows,cols)
    copy <- FALSE
  }
  
  #matrix(.C("qnorm_c", as.double(as.vector(x)), as.integer(rows), as.integer(cols))[[1]], rows, cols)
  
  ##  .Call("R_qnorm_c",x,copy, PACKAGE="preprocessCore");
  .Call("R_qnorm_c_handleNA",x,copy, PACKAGE="preprocessCore");
}

#-----------------------------------------------------#
normalize.quantiles.robust <- 
  function(x,copy=TRUE,weights=NULL,
           remove.extreme=c("variance","mean","both","none"),
           n.remove=1,use.median=FALSE,use.log2=FALSE){
    
    calc.var.ratios <- function(x){
      cols <- dim(x)[2]
      vars <- apply(x,2,var)
      results <- matrix(0,cols,cols)
      for (i in 1:cols-1)
        for (j in (i+1):cols){
          results[i,j] <- vars[i]/vars[j]
          results[j,i] <- vars[j]/vars[i]
        }
      results
    }
    
    calc.mean.dists <- function(x){
      cols <- dim(x)[2]
      means <- colMeans(x)
      results <- matrix(0,cols,cols)
      for (i in 1:cols-1)
        for (j in (i+1):cols){
          results[i,j] <- means[i] - means[j]
          results[j,i] <- means[j] - means[i]
        }
      results
    }
    
    use.huber <- FALSE
    remove.extreme <- match.arg(remove.extreme)
    
    rows <- dim(x)[1]
    cols <- dim(x)[2]
    
    if (is.null(weights)){
      weights <- .Call("R_qnorm_robust_weights",x,remove.extreme,as.integer(n.remove),PACKAGE="preprocessCore")
    } else {
      if (is.numeric(weights)){
        if (length(weights) != cols){
          stop("Weights vector incorrect length\n")
        }
        if (sum(weights > 0) < 1){
          stop("Need at least one non negative weights\n")
        }
        if (any(weights < 0)){
          stop("Can't have negative weights")
        }
      } else {
        if (weights =="huber"){
          use.huber <- TRUE
          weights <- rep(1,cols)
        } else {
          stop("Don't recognise weights argument as valid.")
        }
      }
    }
    
    .Call("R_qnorm_robust_c",x,copy,weights,as.integer(use.median),as.integer(use.log2),as.integer(use.huber),PACKAGE="preprocessCore")
  }

normalize.quantiles(as.matrix(chris$TSH))
#-----------------------------------------------------#

#Adopted from caret package
"normalize2Reference" <- function (data, refData = NULL, ties = TRUE) 
{
  ## adapted from limma's normalizeQuantiles
  if (is.null(dim(data)))
  {
    numSamples <- 1
    numProbes <- length(data)
    data <- as.matrix(data)
  } else {
    numSamples <- dim(data)[2]
    numProbes <- dim(data)[1]   
  }
  
  dataOrder <- sortedData <- array(, dim(data))  
  if (ties) dataRanks <- dataOrder           
  
  nobs <- rep(numProbes, numSamples)
  
  quantProbs <- (0:(numProbes - 1))/(numProbes - 1)  
  
  for (j in 1:numSamples) 
  {
    ## sort each column of data and get indicies and ranks
    sortedSampleI <- sort(data[, j],
                          method = "quick",
                          index.return = TRUE)
    if (ties) dataRanks[, j] <- rank(data[, j])
    nobsj <- length(sortedSampleI$x)
    if (nobsj < numProbes) {
      nobs[j] <- nobsj
      isna <- is.na(data[, j])
      sortedData[, j] <- approx((0:(nobsj - 1))/(nobsj - 1),
                                sortedSampleI$x, 
                                quantProbs,
                                ties = "ordered")$y
      dataOrder[!isna, j] <- ((1:numProbes)[!isna])[sortedSampleI$ix]
    }
    else {
      sortedData[, j] <- sortedSampleI$x
      dataOrder[, j] <- sortedSampleI$ix
    }
  }
  ## at this point, we have the data sorted within each column and a key in dataOrder
  ## to revert back to the original input data: data[,j] == > sortedData[dataOrder[,j],j]
  
  ## use the row means from the data quantiles when no explicit reference is given.    
  if(is.null(refData))
  {
    refQuantile <- rowMeans(sortedData)
  } else refQuantile <- sort(refData)
  
  for (j in 1:numSamples) {
    if (nobs[j] < numProbes)
    {
      isna <- is.na(data[, j])
      if (ties) 
        ## The approx does a linear approx on the line where 
        ## x is the index scaled to [0, 1] and y is the mean
        ## for that index across all samples (or reference quants)
        ## It returns an approximated mean value 
        data[!isna, j] <- approx(quantProbs,
                                 refQuantile, 
                                 (dataRanks[!isna, j] - 1)/(nobs[j] - 1),
                                 ties = "ordered")$y
      else data[dataOrder[!isna, j], j] <- approx(quantProbs,
                                                  refQuantile, 
                                                  (0:(nobs[j] - 1))/(nobs[j] - 1),
                                                  ties = "ordered")$y
    }
    else {
      if (ties) 
        data[, j] <- approx(quantProbs,
                            refQuantile, 
                            (dataRanks[, j] - 1)/(numProbes - 1),
                            ties = "ordered")$y
      else data[dataOrder[, j], j] <- refQuantile
    }
  }
  data
}
#-----------------------------------------------------#

#Cristian's alt for quantile norm

quantile(chris[chris$APTT.Ins == 1, "APTT"], 
         probs = seq(0,1, length.out = length(chris[chris$APTT.Ins == 0, "APTT"])), 
         na.rm = T, names = T, type = 7, digits = 7)

#-----------------------------------------------------#
#Previous used approach for transforming only TSH
#-----------------------------------------------------#

library(caret)

nq <- normalize2Reference(
  chris[chris$TSH.Ins == 0, "TSH"], 
  refData = quantile(chris[chris$TSH.Ins == 1, "TSH"], 
                     probs = seq(0,1, length.out=length(chris[chris$TSH.Ins == 0, "TSH"])), 
                     na.rm = T, names = T, type = 7, digits = 7), ties = T)

#TSH quantile normalized
chris$TSH.q <- chris$TSH
chris[chris$TSH.Ins == 0, "TSH.q"] <- nq
chris$TSH.Std <- chris$TSH.q

summary(cbind(nq, chris[chris$TSH.Ins == 0, "TSH"]))
#-----------------------------------------------------#

#categorized TSH
chris <- 
  chris %>% 
  mutate(TSH_cat = cut(TSH,
                       breaks = c(-Inf, 0.401, 3.799, Inf), 
                       #labels = c("0", "1", "2")
                       labels = c("HyperT", "NormT", "HypoT")))

chris$TSH_cat <- as.character(chris$TSH_cat)
#vcfReg$TSH_cat <- relevel(vcfReg$TSH_cat, ref = 2)

#-----------------------------------------------------#
# Our quantile function iterating normalize2Reference
#-----------------------------------------------------#

Giulia_fun <- function(outcome, instrument) {
  
  subdata_1 <- outcome[instrument==1]
  subdata_0 <- outcome[instrument==0]
  refData   <- quantile(subdata_1,
                        probs = seq(0,1, length.out = length(subdata_0)),
                        na.rm = TRUE, names = TRUE, type = 7, digits = 7)
  
  nq <- normalize2Reference(data = subdata_0, refData = refData, ties = TRUE)
  outcome_std <- outcome
  outcome_std[instrument == 0] <- nq
  return(outcome_std)
}

# Iterating on all quantitative variables

#chris[c("UAlb.Std", "UCr.Std", "SAlb.Std")] <-
map2_df(chris[c("UAlb", "UCr", "SAlb")],
          chris[paste0(c("UAlb", "UCr", "SAlb"), ".Ins")],
          Giulia_fun)

# Quantitative variables
cbind(colnames(chris[quantVars[-c(1:8)]]),
      names(chris[paste0(quantVars[-c(1:8)], ".Ins")]))

varNOTforQuant <- c("Height", "Weight", "BMI", "Body_Fat", "Visceral_Fat",
                    "SBP", "DBP", "Pulse_Rate", "UACR", "APTT")

# saving the quantile-transformed variables inplace
chris[setdiff(quantVars, varNOTforQuant)] <-
  map2_df(chris[setdiff(quantVars, varNOTforQuant)],
          chris[paste0(setdiff(quantVars, varNOTforQuant), ".Ins")],
          Giulia_fun) #%>% View
# rename_with(~paste0(.x, ".qstd")) %>% head()  #View()

#-----------------------------------------------------#

#Checking N. traits with change of the instrument
do.call(rbind, sapply(chris[qualVars], table))
#---------#

# Statistics for the quantiled normalized traits
tbl1 <- 
  map2_df(chris[quantVars[c(12,13,15,16,17,18,19,21:80)]],
          chris[paste0(quantVars[c(12, 13, 15,16,17,18,19,21:80)], ".Ins")],
          Giulia_fun) %>%
  rename_with(~paste0(.x, "_std")) %>% 
  summarise_all(c(mean, sd), na.rm=TRUE) %>%
  round(digits = 2) %>%
  t() #%>% #View()

write.csv(tbl1, "mean&SD.csv")
#---------#

