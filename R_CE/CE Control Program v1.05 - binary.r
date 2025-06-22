###  **********************************************************************;
###  **********************************************************************;
###  *****  CE Control Program.R                                      *****;
###  *****                                                            *****;
###  *****  Control program for R version of CE                       *****;
###  *****  01/31/2013  Created by Nan Flaaten                        *****;
###  *****       *                                                    *****;
###  *****  04/09/2013 v1.05 - added functionality to handle outliers *****;
###  *****                     for continuous dependent variable      *****;
###  **********************************************************************;
###  **********************************************************************;

## Import data from SAS
## If you need to convert your data from SAS to R
## Log in to hqqmgapp90 and run this code, substituting your project path and file names

##  require(lattice)
##  require(Revobase)
##  require(RevoScaleR)
##  test_data <- rxImport("/mnt/projects/public/projects/RTraining/Data/geico/test_data.sas7bdat")
##  save(test_data, file = "/mnt/projects/public/projects/RTraining/Data/geico/test_data.RData")


# install packages in R not in Rstudio
# install.packages('XLConnect')
# install.packages('subselect')
# install.packages('biglm')
# install.packages('corpcor')
# install.packages('lubridate')

require(XLConnect)
require(subselect)
require(biglm)
require(corpcor)
require(lubridate)

## Load functions
source("C:\\Darren\\ziqing\\Merkle\\Model\\CE\\RCE\\code\\CE Functions v1.05.R", echo = F)

## Read in data into dataset
inlib <- "C:\\Darren\\ziqing\\Merkle\\Model\\CE\\RCE\\data"
setwd(inlib)
load("test_data.RData") 

View(head(test_data))

dat0 <- test_data[test_data$model_seg !='outval',]

dat0$wgt<-1

dat0$lrnval <- ifelse(dat0$model_seg=='inval','V','L')

table(dat0$lrnval)
# L     V 
# 54403 23607 
table(dat0$dv_any_purchase_flag_in_attr,dat0$lrnval)
# L     V
# 0 37926 16430
# 1 16477  7177

## Specified Parameters
lab.targ<-'dv_any_purchase_flag_in_attr'   ## Dependent variable
lab.split<- 'lrnval'        ## Value is lrnval if using sampling function. Set to NULL if not using sampling and no split variable is defined.
lab.id<- 'rid'              ## ID variable
lab.wgt<- 'wgt'             ## Value is wgt if using sampling function. Set to NULL if not using sampling and no weight variable is defined.
MissPct <- .9               ## Maximum allowed missing percent
MaxBinPct <- .9             ## Maximum allowed in a single bin (character)
MaxCorr <- .75              ## Maximum correlation allowed between variables
talpha <- .05               ## p(t) to start new bin
samp.flg <- 'Y'             ## Set to Y to do sampling, N to skip
samp.size <- 20000          ## If sampling, number of records in sample file
samp.strat <- 'N'           ## If sampling, should it be done stratified
samp.mult <- 5              ## If sampling and stratified, ratio of 0's to 1's (binary) or above cut to below (continuous)
samp.lv <- .7               ## Proportion of sample in learning 
samp.cut <- 0               ## Value for cutpoint if stratified and continuous
seed1 <- 123456789          ## Seed to use for sampling
seed2 <- 987654321          ## Seed to use for learning / validation split
out <- 'Y'                  ## Indicator to handle outliers. Set to Y to remove, R to recode to out.max, N to do nothing
out.max <- 99               ## Percentile definition of outliers. Anything over this percentile will be removed or recoded down to this level depending on the value of out
dv.type <- 'C'              ## Flag to specify ordinal dependent. Set to 'O' for ordinal and 'C' for continuous.
dv.vals <- c(1,2,3,4,5)     ## If ordinal and stratified, what are the values for the dependent variable
dv.cnts <- c(5000,5000,5000,5000,5000)  ##If ordinal and stratified, how many observations should be selected for each value of the dependent variable
LowerPct <- .05             ## Lower threshold for determining usable numeric variables
UpperPct <- .95             ## Upper threshold for determining usable numeric variables
                            ## If the value at LowerPct = value at UpperPct then the variable is unusable
prof.flg <- 'N'             ## Set to Y to do pre variable reduction profiling, N to skip
prof.breaks <- 5            ## Number of bins to use in profiling
mod.method <- 'backward'     ## Model selection method - forward or backward
mod.maxsteps <- 30          ## Maximum number of variables for forward selection
mod.sls <- .05              ## p value for entry or remove

# Initiate Excel workbook
if (file.exists('EDA.xlsx')) file.remove('EDA.xlsx')
ff <- loadWorkbook('EDA.xlsx', create = TRUE)

# Create formats
csPercent = createCellStyle(ff, name = "pct")
setDataFormat(csPercent, format = "0.00%")

csPercent1 = createCellStyle(ff, name = "pct1")
setDataFormat(csPercent1, format = "0.0%")

csComma = createCellStyle(ff, name = "comma")
setDataFormat(csComma, format = "#,##0")

csDec4 = createCellStyle(ff, name = "dec4")
setDataFormat(csDec4, format = "0.0000")

csDec2 = createCellStyle(ff, name = "dec2")
setDataFormat(csDec2, format = "0.00")

csLine = createCellStyle(ff, name = "line")
setBorder(csLine, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")

cslComma = createCellStyle(ff, name = "lcomma")
setDataFormat(cslComma, format = "#,##0")
setBorder(cslComma, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")

cslPercent1 = createCellStyle(ff, name = "lpct1")
setDataFormat(cslPercent1, format = "0.0%")
setBorder(cslPercent1, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")

cslDec2 = createCellStyle(ff, name = "ldec2")
setDataFormat(cslDec2, format = "#,##0.00")
setBorder(cslDec2, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")

## Call sampling.
## Parameters:  x - dataset to start with
if (samp.flg == 'Y') dat2 <- f.sample(dat0,lab.targ,samp.size,samp.strat,samp.mult,
                                      samp.lv,samp.cut,seed1,seed2,out,out.max)else{
  dat2 <- f.nosample(dat0,lab.targ,samp.lv,seed2)
                                      }

table(dat2$dv_any_purchase_flag_in_attr,dat2$lrnval)

## Remove unnecessary variables
dat2 <- dat2[,names(dat2)[pmax(grepl("name",names(dat2),ignore.case=T),
                               grepl("addr",names(dat2),ignore.case=T),
                               grepl("city",names(dat2),ignore.case=T),
                               grepl("zip",names(dat2),ignore.case=T),
                               grepl("title",names(dat2),ignore.case=T),
                               grepl("nm3",names(dat2),ignore.case=T),
                               grepl("nm4",names(dat2),ignore.case=T),
                               grepl("date",names(dat2),ignore.case=T),
                               grepl("time",names(dat2),ignore.case=T),
                               grepl("latitude",names(dat2),ignore.case=T),
                               grepl("longitude",names(dat2),ignore.case=T),
                               grepl("suffix",names(dat2),ignore.case=T),
                               grepl("_dt",names(dat2),ignore.case=T),
                               grepl("msc",names(dat2),ignore.case=T),
                               grepl("hispfn",names(dat2),ignore.case=T),
                               grepl("hispln",names(dat2),ignore.case=T),
                               grepl("bm019",names(dat2),ignore.case=T),
                               grepl("match",names(dat2),ignore.case=T)
                               # ,
                               # grepl("id",names(dat2),ignore.case=T)
                               ,grepl("dv_dpa_purchase_flag_in_attr",names(dat2),ignore.case=T)
                               ,grepl("cohort_year",names(dat2),ignore.case=T)
                               ,grepl("cohort_month",names(dat2),ignore.case=T)
                               ,grepl("dv_revenue_in_attr",names(dat2),ignore.case=T)
                               ,grepl("dv_margin_in_attr",names(dat2),ignore.case=T)
                               
                               )==FALSE
                          ]]

#####  Apply transform functions  #####

## DS recodes
dat3 <- f.DSVars(dat2)

## Numeric transforms
tmp <- f.numvars(dat3,lab.id,lab.targ,lab.split,lab.wgt,MissPct,LowerPct,UpperPct)
dat4 <- as.data.frame(tmp[1])
dat5 <- as.data.frame(tmp[2])
rm(tmp)

## Character transforms
dat6 <- f.charvars(dat3,lab.id,lab.targ,lab.split,lab.wgt,MissPct,MaxBinPct)

## Combine
if (length(dat6)==1 && is.na(dat6)) { dat7 <- cbind(dat2[,lab.targ],dat5) 
                   names(dat7) <- c(lab.targ,names(dat5))
} else { dat7 <- cbind(dat2[,lab.targ],dat5,dat6)
         names(dat7) <- c(lab.targ,names(dat5),names(dat6))
}

## Profiling
if (prof.flg == 'Y') profile <- f.profile(dat7,lab.targ,setdiff(names(dat7),lab.targ),prof.breaks,ff)

# Save out Excel workbook
saveWorkbook(ff)

## Variable reduction
keep <- f.var_redu(ds=dat7,dv=lab.targ,MaxCorr,tol=.000001)
dat8 <- cbind(dat7[,c(lab.targ,keep)],dat2[,c(lab.id,lab.split,lab.wgt)])


## Model creation
# Initiate Excel workbook
if (file.exists('Model Performance.xlsx')) file.remove('Model Performance.xlsx')
fm <- loadWorkbook('Model Performance.xlsx', create = TRUE)

# Create formats
csPercent1 = createCellStyle(fm, name = "pct1")
setDataFormat(csPercent1, format = "0.0%")

csPercent0 = createCellStyle(fm, name = "pct0")
setDataFormat(csPercent0, format = "0%")

csComma = createCellStyle(fm, name = "comma")
setDataFormat(csComma, format = "#,##0")

csDec2 = createCellStyle(fm, name = "dec2")
setDataFormat(csDec2, format = "#,##0.00")

csDec4 = createCellStyle(fm, name = "dec4")
setDataFormat(csDec4, format = "#,##0.0000")

csDec6 = createCellStyle(fm, name = "dec6")
setDataFormat(csDec6, format = "#,##0.000000")

csLine = createCellStyle(fm, name = "line")
setBorder(csLine, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")

cslComma = createCellStyle(fm, name = "lcomma")
setDataFormat(cslComma, format = "#,##0")
setBorder(cslComma, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")

cslPercent1 = createCellStyle(fm, name = "lpct1")
setDataFormat(cslPercent1, format = "0.0%")
setBorder(cslPercent1, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")

cslDec2 = createCellStyle(fm, name = "ldec2")
setDataFormat(cslDec2, format = "#,##0.00")
setBorder(cslDec2, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")

##Parameters: modtyp - type of model to built (forward or backward) [default='forward']
##            maxsteps - maximum number of variables for forward selection [default=20]
##            sls - p value for entry or remove [defaul=.05]
finmod <- f.model(dat8,dat3,dat4,lab.targ,lab.id,lab.split,lab.wgt,mod.method,mod.sls,mod.maxsteps,names(dat6))

saveWorkbook(fm)

## Save model sample
save(dat2, file = "C:\\Darren\\ziqing\\Merkle\\Model\\CE\\RCE\\output\\model_sample.RData")
