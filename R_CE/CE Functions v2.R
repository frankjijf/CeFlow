# 1. SAMPLING ###################################

## **********************************************
## *****            1.Sampling              *****
## **********************************************

## *** Macro 1: Sampling macro variables ***;
## %let split_if = if lrnval = 'L';      ** your split condition for modeling portion;
## %let include_if=expression(!is.na(resp));     ** Data inclusion if condition applies. e.g. IF &Dep_var != .;
## %let DS_Present= Y;                   ** change to Y only if your data contain standard DS modeling bundle;
## %let path_DS=/mnt/projects/shared/pst_qmgisi/Modeling/CE/;  ** Location of standard DS recodes. Do Not Change;
## %let bootstrap= N;                    ** request oversamping, bootstrap if responders size is small (Y/N). Binary dv only;
## %let oversampled_rr=0.05;             ** your desired oversample response rate on the modeling part of data, 0.1-0.5;
## %let min_num_resp=500;                ** Minimum number of responders if sampling;
## %let seed=123456;                     ** Random selection seed if sampling;

f.sample <- function(dat, 
                     outds = "CE1_Resampled.xdf",
                     dep_var = "resp",
                     binary_dv = "Y",
                     include_if = NULL,
                     split_if = runif(n,0,1)<.7,
                     ds_present = "Y",
                     bootstrap = "N",
                     oversampled_rr = .1,
                     min_num_resp = 500,
                     seed = 123456
                     ){

  # Check dependent variable classification
  n <- names(rxGetVarInfo(dat))
  if (!dep_var %in% n) {
    stop(paste("Could not find dependent variable",dep_var))
  }
  u <- nrow(table(rxReadXdf(dat,varsToKeep=dep_var,reportProgress=0)))
  if (u == 1){ 
    stop("Dependent variable has only 1 value")
  }
  if (toupper(binary_dv)=="Y" & u > 2){
    stop("Dependent variable is classified as binary but there are more than 2 values")
  }
  if (toupper(binary_dv)!="Y" & u == 2){
    stop("Dependent variable is not classified as binary but there are only 2 values")
  }
  # Functions for splitting the data
  if (is.expression(split_if)) {
    f.split <- function(lst){
      lst$mod_val_test <- ifelse(eval(SplitIf, envir=lst), 
                                 ifelse(runif(length(lst[[1]]),0,1)<.5,1,2),
                                 3)
      return(lst)
    }
  } else {
    if (is.numeric(split_if)) {
      f.split <- function(lst){
        lst$mod_val_test <- ifelse(runif(length(lst[[1]])) < SplitIf,
                                   ifelse(runif(length(lst[[1]]),0,1)<.5,1,2),
                                   3)
        return(lst)
      }
    } else {
      cat("Invalid split condition")
      stop
    }
  }
  # Base splitting algorithm
  rxDataStep(inData = dat,
             outFile = outds,
             overwrite = T,
             xdfCompressionLevel = 1,
             rowSelection = include_if,
             transformFunc = f.split,
             transformObjects = list(SplitIf=split_if),
             reportProgress = 0)
  
  # Sample from full file
  if (bootstrap == "Y" & binary_dv == "Y") {
    cat("Perform sampling \n")
    dv <- table(rxDataStep(inData=outds,
                           rowSelection = mod_val_test != 3,
                           varsToKeep = dep_var,
                           reportProgress = 0),useNA="ifany")
    num_resp <- dv["1"]
    num_nresp <- dv["0"]
    if (num_resp >= min_num_resp) {
      spsize <- floor(num_resp / oversampled_rr)
      smp_resp <- num_resp
    } else {
      spsize <- floor(min_num_resp / oversampled_rr)
      smp_resp <- min_num_resp
    }
    smp_nresp <- floor((smp_resp/oversampled_rr)-smp_resp)
    cat("Original training response count: ", num_resp,'\n')
    cat("Original training non-response count: ", num_nresp,'\n')
    cat("New modeling sample size: ", spsize,'\n')
    cat("Sample training response count: ", smp_resp,'\n')
    cat("Sample training non-response count: ", smp_nresp,'\n')
    
    if (num_resp != smp_resp | num_nresp != smp_nresp) {
      rxDataStep(inData = outds,
                 outFile = "CE1_temp.xdf",
                 overwrite = T,
                 transforms = list(rowNum = .rxStartRow : (.rxStartRow + .rxNumRows - 1)),
                 reportProgress = 0)
      tmp <- rxReadXdf(file = "CE1_temp.xdf",
                       varsToKeep = c("rowNum","mod_val_test",dep_var),
                       reportProgress = 0)
      mod <- tmp[tmp$mod_val_test != 3,]
      val <- tmp[tmp$mod_val_test == 3,]
      if (num_resp == smp_resp) {
        rsp <- mod[mod[dep_var]==1,]  # Keep all responders
      } else {
        if (num_resp > smp_resp){
          rsp = mod[mod[dep_var]==1,][sample(num_resp,smp_resp,replace=F),]  # Sample from responders
        } else {
          rsp = mod[mod[dep_var]==1,][sample(num_resp,smp_resp,replace=T),]  # Oversample responders
        }
      }
      if (num_nresp == smp_nresp) {
        nrsp <- mod[mod[dep_var]==0,]  # Keep all non-responders
      } else {
        if (num_nresp > smp_nresp){
          nrsp = mod[mod[dep_var]==0,][sample(num_nresp,smp_nresp,replace=F),]  # Sample from non-responders
        } else {
          nrsp = mod[mod[dep_var]==0,][sample(num_nresp,smp_nresp,replace=T),]  # Oversample non-responders
        }
      }
      # Combine training portions
      tmp2 <- rbind(rsp,nrsp)
      tmp2$mod_val_test <- ifelse(runif(nrow(tmp2),0,1)<.5,1,2)
      # Report counts
      frq1 <- as.data.frame(dv)
      names(frq1) <- c(dep_var,"Count")
      frq1$Percent <- round((frq1$Count / sum(frq1$Count))*100,2)
      frq2 <- as.data.frame(table(tmp2[dep_var],useNA="ifany"))
      names(frq2) <- c(dep_var,"Count")
      frq2$Percent <- round((frq2$Count / sum(frq2$Count))*100,2)
      frq3 <- as.data.frame(table(val[dep_var],useNA="ifany"))
      names(frq3) <- c(dep_var,"Count")
      frq3$Percent <- round((frq3$Count / sum(frq3$Count))*100,2)
      cat('\n', dep_var, " Frequency in Original Dataset \n \n")
      print(frq1)
      cat('\n', dep_var, " Frequency in Oversampled Training Dataset \n \n")
      print(frq2)
      cat('\n', dep_var, " Frequency in Unsampled Validation Dataset \n \n")
      print(frq3)
      # Finish combining parts
      tmp2 <- rbind(tmp2,val)
      tmp2 <- tmp2[order(tmp2$rowNum),]
      rxMerge(inData1 = "CE1_temp.xdf",
              inData2 = tmp2,
              outFile = "CE1_temp.xdf",
              overwrite = T,
              varsToDrop1 = "mod_val_test",
              varsToDrop2 = dep_var,
              matchVars = "rowNum",
              reportProgress=0)
      rxDataStep(inData = "CE1_temp.xdf",
                 outFile = outds,
                 overwrite = T,
                 varsToDrop = "rowNum",
                 reportProgress = 0)
      file.remove('CE1_temp.xdf')
    } 
  }
    
  
  # Apply standard DataSource recodes
  if (ds_present == "Y") {
    cat("Apply standard DataSource recodes \n")
    source("/mnt/projects/shared/pst_qmgisi/Modeling/RCE/ds recodes.R", echo = F)
    DSvars <- c('AU003','XB005','XM915',
                'FD027','EC001','EC002','EC003','EC004','EC005','EC006','EC007','MS456','XM008',
                'XB044','XB107','XB108','XB111','XB112','XB113','XB114','XB115','XB116','XB117',
                'XB118','XB119','XB120','XB121','XB123','XB124','XB125','XB126','XB127','XB128',
                'VW526','MS496','MS112','MS127','MV012','MS052',
                'XB040','XB043','XB101','XB102','XB103','XB105',
                'DS915','DS922','PD017','MS459','MS525','AU009','MS074','AU016','MS919',
                "ET012","ET022","ET032","ET042",
                "ET014","ET024","ET034","ET044",
                'FT050','FT024','XB036','XB037','VW052','VW107','PD012','PD013','PD014','PD015',
                'MV002','MV010',
                'IN057','IN077','IN097','IN117','IN137',
                'XB151','XB152','XB153','XB154','XB155','XB156','XB157','XB158','XB159',
                'XB160','XB161','XB162','XB163','XB164','XB165','XB166','XB167','XB168',
                'XB169','XB170','XB171','XB172','XB173','XB174','XB175','XB176')
    DvarsIn <- c('AU003N','XB005N','XM915N',
                 'FD027r','EC001r','EC002r','EC003r','EC004r','EC005r','EC006r','EC007r','MS456r','XM008r',
                 'XB044r','XB107r','XB108r','XB111r','XB112r','XB113r','XB114r','XB115r','XB116r','XB117r',
                 'XB118r','XB119r','XB120r','XB121r','XB123r','XB124r','XB125r','XB126r','XB127r','XB128r',
                 'VW526r','MS496r','MS112r','MS127r','MV012r','MS052r',
                 'XB040r','XB043r','XB101r','XB102r','XB103r','XB105r',
                 'DS915N','DS922N','PD017N','MS459C','MS525C','MS074_INF','MS074_INF_FLAG',
                 'MS919_INF','MS919_INF_FLAG',
                 "DV001","DV002","DV003",
                 'FT050_AGE','FT024_AGE','XB036_AGE','XB037_AGE','VW052r','VW107r','PD012r','PD013r',
                 'PD014r','PD015r','MV002r','MV010r',
                 'IN057r','IN077r','IN097r','IN117r','IN137r',
                 'XB151r','XB152r','XB153r','XB154r','XB155r','XB156r','XB157r','XB158r','XB159r',
                 'XB160r','XB161r','XB162r','XB163r','XB164r','XB165r','XB166r','XB167r','XB168r',
                 'XB169r','XB170r','XB171r','XB172r','XB173r','XB174r','XB175r','XB176r')
    DvarsOut <- c('AU003','XB005','XM915',
                  'FD027','EC001','EC002','EC003','EC004','EC005','EC006','EC007','MS456','XM008',
                  'XB044','XB107','XB108','XB111','XB112','XB113','XB114','XB115','XB116','XB117',
                  'XB118','XB119','XB120','XB121','XB123','XB124','XB125','XB126','XB127','XB128',
                  'VW526','MS496','MS112','MS127','MV012','MS052',
                  'XB040','XB043','XB101','XB102','XB103','XB105',
                  'DS915','DS922','PD017','MS459','MS525',
                  "ET022","ET032","ET042","ET024","ET034","ET044",
                  'FT050','FT024','XB036','XB037','VW052','VW107','PD012','PD013','PD014','PD015',
                  'MV002','MV010',
                  'IN057','IN077','IN097','IN117','IN137',
                  'XB151','XB152','XB153','XB154','XB155','XB156','XB157','XB158','XB159',
                  'XB160','XB161','XB162','XB163','XB164','XB165','XB166','XB167','XB168',
                  'XB169','XB170','XB171','XB172','XB173','XB174','XB175','XB176')
    vn <- rxGetVarNames(outds)
    prcvars <-DSvars[DSvars %in% vn]
    drvarsin <- DvarsIn[DvarsIn %in% vn]
    drvarsout <- DvarsOut[DvarsOut %in% vn]
    if (length(drvarsin) < 1){drvarsin <- NULL}
    if (length(drvarsout) < 1){drvarsout <- NULL}
    if (length(prcvars) > 0) {
      rxDataStep(inData = outds,
                 outFile = "CE1_temp.xdf",
                 overwrite = T,
                 varsToDrop = drvarsin,
                 transformFunc = f.DSVars,
                 transformVars = prcvars,
                 transformPackages = "car",
                 reportProgress = 0)
      rxDataStep(inData = "CE1_temp.xdf",
                 outFile = outds,
                 overwrite = T,
                 varsToDrop = drvarsout,
                 reportProgress = 0)
      file.remove('CE1_temp.xdf')
    }
  }
  
  # Summary report by mod_val_test
  tmp <- rxDataStep(inData=outds,
                    varsToKeep = c("mod_val_test",dep_var),
                    reportProgress = 0)
  names(tmp) <- c("mod_val_test","dv")
  if (tolower(binary_dv)=="y") {
    frq4 <- ddply(tmp,"mod_val_test",summarize,
                  Count=length(dv),
                  Mean=round(mean(dv),4),
                  Std_Dev=round(sd(dv),4),
                  Min=min(dv),
                  Max=max(dv))
  } else {
    frq4 <- ddply(tmp,"mod_val_test",summarize,
                  Count=length(dv),
                  Mean=round(mean(dv),2),
                  Std_Dev=round(sd(dv),2),
                  Min=round(min(dv),2),
                  Max=round(max(dv),2))
  }
  cat('\n File Summary \n \n')
  print(frq4)
  
  # Clean up
  rm(list=ls())
}

# 2. EDA, PROFILING AND RECODE ##################

## **********************************************
## *****     2.EDA,profiling and Recode     *****
## **********************************************

## *** Macro 2: Recoding macro variables ***;
## %let Profiling = Y;                   ** Request profiling report on all variables (Y/N);
## %let missrate = .75;                  ** maximum missing rate allowed;
## * Binary and Nominal variable types;
## %let concrate = .9;                   ** Maximum amount of file that can be in a single value of an independent variable;
## * Nominal variable types;
## %let valcnt = 50;                     ** Maximum number of unique values allowed. Set to 0 to allow any number;
## %let minbinnc = 500;                  ** Minimum count in a bin to be usable. set to 0 to use minbinpct only;
## %let minbinnp = .05;                  ** Minimum percent of file in a bin to be usable. set to 0 to use minbincnt only;
## %let talpha  = 0.05;                  ** T Test significance level for collapse of bins;
## %let bonfer = N;                      ** Do Bonferoni adjustment for talpha? (Y/N);
## %let nom_method = INDEX;              ** Recoding method for nominal variables: Binary, Index or Mean;
## * Continuous and Ordinal variable types;
## %let pvalue = .05;                    ** Pvalue threshold to include variable;
## %let min_size = 500;                  ** Minimum missing group size for Equal Response imputation;
## %let num_category = 10;               ** Maximum number of categories for profiling variables;
## %let equal_dist = N;                  ** Use equal distance for dividing variables into groups for profiling (Y/N);
## * Continuous variable types;
## %let p_lo =  1;                       ** Lower percentile for checking constant value;
## %let p_hi = 99;                       ** Upper percentile for checking constant value;
## %let impmethodC = median;               ** What method to use for missing imputation for continuous variables?;
## ** Options are ER for Equal Reponse or any proc stdize method;
## ** The location associated with the method will be used for missing values. Be sure to include c or p value if needed;
## **    Method      Location 
## **    ______      ________
## **    MEAN        Mean
## **    MEDIAN      Median
## **    SUM         0
## **    EUCLEN      0
## **    USTD        0
## **    STD         Mean
## **    RANGE       Minimum
## **    MIDRANGE    Midrange
## **    MAXABS      0
## **    IQR         Median
## **    MAD         Median
## **    ABW(c)      Biweight one-step M-estimate
## **    AHUBER(c)   Huber one-step M-estimate
## **    AWAVE(c)    Wave one-step M-estimate
## **    AGK(p)      Mean
## **    SPACING(p)  Mid-minimum spacing
## **    L(p)        L(p)  ;
## %let stdmethodC = STD;                ** Standardization options: any method allowed in proc stdize or NO to skip;
## %let cap_flrC = Y;                    ** Do you want to do capping / flooring to handle outliers? (Y/N);
## %let transformationC = Y;             ** Include transformed variables in evaluation (Y/N);
## * Ordinal variable types;
## %let impmethodO = mean;               ** What method to use for missing imputation for continuous variables?;
## ** Options are ER for Equal Reponse or any proc stdize method;
## %let stdmethodO = No;                 ** Standardization options: any method allowed in proc stdize or NO to skip;
## %let cap_flrO = N;                    ** Do you want to do capping / flooring to handle outliers? (Y/N);
## %let transformationO = N;             ** Include transformed variables in evaluation? (Y/N);

f.recode <- function(dat, keep_list,
                     dep_var = "resp",
                     binary_dv = "Y",
                     prefix = "R1_",
                     missrate = .75,
                     concrate = .9,
                     valcnt = 50,
                     minbinnc = 500,
                     minbinnp = .05,
                     talpha  = 0.05,
                     nom_method = "INDEX",
                     min_size = 500,
                     cap_flrO = "Y",
                     transformationO = "Y",
                     impmethodO = "MEDIAN",
                     stdmethodO = "NO",
                     p_lo = .01,
                     p_hi = .99,
                     cap_flrC = "Y",
                     transformationC = "Y",
                     impmethodC = "MEAN",
                     stdmethodC = "STD",
                     profiling = 'N',
                     equal_dist = 'N',
                     num_category = 10)
{
  
  #************************************************
  #****    Binary variable recode function    *****
  #************************************************
  f.binary <- function(x,binvar,prefix="R1_",missrate=.75,concrate=.9){
    f.catwrite<-function(...,f='Binary Recodes.txt',app=T) cat(...,file=f,append=app)
    
    ##  Binary recodes
    cat("Starting number of binary variables is",length(binvar),"\n")
    nmiss <- apply(is.na(x[binvar]),2,sum)  # Get missing count for each variable
    keep <- names(nmiss[nmiss <= (nrow(x)*missrate)])  # Names of variables with acceptable missing rate
    nvars <- binvar[sapply(x[binvar],is.numeric)]  # Names of numeric variables
    cvars <- binvar[sapply(x[binvar],is.character)]  # Names of character variables
    
    f.nrec <- function(y) {as.numeric(ifelse(is.na(y),0,ifelse(y==1,1,0)))}  # Function to recode numeric variables
    f.crec <- function(y) {as.numeric(ifelse(is.na(y),0,ifelse(y %in% c("Y","y","1"),1,0)))}  # Function to recode character variables
    
    # Start statistics file
    if (length(nvars)>0 & length(cvars)>0){
      bstats <- merge(data.frame(Name=binvar,miss_cnt=nmiss,stringsAsFactors=F),
                      rbind(data.frame(Name=nvars,Type=2,stringsAsFactors=F),
                            data.frame(Name=cvars,Type=1,stringsAsFactors=F)),
                      by="Name", all=TRUE)  # Start statistics file
    } else {
      if (length(nvars)>0){
        bstats <- data.frame(Name=binvar,miss_cnt=nmiss,Type=2,stringsAsFactors=F)
      } else {
        bstats <- data.frame(Name=binvar,miss_cnt=nmiss,Type=1,stringsAsFactors=F)
      }
    }
    
    nvars <- keep[sapply(x[keep],is.numeric)]  # Names of numeric variables
    cvars <- keep[sapply(x[keep],is.character)]  # Names of character variables
    
    if (length(nvars)>0 | length(cvars)>0) {
      # Create new dataset with recoded variables
      if (length(nvars)>0 & length(cvars)>0) {
        x2 <- data.frame(apply(x[nvars],2,f.nrec),
                         apply(x[cvars],2,f.crec))
      } else {
        if (length(nvars)>0) {
          x2 <- data.frame(apply(x[nvars],2,f.nrec))
        } else {
          x2 <- data.frame(apply(x[cvars],2,f.crec))
        }
      }
      
      ocnt <- apply(x2,2,sum)  # Count ones for each variable
      keep2 <- names(ocnt[ocnt >= (nrow(x)*(1-concrate)) & ocnt <= (nrow(x)*concrate)])  # Names of variables with acceptable concentration rate
      x3 <- x2[keep2]  # Dataset with only usable variables
      if (length(keep2) > 0) {
        nn <- paste(prefix, keep2,sep="")  # New names for variables
        names(x3) <-  nn  # Apply new names  
        bstats <- merge(merge(bstats,
                              data.frame(Name=names(ocnt),ones_cnt=ocnt,stringsAsFactors=F),
                              by="Name", all=TRUE),
                        data.frame(Name=keep2,Usable="Yes",stringsAsFactors=F),
                        by="Name", all=TRUE)  # Finish stats file
        rm(x2, ocnt, keep2, nn)  # Clean up
      } else {
        bstats$ones_cnt <- NA
        bstats$Usable <- NA
      }
    } else {bstats$Usable <- NA}
    bstats <- bstats[,c("Name","Type","miss_cnt","ones_cnt","Usable")]
    rm(nmiss, keep, nvars, cvars, binvar)  # Clean up
    bstats$New_Var <- ifelse(bstats$Usable=="Yes",paste0(prefix,bstats$Name),NA)
    save(bstats, file = "CE2_binary_stats.RData",compress=T)  # Save stats file
    cat("Usable binary variable count is",nrow( bstats[!is.na(bstats$Usable),]),"\n")
    
    # Output statistics
    allRows <- seq(length = nrow(bstats)) + 1
    createSheet(ff, name = "Binary Variables")
    writeWorksheet(ff, bstats, sheet = "Binary Variables", startRow = 1, startCol = 1)
    setColumnWidth(ff, sheet = "Binary Variables", column = 1, width = 8192)
    setColumnWidth(ff, sheet = "Binary Variables", column = 6, width = 8192)
    setCellStyle(ff, sheet = "Binary Variables", row = allRows, col = 3, cellstyle = csComma)
    setCellStyle(ff, sheet = "Binary Variables", row = allRows, col = 4, cellstyle = csComma)
    
    #############################################################################
    # Create code
    f.write <- function(y){
      f.catwrite(paste("### Recode: ",y["Name"],'\n',sep=""))
      if (y["Type"]==2) {
        f.catwrite(paste("x$",prefix,y["Name"]," <- ifelse(is.na(x$",y["Name"],"),0,ifelse(x$",y["Name"],"== 1,1,0))",'\n',sep=""))
      } else {
        f.catwrite(paste("x$",prefix,y["Name"]," <- as.numeric(x$",y["Name"]," %in% c('Y','y','1'))",'\n',sep=""))
      }
      f.catwrite('\n')
    }
    if (nrow(bstats[!is.na(bstats$Usable),]) > 0) {
      f.catwrite('######################################################################\n',app=F)
      f.catwrite('#####                 Binary Variable Recodes                    #####\n')
      f.catwrite('######################################################################\n')
      f.catwrite('\n')
      f.catwrite('f.bin_rec <- function(x){ \n')
      
      apply(bstats[!is.na(bstats$Usable),], 1, f.write)
      
      f.catwrite('return(x) \n')
      f.catwrite('} \n')
      
      # Write out list of variables to keep
      keep6 <- paste(prefix,bstats[!is.na(bstats$Usable),"Name"],sep="")
      if (length(keep6) > 1) {
        temp <- data.frame(name=keep6,len=nchar(keep6),splt=NA,outvar=NA)
        temp[1,"splt"] <- 1
        temp[1,"outvar"] <- paste("keep_b <- c('",temp[1,"name"],"'",sep="")
        cnt <- temp[1,"len"]
        ov <- temp[1,"outvar"]
        splt <- 1
        
        for (i in 2:nrow(temp)) {
          if (is.na(temp[i,"splt"])) {
            if (cnt + temp[i,"len"] > 250) {
              splt <- splt + 1
              cnt <- temp[i,"len"]
              ov <- paste("'",temp[i,"name"],"'",sep="")
            } else
            { ov <- paste(ov,",'",temp[i,"name"],"'",sep="")
              cnt <- cnt + temp[i,"len"] + 3
            }
            temp[i,"splt"] <- splt
            temp[i,"outvar"] <- ov
          }
        }
        # Mark first and last record for strings
        temp$slast <- !duplicated( temp[, "splt" ], fromLast=T)
        # Keep lest record for each string
        temp2 <- temp[temp$slast==T,]
        temp2$outvar <- ifelse(temp2$splt==splt,
                               paste(temp2$outvar,")",sep=""),
                               paste(temp2$outvar,",",sep=""))
        
        f.catwrite(paste(temp2$outvar,'\n',sep=""))
      } else {
        if (length(keep6)==1){
          f.catwrite(paste("keep_b <- c('",keep6,"')",'\n',sep=""))
        }
      }      
    }
    # Return recoded data
    return(x3)
  } # End binary recode function
  
  #************************************************
  #****    Nominal variable recode function   *****
  #************************************************
  f.nominal <- function(x,dep_var = "resp",
                        nomvar,
                        prefix="R1_",
                        missrate = .75,
                        concrate = .9,
                        valcnt = 50,
                        minbinnc = 500,
                        minbinnp = .05,
                        talpha  = 0.05,
                        nom_method = "INDEX") {
    
    ##  Nominal recodes
    f.catwrite<-function(...,f='Nominal Recodes.txt',app=T) cat(...,file=f,append=app)
    f.prfwrite<-function(...,f='Nominal Profile Recodes.txt',app=T) cat(...,file=f,append=app)
    
    mincount <- max(nrow(x)*minbinnp,minbinnc)  # Determine minimum bin size
    cat("Starting number of nominal variables is",length(nomvar),"\n")
    
    nvars <- nomvar[sapply(x[nomvar],is.numeric)]  # Names of numeric variables
    cvars <- nomvar[sapply(x[nomvar],is.character)]  # Names of character variables
    
    ## Statistics to clean character variables
    f.max <- function(z) {max(table(z,useNA="ifany"))}
    f.cnt <- function(z) {nrow(table(z,useNA="ifany"))}
    f.test <- function(z) {sum(x[!is.na(z),dep_var])}
    
    # Calculate nominal variable stats
    if (length(nvars) > 0 & length(cvars) > 0) {
      nstats <- data.frame(
        merge(data.frame(Name=nomvar,stringsAsFactors=F),
              rbind(data.frame(Name=nvars,Type=1,stringsAsFactors=F),
                    data.frame(Name=cvars,Type=2,stringsAsFactors=F)),
              by="Name", all=TRUE),
        Obs = apply(!is.na(x[nomvar]), 2, sum),
        MissingObs = apply(is.na(x[nomvar]), 2, sum),
        SumDV =  apply(x[nomvar],2,f.test),
        UniqueCnt = apply(x[nomvar],2,f.cnt),
        MaxCnt = apply(x[nomvar],2,f.max),
        stringsAsFactors=F)
    } else {
      if (length(nvars) > 1) {
        nstats <- data.frame(Name=nomvar,
                             Type=1,
                             Obs = apply(!is.na(x[nomvar]), 2, sum),
                             MissingObs = apply(is.na(x[nomvar]), 2, sum),
                             SumDV =  apply(x[nomvar],2,f.test),
                             UniqueCnt = apply(x[nomvar],2,f.cnt),
                             MaxCnt = apply(x[nomvar],2,f.max),
                             stringsAsFactors=F)
      } else {
        if (length(nvars) == 1) {
          nstats <- data.frame(Name=nomvar,
                               Type=1,
                               Obs = sum(!is.na(x[nomvar])),
                               MissingObs = sum(is.na(x[nomvar])),
                               SumDV =  f.test(x[nomvar]),
                               UniqueCnt = f.cnt(x[nomvar]),
                               MaxCnt = f.max(x[nomvar]),
                               stringsAsFactors=F)
        } else {
          if (length(cvars) > 1) {
            nstats <- data.frame(Name=nomvar,
                                 Type=2,
                                 Obs = apply(!is.na(x[nomvar]), 2, sum),
                                 MissingObs = apply(is.na(x[nomvar]), 2, sum),
                                 SumDV =  apply(x[nomvar],2,f.test),
                                 UniqueCnt = apply(x[nomvar],2,f.cnt),
                                 MaxCnt = apply(x[nomvar],2,f.max),
                                 stringsAsFactors=F)
          } else {
            if (length(cvars) == 1) {
              nstats <- data.frame(Name=nomvar,
                                   Type=2,
                                   Obs = sum(!is.na(x[nomvar])),
                                   MissingObs = sum(is.na(x[nomvar])),
                                   SumDV =  f.test(x[nomvar]),
                                   UniqueCnt = f.cnt(x[nomvar]),
                                   MaxCnt = f.max(x[nomvar]),
                                   stringsAsFactors=F)
            }
          }
        }
      }
    }
    
    # Determine usability
    nstats$Usable <- ifelse((nstats$MissingObs / (nstats$Obs+nstats$MissingObs) < missrate & 
                               nstats$MaxCnt / (nstats$Obs+nstats$MissingObs) < concrate & 
                               nstats$SumDV != 0 & 
                               nstats$SumDV != nstats$Obs &
                               (nstats$UniqueCnt <= valcnt | valcnt == 0)),"Yes",NA)
    keep2 <- as.character(nstats[!is.na(nstats$Usable),"Name"])  # Names of usable variables
    typ2 <- as.character(nstats[!is.na(nstats$Usable),"Type"])  # Types for usable variables
    nstats$SumDV <- NULL
    
    if (length(keep2) > 0) {
      f.catwrite('######################################################################\n',app=F)
      f.catwrite('#####                 Nominal Variable Recodes                   #####\n')
      f.catwrite('######################################################################\n')
      f.catwrite('\n')
      f.catwrite('f.nom_rec <- function(x){ \n')
      f.prfwrite('f.pnom_rec <- function(x){ \n',app=F)
      
      if (toupper(nom_method)=="INDEX") {oa <- mean(x[,dep_var])}  # Get overall average if using method Index
      # Loop through each variable
      for (j in 1:length(keep2)) {
        # Summarize and initialize
        tmp <- data.frame(value=x[,keep2[j]],dv=x[,dep_var],stringsAsFactors=F)
        temp1 <- ddply(tmp,"value",summarize,
                       freq=length(dv),
                       mean=mean(dv),
                       var=var(dv))
        rm(tmp)
        temp1$var <- ifelse(is.na(temp1$var),0,temp1$var)
        temp1 <- temp1[order(temp1$mean, decreasing=T),]
        temp1$group <- NA
        temp1[1,5] <- 1
        
        grpfreq <- temp1[1,2]
        grpmean <- temp1[1,3]
        grpvar <- temp1[1,4]
        grp <- 1
        
        # Loop to collapse bins
        for (i in 2:nrow(temp1)) {
          if (is.na(temp1[i,5])) {
            if (grpfreq < mincount) {
              temp1[i,5] <- grp
              grpmean <- ((grpmean * grpfreq) + (temp1[i,3] * temp1[i,2])) / (grpfreq + temp1[i,2])
              if (grpfreq + temp1[i,2] > 2) grpvar <- (((grpfreq - 1) * grpvar) + ((temp1[i,2] - 1) * temp1[i,4])) / (grpfreq + temp1[i,2] - 2) else grpvar <- 0
              grpfreq <- grpfreq + temp1[i,2]
            } else
            { pvar <- (((grpfreq - 1) * grpvar) + ((temp1[i,2] - 1) * temp1[i,4]))/ (grpfreq + temp1[i,2] - 2)
              t <- (grpmean - temp1[i,3])/sqrt(pvar * ((1/grpfreq) + (1/temp1[i,2])))
              prob = 1-pt(abs(t),df=(grpfreq + temp1[i,2] - 2))
              if (prob <= talpha) {
                grp <- grp + 1
                grpfreq <- temp1[i,2]
                grpmean <- temp1[i,3]
                grpvar <- temp1[i,4]
                temp1[i,5] <- grp
              } else
              { temp1[i,5] <-grp
                grpmean <- ((grpmean * grpfreq) + (temp1[i,3] * temp1[i,2])) / (grpfreq + temp1[i,2])
                grpvar <- (((grpfreq - 1) * grpvar) + ((temp1[i,2] - 1) * temp1[i,4])) / (grpfreq + temp1[i,2] - 2)
                grpfreq <- grpfreq + temp1[i,2]
              }
            }
          }
        }
        # Check for multiple buckets
        if (temp1[1,"group"]!=grp) {f.catwrite("### Recode: ",keep2[j],"\n")} else {
          nstats$Usable[nstats$Name == keep2[j]] <- NA
          next
        }
        # Calculate bin value to use for INDEX and MEAN methods
        if (toupper(nom_method)=="INDEX" | toupper(nom_method)=="MEAN") {
          tst <- ddply(temp1,"group",summarize,
                       frq=sum(freq),
                       tot=sum(freq*mean))
          if (toupper(nom_method)=="INDEX") {tst$val <- ((tst$tot/tst$frq)/oa)*100 } else 
          {tst$val <- tst$tot/tst$frq} 
          temp1 <- merge(temp1,tst[c("group","val")],by="group", all=TRUE)
          temp1 <- temp1[order(temp1$group, decreasing=T),]
        }
        # Copy variable value to new variable with formatting
        if (typ2[j]==2){
          temp1$value2 <- ifelse(is.na(temp1$value),NA,paste("'",temp1$value,"'",sep=""))
        } else {
          temp1$value2 <- temp1$value
        }
        # Create strings to use in new code
        temp1$len <- nchar(temp1$value2)
        temp1$splt <- NA
        temp1$outvar <- NA
        temp1[1,"splt"] <- 1
        temp1[1,"outvar"] <- as.character(temp1[1,"value2"])
        cnt <- temp1[1,"len"]
        ov <- temp1[1,"outvar"]
        splt <- 1
        
        for (i in 2:nrow(temp1)) {
          if (is.na(temp1[i,"splt"])) {
            if (cnt + temp1[i,"len"] > 130 | temp1[i,"group"] != temp1[i-1,"group"]) {
              splt <- splt + 1
              cnt <- temp1[i,"len"]
              ov <- as.character(temp1[i,"value2"])
            } else
            { ov <- paste(ov,temp1[i,"value2"],sep=",",collapse=",")
              cnt <- cnt + temp1[i,"len"] + 1
            }
            temp1[i,"splt"] <- splt
            temp1[i,"outvar"] <- ov
          }
        }
        # Mark first and last record for strings
        temp1$sfirst <- !duplicated( temp1[, "splt" ]) 
        temp1$slast <- !duplicated( temp1[, "splt" ], fromLast=T)
        # Keep lest record for each string
        temp2 <- temp1[temp1$slast==T,]
        # Mark first and last record for each group
        temp2$first <- !duplicated( temp2[, "group" ]) 
        temp2$last <- !duplicated( temp2[, "group" ], fromLast=T)
        
        # Create code lines for writing out
        if (toupper(nom_method)=="BINARY"){
          #temp2 <- temp2[temp2$group != grp,]  # Remove last group as not needed
          temp2$wvar <-
            ifelse(temp2$first & temp2$last, paste("x$",prefix,keep2[j],".",temp2$group," <- as.numeric(x$",keep2[j]," %in% c(",temp2$outvar,"))\n",sep=""),
                   ifelse(temp2$first, paste("x$",prefix,keep2[j],".",temp2$group," <- as.numeric(x$",keep2[j]," %in% c(",temp2$outvar,",\n",sep=""),
                          ifelse(temp2$last, paste(temp2$outvar,"))\n",sep=""),paste(temp2$outvar,",\n",sep=""))))
        } else {
          temp2 <- temp2[(temp2$first & temp2$group==1)|temp2$group > 1,]
          temp2$wvar <- ifelse(temp2$group==1,paste(temp2$val,paste(rep(")",grp-1),sep="",collapse=""),"\n",sep="") ,NA)
          temp2[1,"wvar"] <- ifelse(temp2[1,"last"], paste("x$",prefix,keep2[j]," <- ifelse(x$",keep2[j]," %in% c(",temp2[1,"outvar"],"),",temp2[1,"val"],",\n",sep=""),
                                    paste("x$",prefix,keep2[j]," <- ifelse(x$",keep2[j]," %in% c(",temp2[1,"outvar"],",\n",sep=""))
          temp2$wvar <- ifelse(is.na(temp2$wvar) & temp2$first & temp2$last,paste("ifelse(x$",keep2[j]," %in% c(",temp2$outvar,"),",temp2$val,",\n",sep=""),
                               ifelse(is.na(temp2$wvar) & temp2$first ,paste("ifelse(x$",keep2[j]," %in% c(",temp2$outvar,",\n",sep=""),
                                      ifelse(is.na(temp2$wvar) & temp2$last ,paste(temp2$outvar,"),",temp2$val,",\n",sep=""),
                                             ifelse(is.na(temp2$wvar),paste(temp2$outvar,",\n",sep=""),temp2$wvar))))
        }
        # Create lookup table
        if (!exists("lkup", mode="list")){
          if (toupper(nom_method) %in% c("MEAN","INDEX")){
            lkup <- data.frame(variable=keep2[j],
                               ddply(temp1,"val",function(x){c(Category=paste(x$value,collapse=","))}),
                               stringsAsFactors = F)} else {
                                 lkup <- data.frame(variable=keep2[j],
                                                    ddply(temp1,"group",function(x){c(Category=paste(x$value,collapse=","))}),
                                                    stringsAsFactors = F)}
        } else {
          if (toupper(nom_method) %in% c("MEAN","INDEX")){
            lkup <- rbind(lkup,data.frame(variable=keep2[j],
                                          ddply(temp1,"val",function(x){c(Category=paste(x$value,collapse=","))}),
                                          stringsAsFactors = F))} else {
                                            lkup <- rbind(lkup,data.frame(variable=keep2[j],
                                                                          ddply(temp1,"group",function(x){c(Category=paste(x$value,collapse=","))}),
                                                                          stringsAsFactors = F))}
        }
        # Write out recodes
        f.prfwrite(temp2$wvar)
        f.prfwrite('\n')
        if (toupper(nom_method)=="BINARY"){temp2 <- temp2[temp2$group != grp,]}  # Remove last group as not needed
        f.catwrite(temp2$wvar)
        f.catwrite('\n')
      }
    }
    f.catwrite('return(x) \n')
    f.catwrite('} \n')
    f.prfwrite('return(x) \n')
    f.prfwrite('} \n')
    
    # Final list of usable variables
    keep3 <- as.character(nstats[!is.na(nstats$Usable),"Name"])  # Names of usable variables
    
    # Save stats file
    save(nstats, file = "CE2_nominal_stats.RData",compress=T)
    cat("Usable nominal variable count is",length(keep3),"\n")
    
    # Save lookup table
    if (exists("lkup", mode="list")){save(lkup, file = "nominal_lkup.RData",compress=T)}
    
    # Output statistics
    allRows <- seq(length = nrow(nstats)) + 1
    createSheet(ff, name = "Nominal Variables")
    writeWorksheet(ff, nstats, sheet = "Nominal Variables", startRow = 1, startCol = 1)
    setColumnWidth(ff, sheet = "Nominal Variables", column = 1, width = 8192)
    setCellStyle(ff, sheet = "Nominal Variables", row = allRows, col = 3, cellstyle = csComma)
    setCellStyle(ff, sheet = "Nominal Variables", row = allRows, col = 4, cellstyle = csComma)
    setCellStyle(ff, sheet = "Nominal Variables", row = allRows, col = 5, cellstyle = csComma)
    setCellStyle(ff, sheet = "Nominal Variables", row = allRows, col = 6, cellstyle = csComma)
    
    # Write out list of variables to keep
    if (length(keep3) > 1){
      # Recoded data
      source("Nominal Recodes.txt", local=T,echo = F)
      x <- rxDataStep(inData = dat,
                      rowSelection = mod_val_test != 3,
                      varsToKeep = keep3,
                      transformFunc = f.nom_rec,
                      transformVars = keep3,
                      maxRowsByCols = mx,
                      reportProgress = 0)
      x2 <- x[,setdiff(names(x),c("mod_val_test",keep3))]
      keep4 <- names(x2)
      
      temp <- data.frame(name=keep4,len=nchar(keep4),splt=NA,outvar=NA,stringsAsFactors=F)
      temp[1,"splt"] <- 1
      temp[1,"outvar"] <- paste("keep_n <- c('",temp[1,"name"],"'",sep="")
      cnt <- temp[1,"len"]
      ov <- temp[1,"outvar"]
      splt <- 1
      
      for (i in 2:nrow(temp)) {
        if (is.na(temp[i,"splt"])) {
          if (cnt + temp[i,"len"] > 250) {
            splt <- splt + 1
            cnt <- temp[i,"len"]
            ov <- paste("'",temp[i,"name"],"'",sep="")
          } else
          { ov <- paste(ov,",'",temp[i,"name"],"'",sep="")
            cnt <- cnt + temp[i,"len"] + 3
          }
          temp[i,"splt"] <- splt
          temp[i,"outvar"] <- ov
        }
      }
      # Mark first and last record for strings
      temp$slast <- !duplicated( temp[, "splt" ], fromLast=T)
      # Keep lest record for each string
      temp2 <- temp[temp$slast==T,]
      temp2$outvar <- ifelse(temp2$splt==splt,
                             paste(temp2$outvar,")",sep=""),
                             paste(temp2$outvar,",",sep=""))
      
      f.catwrite(paste(temp2$outvar,'\n',sep=""))
      
      return(x2)
    } else {
      if (length(keep3)==1){
        # Recoded data
        source("Nominal Recodes.txt", local=T,echo = F)
        x <- rxDataStep(inData = dat,
                        rowSelection = mod_val_test != 3,
                        varsToKeep = keep3,
                        transformFunc = f.nom_rec,
                        transformVars = keep3,
                        maxRowsByCols = mx,
                        reportProgress = 0)
        x2 <- x[,setdiff(names(x),c("mod_val_test",keep3))]
        keep4 <- names(x2)
        
        f.catwrite(paste("keep_n <- c('",keep4,"')",'\n',sep=""))
        
        return(x2)
      } else {
        f.catwrite("keep_n <- NA",'\n',sep="")
        return(NULL)
      }
    }
  } # End nominal recode function
  
  #************************************************
  #****    Ordinal variable recode function   *****
  #************************************************
  
  f.ordinal <- function(x, dep_var = "resp",
                        ordvar,
                        prefix = "R1_",
                        binary_dv = "Y",
                        missrate = .75,
                        concrate = .9,
                        cap_flrO = "Y",
                        transformationO = "Y",
                        impmethodO = "MEAN",
                        stdmethodO = "STD",
                        min_size = 500) {
    f.catwrite<-function(...,f='Ordinal Recodes.txt',app=T) cat(...,file=f,append=app)
    
    # Ordinal Recodes
    cat("Starting number of ordinal variables is",length(ordvar),"\n")
    ## Basic statistics needed for transformations
    f.ckn <- function(z) {min(x[!is.na(z),dep_var])}
    f.ckx <- function(z) {max(x[!is.na(z),dep_var])}
    f.max <- function(z) {max(table(z,useNA="ifany"))}
    f.cnt <- function(z) {nrow(table(z,useNA="ifany"))}
    f.stats1 <- function(y){
      tmp <- data.frame(
        Name = ordvar,
        Obs = apply(!is.na(y), 2, sum),
        MissingObs = apply(is.na(y), 2, sum),
        stringsAsFactors=F)
    }
    f.stats2 <- function(y){
      tmp <- data.frame(
        Name = ordvar2,
        MinCk = apply(y,2,f.ckn),
        MaxCk = apply(y,2,f.ckx),
        UniqueCnt = apply(y,2,f.cnt),
        MaxCnt = apply(y,2,f.max),
        Median = apply(y, 2,quantile, probs=0.50, na.rm=TRUE),
        stringsAsFactors=F)
    }
    ostats <- f.stats1(x[ordvar])
    ordvar2 <- as.character(ostats[ostats$Obs>0,"Name"])
    ostats <- merge(ostats,f.stats2(x[ordvar2]),by="Name", all=TRUE)
    ostats$Usable <- ifelse(ostats$MissingObs / (ostats$Obs+ostats$MissingObs) < missrate &
                              ostats$MaxCnt / (ostats$Obs+ostats$MissingObs) < concrate & 
                              ostats$MinCk != ostats$MaxCk,"Yes",NA)
    ostats$MinCk <- NULL
    ostats$MaxCk <- NULL
    keep2 <- as.character(ostats[!is.na(ostats$Usable),"Name"])  # Names of usable variables
    cat("Usable ordinal variable count is",length(keep2),"\n")
    if (length(keep2) > 0) {
      if (toupper(cap_flrO)=="Y") {
        f.bounds<-function(y){
          quants<-quantile(y,c(0,.01,.25,.5,.75,.99,1),na.rm=T)
          iqr <- pmax((quants[5]-quants[3]),(quants[6]-quants[5]),(quants[3]-quants[2]))
          lb <- pmin(pmax((quants[3] - (1.5*iqr)), quants[1]), quants[2])
          ub <- pmin(pmax((quants[5] + (1.5*iqr)), quants[7]), quants[6])
          lb <- ifelse(lb == ub, quants[1], lb)
          ub <- ifelse(lb == ub, quants[7], ub)
          bnds <- cbind(lb,ub)
          return(bnds)
        }
        tmp <- as.data.frame(apply(x[keep2],2,f.bounds),stringsAsFactors=F)
        ostats <- merge(ostats,
                        data.frame("Name"=keep2,
                                   lb=as.numeric(tmp[1,]),
                                   ub=as.numeric(tmp[2,]),
                                   stringsAsFactors=F),
                        by="Name", all=TRUE)
      }
      
      # Determine missing imputation value: options ER, MEAN, MEDIAN, ZERO, MINIMUM, MIDRANGE
      if (toupper(impmethodO)=="ER"){
        lst1 <- ostats[ostats$MissingObs >= min_size & !is.na(ostats$Usable),"Name"]  # Use regression
        lst2 <- ostats[ostats$MissingObs < min_size & !is.na(ostats$Usable),"Name"]  # Use Median
        
        if (length(lst1) > 1) {
          f.miss_rr <- function(z){mean(x[is.na(z),dep_var])}
          miss_rr <- apply(x[,lst1],2,f.miss_rr)
        } else {
          if (length(lst1) == 1) {
            miss_rr <- mean(x[is.na(x[,lst1]),dep_var])
          } 
        }
        
        if (binary_dv == "Y" & length(lst1) > 0) {
          miss_rr <- ifelse(miss_rr==0,.0001,ifelse(miss_rr==1,.9999,miss_rr))
        }
        
        f.eval <- function(y){
          tmp <- data.frame(dv=dv,iv=y,stringsAsFactors=F)
          if (binary_dv == "Y") {
            mod1 <- glm(dv~iv, binomial,data=tmp)
          } else {
            mod1 <- lm(dv~iv, data=tmp)
          }
          ck <- data.frame(summary(mod1)["coefficients"],stringsAsFactors=F)
          val <- data.frame(Intercept=ck[1,1],Estimate=ck[2,1],Prob=ck[2,4],stringsAsFactors=F)
          return(val)
        }
        dv <- x[,dep_var]
        if (length(lst1) > 1) {
          coef <- data.frame(do.call(rbind,apply(x[,lst1],2,f.eval)),miss_rr=miss_rr,stringsAsFactors=F)
          coef$Name = row.names(coef)
        } else {
          if (length(lst1) == 1) {
            tmp <- data.frame(dv=dv,iv=x[,lst1],stringsAsFactors=F)
            if (binary_dv == "Y") {
              mod1 <- glm(dv~iv, binomial,data=tmp)
            } else {
              mod1 <- lm(dv~iv, data=tmp)
            }
            ck <- data.frame(summary(mod1)["coefficients"],stringsAsFactors=F)
            val <- data.frame(Intercept=ck[1,1],Estimate=ck[2,1],Prob=ck[2,4],stringsAsFactors=F)
            coef <- data.frame(Intercept=ck[1,1],Estimate=ck[2,1],Prob=ck[2,4],miss_rr=miss_rr,stringsAsFactors=F)
            coef$Name = lst1
          }
        }
        if (length(lst1) > 0) {
          er <- merge(coef,ostats[,c("Name","Median","lb","ub")],by="Name")
          
          if (binary_dv == "Y"){
            er$Miss_Impute <- ifelse(er$Prob > .05, er$Median,
                                     (log(er$miss_rr/(1-er$miss_rr))-er$Intercept)/er$Estimate)
          } else {
            er$Miss_Impute <- ifelse(er$Prob > .05, er$Median,
                                     (er$miss_rr-er$Intercept)/er$Estimate)
          }
          er$Miss_Impute <- ifelse(er$Miss_Impute < er$lb,er$lb,
                                   ifelse(er$Miss_Impute > er$ub,er$ub,er$Miss_Impute))
        }
        
        
        if (length(lst1) > 0 & length(lst2) > 0) {
          tmp <- ostats[ostats$Name %in% lst2, c("Name","Median")]
          names(tmp) <- c("Name","Miss_Impute")
          tmp2 <- rbind(tmp,er[,c("Name","Miss_Impute")])
          ostats <- merge(ostats, tmp2, by="Name", all=TRUE) 
        } else {
          if (length(lst1) > 0) {
            ostats <- merge(ostats, er[,c("Name","Miss_Impute")], by="Name", all=TRUE) 
            missval <- er$Miss_Impute
          } else {
            tmp <- data.frame(ostats[ostats$Name %in% lst2, c("Name","Median")])
            names(tmp) <- c("Name","Miss_Impute")
            ostats <- merge(ostats, tmp, by="Name", all=TRUE)
          }
        }
        missval <- ostats[ostats$Name %in% keep2,'Miss_Impute']
      } else {
        if (toupper(impmethodO)=="MEDIAN") {missval <- ostats[ostats$Name %in% keep2,"Median"]} else {
          if (toupper(impmethodO)=="ZERO") {missval <- rep.int(0,length(keep2))} else {
            if (toupper(impmethodO)=="MINIMUM") {missval <- apply(x[keep2],2,min,na.rm=T)} else {
              if (toupper(impmethodO)=="MIDRANGE") {missval <- (apply(x[keep2],2,min,na.rm=T)+apply(x[keep2],2,max,na.rm=T))/2} else {
                missval <- apply(x[keep2],2,mean,na.rm=T)  # Default is mean
              } } } }
        ostats <- merge(ostats,
                        data.frame("Name"=keep2,
                                   Miss_Impute=missval,
                                   stringsAsFactors=F),
                        by="Name", all=TRUE)
      }
      
      # Apply capping / flooring if option is chosen
      if (toupper(cap_flrO)=="Y") {
        f.cap <- function(y){
          quants<-quantile(y,c(0,.01,.25,.5,.75,.99,1),na.rm=T)
          iqr <- pmax((quants[5]-quants[3]),(quants[6]-quants[5]),(quants[3]-quants[2]))
          lb <- pmin(pmax((quants[3] - (1.5*iqr)), quants[1]), quants[2])
          ub <- pmin(pmax((quants[5] + (1.5*iqr)), quants[7]), quants[6])
          lb <- ifelse(lb == ub, quants[1], lb)
          ub <- ifelse(lb == ub, quants[7], ub)
          RW = pmin(pmax(y,lb),ub)
          return(RW)
        }
        x2<-data.frame(do.call(cbind,lapply(x[keep2],f.cap)),stringsAsFactors=F)
      } else {x2 <- x[,keep2]}
      
      # Replace missings
      if (length(keep2) > 1) {
        x2[is.na(x2)] <- matrix(missval,nrow(x2),length(keep2),byrow=T)[is.na(x2)]  
      } else {
        x2[is.na(x2[keep2]),keep2] <- missval
      }
      
      # Apply transforms if option is chosen
      if (toupper(transformationO)=="Y") {
        f.transform<-function(y,sqrt_lo=0,log_lo=0.0001){
          cbind(
            RW = y,
            SR = sqrt(pmax(y,sqrt_lo,na.rm=T)),
            SQ = y^2,
            LN = log(pmax(y,log_lo,na.rm=T)),
            EP = exp(pmin(y,100,na.rm=T)),
            IV = ifelse(y==0, 0, 1/y))
        }
        x3<-data.frame(do.call(cbind,lapply(x2,f.transform)),stringsAsFactors=F)
        keep3<-apply(expand.grid(colnames(f.transform(NA)),keep2)[,2:1],1,paste,collapse='.')
        keep3 <- paste(prefix,keep3,sep="")
        names(x3)<-keep3  
      } else {
        x3 <- x2
        keep3 <- paste(prefix,keep2,".RW",sep="")
        names(x3) <- keep3
      }
      
      # Apply standardization: MEAN,MEDIAN,SUM,EUCLEN,USTD,STD,RANGE,MIDRANGE,MAXABS,IQR,MAD,NO to skip
      if (toupper(stdmethodO) %in% c("MEAN","MEDIAN","SUM","EUCLEN","USTD","STD",
                                     "RANGE","MIDRANGE","MAXABS","IQR","MAD")){
        if (toupper(stdmethodO) == "MEAN"){
          location <- apply(x3,2,mean,na.rm=T)
          scale <- rep.int(1,length(keep3))
        }
        if (toupper(stdmethodO) == "MEDIAN"){
          location <- apply(x3,2,median,na.rm=T)
          scale <- rep.int(1,length(keep3))
        }
        if (toupper(stdmethodO) == "SUM"){
          location <- rep.int(0,length(keep3))
          scale <- apply(x3,2,sum,na.rm=T)
        }
        if (toupper(stdmethodO) == "EUCLEN"){
          f.euclen <- function(y){
            y2 <- sqrt(sum(y^2,na.rm=T))
          }
          location <- rep.int(0,length(keep3))
          scale <- apply(x3,2,f.euclen)
        }
        if (toupper(stdmethodO) == "USTD"){
          f.ustd <- function(y){
            y2 <- sqrt(sum((y-0)^2,na.rm=T)/sum(!is.na(y)))
          }
          location <- rep.int(0,length(keep3))
          scale <- apply(x3,2,f.ustd)
        }
        if (toupper(stdmethodO) == "STD"){
          location <- apply(x3,2,mean,na.rm=T)
          scale <- apply(x3,2,sd,na.rm=T)
        }
        if (toupper(stdmethodO) == "RANGE"){
          f.range <- function(y){
            y2 <- max(y,na.rm=T)-min(y,na.rm=T)
          }
          location <- apply(x3,2,min,na.rm=T)
          scale <- apply(x3,2,f.range)
        }
        if (toupper(stdmethodO) == "MIDRANGE"){
          f.midrangeL <- function(y){
            y2 <- (max(y,na.rm=T)+min(y,na.rm=T))/2
          }
          f.midrangeS <- function(y){
            y2 <- (max(y,na.rm=T)-min(y,na.rm=T))/2
          }
          location <- apply(x3,2,f.midrangeL)
          scale <- apply(x3,2,f.midrangeS)
        }
        if (toupper(stdmethodO) == "MAXABS"){
          f.maxabs <- function(y){
            y2 <- max(abs(y),na.rm=T)
          }
          location <- rep.int(0,length(keep3))
          scale <- apply(x3,2,f.maxabs)
        }
        if (toupper(stdmethodO) == "IQR"){
          f.iqr <- function(y){
            y2 <- quantile(y,probs=.75,na.rm=T)-quantile(y,probs=.25,na.rm=T)
          }
          location <- apply(x3,2,median,na.rm=T)
          scale <- apply(x3,2,f.iqr)
        }
        if (toupper(stdmethodO) == "MAD"){
          f.mad <- function(y){
            y2 <- median(abs(y-median(y,na.rm=T)),na.rm=T)
          }
          location <- apply(x3,2,median,na.rm=T)
          scale <- apply(x3,2,f.mad)
        }
        x4<-data.frame(scale(x3,center=location,scale=scale),stringsAsFactors=F)
      } else {x4 <- x3}
      
      ## Eliminate variables with 0 standard deviation
      ck <- round(apply(x4,2,sd),4)
      keep4 <- names(ck[ck>0])
      keep4 <- keep4[!is.na(keep4)]
      x4 <- x4[,keep4]
      
      ## Choose best form for each variable
      ck <- as.data.frame(cor(x4,x[,dep_var]),stringsAsFactors=F)
      names(ck) <- "Corr"
      ck$New_Var <- row.names(ck)
      ck$Name <- substr(ck$New_Var,nchar(prefix)+1,nchar(ck$New_Var)-3)
      ck$abscorr <- abs(ck$Corr)
      ck <- ck[order(ck$Name,ck$abscorr),]
      ck <- ck[!duplicated( ck[, "Name" ], fromLast=T),]
      ck$abscorr <- NULL
      keep5 <- ck$New_Var
      
      ## Finalize stats file
      if (toupper(stdmethodO) %in% c("MEAN","MEDIAN","SUM","EUCLEN","USTD","STD",
                                     "RANGE","MIDRANGE","MAXABS","IQR","MAD")){
        ostats <- merge(ostats,
                        merge(ck,
                              data.frame("New_Var"=keep3,
                                         Location=location,
                                         Scale=scale,
                                         stringsAsFactors=F),
                              by="New_Var"),
                        by = "Name", all=T)
      } else {
        ostats <- merge(ostats,ck, by="Name", all=T)
      }
      keep6 <- ostats$New_Var[!is.na(ostats$New_Var)]
    } else {
      ostats$Usable <- NA
      keep6 <- NULL
    }
    
    save(ostats, file = "CE2_ordinal_stats.RData",compress=T)  # Save stats file
    
    # Output statistics
    allRows <- seq(length = nrow(ostats)) + 1
    createSheet(ff, name = "Ordinal Variables")
    writeWorksheet(ff, ostats, sheet = "Ordinal Variables", startRow = 1, startCol = 1)
    setColumnWidth(ff, sheet = "Ordinal Variables", column = 1, width = 8192)
    setColumnWidth(ff, sheet = "Ordinal Variables", column = c(2,3,4,5,6,7,8,9,10,11,12), width = 3072)
    setCellStyle(ff, sheet = "Ordinal Variables", row = allRows, col = 2, cellstyle = csComma)
    setCellStyle(ff, sheet = "Ordinal Variables", row = allRows, col = 3, cellstyle = csComma)
    setCellStyle(ff, sheet = "Ordinal Variables", row = allRows, col = 4, cellstyle = csComma)
    setCellStyle(ff, sheet = "Ordinal Variables", row = allRows, col = 5, cellstyle = csComma)
    setCellStyle(ff, sheet = "Ordinal Variables", row = allRows, col = 6, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Ordinal Variables", row = allRows, col = 8, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Ordinal Variables", row = allRows, col = 9, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Ordinal Variables", row = allRows, col = 10, cellstyle = csDec4)
    setCellStyle(ff, sheet = "Ordinal Variables", row = allRows, col = 11, cellstyle = csDec4)
    setCellStyle(ff, sheet = "Ordinal Variables", row = allRows, col = 12, cellstyle = csDec4)
    
    #************************************************
    # Create code
    f.trans <- function(z) {
      tmp <- ifelse (z["type"] == 'SR', paste("x$",z["New_Var"]," <- sqrt(pmax(x$",prefix,z["Name"],".RW,0,na.rm=T))",sep=""),
                     ifelse (z["type"] == 'SQ', paste("x$",z["New_Var"]," <- x$",prefix,z["Name"],".RW^2",sep=""),
                             ifelse (z["type"] == 'LN', paste("x$",z["New_Var"]," <- log(pmax(x$",prefix,z["Name"],".RW,.0001,na.rm=T))",sep=""),
                                     ifelse(z["type"] == 'EP', paste("x$",z["New_Var"]," <- exp(pmin(x$",prefix,z["Name"],".RW,100,na.rm=T))",sep=""),
                                            ifelse(z["type"] == 'IV', paste("x$",z["New_Var"]," <- ifelse(x$",prefix,z["Name"],".RW==0,0,1/x$",prefix,z["Name"],".RW)",sep=""),
                                                   " ")))))
    }
    
    f.write <- function(y){
      f.catwrite(paste("### Recode: ",y["Name"],'\n',sep="")) 
      f.catwrite(paste("x$",prefix,y["Name"],".RW <- x$",y["Name"],'\n',sep=""))                                           # Create new variable
      f.catwrite(paste("x$",prefix,y["Name"],".RW[is.na(x$",prefix,y["Name"],".RW)] <- ",y["Miss_Impute"],'\n',sep=""))                # Set missing values to median
      if (toupper(cap_flrO)=="Y") {
        f.catwrite(paste("x$",prefix,y["Name"],".RW <- pmin(pmax(x$",prefix,y["Name"],".RW,",y["lb"],"),",y["ub"],")",'\n',sep="")) # Set variable bounds  
      }
      if (toupper(transformationO)=="Y") {
        trans <- f.trans(y)
        f.catwrite(paste(trans,'\n',sep=""))
      }
      if (toupper(stdmethodO) != "NO") {
        f.catwrite(paste("x$",y["New_Var"]," <- (x$",y["New_Var"]," - ",y["Location"],") / ",y["Scale"],'\n',sep=""))         # transform  
      }
      f.catwrite('\n')
    }
    if (length(keep6) > 0) {
      ostats$type <- substr(ostats$New_Var,nchar(ostats$New_Var)-1,nchar(ostats$New_Var))
      f.catwrite('######################################################################\n',app=F)
      f.catwrite('#####               Ordinal Variable Recodes                     #####\n')
      f.catwrite('######################################################################\n')
      f.catwrite('\n')
      f.catwrite('f.ord_rec <- function(x){ \n')
      
      apply(ostats[ostats$New_Var %in% keep6,], 1, f.write)
      
      f.catwrite('return(x) \n')
      f.catwrite('} \n')
      
      # Write out list of variables to keep
      if (length(keep6 ) > 1){
        temp <- data.frame(name=keep6,len=nchar(keep6),splt=NA,outvar=NA,stringsAsFactors=F)
        temp[1,"splt"] <- 1
        temp[1,"outvar"] <- paste("keep_o <- c('",temp[1,"name"],"'",sep="")
        cnt <- temp[1,"len"]
        ov <- temp[1,"outvar"]
        splt <- 1
        
        for (i in 2:nrow(temp)) {
          if (is.na(temp[i,"splt"])) {
            if (cnt + temp[i,"len"] > 250) {
              splt <- splt + 1
              cnt <- temp[i,"len"]
              ov <- paste("'",temp[i,"name"],"'",sep="")
            } else
            { ov <- paste(ov,",'",temp[i,"name"],"'",sep="")
              cnt <- cnt + temp[i,"len"] + 3
            }
            temp[i,"splt"] <- splt
            temp[i,"outvar"] <- ov
          }
        }
        # Mark first and last record for strings
        temp$slast <- !duplicated( temp[, "splt" ], fromLast=T)
        # Keep lest record for each string
        temp2 <- temp[temp$slast==T,]
        temp2$outvar <- ifelse(temp2$splt==splt,
                               paste(temp2$outvar,")",sep=""),
                               paste(temp2$outvar,",",sep=""))
        
        f.catwrite(paste(temp2$outvar,'\n',sep=""))
      } else {
        if (length(keep6)==1){
          f.catwrite(paste("keep_o <- c('",keep6,"')",'\n',sep=""))
        }
      }     
    }
    
    # Return recoded data
    if (length(keep6) > 0){
      x5 <- data.frame(x4[,keep6],stringsAsFactors=F)
      names(x5) <- keep6
      return(x5)
    } else {return(NULL)}
  } # End ordinal recode function
  
  #************************************************
  #****  Continuous variable recode function  *****
  #************************************************
  f.cont <- function(x, dep_var = "resp",
                     contvar,
                     prefix = "R1_",
                     binary_dv = "Y",
                     missrate = .75,
                     p_lo = .01,
                     p_hi = .99,
                     cap_flrC = "Y",
                     transformationC = "Y",
                     impmethodC = "MEAN",
                     stdmethodC = "STD",
                     min_size = 500) {
    f.catwrite<-function(...,f='Continuous Recodes.txt',app=T) cat(...,file=f,append=app)
    
    # Continuous Recodes
    cat("Starting number of continuous variables is",length(contvar),"\n")
    ## Basic statistics needed for transformations
    f.ckn <- function(z) {min(x[!is.na(z),dep_var])}
    f.ckx <- function(z) {max(x[!is.na(z),dep_var])}
    f.stats1 <- function(y){
      tmp <- data.frame(
        Name = contvar,
        Obs = apply(!is.na(y), 2, sum),
        MissingObs = apply(is.na(y), 2, sum),
        stringsAsFactors=F)
    }
    f.stats2 <- function(y){
      tmp <- data.frame(
        Name = contvar2,
        MinCk = apply(y,2,f.ckn),
        MaxCk = apply(y,2,f.ckx),
        LowerPct = apply(y, 2,quantile, probs=p_lo, na.rm=TRUE),
        UpperPct = apply(y, 2,quantile, probs=p_hi, na.rm=TRUE),
        Median = apply(y, 2,quantile, probs=0.50, na.rm=TRUE),
        stringsAsFactors=F)
    }
    cstats <- f.stats1(x[contvar])
    contvar2 <- as.character(cstats[cstats$Obs>0,"Name"])
    cstats <- merge(cstats,f.stats2(x[contvar2]),by="Name", all=TRUE)
    cstats$Usable <- ifelse(cstats$MissingObs / (cstats$Obs+cstats$MissingObs) < missrate &
                              cstats$LowerPct != cstats$UpperPct & 
                              cstats$MinCk != cstats$MaxCk,"Yes",NA)
    cstats$MinCk <- NULL
    cstats$MaxCk <- NULL
    keep2 <- as.character(cstats[!is.na(cstats$Usable),"Name"])  # Names of usable variables
    cat("Usable continuous variable count is",length(keep2),"\n")
    
    if (length(keep2) > 0) {
      if (toupper(cap_flrC)=="Y") {
        f.bounds<-function(y){
          quants<-quantile(y,c(0,.01,.25,.5,.75,.99,1),na.rm=T)
          iqr <- pmax((quants[5]-quants[3]),(quants[6]-quants[5]),(quants[3]-quants[2]))
          lb <- pmin(pmax((quants[3] - (1.5*iqr)), quants[1]), quants[2])
          ub <- pmin(pmax((quants[5] + (1.5*iqr)), quants[7]), quants[6])
          lb <- ifelse(lb == ub, quants[1], lb)
          ub <- ifelse(lb == ub, quants[7], ub)
          bnds <- cbind(lb,ub)
          return(bnds)
        }
        tmp <- as.data.frame(apply(x[keep2],2,f.bounds),stringsAsFactors=F)
        cstats <- merge(cstats,
                        data.frame("Name"=keep2,
                                   lb=as.numeric(tmp[1,]),
                                   ub=as.numeric(tmp[2,]),
                                   stringsAsFactors=F),
                        by="Name", all=TRUE)
      }
      
      # Determine missing imputation value: options ER, MEAN, MEDIAN, ZERO, MINIMUM, MIDRANGE
      if (toupper(impmethodC)=="ER"){
        lst1 <- cstats[cstats$MissingObs >= min_size & !is.na(cstats$Usable),"Name"]  # Use regression
        lst2 <- cstats[cstats$MissingObs < min_size & !is.na(cstats$Usable),"Name"]  # Use Median
        
        if (length(lst1) > 1) {
          f.miss_rr <- function(z){mean(x[is.na(z),dep_var])}
          miss_rr <- apply(x[,lst1],2,f.miss_rr)
        } else {
          if (length(lst1) == 1) {
            miss_rr <- mean(x[is.na(x[,lst1]),dep_var])
          } 
        }
        
        if (binary_dv == "Y" & length(lst1) > 0) {
          miss_rr <- ifelse(miss_rr==0,.0001,ifelse(miss_rr==1,.9999,miss_rr))
        }
        
        f.eval <- function(y){
          tmp <- data.frame(dv=dv,iv=y,stringsAsFactors=F)
          if (binary_dv == "Y") {
            mod1 <- glm(dv~iv, binomial,data=tmp)
          } else {
            mod1 <- lm(dv~iv, data=tmp)
          }
          ck <- data.frame(summary(mod1)["coefficients"],stringsAsFactors=F)
          val <- data.frame(Intercept=ck[1,1],Estimate=ck[2,1],Prob=ck[2,4],stringsAsFactors=F)
          return(val)
        }
        dv <- x[,dep_var]
        if (length(lst1) > 1) {
          coef <- data.frame(do.call(rbind,apply(x[,lst1],2,f.eval)),miss_rr=miss_rr,stringsAsFactors=F)
          coef$Name = row.names(coef)
        } else {
          if (length(lst1) == 1) {
            tmp <- data.frame(dv=dv,iv=x[,lst1],stringsAsFactors=F)
            if (binary_dv == "Y") {
              mod1 <- glm(dv~iv, binomial,data=tmp)
            } else {
              mod1 <- lm(dv~iv, data=tmp)
            }
            ck <- data.frame(summary(mod1)["coefficients"],stringsAsFactors=F)
            val <- data.frame(Intercept=ck[1,1],Estimate=ck[2,1],Prob=ck[2,4],stringsAsFactors=F)
            coef <- data.frame(Intercept=ck[1,1],Estimate=ck[2,1],Prob=ck[2,4],miss_rr=miss_rr,stringsAsFactors=F)
            coef$Name = lst1
          }
        }
        if (length(lst1) > 0) {
          er <- merge(coef,cstats[,c("Name","Median","lb","ub")],by="Name")
          
          if (binary_dv == "Y"){
            er$Miss_Impute <- ifelse(er$Prob > .05, er$Median,
                                     (log(er$miss_rr/(1-er$miss_rr))-er$Intercept)/er$Estimate)
          } else {
            er$Miss_Impute <- ifelse(er$Prob > .05, er$Median,
                                     (er$miss_rr-er$Intercept)/er$Estimate)
          }
          er$Miss_Impute <- ifelse(er$Miss_Impute < er$lb,er$lb,
                                   ifelse(er$Miss_Impute > er$ub,er$ub,er$Miss_Impute))
        }
        
        if (length(lst1) > 0 & length(lst2) > 0) {
          tmp <- cstats[cstats$Name %in% lst2, c("Name","Median")]
          names(tmp) <- c("Name","Miss_Impute")
          tmp2 <- rbind(tmp,er[,c("Name","Miss_Impute")])
          cstats <- merge(cstats, tmp2, by="Name", all=TRUE) 
        } else {
          if (length(lst1) > 0) {
            cstats <- merge(cstats, er[,c("Name","Miss_Impute")], by="Name", all=TRUE) 
            missval <- er$Miss_Impute
          } else {
            tmp <- data.frame(cstats[cstats$Name %in% lst2, c("Name","Median")])
            names(tmp) <- c("Name","Miss_Impute")
            cstats <- merge(cstats, tmp, by="Name", all=TRUE)
          }
        }
        missval <- cstats[cstats$Name %in% keep2,'Miss_Impute']
      } else {
        if (toupper(impmethodC)=="MEDIAN") {missval <- cstats[cstats$Name %in% keep2,"Median"]} else {
          if (toupper(impmethodC)=="ZERO") {missval <- rep.int(0,length(keep2))} else {
            if (toupper(impmethodC)=="MINIMUM") {missval <- apply(x[keep2],2,min,na.rm=T)} else {
              if (toupper(impmethodC)=="MIDRANGE") {missval <- (apply(x[keep2],2,min,na.rm=T)+apply(x[keep2],2,max,na.rm=T))/2} else {
                missval <- apply(x[keep2],2,mean,na.rm=T)  # Default is mean
              } } } }
        cstats <- merge(cstats,
                        data.frame("Name"=keep2,
                                   Miss_Impute=missval,
                                   stringsAsFactors=F),
                        by="Name", all=TRUE)
      }
      
      # Apply capping / flooring if option is chosen
      if (toupper(cap_flrC)=="Y") {
        f.cap <- function(y){
          quants<-quantile(y,c(0,.01,.25,.5,.75,.99,1),na.rm=T)
          iqr <- pmax((quants[5]-quants[3]),(quants[6]-quants[5]),(quants[3]-quants[2]))
          lb <- pmin(pmax((quants[3] - (1.5*iqr)), quants[1]), quants[2])
          ub <- pmin(pmax((quants[5] + (1.5*iqr)), quants[7]), quants[6])
          lb <- ifelse(lb == ub, quants[1], lb)
          ub <- ifelse(lb == ub, quants[7], ub)
          RW = pmin(pmax(y,lb),ub)
          return(RW)
        }
        x2<-data.frame(do.call(cbind,lapply(x[keep2],f.cap)),stringsAsFactors=F)
      } else {x2 <- x[,keep2]}
      
      # Replace missings
      if (length(keep2) > 1) {
        x2[is.na(x2)] <- matrix(missval,nrow(x2),length(keep2),byrow=T)[is.na(x2)]  
      } else {
        x2[is.na(x2[keep2]),keep2] <- missval
      }
      
      # Apply transforms if option is chosen
      if (toupper(transformationC)=="Y") {
        f.transform<-function(y,sqrt_lo=0,log_lo=0.0001){
          cbind(
            RW = y,
            SR = sqrt(pmax(y,sqrt_lo,na.rm=T)),
            SQ = y^2,
            LN = log(pmax(y,log_lo,na.rm=T)),
            EP = exp(pmin(y,100,na.rm=T)),
            IV = ifelse(y==0, 0, 1/y))
        }
        x3<-data.frame(do.call(cbind,lapply(x2,f.transform)),stringsAsFactors=F)
        keep3<-apply(expand.grid(colnames(f.transform(NA)),keep2)[,2:1],1,paste,collapse='.')
        keep3 <- paste(prefix,keep3,sep="")
        names(x3)<-keep3  
      } else {
        x3 <- x2
        keep3 <- paste(prefix,keep2,".RW",sep="")
        names(x3) <- keep3
      }
      
      # Apply standardization: MEAN,MEDIAN,SUM,EUCLEN,USTD,STD,RANGE,MIDRANGE,MAXABS,IQR,MAD,NO to skip
      if (toupper(stdmethodC) %in% c("MEAN","MEDIAN","SUM","EUCLEN","USTD","STD",
                                     "RANGE","MIDRANGE","MAXABS","IQR","MAD")){
        if (toupper(stdmethodC) == "MEAN"){
          location <- apply(x3,2,mean,na.rm=T)
          scale <- rep.int(1,length(keep3))
        }
        if (toupper(stdmethodC) == "MEDIAN"){
          location <- apply(x3,2,median,na.rm=T)
          scale <- rep.int(1,length(keep3))
        }
        if (toupper(stdmethodC) == "SUM"){
          location <- rep.int(0,length(keep3))
          scale <- apply(x3,2,sum,na.rm=T)
        }
        if (toupper(stdmethodC) == "EUCLEN"){
          f.euclen <- function(y){
            y2 <- sqrt(sum(y^2,na.rm=T))
          }
          location <- rep.int(0,length(keep3))
          scale <- apply(x3,2,f.euclen)
        }
        if (toupper(stdmethodC) == "USTD"){
          f.ustd <- function(y){
            y2 <- sqrt(sum((y-0)^2,na.rm=T)/sum(!is.na(y)))
          }
          location <- rep.int(0,length(keep3))
          scale <- apply(x3,2,f.ustd)
        }
        if (toupper(stdmethodC) == "STD"){
          location <- apply(x3,2,mean,na.rm=T)
          scale <- apply(x3,2,sd,na.rm=T)
        }
        if (toupper(stdmethodC) == "RANGE"){
          f.range <- function(y){
            y2 <- max(y,na.rm=T)-min(y,na.rm=T)
          }
          location <- apply(x3,2,min,na.rm=T)
          scale <- apply(x3,2,f.range)
        }
        if (toupper(stdmethodC) == "MIDRANGE"){
          f.midrangeL <- function(y){
            y2 <- (max(y,na.rm=T)+min(y,na.rm=T))/2
          }
          f.midrangeS <- function(y){
            y2 <- (max(y,na.rm=T)-min(y,na.rm=T))/2
          }
          location <- apply(x3,2,f.midrangeL)
          scale <- apply(x3,2,f.midrangeS)
        }
        if (toupper(stdmethodC) == "MAXABS"){
          f.maxabs <- function(y){
            y2 <- max(abs(y),na.rm=T)
          }
          location <- rep.int(0,length(keep3))
          scale <- apply(x3,2,f.maxabs)
        }
        if (toupper(stdmethodC) == "IQR"){
          f.iqr <- function(y){
            y2 <- quantile(y,probs=.75,na.rm=T)-quantile(y,probs=.25,na.rm=T)
          }
          location <- apply(x3,2,median,na.rm=T)
          scale <- apply(x3,2,f.iqr)
        }
        if (toupper(stdmethodC) == "MAD"){
          f.mad <- function(y){
            y2 <- median(abs(y-median(y,na.rm=T)),na.rm=T)
          }
          location <- apply(x3,2,median,na.rm=T)
          scale <- apply(x3,2,f.mad)
        }
        x4<-data.frame(scale(x3,center=location,scale=scale),stringsAsFactors=F)
      } else {x4 <- x3}
      
      ## Eliminate variables with 0 standard deviation
      ck <- round(apply(x4,2,sd),4)
      keep4 <- names(ck[ck>0])
      keep4 <- keep4[!is.na(keep4)]
      x4 <- x4[,keep4]
      
      ## Choose best form for each variable
      ck <- as.data.frame(cor(x4,x[,dep_var]),stringsAsFactors=F)
      names(ck) <- "Corr"
      ck$New_Var <- row.names(ck)
      ck$Name <- substr(ck$New_Var,nchar(prefix)+1,nchar(ck$New_Var)-3)
      ck$abscorr <- abs(ck$Corr)
      ck <- ck[order(ck$Name,ck$abscorr),]
      ck <- ck[!duplicated( ck[, "Name" ], fromLast=T),]
      ck$abscorr <- NULL
      keep5 <- ck$New_Var
      
      ## Finalize stats file
      if (toupper(stdmethodC) %in% c("MEAN","MEDIAN","SUM","EUCLEN","USTD","STD",
                                     "RANGE","MIDRANGE","MAXABS","IQR","MAD")){
        cstats <- merge(cstats,
                        merge(ck,
                              data.frame("New_Var"=keep3,
                                         Location=location,
                                         Scale=scale,
                                         stringsAsFactors=F),
                              by="New_Var"),
                        by = "Name", all=T)
      } else {
        cstats <- merge(cstats,ck, by="Name", all=T)
      }
      keep6 <- cstats$New_Var[!is.na(cstats$New_Var)]  
    } else {
      cstats$Usable <- NA
      keep6 <- NULL
    }
    
    save(cstats, file = "CE2_continuous_stats.RData",compress=T)  # Save stats file
    
    # Output statistics
    allRows <- seq(length = nrow(cstats)) + 1
    createSheet(ff, name = "Continuous Variables")
    writeWorksheet(ff, cstats, sheet = "Continuous Variables", startRow = 1, startCol = 1)
    setColumnWidth(ff, sheet = "Continuous Variables", column = 1, width = 8192)
    setColumnWidth(ff, sheet = "Continuous Variables", column = c(2,3,4,5,6,7,8,9,10,11,12), width = 3072)
    setCellStyle(ff, sheet = "Continuous Variables", row = allRows, col = 2, cellstyle = csComma)
    setCellStyle(ff, sheet = "Continuous Variables", row = allRows, col = 3, cellstyle = csComma)
    setCellStyle(ff, sheet = "Continuous Variables", row = allRows, col = 4, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Continuous Variables", row = allRows, col = 5, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Continuous Variables", row = allRows, col = 6, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Continuous Variables", row = allRows, col = 8, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Continuous Variables", row = allRows, col = 9, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Continuous Variables", row = allRows, col = 10, cellstyle = csDec4)
    setCellStyle(ff, sheet = "Continuous Variables", row = allRows, col = 11, cellstyle = csDec4)
    setCellStyle(ff, sheet = "Continuous Variables", row = allRows, col = 12, cellstyle = csDec4)
    
    #************************************************
    # Create code
    f.trans <- function(y) {
      tmp <- ifelse (y["type"] == 'SR', paste("x$",y["New_Var"]," <- sqrt(pmax(x$",prefix,y["Name"],".RW,0,na.rm=T))",sep=""),
                     ifelse (y["type"] == 'SQ', paste("x$",y["New_Var"]," <- x$",prefix,y["Name"],".RW^2",sep=""),
                             ifelse (y["type"] == 'LN', paste("x$",y["New_Var"]," <- log(pmax(x$",prefix,y["Name"],".RW,.0001,na.rm=T))",sep=""),
                                     ifelse(y["type"] == 'EP', paste("x$",y["New_Var"]," <- exp(pmin(x$",prefix,y["Name"],".RW,100,na.rm=T))",sep=""),
                                            ifelse(y["type"] == 'IV', paste("x$",y["New_Var"]," <- ifelse(x$",prefix,y["Name"],".RW==0,0,1/x$",prefix,y["Name"],".RW)",sep=""),
                                                   "")))))
    }
    
    f.write <- function(y){
      f.catwrite(paste("### Recode: ",y["Name"],'\n',sep="")) 
      f.catwrite(paste("x$",prefix,y["Name"],".RW <- x$",y["Name"],'\n',sep=""))                                           # Create new variable
      f.catwrite(paste("x$",prefix,y["Name"],".RW[is.na(x$",prefix,y["Name"],".RW)] <- ",y["Miss_Impute"],'\n',sep=""))                # Set missing values to median
      if (toupper(cap_flrC)=="Y") {
        f.catwrite(paste("x$",prefix,y["Name"],".RW <- pmin(pmax(x$",prefix,y["Name"],".RW,",y["lb"],"),",y["ub"],")",'\n',sep="")) # Set variable bounds  
      }
      if (toupper(transformationC)=="Y") {
        trans <- f.trans(y)
        f.catwrite(paste(trans,'\n',sep=""))
      }
      if (toupper(stdmethodC) != "NO") {
        f.catwrite(paste("x$",y["New_Var"]," <- (x$",y["New_Var"]," - ",y["Location"],") / ",y["Scale"],'\n',sep=""))         # transform  
      }
      f.catwrite('\n')
    }
    if (length(keep6) > 0) {
      cstats$type <- substr(cstats$New_Var,nchar(cstats$New_Var)-1,nchar(cstats$New_Var))
      f.catwrite('######################################################################\n',app=F)
      f.catwrite('#####             Continuous Variable Recodes                    #####\n')
      f.catwrite('######################################################################\n')
      f.catwrite('\n')
      f.catwrite('f.con_rec <- function(x){ \n')
      
      apply(cstats[cstats$New_Var %in% keep6,], 1, f.write)
      
      f.catwrite('return(x) \n')
      f.catwrite('} \n')
      
      # Write out list of variables to keep
      if (length(keep6) > 1) {
        temp <- data.frame(name=keep6,len=nchar(keep6),splt=NA,outvar=NA,stringsAsFactors=F)
        temp[1,"splt"] <- 1
        temp[1,"outvar"] <- paste("keep_c <- c('",temp[1,"name"],"'",sep="")
        cnt <- temp[1,"len"]
        ov <- temp[1,"outvar"]
        splt <- 1
        
        for (i in 2:nrow(temp)) {
          if (is.na(temp[i,"splt"])) {
            if (cnt + temp[i,"len"] > 250) {
              splt <- splt + 1
              cnt <- temp[i,"len"]
              ov <- paste("'",temp[i,"name"],"'",sep="")
            } else
            { ov <- paste(ov,",'",temp[i,"name"],"'",sep="")
              cnt <- cnt + temp[i,"len"] + 3
            }
            temp[i,"splt"] <- splt
            temp[i,"outvar"] <- ov
          }
        }
        # Mark first and last record for strings
        temp$slast <- !duplicated( temp[, "splt" ], fromLast=T)
        # Keep lest record for each string
        temp2 <- temp[temp$slast==T,]
        temp2$outvar <- ifelse(temp2$splt==splt,
                               paste(temp2$outvar,")",sep=""),
                               paste(temp2$outvar,",",sep=""))
        
        f.catwrite(paste(temp2$outvar,'\n',sep=""))
      } else {
        if (length(keep6)==1){
          f.catwrite(paste("keep_c <- c('",keep6,"')",'\n',sep=""))
        }
      }     
    }
    
    # Return recoded data
    if (length(keep6) > 0){
      x5 <- data.frame(x4[,keep6],stringsAsFactors=F)
      names(x5) <- keep6
      return(x5)
    } else {return(NULL)}
  }  # End continuous recode function
  
  #************************************************
  #************************************************
  
  # Initiate Excel workbook
  if (file.exists('CE2_EDA_report.xlsx')) file.remove('CE2_EDA_report.xlsx')
  ff <- loadWorkbook('CE2_EDA_report.xlsx', create = TRUE)
  
  # Create formats
  csComma = createCellStyle(ff, name = "comma")
  setDataFormat(csComma, format = "#,##0")
  
  csDec4 = createCellStyle(ff, name = "dec4")
  setDataFormat(csDec4, format = "0.0000")
  
  csDec2 = createCellStyle(ff, name = "dec2")
  setDataFormat(csDec2, format = "0.00")
  
  ## Clean up before starting
  if (file.exists('Binary Recodes.txt')) file.remove('Binary Recodes.txt')
  if (file.exists('Nominal Recodes.txt')) file.remove('Nominal Recodes.txt')
  if (file.exists('Ordinal Recodes.txt')) file.remove('Ordinal Recodes.txt')
  if (file.exists('Continuous Recodes.txt')) file.remove('Continuous Recodes.txt')
  if (file.exists('CE2_binary_stats.RData')) file.remove('CE2_binary_stats.RData')
  if (file.exists('CE2_nominal_stats.RData')) file.remove('CE2_nominal_stats.RData')
  if (file.exists('nominal_lkup.RData')) file.remove('nominal_lkup.RData')
  if (file.exists('CE2_ordinal_stats.RData')) file.remove('CE2_ordinal_stats.RData')
  if (file.exists('CE2_continuous_stats.RData')) file.remove('CE2_continuous_stats.RData')
  if (file.exists('CE2_vars.RData')) file.remove('CE2_vars.RData')
  
  # Determine maximum file size
  tst <- rxGetInfo(data=dat)
  mx <- tst$numRows*tst$numVars
  rm(tst)
  
  if (profiling =="Y"){
    # Global vectors
    dv <- rxDataStep(inData=dat,
                     rowSelection = mod_val_test != 3,
                     varsToKeep = dep_var,
                     maxRowsByCols = mx,
                     reportProgress=0)[,1]
    cnt <- length(dv)
    mn <- mean(dv)
    
    ## Load necessary functions
    # Summary functions
    f.f1 <- function(y){
      f1 <- as.data.frame(table(y))
      f1$Percent <- f1$Freq / cnt
      f2 <- as.data.frame(as.table(by(dv,y,mean)))
      f3 <- cbind(f1,f2)
      f3 <- f3[,-4]
      names(f3) <- c("Values","Count","Percent","MeanDV")
      f3$Index <- (f3[,4] / mn)*100
      ov <- data.frame(Values="Overall",Count=cnt,Percent=1,MeanDV=mn,Index=100,stringsAsFactors=F)
      f4 <- rbind(ov,f3)
      return(f4)
    }
    f.f2 <- function(y,num_category=10){
      f1 <- as.data.frame(table(y))
      f1$Percent <- f1$Freq / cnt
      f2 <- as.data.frame(as.table(by(dv,y,mean)))
      f3 <- cbind(f1,f2)
      f3 <- f3[,-4]
      names(f3) <- c("Category","Count","Percent","MeanDV")
      f3$Index <- (f3[,4] / mn)*100
      ov <- data.frame(Category="Overall",Count=cnt,Percent=1,MeanDV=mn,Index=100,stringsAsFactors=F)
      if (length(y[is.na(y)])>0){
        ms <- data.frame(Category="Missing",
                         Count=length(y[is.na(y)]),
                         MeanDV=mean(dv[is.na(y)]),
                         stringsAsFactors=F)
        ms$Percent <- ms$Count / cnt
        ms$Index <- (ms$MeanDV / mn)*100
        ms <- ms[,c("Category","Count","Percent","MeanDV","Index")]
        f4 <- rbind(ov,ms,f3)
      } else f4 <- rbind(ov,f3)
      
      return(f4)
    }
    f.f3 <- function(y,num_category=10){
      if (equal_dist == 'Y') {
        f1 <- as.data.frame(table(cut(y,b=num_category,include.lowest=T)))
        f2 <- as.data.frame(as.table(by(dv,cut(y,b=num_category,include.lowest=T),mean)))
      } else {
        f1 <- as.data.frame(table(cut(y,b=unique(quantile(y, seq(0, 1, 1/num_category),na.rm=T)),include.lowest=T)))
        f2 <- as.data.frame(as.table(by(dv,cut(y,b=unique(quantile(y, seq(0, 1, 1/num_category),na.rm=T)),include.lowest=T),mean)))
      }
      f1$Percent <- f1$Freq / cnt
      f3 <- cbind(f1,f2)
      f3 <- f3[,-4]
      names(f3) <- c("Category","Count","Percent","MeanDV")
      f3$Index <- (f3[,4] / mn)*100
      ov <- data.frame(Category="Overall",Count=cnt,Percent=1,MeanDV=mn,Index=100,stringsAsFactors=F)
      if (length(y[is.na(y)])>0){
        ms <- data.frame(Category="Missing",
                         Count=length(y[is.na(y)]),
                         MeanDV=mean(dv[is.na(y)]),
                         stringsAsFactors=F)
        ms$Percent <- ms$Count / cnt
        ms$Index <- (ms$MeanDV / mn)*100
        ms <- ms[,c("Category","Count","Percent","MeanDV","Index")]
        f4 <- rbind(ov,ms,f3)
      } else f4 <- rbind(ov,f3)
      return(f4)
    }
    # Binary
    f.prof1 <- function(x){
      tmp<-data.frame(do.call(rbind,lapply(x,f.f1)))
      rn <- row.names(tmp)
      tmp$Variable <- ifelse(substr(rn,nchar(rn)-1,nchar(rn)-1)==".",
                             substr(rn,1,nchar(rn)-2),
                             substr(rn,1,nchar(rn)-3))
      tmp$Name <- substr(tmp$Variable,nchar(prefix)+1,nchar(tmp$Variable))
      
      tmp <- merge(tmp,bstats[,c("Name","Type")],by="Name")
      tmp$Category <- ifelse(tmp$Values=="Overall","Overall",
                             ifelse(tmp$Type=="1",
                                    ifelse(tmp$Values==0,"Missing,0,N","1,Y"),
                                    ifelse(tmp$Values==0,"Missing,0","1")))
      
      tmp2 <- tmp[,c("Variable","Category","Count","Percent","MeanDV","Index")]
      tmp2$Star <- ifelse(tmp2$Index >= 110,'* (+)',
                          ifelse(tmp2$Index > 100,'  (+)',
                                 ifelse(tmp2$Index <= 90,'* (-)',
                                        ifelse(tmp2$Index < 100,'  (-)','  (0)'))))
      
      return(tmp2)
    }
    
    # Nominal
    f.prof2 <- function(x,meth){
      # Just keep recoded variables
      ck <- names(x)
      ck2 <- substr(ck,1,nchar(prefix))
      ck3 <- ck[ck2 == prefix]
      x <- x[,ck3]
      rm(ck,ck2,ck3)
      # Proceed with profiling
      if (length(names(x))> 1) {
        tmp<-data.frame(do.call(rbind,lapply(x,f.f1)))  
      } else {
        tmp <- f.f1(x)
      }
      
      rn <- row.names(tmp)
      tmp$Variable <- ifelse(substr(rn,nchar(rn)-1,nchar(rn)-1)==".",
                             substr(rn,1,nchar(rn)-2),
                             substr(rn,1,nchar(rn)-3))
      tmp$Name <- substr(tmp$Variable,nchar(prefix)+1,nchar(tmp$Variable))
      tmp$sord <- seq(1:nrow(tmp))
      
      load("nominal_lkup.RData")
      if (toupper(meth)=="BINARY") {
        lkup$Name <- paste0(lkup$variable,".",lkup$group)
        tmp$Variable <- ifelse(substr(tmp$Variable,nchar(tmp$Variable)-1,nchar(tmp$Variable)-1)==".",
                               substr(tmp$Variable,1,nchar(tmp$Variable)-2),
                               substr(tmp$Variable,1,nchar(tmp$Variable)-3))
        tmp$first <- !duplicated(tmp[, "Variable" ], fromLast=F)
        tmp1 <- tmp[tmp$Values==1 | tmp$first,]
        tmp2 <- merge(tmp1,lkup,by="Name",all.x=T)
      } else {
        tmp2 <- merge(tmp,lkup,by.x=c("Name","Values"),by.y=c("variable","val"),all.x=T)
        tmp2 <- tmp2[order(tmp2$sord),]
      }
      tmp2$Category <- ifelse(tmp2$Values=="Overall","Overall",tmp2$Category)
      tmp2 <- tmp2[,c("Variable","Category","Count","Percent","MeanDV","Index")]
      tmp2$Star <- ifelse(tmp2$Index >= 110,'* (+)',
                          ifelse(tmp2$Index > 100,'  (+)',
                                 ifelse(tmp2$Index <= 90,'* (-)',
                                        ifelse(tmp2$Index < 100,'  (-)','  (0)'))))
      
      return(tmp2)
    }
    
    # Ordinal
    f.prof3 <- function(x){
      if (length(ostats[!is.na(ostats$Usable)&ostats$UniqueCnt>num_category,"New_Var"]) > 1) {
        v <- ostats[!is.na(ostats$Usable)&ostats$UniqueCnt>num_category,"Name"]
        tst <- rxDataStep(inData=dat,
                          rowSelection = mod_val_test != 3,
                          varsToKeep = v,
                          maxRowsByCols = mx,
                          reportProgress=0)
        tst[!is.na(tst)] <- matrix(1,nrow(tst),length(tst),byrow=T)[!is.na(tst)]
        tst2 <- x[,ostats[!is.na(ostats$Usable)&ostats$UniqueCnt>num_category,"New_Var"]]*tst
        t1<-data.frame(do.call(rbind,lapply(tst2[ostats[!is.na(ostats$Usable)&ostats$UniqueCnt>num_category,"New_Var"]],f.f3,num_category)))
        rn <- row.names(t1)
        t1$Variable <- ifelse(substr(rn,nchar(rn)-1,nchar(rn)-1)==".",
                              substr(rn,nchar(prefix)+1,nchar(rn)-2),
                              substr(rn,nchar(prefix)+1,nchar(rn)-3))
        
      } else {
        if (length(ostats[!is.na(ostats$Usable)&ostats$UniqueCnt>num_category,"New_Var"]) == 1) {
          v <- ostats[!is.na(ostats$Usable)&ostats$UniqueCnt>num_category,"Name"]
          tst <- rxDataStep(inData=dat,
                            rowSelection = mod_val_test != 3,
                            varsToKeep = v,
                            maxRowsByCols = mx,
                            reportProgress=0)
          tst[!is.na(tst)] <- 1
          tst2 <- x[,ostats[!is.na(ostats$Usable)&ostats$UniqueCnt>num_category,"New_Var"]]*tst
          t1 <- data.frame(f.f3(tst2[,1],num_category))
          rn <- ostats[!is.na(ostats$Usable)&ostats$UniqueCnt>num_category,"New_Var"]
          t1$Variable <- ifelse(substr(rn,nchar(rn)-1,nchar(rn)-1)==".",
                                substr(rn,nchar(prefix)+1,nchar(rn)-2),
                                substr(rn,nchar(prefix)+1,nchar(rn)-3))
        }
      }
      
      if (length(ostats[!is.na(ostats$Usable)&ostats$UniqueCnt<=num_category,"New_Var"]) > 1) {
        v <- ostats[!is.na(ostats$Usable)&ostats$UniqueCnt<=num_category,"Name"]
        tst <- rxDataStep(inData=dat,
                          rowSelection = mod_val_test != 3,
                          varsToKeep = v,
                          maxRowsByCols = mx,
                          reportProgress=0)
        tst[!is.na(tst)] <- matrix(1,nrow(tst),length(tst),byrow=T)[!is.na(tst)]
        tst2 <- x[,ostats[!is.na(ostats$Usable)&ostats$UniqueCnt<=num_category,"New_Var"]]*tst
        t2<-data.frame(do.call(rbind,lapply(tst2[ostats[!is.na(ostats$Usable)&ostats$UniqueCnt<=num_category,"New_Var"]],f.f2,num_category)))
        rn <- row.names(t2)
        t2$Variable <- ifelse(substr(rn,nchar(rn)-1,nchar(rn)-1)==".",
                              substr(rn,nchar(prefix)+1,nchar(rn)-2),
                              substr(rn,nchar(prefix)+1,nchar(rn)-3))
        
      } else {
        if (length(ostats[!is.na(ostats$Usable)&ostats$UniqueCnt<=num_category,"New_Var"]) == 1) {
          v <- ostats[!is.na(ostats$Usable)&ostats$UniqueCnt<=num_category,"Name"]
          tst <- rxDataStep(inData=dat,
                            rowSelection = mod_val_test != 3,
                            varsToKeep = v,
                            maxRowsByCols = mx,
                            reportProgress=0)
          tst[!is.na(tst)] <- 1
          tst2 <- x[,ostats[!is.na(ostats$Usable)&ostats$UniqueCnt<=num_category,"New_Var"]]*tst
          t2 <- data.frame(f.f2(tst2,num_category))
          rn <- ostats[!is.na(ostats$Usable)&ostats$UniqueCnt<=num_category,"New_Var"]
          t2$Variable <- ifelse(substr(rn,nchar(rn)-1,nchar(rn)-1)==".",
                                substr(rn,nchar(prefix)+1,nchar(rn)-2),
                                substr(rn,nchar(prefix)+1,nchar(rn)-3))
        }
      }
      
      if (length(ostats[!is.na(ostats$Usable)&ostats$UniqueCnt<=num_category,"New_Var"]) > 0 &
            length(ostats[!is.na(ostats$Usable)&ostats$UniqueCnt>num_category,"New_Var"])) {
        tmp<-rbind(t1,t2)  
      } else {
        if (length(ostats[!is.na(ostats$Usable)&ostats$UniqueCnt<=num_category,"New_Var"]) > 0) {
          tmp<-t2
        } else {
          tmp<-t1
        }
      }
      
      tmp2 <- tmp[,c("Variable","Category","Count","Percent","MeanDV","Index")]
      tmp2$Star <- ifelse(tmp2$Index >= 110,'* (+)',
                          ifelse(tmp2$Index > 100,'  (+)',
                                 ifelse(tmp2$Index <= 90,'* (-)',
                                        ifelse(tmp2$Index < 100,'  (-)','  (0)'))))
      
      return(tmp2)
    }
    
    # Continuous
    f.prof4 <- function(x){
      v <- cstats[!is.na(cstats$Usable),"Name"]
      tst <- rxDataStep(inData=dat,
                        rowSelection = mod_val_test != 3,
                        varsToKeep = v,
                        maxRowsByCols = mx,
                        reportProgress=0)
      if (length(cstats[!is.na(cstats$Usable),"Name"]) > 1) {
        tst[!is.na(tst)] <- matrix(1,nrow(tst),length(tst),byrow=T)[!is.na(tst)]
        tst2 <- x*tst
        tmp<-data.frame(do.call(rbind,lapply(tst2[cstats[!is.na(cstats$Usable),"New_Var"]],f.f3,num_category)))
        rn <- row.names(tmp)
      } else {
        tst[!is.na(tst)] <- 1
        tst2 <- x[,cstats[!is.na(cstats$Usable),"New_Var"]]*tst
        tmp<-data.frame(f.f3(tst2[,1],num_category))
        rn <- cstats[!is.na(cstats$Usable),"New_Var"]
      }
      
      tmp$Variable <- ifelse(substr(rn,nchar(rn)-1,nchar(rn)-1)==".",
                             substr(rn,nchar(prefix)+1,nchar(rn)-2),
                             substr(rn,nchar(prefix)+1,nchar(rn)-3))
      tmp2 <- tmp[,c("Variable","Category","Count","Percent","MeanDV","Index")]
      tmp2$Star <- ifelse(tmp2$Index >= 110,'* (+)',
                          ifelse(tmp2$Index > 100,'  (+)',
                                 ifelse(tmp2$Index <= 90,'* (-)',
                                        ifelse(tmp2$Index < 100,'  (-)','  (0)'))))
      
      return(tmp2)
    }
  }
  
  # Initialize variable name file
  CE2_vars <- data.frame(Name=NA,New_Var=NA,stringsAsFactors=F)
  
  # Apply recode logic
  if (file.exists("binvar_list.txt") ){
    v <- readLines("binvar_list.txt")  # Read in list of binary variables
    v <- v[!is.na(v) & v!=""] # Make sure there are no extraneous items
    if (length(v) > 0) {
      bin_dat <- f.binary(rxDataStep(inData=dat,
                                     rowSelection = mod_val_test != 3,
                                     varsToKeep = v,
                                     maxRowsByCols = mx,
                                     reportProgress=0),
                          v,
                          prefix,
                          missrate,
                          concrate)
      
      load("CE2_binary_stats.RData")
      CE2_vars <- rbind(CE2_vars,bstats[!is.na(bstats$New_Var),c("Name","New_Var")])
      # Profiling
      if (profiling =="Y" & file.exists("Binary Recodes.txt") ){
        bin_prof <- f.prof1(bin_dat)
      }
      rm(v,bin_dat)
    } else {rm(v)}
  }
  
  if (file.exists("nomvar_list.txt") ){
    v <- readLines("nomvar_list.txt")  # Read in list of nominal variables
    v <- v[!is.na(v) & v!=""] # Make sure there are no extraneous items
    if (length(v) > 0) {
      nom_dat <- f.nominal(rxDataStep(inData=dat,
                                      rowSelection = mod_val_test != 3,
                                      varsToKeep = c(dep_var,v),
                                      maxRowsByCols = mx,
                                      reportProgress=0),
                           dep_var,
                           v,
                           prefix,
                           missrate,
                           concrate,
                           valcnt,
                           minbinnc,
                           minbinnp,
                           talpha,
                           nom_method)
      
      load("CE2_nominal_stats.RData")
      if (length(nstats[!is.na(nstats$Usable),"Name"])>0){
        load("nominal_lkup.RData")
        if (toupper(nom_method) == "BINARY") {lkup$New_Var <- paste0(prefix,lkup$variable,".",lkup$group)} else {
          lkup$New_Var <- paste0(prefix,lkup$variable)
        }
        nv <- lkup[,c("variable","New_Var")]
        names(nv) <- c("Name","New_Var")
        nv <- unique(nv)
        CE2_vars <- rbind(CE2_vars,nv)
        rm(lkup,nv)
      }
      
      # Profiling
      if (profiling =="Y" & file.exists("Nominal Recodes.txt") ){
        load("CE2_nominal_stats.RData")
        v <- nstats[!is.na(nstats$Usable),"Name"]  # Variables to keep
        if (length(v)>0){
          source("Nominal Profile Recodes.txt", local=T,echo = F)  # Load recode function
          nom_prof <- f.prof2(rxDataStep(inData=dat,
                                         rowSelection = mod_val_test != 3,
                                         varsToKeep = c(v),
                                         transformFunc = f.pnom_rec,
                                         transformVars = v,
                                         maxRowsByCols = mx,
                                         reportProgress=0),nom_method)
        }
      }
      rm(v,nom_dat)
    } else {rm(v)}
  }
  
  if (file.exists("ordvar_list.txt") ){
    v <- readLines("ordvar_list.txt")  # Read in list of ordinal variables
    v <- v[!is.na(v) & v!=""] # Make sure there are no extraneous items
    # Check for non-numeric variables
    tp <- sapply(rxGetVarInfo(dat,varsToKeep=v),getElement,"varType")
    bd <- v[tp!="numeric"]
    if (length(bd) > 0) {
      cat("The following non-numeric variables were dropped from the ordinal variable list:\n")
      cat(paste(bd, collapse = ","),"\n")
      v <- v[tp=="numeric"]
    }
    if (length(v) > 0) {
      ord_dat <- f.ordinal(rxDataStep(inData=dat,
                                      rowSelection = mod_val_test != 3,
                                      varsToKeep = c(dep_var,v),
                                      maxRowsByCols = mx,
                                      reportProgress=0),
                           dep_var,
                           v,
                           prefix,
                           binary_dv,
                           missrate,
                           concrate,
                           cap_flrO,
                           transformationO,
                           impmethodO,
                           stdmethodO,
                           min_size)
      
      load("CE2_ordinal_stats.RData")
      CE2_vars <- rbind(CE2_vars,ostats[!is.na(ostats$New_Var),c("Name","New_Var")])
      # Profiling
      if (profiling =="Y" & file.exists("Ordinal Recodes.txt") ){
        ord_prof <- f.prof3(ord_dat)
      }
      rm(v,ord_dat)
    } else {rm(v)}
  }
  
  if (file.exists("contvar_list.txt") ){
    v <- readLines("contvar_list.txt")  # Read in list of continuous variables
    v <- v[!is.na(v) & v!=""] # Make sure there are no extraneous items
    # Check for non-numeric variables
    tp <- sapply(rxGetVarInfo(dat,varsToKeep=v),getElement,"varType")
    bd <- v[tp!="numeric"]
    if (length(bd) > 0) {
      cat("The following non-numeric variables were dropped from the continuous variable list:\n")
      cat(paste(bd, collapse = ","),"\n")
      v <- v[tp=="numeric"]
    }
    if (length(v) > 0) {
      cont_dat <- f.cont(rxDataStep(inData=dat,
                                    rowSelection = mod_val_test != 3,
                                    varsToKeep = c(dep_var,v),
                                    maxRowsByCols = mx,
                                    reportProgress=0),
                         dep_var,
                         v,
                         prefix,
                         binary_dv,
                         missrate,
                         p_lo,
                         p_hi,
                         cap_flrC,
                         transformationC,
                         impmethodC,
                         stdmethodC,
                         min_size)
      
      load("CE2_continuous_stats.RData")
      CE2_vars <- rbind(CE2_vars,cstats[!is.na(cstats$New_Var),c("Name","New_Var")])
      # Profiling
      if (profiling =="Y" & file.exists("Continuous Recodes.txt") ){
        cont_prof <- f.prof4(cont_dat)
      }
      rm(v,cont_dat)
    } else {rm(v)}
  }
  if (profiling == "Y") {rm(dv,cnt,mn)}
  
  # Finalize variable file
  CE2_vars <- CE2_vars[!is.na(CE2_vars$Name),]    # Dedupe
  
  # Get variable descriptions
  tst <- rxGetVarInfo(data=dat)
  ds <- sapply(tst,getElement,"description")
  ds2 <- character(length=length(ds))
  for (i in 1:length(ds)) {
    ds2[i] <- if(is.null(ds[[i]])) {NA} else {ds[[i]]}
  }
  contents <- data.frame(Name=names(tst),
                         orig_label=ds2,
                         stringsAsFactors=F)
  
  CE2_vars <- merge(CE2_vars,contents,by="Name",all.x=T)
  names(CE2_vars) <- c("orig_var","Variable","orig_label")
  save(CE2_vars, file = "CE2_vars.RData",compress=T)  # Save stats file
  rm(CE2_vars,tst,ds,ds2,i,contents)
  
  # Save out Excel workbook
  saveWorkbook(ff)
  rm(ff,csComma,csDec2,csDec4)
  
  # Apply recodes to full file
  #CE2_Recoded <- dat[,keep_list]
  rxDataStep(inData = dat,
             outFile = "CE2_temp.xdf",
             overwrite=T,
             reportProgress=0)
  kvars <- keep_list
  
  if (file.exists("Binary Recodes.txt") ){
    cat("Apply binary \n")
    load("CE2_binary_stats.RData")
    v <- bstats[!is.na(bstats$Usable),"Name"]  # Variables to keep
    if (length(v)>0){
      source("Binary Recodes.txt", local=T,echo = F)  # Load recode function
      rxDataStep(inData = "CE2_temp.xdf",
                 outFile = "CE2_temp.xdf",
                 overwrite = T,
                 transformFunc = f.bin_rec,
                 transformVars = v,
                 reportProgress = 0)
      kvars <- c(kvars, keep_b)
    }
  }
  
  if (file.exists("Nominal Recodes.txt") ){
    cat("Apply nominal \n")
    load("CE2_nominal_stats.RData")
    v <- nstats[!is.na(nstats$Usable),"Name"]  # Variables to keep
    if (length(v)>0){
      source("Nominal Recodes.txt", local=T,echo = F)  # Load recode function
      rxDataStep(inData = "CE2_temp.xdf",
                 outFile = "CE2_temp.xdf",
                 overwrite = T,
                 transformFunc = f.nom_rec,
                 transformVars = v,
                 reportProgress = 0)
      kvars <- c(kvars, keep_n)
    }
  }
  
  if (file.exists("Ordinal Recodes.txt") ){
    cat("Apply ordinal \n")
    load("CE2_ordinal_stats.RData")
    v <- ostats[!is.na(ostats$Usable),"Name"]  # Variables to keep
    if (length(v)>0){
      source("Ordinal Recodes.txt", local=T,echo = F)  # Load recode function
      rxDataStep(inData = "CE2_temp.xdf",
                 outFile = "CE2_temp.xdf",
                 overwrite = T,
                 transformFunc = f.ord_rec,
                 transformVars = v,
                 reportProgress = 0)
      kvars <- c(kvars, keep_o)
    }
  }
  
  if (file.exists("Continuous Recodes.txt") ){
    cat("Apply continuous \n")
    load("CE2_continuous_stats.RData")
    v <- cstats[!is.na(cstats$Usable),"Name"]  # Variables to keep
    if (length(v)>0){
      source("Continuous Recodes.txt", local=T,echo = F)  # Load recode function
      rxDataStep(inData = "CE2_temp.xdf",
                 outFile = "CE2_temp.xdf",
                 overwrite = T,
                 transformFunc = f.con_rec,
                 transformVars = v,
                 reportProgress = 0)
      kvars <- c(kvars, keep_c)
    }
  }
  # Final dataset
  rxDataStep(inData = "CE2_temp.xdf",
             outFile = "CE2_Recoded.xdf",
             overwrite = T,
             #varsToKeep = kvars,
             xdfCompressionLevel = 1,
             reportProgress = 0)
  file.remove('CE2_temp.xdf')
  
  # Write out list of variables to keep
  f.catwrite<-function(...,f='CE2 Variables.txt',app=F) cat(...,file=f,append=app)
  keep_vars <- setdiff(kvars,keep_list)
  temp <- data.frame(name=keep_vars,len=nchar(keep_vars),splt=NA,outvar=NA,stringsAsFactors=F)
  temp[1,"splt"] <- 1
  temp[1,"outvar"] <- paste("keep_vars <- c('",temp[1,"name"],"'",sep="")
  cnt <- temp[1,"len"]
  ov <- temp[1,"outvar"]
  splt <- 1
  
  for (i in 2:nrow(temp)) {
    if (is.na(temp[i,"splt"])) {
      if (cnt + temp[i,"len"] > 250) {
        splt <- splt + 1
        cnt <- temp[i,"len"]
        ov <- paste("'",temp[i,"name"],"'",sep="")
      } else
      { ov <- paste(ov,",'",temp[i,"name"],"'",sep="")
        cnt <- cnt + temp[i,"len"] + 3
      }
      temp[i,"splt"] <- splt
      temp[i,"outvar"] <- ov
    }
  }
  # Mark first and last record for strings
  temp$slast <- !duplicated( temp[, "splt" ], fromLast=T)
  # Keep lest record for each string
  temp2 <- temp[temp$slast==T,]
  temp2$outvar <- ifelse(temp2$splt==splt,
                         paste(temp2$outvar,")",sep=""),
                         paste(temp2$outvar,",",sep=""))
  
  f.catwrite(paste(temp2$outvar,'\n',sep=""))
  rm(f.catwrite,temp,cnt,ov,splt,temp2,i)
  
  ## Correlation report
  cat("Run correlation \n")
  form <- as.formula(paste("~", paste(c(dep_var,keep_vars), collapse = "+")))
  cor <- data.frame(rxCor(form, data="CE2_Recoded.xdf", 
                          rowSelection=mod_val_test != 3,
                          reportProgress=0)[-1,1],
                    stringsAsFactors=F)
  names(cor) <- "Correlation"
  cor$Variable <- row.names(cor)
  cor <- cor[order(abs(cor$Correlation),decreasing=T),]
  cor <- cor[,c("Variable","Correlation")]
  rxDataStep(inData=cor,
             outFile="CE2_cor.xdf",
             overwrite = T,
             xdfCompressionLevel=1,
             reportProgress=0)
  rm(form,keep_vars)
  
  if (file.exists('CE2_Correlation_report.xlsx')) file.remove('CE2_Correlation_report.xlsx')
  cf <- loadWorkbook('CE2_Correlation_report.xlsx', create = TRUE)
  
  csDec4 = createCellStyle(cf, name = "dec4")
  setDataFormat(csDec4, format = "0.0000")
  
  allRows <- seq(length = nrow(cor)) + 1
  createSheet(cf, name = "Correlations")
  writeWorksheet(cf, cor, sheet = "Correlations", startRow = 1, startCol = 1)
  setColumnWidth(cf, sheet = "Correlations", column = 1, width = 8192)
  setColumnWidth(cf, sheet = "Correlations", column = 2, width = 3072)
  setCellStyle(cf, sheet = "Correlations", row = allRows, col = 2, cellstyle = csDec4)
  saveWorkbook(cf)
  rm(cf,csDec4,allRows,cor)
  
  
  #************************************************
  #****          Profiling function           *****
  #************************************************
  if (profiling =="Y"){
    
    cat("Initiate profiling report \n")
    # Initiate Excel workbook
    if (file.exists('CE2_Profile_report.xlsx')) file.remove('CE2_Profile_report.xlsx')
    pf <- loadWorkbook('CE2_Profile_report.xlsx', create = TRUE)
    
    # Create formats
    csPercent = createCellStyle(pf, name = "pct")
    setDataFormat(csPercent, format = "0.00%")
    
    csPercent1 = createCellStyle(pf, name = "pct1")
    setDataFormat(csPercent1, format = "0.0%")
    
    csComma = createCellStyle(pf, name = "comma")
    setDataFormat(csComma, format = "#,##0")
    
    csDec4 = createCellStyle(pf, name = "dec4")
    setDataFormat(csDec4, format = "0.0000")
    
    csDec2 = createCellStyle(pf, name = "dec2")
    setDataFormat(csDec2, format = "0.00")
    
    csLine = createCellStyle(pf, name = "line")
    setBorder(csLine, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")
    
    cslComma = createCellStyle(pf, name = "lcomma")
    setDataFormat(cslComma, format = "#,##0")
    setBorder(cslComma, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")
    
    cslPercent1 = createCellStyle(pf, name = "lpct1")
    setDataFormat(cslPercent1, format = "0.0%")
    setBorder(cslPercent1, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")
    
    cslDec2 = createCellStyle(pf, name = "ldec2")
    setDataFormat(cslDec2, format = "#,##0.00")
    setBorder(cslDec2, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")
    
    # Binary
    f.prt_prof1 <- function(tmp2){
      # Order file
      ord <- ddply(tmp2,"Variable",summarize,diff=max(abs(Index),na.rm=T)-min(abs(Index),na.rm=T))
      ord <- ord[order(ord$diff,decreasing=T),]
      ord$sord <- seq(1,nrow(ord))
      tmp2$s2 <- seq(1,nrow(tmp2))
      tmp2 <- merge(tmp2,ord[,c("Variable","sord")],by="Variable",all=T)
      tmp2 <- tmp2[order(tmp2$sord,tmp2$s2),]
      tmp2$sord <- NULL
      tmp2$s2 <- NULL
      ## Write out report
      allRows <- seq(length = nrow(tmp2)) + 1
      rows <- as.character(tmp2[,1])
      lrows <- c(NA,rows[1:(length(rows)-1)])
      chg <- which(rows != lrows)+1
      createSheet(pf, name = "Binary Variables")
      writeWorksheet(pf, tmp2, sheet = "Binary Variables", startRow = 1, startCol = 1)
      setColumnWidth(pf, sheet = "Binary Variables", column = c(1,2), width = 8192)
      setColumnWidth(pf, sheet = "Binary Variables", column = c(3,4,5,6,7), width = 3072)
      if (length(chg) > 0) {
        setCellStyle(pf, sheet = "Binary Variables", row = chg, col = 1, cellstyle = csLine)
        setCellStyle(pf, sheet = "Binary Variables", row = chg, col = 2, cellstyle = csLine)
      }
      setCellStyle(pf, sheet = "Binary Variables", row = allRows, col = 3, cellstyle = csComma)
      setCellStyle(pf, sheet = "Binary Variables", row = allRows, col = 4, cellstyle = csPercent1)
      setCellStyle(pf, sheet = "Binary Variables", row = allRows, col = 5, cellstyle = csDec2)
      setCellStyle(pf, sheet = "Binary Variables", row = allRows, col = 6, cellstyle = csDec2)
      if (length(chg) > 0) {
        setCellStyle(pf, sheet = "Binary Variables", row = chg, col = 3, cellstyle = cslComma)
        setCellStyle(pf, sheet = "Binary Variables", row = chg, col = 4, cellstyle = cslPercent1)
        setCellStyle(pf, sheet = "Binary Variables", row = chg, col = 5, cellstyle = cslDec2)
        setCellStyle(pf, sheet = "Binary Variables", row = chg, col = 6, cellstyle = cslDec2)
        setCellStyle(pf, sheet = "Binary Variables", row = chg, col = 7, cellstyle = csLine)  
      }
    }
    
    # Nominal
    f.prt_prof2 <- function(tmp2){
      # Order file
      ord <- ddply(tmp2,"Variable",summarize,diff=max(abs(Index),na.rm=T)-min(abs(Index),na.rm=T))
      ord <- ord[order(ord$diff,decreasing=T),]
      ord$sord <- seq(1,nrow(ord))
      tmp2$s2 <- seq(1,nrow(tmp2))
      tmp2 <- merge(tmp2,ord[,c("Variable","sord")],by="Variable",all=T)
      tmp2 <- tmp2[order(tmp2$sord,tmp2$s2),]
      tmp2$sord <- NULL
      tmp2$s2 <- NULL
      ## Write out report
      allRows <- seq(length = nrow(tmp2)) + 1
      rows <- as.character(tmp2[,1])
      lrows <- c(NA,rows[1:(length(rows)-1)])
      chg <- which(rows != lrows)+1
      createSheet(pf, name = "Nominal Variables")
      writeWorksheet(pf, tmp2, sheet = "Nominal Variables", startRow = 1, startCol = 1)
      setColumnWidth(pf, sheet = "Nominal Variables", column = c(1,2), width = 8192)
      setColumnWidth(pf, sheet = "Nominal Variables", column = c(3,4,5,6,7), width = 3072)
      if (length(chg) > 0) {
        setCellStyle(pf, sheet = "Nominal Variables", row = chg, col = 1, cellstyle = csLine)
        setCellStyle(pf, sheet = "Nominal Variables", row = chg, col = 2, cellstyle = csLine)  
      }
      setCellStyle(pf, sheet = "Nominal Variables", row = allRows, col = 3, cellstyle = csComma)
      setCellStyle(pf, sheet = "Nominal Variables", row = allRows, col = 4, cellstyle = csPercent1)
      setCellStyle(pf, sheet = "Nominal Variables", row = allRows, col = 5, cellstyle = csDec2)
      setCellStyle(pf, sheet = "Nominal Variables", row = allRows, col = 6, cellstyle = csDec2)
      if (length(chg) > 0) {
        setCellStyle(pf, sheet = "Nominal Variables", row = chg, col = 3, cellstyle = cslComma)
        setCellStyle(pf, sheet = "Nominal Variables", row = chg, col = 4, cellstyle = cslPercent1)
        setCellStyle(pf, sheet = "Nominal Variables", row = chg, col = 5, cellstyle = cslDec2)
        setCellStyle(pf, sheet = "Nominal Variables", row = chg, col = 6, cellstyle = cslDec2)
        setCellStyle(pf, sheet = "Nominal Variables", row = chg, col = 7, cellstyle = csLine)  
      }
    }
    
    # Ordinal
    f.prt_prof3 <- function(tmp2){
      # Order file
      ord <- ddply(tmp2,"Variable",summarize,diff=max(abs(Index),na.rm=T)-min(abs(Index),na.rm=T))
      ord <- ord[order(ord$diff,decreasing=T),]
      ord$sord <- seq(1,nrow(ord))
      tmp2$s2 <- seq(1,nrow(tmp2))
      tmp2 <- merge(tmp2,ord[,c("Variable","sord")],by="Variable",all=T)
      tmp2 <- tmp2[order(tmp2$sord,tmp2$s2),]
      tmp2$sord <- NULL
      tmp2$s2 <- NULL
      ## Write out report
      allRows <- seq(length = nrow(tmp2)) + 1
      rows <- as.character(tmp2[,1])
      lrows <- c(NA,rows[1:(length(rows)-1)])
      chg <- which(rows != lrows)+1
      createSheet(pf, name = "Ordinal Variables")
      writeWorksheet(pf, tmp2, sheet = "Ordinal Variables", startRow = 1, startCol = 1)
      setColumnWidth(pf, sheet = "Ordinal Variables", column = c(1,2), width = 8192)
      setColumnWidth(pf, sheet = "Ordinal Variables", column = c(3,4,5,6,7), width = 3072)
      if (length(chg) > 0) {
        setCellStyle(pf, sheet = "Ordinal Variables", row = chg, col = 1, cellstyle = csLine)
        setCellStyle(pf, sheet = "Ordinal Variables", row = chg, col = 2, cellstyle = csLine)  
      }
      setCellStyle(pf, sheet = "Ordinal Variables", row = allRows, col = 3, cellstyle = csComma)
      setCellStyle(pf, sheet = "Ordinal Variables", row = allRows, col = 4, cellstyle = csPercent1)
      setCellStyle(pf, sheet = "Ordinal Variables", row = allRows, col = 5, cellstyle = csDec2)
      setCellStyle(pf, sheet = "Ordinal Variables", row = allRows, col = 6, cellstyle = csDec2)
      if (length(chg) > 0) {
        setCellStyle(pf, sheet = "Ordinal Variables", row = chg, col = 3, cellstyle = cslComma)
        setCellStyle(pf, sheet = "Ordinal Variables", row = chg, col = 4, cellstyle = cslPercent1)
        setCellStyle(pf, sheet = "Ordinal Variables", row = chg, col = 5, cellstyle = cslDec2)
        setCellStyle(pf, sheet = "Ordinal Variables", row = chg, col = 6, cellstyle = cslDec2)
        setCellStyle(pf, sheet = "Ordinal Variables", row = chg, col = 7, cellstyle = csLine)
      }
    }
    
    # Continuous
    f.prt_prof4 <- function(tmp2){
      # Order file
      ord <- ddply(tmp2,"Variable",summarize,diff=max(abs(Index),na.rm=T)-min(abs(Index),na.rm=T))
      ord <- ord[order(ord$diff,decreasing=T),]
      ord$sord <- seq(1,nrow(ord))
      tmp2$s2 <- seq(1,nrow(tmp2))
      tmp2 <- merge(tmp2,ord[,c("Variable","sord")],by="Variable",all=T)
      tmp2 <- tmp2[order(tmp2$sord,tmp2$s2),]
      tmp2$sord <- NULL
      tmp2$s2 <- NULL
      ## Write out report
      allRows <- seq(length = nrow(tmp2)) + 1
      rows <- as.character(tmp2[,1])
      lrows <- c(NA,rows[1:(length(rows)-1)])
      chg <- which(rows != lrows)+1
      createSheet(pf, name = "Continuous Variables")
      writeWorksheet(pf, tmp2, sheet = "Continuous Variables", startRow = 1, startCol = 1)
      setColumnWidth(pf, sheet = "Continuous Variables", column = c(1,2), width = 8192)
      setColumnWidth(pf, sheet = "Continuous Variables", column = c(3,4,5,6,7), width = 3072)
      if (length(chg) > 0) {
        setCellStyle(pf, sheet = "Continuous Variables", row = chg, col = 1, cellstyle = csLine)
        setCellStyle(pf, sheet = "Continuous Variables", row = chg, col = 2, cellstyle = csLine)  
      }
      setCellStyle(pf, sheet = "Continuous Variables", row = allRows, col = 3, cellstyle = csComma)
      setCellStyle(pf, sheet = "Continuous Variables", row = allRows, col = 4, cellstyle = csPercent1)
      setCellStyle(pf, sheet = "Continuous Variables", row = allRows, col = 5, cellstyle = csDec2)
      setCellStyle(pf, sheet = "Continuous Variables", row = allRows, col = 6, cellstyle = csDec2)
      if (length(chg) > 0) {
        setCellStyle(pf, sheet = "Continuous Variables", row = chg, col = 3, cellstyle = cslComma)
        setCellStyle(pf, sheet = "Continuous Variables", row = chg, col = 4, cellstyle = cslPercent1)
        setCellStyle(pf, sheet = "Continuous Variables", row = chg, col = 5, cellstyle = cslDec2)
        setCellStyle(pf, sheet = "Continuous Variables", row = chg, col = 6, cellstyle = cslDec2)
        setCellStyle(pf, sheet = "Continuous Variables", row = chg, col = 7, cellstyle = csLine)  
      }
    }
    
    # Binary
    if (exists("bin_prof") ){
      f.prt_prof1(bin_prof)
    }
    
    # Nominal
    if (exists("nom_prof") ){
      f.prt_prof2(nom_prof)
    }
    
    # Ordinal
    if (exists("ord_prof") ){
      f.prt_prof3(ord_prof)
    }
    
    # Continuous
    if (exists("cont_prof") ){
      f.prt_prof4(cont_prof)
    }
    
    # Save out Excel workbook
    saveWorkbook(pf)
    rm(pf,csPercent,csPercent1,csComma,csDec4,csDec2,csLine,cslComma,cslPercent1,cslDec2)
  }
  # Clean up
  rm(list=ls())
}

# 3. VARIABLE REDUCTION #########################

## ****************************************************
## *****     3.Variable reduction and ranking     *****
## ****************************************************

## *** Macro 3: Variable Reduction macro variables ***;
## %let samplesize = 50000;              ** Sample size to use for variable reduction;
## %let redu_weight = N;                 ** Use weights in variables reduction? (Y/N);
## %let sources = 3;                     ** Minimum number of sources to be selected;
## * Univariate regression option;
## %let univ_reg = Y;                    ** Use univariate regression to choose variables? (Y/N);
## %let maxpuni = .0001;                 ** Maximum p value correlation for selecting via univariate regression;
## * Correlation option;
## %let correlation = Y;                 ** Use correlation to choose variables? (Y/N);
## %let corrcut = .01;                   ** Minimum correlation between independent variable and dependent variable;
## * Principal components option;
## %let principal = Y;                   ** Use principal components to choose variables? (Y/N);
## %let nprin = 10;                      ** Number of principal components desired;
## %let minprin = .5;                    ** Minimum factor correlation for selecting via principal component;
## * Cluster option;
## %let cluster = Y;                     ** Use cluster analysis to choose variables? (Y/N);
## %let maxc = 20;                       ** Number of clusters desired;
## %let maxratio = .5;                   ** Maximum R Squared ratio for selecting variables via clustering;
## * Linear regression option;
## %let regression = Y;                  ** Use linear regression to choose variables? (Y/N);
## %let alphareg = .05;                  ** Alpha level for forward selection in linear regression;
## * Logistic regression option - only applicable if binary dependent variable;
## %let logistic = Y;                    ** Use logistic regression to choose variables? (Y/N);
## %let alphalog = .05;                  ** Alpha level for forward selection in logistic regression;
## * Information value option;
## %let information = Y;                 ** Use information value to choose variables? (Y/N);
## %let decile = 20;                     ** Number of groups to use when calculate information values;
## %let infvcut = .01;                   ** Minimum information value for selecting via information value;
## * Maximum correlation between independent variables option;
## %let ind_correlation = Y;             ** Exclude variables with high correlation to others? (Y/N);
## %let maxcorr = .7;                    ** Maximum correlation allowed between independent variables;
## * Maximum correlation to dependent variable option;
## %let ind_dv_corr=N;                   ** Exclude variables with high correlation to dependent variable? (Y/N);
## %let max_dv_corr=.7;                  ** Maximum correlation allowed to dependent variable;

f.var_redu <- function(dat,
                       dep_var = "resp",
                       binary_dv = "Y",
                       samplesize = 50000,
                       redu_weight = "N",
                       weight = NULL,
                       sources = 3,
                       maxnum = 1000,
                       maxcorr = .7,
                       ind_dv_corr="N",
                       max_dv_corr=.7,
                       univ_reg = "Y",
                       maxpuni = .05,
                       correlation = "Y",
                       corrcut = .01,
                       factor = "Y",
                       nfact = 10,
                       minfact = .5,
                       regression = "Y",
                       alphareg = .05,
                       logistic = "Y",
                       alphalog = .05,
                       information = "Y",
                       decile = 20,
                       infvcut = .01) {
  
  # Load independent variable list
  source("CE2 Variables.txt", local=T,echo = F)
  keep_vars <- keep_vars[!is.na(keep_vars) & keep_vars!=""] # Make sure there are no extraneous items
  
  # Get basic metrics
  form <- as.formula(paste("~", paste(keep_vars,collapse="+"), collapse = ""))        
  sm <- rxSummary(form,
                  data=dat,
                  summaryStats = c("ValidObs","Min","Max","Mean"),
                  rowSelection = mod_val_test != 3,
                  reportProgress = 0)[["sDataFrame"]]
  names(sm) <- c("Variable","Mean","Minimum","Maximum","Count")
  rm(form)
  sm <- sm[,c("Variable","Count","Mean","Minimum","Maximum")]
  
  # Get number of records
  rws <- max(sm$Count)
  
  cat("Number of rows: ", rws, '\n')
  if (rws > samplesize) cat("Sample rate: ", round(samplesize / rws,4), '\n') else cat("No need to sample \n")
  
  # Need to sample?
  if (rws > samplesize) {
    smp <- samplesize / rws
    f.tmp <- function(lst) {
      lst[["sel"]] <- runif(.rxNumRows)<=fsmp
      return(lst)
    }
    rxDataStep(inData=dat,
               outFile="CE3_Sample.xdf",
               varsToKeep = c("mod_val_test",dep_var,weight,keep_vars),
               rowSelection = mod_val_test != 3 & sel,
               transformFunc = f.tmp,
               transformObjects = list(fsmp = smp),
               xdfCompressionLevel = 1,
               overwrite = T,
               reportProgress=0)
  } else {
    smp = 1
    rxDataStep(inData=dat,
               outFile="CE3_Sample.xdf",
               varsToKeep = c("mod_val_test",dep_var,weight,keep_vars),
               rowSelection = mod_val_test != 3,
               xdfCompressionLevel = 1,
               overwrite = T,
               reportProgress=0)
  }
  
  # Clean up correlation issues
  f.cleancorr <- function(dv,redu_weight,weight,maxcorr=.75,tol=.000001) {
    form <- as.formula(paste("~", paste(c(dv,keep_vars), collapse = "+")))
    if (toupper(redu_weight)=="Y") {
      cor1 <- rxCor(form, "CE3_Sample.xdf", pweights=weight, reportProgress=0)
    } else {
      cor1 <- rxCor(form, "CE3_Sample.xdf", reportProgress=0)  
    }
    cor2 <- cor1[order(abs(cor1[dv,]),decreasing=T),]  ## Sort by descending correlation to dependent
    cor3 <- cor2[-which(rownames(cor2)==dv),-which(colnames(cor2)==dv)]  ## Remove dependent variable
    rsums <- apply((cor3 >= maxcorr & cor3 < 1),1,sum)  ## Count by row cells with high correlations
    csums <- apply((cor3 >= maxcorr & cor3 < 1),2,sum)  ## Count by columns cells with high correlations
    cor4 <- cor3[rsums>0,csums>0]   ## Keep problems
    d <- dim(cor4)[2]               ## Number of columns in problem matrix
    if (d > 0) {
      stored.names <- colnames(cor4)  ## Capture names for later use
      colnames(cor4) <- 1:d           ## Replace column names with numebers
      cnt <- dim(cor4)[1]             ## Count of rows to control loop
      if (is.null(cnt)) cnt<-0
      drop <- rep(FALSE, d)           ## Set up list to track variable disposition
      while (cnt > 0){
        ds <- data.frame(colnum=as.numeric(colnames(cor4)),colval=cor4[1,])  ## get variable numbers and correlations for current variable
        dropvar <- ds[ds[,"colval"] >= maxcorr & ds[,"colval"] < 1,1]  ## get variable numbers for bad vars
        keepvar <- ds[ds[,"colval"] == 1,1]  ## get variable number for current variable
        drop[dropvar] <- TRUE  ## Set flag for bad variables
        drop[keepvar] <- NA    ## set current variable to NA
        cor4 <- cor4[(row.names(cor4) %in% stored.names[drop==F]),(!colnames(cor4) %in% c(dropvar,keepvar))]  ## Remove rows and columns tested bad
        cnt <- dim(cor4)[1]  ## Count of rows to control loop
        if (is.null(cnt)) cnt<-0
      }
      keep <- c(colnames(cor3[rsums==0,csums==0]),stored.names[is.na(drop)])  ## final list of variables
    } else {keep <- colnames(cor2)}
    trim <- trim.matrix(cor1[rownames(cor1) %in% keep,colnames(cor1) %in% keep],tol)
    keep2 <- row.names(as.data.frame(trim["trimmedmat"]))
    return(keep2)
  }
  cat("Clean up correlation problems \n")
  keep_vars2 <- setdiff(f.cleancorr(dv=dep_var,redu_weight,weight,maxcorr=maxcorr,tol=.000001),dep_var)
  sm$drop_corr <- ifelse(sm$Variable %in% keep_vars2, NA,"Y")
  
  cat("Initial number of variables: ", length(keep_vars), '\n')
  cat("Usable number of variables: ", length(keep_vars2), '\n')
  
  if (toupper(ind_dv_corr) == "Y") {
    d_cor <- rxDataStep(inData="CE2_cor.xdf", reportProgress=0)
    d_cor <- d_cor[d_cor$Variable %in% keep_vars2,]
    names(d_cor) <- c("Variable","dv_corr")
    d_cor$drop_dv_corr <- ifelse(abs(d_cor$dv_corr) > max_dv_corr, "Y",NA)
    keep_vars2 <- keep_vars2[keep_vars2 %in% d_cor[is.na(d_cor$drop_dv_corr),"Variable"]]
    sm <- merge(sm, d_cor, by="Variable", all=T)
  }
  
  vrs <- length(keep_vars2)
  set.seed <- 274923;
  keep_vars2 <- keep_vars2[order(runif(vrs))]
  if (vrs > 100 & (toupper(regression) =="Y" | toupper(logistic) == "Y")) {
    reg_vars <- data.frame(variable=keep_vars2,
                           grp=floor(cumsum(rep(1,times=vrs))/ceiling(vrs/ceiling(vrs/100)))+1,
                           stringsAsFactors=F)
    grp_num <- max(reg_vars$grp)
  }
  
  # Maximum dataframe size
  mx <- rws * vrs
  
  # Univariate regression function
  f.univ_reg <- function(var,x,dv,binary_dv,wgt,redu_weight=N){
    form <- as.formula(paste(dv, "~", var, collapse = ""))
    if (toupper(binary_dv) =="Y") {         # Binary dependent
      if (toupper(redu_weight) == "Y") {    # Use weights
        md <- rxLogit(form,
                      data = x,
                      pweights = wgt,
                      maxIterations = 5,
                      coeffTolerance = .1,
                      objectiveFunctionTolerance = 1,
                      reportProgress=0)
      } else {                              # Don't use weights
        md <- rxLogit(form,
                      data = x,
                      coeffTolerance = .1,
                      objectiveFunctionTolerance = 1,
                      reportProgress=0)
      }
    } else {                                # Continuous dependent
      if (toupper(redu_weight) == "Y") {    # Use weights
        md <- rxLinMod(form,
                       data = x,
                       pweights = wgt,
                       reportProgress=0)
        
      } else {                              # Don't use weights
        md <- rxLinMod(form,
                       data = x,
                       reportProgress=0)
      }
    }
    # Returned items
    rt <- md$coef.p.value[2]
    
    return(rt)
  } 
  
  # Multivariate linear regression function
  f.reg <- function(vars,form,redu_weight,weight){
    scp <- as.formula(paste("~", paste(vars, collapse = "+")))
    sel <- rxStepControl(method="forward",scope = scp,stepCriterion = "SigLevel",maxSigLevelToAdd = alphareg)
    if (toupper(redu_weight) == "Y") {
      md <- rxLinMod(form,
                     data = "CE3_Sample.xdf",
                     pweights = weight,
                     variableSelection = sel,
                     reportProgress=0)
    } else {
      md <- rxLinMod(form,
                     data = "CE3_Sample.xdf",
                     variableSelection = sel,
                     reportProgress=0)
    }
    rt <- row.names(md$coefficients)
    rt <- rt[rt != "(Intercept)"]
    return(rt)
  }
  
  # Multivariate logistic regression function
  f.log <- function(vars,form,redu_weight,weight){
    scp <- as.formula(paste("~", paste(vars, collapse = "+")))
    sel <- rxStepControl(method="forward",scope = scp,stepCriterion = "SigLevel",
                         maxSigLevelToAdd = alphalog,test="Chisq", refitEachStep = F)
    if (toupper(redu_weight) == "Y") {
      md <- rxLogit(form,
                    data = "CE3_Sample.xdf",
                    pweights = weight,
                    variableSelection = sel,
                    initialValues = NA,
                    coeffTolerance = .1,
                    reportProgress=0)
    } else {
      md <- rxLogit(form,
                    data = "CE3_Sample.xdf",
                    variableSelection = sel,
                    initialValues = NA,
                    coeffTolerance = .1,
                    reportProgress=0)
    }
    rt <- names(md$coefficients)
    rt <- rt[rt != "(Intercept)"]
    return(rt)
  }
  
  
  # Information value calculation function
  f.info_val <- function(var,x,dv,wgt,redu_weight=N,decile=10,set_size){
    if (toupper(redu_weight) == "Y") {
      tmp <- rxSort(inData=x,
                    sortByVars=var,
                    varsToKeep=c(dv,var,wgt),
                    maxRowsByCols=mx,
                    reportProgress=0)
      names(tmp) <- c("dv","iv","wgt")
      unq <- length(unique(tmp$iv))
      if (unq > decile) {
        tmp$bin <- floor(cumsum(tmp$wgt)*decile/(sum(tmp$wgt)+1))
        f1 <- ddply(tmp,"bin",summarize,
                    cnt=sum(wgt),
                    mean_dv=weighted.mean(dv,wgt),
                    max_dv=max(dv),
                    mean_var=weighted.mean(iv,wgt),
                    min_var=min(iv),
                    max_var=max(iv))
      } else {
        f1 <- ddply(tmp,"iv",summarize,
                    cnt=sum(wgt),
                    mean_dv=weighted.mean(dv,wgt),
                    max_dv=max(dv),
                    mean_var=weighted.mean(iv,wgt),
                    min_var=min(iv),
                    max_var=max(iv))
      }
    } else {
      tmp <- rxSort(inData=x,
                    sortByVars=var,
                    varsToKeep=c(dv,var),
                    maxRowsByCols=mx,
                    reportProgress=0)
      names(tmp) <- c("dv","iv")
      unq <- length(unique(tmp$iv))
      if (unq > decile) {
        trws <- nrow(tmp)
        tmp$bin <- floor(cumsum(rep(1,times=trws))*decile/(trws+1))
        f1 <- ddply(tmp,"bin",summarize,
                    cnt=length(iv),
                    mean_dv=mean(dv),
                    max_dv=max(dv),
                    mean_var=mean(iv),
                    min_var=min(iv),
                    max_var=max(iv))
      } else {
        f1 <- ddply(tmp,"iv",summarize,
                    cnt=length(iv),
                    mean_dv=mean(dv),
                    max_dv=max(dv),
                    mean_var=mean(iv),
                    min_var=min(iv),
                    max_var=max(iv))
      }
    }
    rm(tmp)
    totalresp <- sum(f1$cnt*f1$mean_dv)
    norm <- totalresp / sum(f1$cnt)
    ntotal <- sum(f1$cnt)
    maxdv <- max(f1$max_dv)
    
    cumpct_resp <- cumsum(f1$mean_dv*f1$cnt)/totalresp
    cumpct_freq <- cumsum(f1$cnt)/ntotal
    prev_dv = c(NA,cumpct_resp[1:length(cumpct_resp)-1])
    prev_pct = c(NA,cumpct_freq[1:length(cumpct_freq)-1])
    gini = 2*((cumpct_freq+prev_pct)/2-(cumpct_resp+prev_dv)/2)*(cumpct_freq-prev_pct)
    gini[1] <- 0
    gini <- cumsum(gini)
    
    badratio <- (f1$cnt * f1$mean_dv)/totalresp
    goodratio <- (f1$cnt - ((f1$mean_dv/maxdv) * f1$cnt ))/(ntotal-totalresp/maxdv)
    cumbad <- cumsum(badratio)
    cumgood <- cumsum(goodratio)
    
    infv  <- cumsum(ifelse(badratio+goodratio == 0,0,
                           ifelse(badratio*goodratio==0,f1$cnt*1.0/set_size*pmax(badratio,goodratio),
                                  (badratio-goodratio)*log(pmax(badratio/goodratio,0.00001)))))
    
    rt <- data.frame(infv = infv[length(infv)],
                     gini = gini[length(gini)],
                     ks = max(abs(cumbad-cumgood)),
                     stringsAsFactors = F)
    
    return(rt)
  }
  
  if (toupper(univ_reg) == "Y") {
    cat("Start univariate regression testing \n")
    PValue <- do.call(cbind,data.frame(sapply(keep_vars2,f.univ_reg,"CE3_Sample.xdf",
                                              dep_var,binary_dv,weight,redu_weight),
                                       stringsAsFactors=F))
    
    f1 <- data.frame(Variable = keep_vars2,
                     sign = ifelse(PValue[,1] < 0, "(-)", "(+)"),
                     PValue = PValue[,1],
                     stringsAsFactors = F)
    
    sm <- merge(sm, f1, by="Variable", all=T)
    cat("Number of variables passing univariate regression testing:", nrow(f1[f1$PValue <= maxpuni,]),'\n')
    rm(f1,"PValue")
  }
  
  if (toupper(correlation) == "Y") {
    cat("Start correlation testing \n")
    form <- as.formula(paste("~", paste(c(dep_var,keep_vars2), collapse = "+")))
    if (toupper(redu_weight) == "Y") {
      f1 <- data.frame(rxCor(form, data="CE3_Sample.xdf", 
                             pweights = weight,
                             reportProgress=0)[-1,1],
                       stringsAsFactors=F)
    } else {
      f1 <- data.frame(rxCor(form, data="CE3_Sample.xdf", 
                             reportProgress=0)[-1,1],
                       stringsAsFactors=F)
    }
    names(f1) <- "Corr"
    f1$Variable <- row.names(f1)
    sm <- merge(sm, f1, by="Variable", all=T)
    cat("Number of variables passing correlation testing:", nrow(f1[abs(f1$Corr) >= corrcut,]),'\n')
    rm(form,f1)
  }
  
  if (toupper(regression) == "Y"){
    cat("Start regression testing \n")
    form <- as.formula(paste(dep_var, "~ 1"))
    
    # Work with subsets of variables first if necessary
    if (vrs > 100){
      tvars <- f.reg(reg_vars[reg_vars$grp==1,"variable"],form,redu_weight,weight)
      for (i in 2:grp_num){
        tvars <- c(tvars,f.reg(reg_vars[reg_vars$grp==i,"variable"],form,redu_weight,weight))
      }
    } else {
      tvars <- keep_vars2
    }
    
    # Final model build
    scp <- as.formula(paste("~", paste(tvars, collapse = "+")))
    sel <- rxStepControl(method="forward",scope = scp,stepCriterion = "SigLevel",maxSigLevelToAdd = alphareg)
    if (toupper(redu_weight) == "Y") {
      md <- rxLinMod(form,
                     data = "CE3_Sample.xdf",
                     pweights = weight,
                     variableSelection = sel,
                     reportProgress=0)
    } else {
      md <- rxLinMod(form,
                     data = "CE3_Sample.xdf",
                     variableSelection = sel,
                     reportProgress=0)
    }
    f1 <- data.frame(md$coef.p.value,stringsAsFactors=F)
    names(f1) <- "RegPValue"
    f1$Variable <- row.names(f1)
    f1 <- f1[f1$Variable != "(Intercept)",]
    sm <- merge(sm, f1, by="Variable", all=T)
    cat("Number of variables passing regression testing:", nrow(f1),'\n')
    rm(form,scp,sel,tvars,md,f1)
  }
  
  if (toupper(logistic) == "Y" & toupper(binary_dv) == "Y"){
    cat("Start logistic testing \n")
    form <- as.formula(paste(dep_var, "~ 1"))
    
    # Work with subsets of variables first if necessary
    if (vrs > 100){
      tvars <- f.log(reg_vars[reg_vars$grp==1,"variable"],form,redu_weight,weight)
      for (i in 2:grp_num){
        tvars <- c(tvars,f.log(reg_vars[reg_vars$grp==i,"variable"],form,redu_weight,weight))
      }
    } else {
      tvars <- keep_vars2
    }
    
    # Final model build
    scp <- as.formula(paste("~", paste(tvars, collapse = "+")))
    sel <- rxStepControl(method="forward",scope = scp,stepCriterion = "SigLevel",
                         maxSigLevelToAdd = alphalog, test="Chisq", refitEachStep = F)
    
    if (toupper(redu_weight) == "Y") {
      md <- rxLogit(form,
                    data = "CE3_Sample.xdf",
                    pweights = weight,
                    variableSelection = sel,
                    initialValues = NA,
                    coeffTolerance = .1,
                    reportProgress=0)
    } else {
      md <- rxLogit(form,
                    data = "CE3_Sample.xdf",
                    variableSelection = sel,
                    initialValues = NA,
                    coeffTolerance = .1,
                    reportProgress=0)
    }
    f1 <- data.frame(md$coef.p.value,stringsAsFactors=F)
    names(f1) <- "LogPValue"
    f1$Variable <- names(md$coefficients)
    f1 <- f1[f1$Variable != "(Intercept)",]
    sm <- merge(sm, f1, by="Variable", all=T)
    cat("Number of variables passing logistic testing:", nrow(f1),'\n')
    rm(form,scp,sel,md,f1,tvars)
  }
  
  if (toupper(factor)=="Y"){
    cat("Start factor testing \n")
    form <- as.formula(paste("~", paste(keep_vars2, collapse = "+")))
    if (toupper(redu_weight) == "Y") {
      cv<-rxCor(form,
                "CE3_Sample.xdf",
                pweights=weight,
                reportProgress=0)
      nobs <- sum(rxDataStep(inData="CE3_Sample.xdf",
                             varsToKeep = "wgt",
                             reportProgress=0))
      f <- factanal(covmat=cv,factors=nfact,n.obs=nobs,start=rep(0,vrs))
    } else {
      cv<-rxCor(form,
                "CE3_Sample.xdf",
                reportProgress=0)
      f <- factanal(covmat=cv,factors=nfact,n.obs=rws,start=rep(0,vrs))
    }
    f1 <- data.frame(Variable=row.names(f$loadings),
                     factor=apply(f$loadings,1,function(x){max(abs(x))}),
                     stringsAsFactors=F)
    sm <- merge(sm, f1, by="Variable", all=T)
    cat("Number of variables passing factor testing:", nrow(f1[f1$factor>=minfact,]),'\n')
    rm(form,cv,f,f1)
  }
  
  if (toupper(information)=="Y") {
    cat("Start information value testing \n")
    f1 <- data.frame(sapply(keep_vars2,f.info_val,"CE3_Sample.xdf",
                            dep_var,weight,redu_weight,decile,length(keep_vars2)),
                     stringsAsFactors=F)
    f2 <- data.frame(Variable = names(f1),
                     infv = sapply(f1,getElement,"infv"),
                     gini = sapply(f1,getElement,"gini"),
                     ks = sapply(f1,getElement,"ks"),
                     stringsAsFactors = F)
    sm <- merge(sm, f2, by="Variable", all=T)
    cat("Number of variables passing information value testing:", nrow(f2[f2$infv>=infvcut,]),'\n')
    rm(f1,f2)
  }
  
  # Create flags and source count
  cat("Create flags \n")
  sm$num_sources <- 0
  if (toupper(univ_reg) == "Y") {
    sm$univsource <- as.numeric(sm$PValue <= maxpuni)
    sm$num_sources <- sm$num_sources + sm$univsource
  }
  if (toupper(correlation) == "Y") {
    sm$corrsource <- as.numeric(abs(sm$Corr) >= corrcut)
    sm$num_sources <- sm$num_sources + sm$corrsource
  }
  if (toupper(regression) == "Y"){
    sm$regsource <- as.numeric(!(is.na(sm$RegPValue)))
    sm$num_sources <- sm$num_sources + sm$regsource
  }
  if (toupper(logistic) == "Y" & toupper(binary_dv) == "Y"){
    sm$logsource <- as.numeric(!(is.na(sm$LogPValue)))
    sm$num_sources <- sm$num_sources + sm$logsource
  }
  if (toupper(factor) == "Y") {
    sm$factsource = as.numeric(sm$factor>=minfact)
    sm$num_sources <- sm$num_sources + sm$factsource
  }
  if (toupper(information) == "Y") {
    sm$infvsource = as.numeric(sm$infv>=infvcut )
    sm$num_sources <- sm$num_sources + sm$infvsource
  }
  
  # Sort dataset
  cat("Sort dataset \n")
  if (toupper(information)=="Y"){
    sm <- sm[order(sm$num_sources,sm$infv,decreasing=T),]
  } else {
    sm <- sm[order(sm$num_sources,decreasing=T),]
  }
  
  # Save dataset
  cat("Save dataset \n")
  CE3_var_redu <- sm
  save(CE3_var_redu, file = "CE3_var_redu.RData",compress=T)
  rm(CE3_var_redu)
  
  # Initiate Excel workbook
  cat("Create Excel report \n")
  if (file.exists('CE3_Var_Redu_Results.xlsx')) file.remove('CE3_Var_Redu_Results.xlsx')
  vr <- loadWorkbook('CE3_Var_Redu_Results.xlsx', create = TRUE)
  
  ## Write out report
  #allRows <- seq(length = nrow(sm)) + 1
  createSheet(vr, name = "Variables")
  writeWorksheet(vr, sm, sheet = "Variables", startRow = 1, startCol = 1)
  setColumnWidth(vr, sheet = "Variables", column = c(1), width = 8192)
  
  saveWorkbook(vr)
  
  # Write out list of variables to keep
  cat("Write out variables to keep \n")
  f.catwrite<-function(...,f='CE3_Varlist_redu.txt',app=F) cat(...,file=f,append=app)
  
  if (toupper(ind_dv_corr)=="Y") {
    keep_vars3 <- sm[is.na(sm$drop_dv_corr) & !is.na(sm$num_sources) & sm$num_sources >= sources,"Variable"]
  } else {
    keep_vars3 <- sm[!is.na(sm$num_sources) & sm$num_sources >= sources,"Variable"]
  }
  #if (toupper(logistic)=="Y" & toupper(binary_dv)=="Y") {
  #  keep_vars3 <- unique(keep_vars3,sm[sm$logsource == 1,"Variable"])
  #} 
  #if (toupper(regression)=="Y") {
  #  keep_vars3 <- unique(keep_vars3,sm[sm$regsource == 1,"Variable"])
  #} 
  # Restrict to maximum number of variables
  if (length(keep_vars3)>maxnum) {keep_vars3 <- keep_vars3[1:maxnum]}
  
  if (length(keep_vars3)> 1){
    temp <- data.frame(name=keep_vars3,len=nchar(keep_vars3),splt=NA,outvar=NA,stringsAsFactors=F)
    temp[1,"splt"] <- 1
    temp[1,"outvar"] <- paste("varlist_redu <- c('",temp[1,"name"],"'",sep="")
    cnt <- temp[1,"len"]
    ov <- temp[1,"outvar"]
    splt <- 1
    
    for (i in 2:nrow(temp)) {
      if (is.na(temp[i,"splt"])) {
        if (cnt + temp[i,"len"] > 250) {
          splt <- splt + 1
          cnt <- temp[i,"len"]
          ov <- paste("'",temp[i,"name"],"'",sep="")
        } else
        { ov <- paste(ov,",'",temp[i,"name"],"'",sep="")
          cnt <- cnt + temp[i,"len"] + 3
        }
        temp[i,"splt"] <- splt
        temp[i,"outvar"] <- ov
      }
    }
    # Mark first and last record for strings
    temp$slast <- !duplicated( temp[, "splt" ], fromLast=T)
    # Keep lest record for each string
    temp2 <- temp[temp$slast==T,]
    temp2$outvar <- ifelse(temp2$splt==splt,
                           paste(temp2$outvar,")",sep=""),
                           paste(temp2$outvar,",",sep=""))
    
    f.catwrite(paste(temp2$outvar,'\n',sep=""))
    rm(f.catwrite,temp,cnt,ov,splt,temp2,i)
    cat(length(keep_vars3)," variables selected \n")
  } else {
    if (length(keep_vars3)==1){
      f.catwrite(paste0("varlist_redu <- c('",keep_vars3,"')\n"))  
      cat("1 variable selected \n")
    } else {
      cat("No variables selected \n")
    }
  }
  # Clean up
  rm(list=ls())
}

# 4. MODEL SELECTION AND TUNING #################

## ************************************************
## *****     4.Model selection and tuning     *****
## ************************************************

## *** Macro 4: Model selection and tuning ***;
## %let sel_alpha = .05;                 ** Alpha level when selecting/removing variable into the model;
## %let includelist = ;                  ** List of variables that has to be in the final model;
## %let startlist = ;                    ** List of variables that need to be in the first step;
## %let excludelist = ;                  ** List of variables to be excluded from modeling;
## %let criteria = c;                    ** Metric to use during variable tuning. Default is c for binary, AdjRsq for other;
## %let threshold = 0;                   ** Minimum change in evaluation metric to include variable;
## %let SQL_join = union;                ** Type of join to use between file portions during variable tuning;
## %let minimp = .01;                    ** Minimum relative importance to keep variable;
## %let graph_plot = Y;                  ** Include graphing? (Y/N);

f.model_val <- function(insdn, 
                        varlist,
                        dep_var,
                        binary_dv,
                        weight = NULL,
                        sel_alpha = .05,
                        refit = F,
                        includelist = NA,
                        startlist = NA,
                        excludelist = NA,
                        criteria = "c",
                        threshold = 0,
                        SQL_join = "union",
                        minimp = .01,
                        graph_plot = "Y") {
  
  varlist <- varlist[!is.na(varlist) & varlist!=""] # Make sure there are no extraneous items
  
  # Master function for binary dependent variable
  f.Model_Val_Logistic <- function(insdn,
                                   regvlist, 
                                   includelist = NULL,
                                   startlist = NULL,
                                   excludelist = NULL,
                                   sle = 0.05,
                                   sls = 0.05){
    
    cat('Get variable list from forward & backward select \n')
    origlst <- regvlist
    
    # Remove excluded variables from lists
    if (length(excludelist) > 0) {
      regvlist <- setdiff(regvlist,excludelist)
      if (length(startlist) > 0) {startlist <- setdiff(startlist,excludelist)}
      if (length(includelist) > 0) {includelist <- setdiff(includelist,excludelist)}
    }
    
    # Build base variable list
    if (length(startlist) > 0 & length(includelist) > 0) {basevars <- unique(startlist,includelist)}
    if (length(startlist) > 0 & length(includelist) == 0) {basevars <- unique(startlist)}
    if (length(startlist) == 0 & length(includelist) > 0) {basevars <- unique(includelist)}
    if (length(startlist) == 0 & length(includelist) == 0) {basevars <- 1}
    
    # Forward selection
    form <- as.formula(paste(dep_var, "~", paste(basevars, collapse = "+")))
    if (length(regvlist) > 1){
      scp <- as.formula(paste("~", paste(regvlist, collapse = "+")))  
    } else {
      scp <- as.formula(paste("~", regvlist))
    }
    
    stp <- rxStepControl(method="forward",scope = scp,stepCriterion = "SigLevel",
                         maxSigLevelToAdd = sle,test="Chisq", refitEachStep = refit)
    # Model
    cat("Forward - Model \n")
    md <- rxLogit(form,
                  data = insdn,
                  rowSelection = mod_val_test == 1,
                  variableSelection = stp,
                  initialValues = NA,
                  coeffTolerance = .1,
                  reportProgress=0)
    Tforward <- names(md$coefficients)
    Tforward <- Tforward[Tforward != "(Intercept)"]
    rm(md)
    
    # Validation
    cat("Forward - Validation \n")
    md <- rxLogit(form,
                  data = insdn,
                  rowSelection = mod_val_test == 2,
                  variableSelection = stp,
                  initialValues = NA,
                  coeffTolerance = .1,
                  reportProgress=0)
    Vforward <- names(md$coefficients)
    Vforward <- Vforward[Vforward != "(Intercept)"]
    rm(md)
    
    # Backward selection
    
    if (length(regvlist) > 1){
      form <- as.formula(paste(dep_var, "~", paste(regvlist, collapse = "+")))
    } else {
      form <- as.formula(paste(dep_var, "~", regvlist))
    }
    scp <- as.formula(paste("~", paste(basevars, collapse = "+")))
    stp <- rxStepControl(method="backward",scope = list(lower=scp),stepCriterion = "SigLevel",
                         minSigLevelToDrop = sls,test="Chisq", refitEachStep = refit)
    # Model
    cat("Backward - Model \n")
    md <- rxLogit(form,
                  data = insdn,
                  rowSelection = mod_val_test == 1,
                  variableSelection = stp,
                  initialValues = NA,
                  coeffTolerance = .1,
                  reportProgress=0)
    Tbackward <- names(md$coefficients)
    Tbackward <- Tbackward[Tbackward != "(Intercept)"]
    rm(md)
    
    # Validation
    cat("Backward - Validation \n")
    md <- rxLogit(form,
                  data = insdn,
                  rowSelection = mod_val_test == 2,
                  variableSelection = stp,
                  initialValues = NA,
                  coeffTolerance = .1,
                  reportProgress=0)
    Vbackward <- names(md$coefficients)
    Vbackward <- Vbackward[Vbackward != "(Intercept)"]
    rm(md) 
    
    if (max(length(Tforward),length(Tbackward),length(Vforward),length(Vbackward)) == 0) {
      stop("No variables selected")
    } else {
      # Get statistics - model
      CE4_Model_Metric_Mod <- rbind(f.ctl_Stats_Log(insdn, 1, dep_var, union(Tforward,Tbackward),"Dev"),
                                    f.ctl_Stats_Log(insdn, 2, dep_var, union(Tforward,Tbackward),"Val"))
      save(CE4_Model_Metric_Mod, file = "CE4_Model_Metric_Mod.RData",compress=T)
      
      # Get statistics - validation
      CE4_Model_Metric_Val <- rbind(f.ctl_Stats_Log(insdn, 2, dep_var, union(Vforward,Vbackward),"Dev"),
                                    f.ctl_Stats_Log(insdn, 1, dep_var, union(Vforward,Vbackward),"Val"))
      save(CE4_Model_Metric_Val, file = "CE4_Model_Metric_Val.RData",compress=T)
      
      CE4_Variables <- data.frame(variable = origlst,
                                  Tforward = ifelse(origlst %in% Tforward,"Y",NA),
                                  Tbackward = ifelse(origlst %in% Tbackward,"Y",NA),
                                  Vforward = ifelse(origlst %in% Vforward,"Y",NA),
                                  Vbackward = ifelse(origlst %in% Vbackward,"Y",NA),
                                  stringsAsFactors=F)
      save(CE4_Variables, file = "CE4_Variables.RData",compress=T)
    }
    
  }
  
  # Binary model control
  f.ctl_Stats_Log <- function(insdn, part, dep_var, varlist, file) {
    
    cat('Get statistics for mod_val_test = ',part,'\n')
    
    # Separate data needed into dataframe
    f.tmp <- function(lst) {
      lst$sel <- lst$mod_val_test == ft
      return(lst)
    }
    rxDataStep(inData = insdn,
               outFile = "CE4_temp.xdf",
               overwrite = T,
               rowSelection = sel,
               varsToKeep  = c(dep_var,weight,varlist),
               transformFunc = f.tmp,
               transformObjects = list(ft = part),
               transformVars = "mod_val_test",
               reportProgress=0)
    
    # Full model
    cat("Full model \n")
    null.deviance <- rxLogit(formula = as.formula(paste(dep_var,"~1")),
                             data = "CE4_temp.xdf",
                             pweights = weight,
                             initialValues = NA,
                             coeffTolerance = .1,
                             reportProgress=0)$deviance
    frm <- as.formula(paste(dep_var, "~", paste(varlist, collapse = "+")))
    md <- rxLogit(formula = frm, data = "CE4_temp.xdf", pweights = weight, 
                  initialValues = NA, coeffTolerance = .1, reportProgress=0)
    # Number of valid observations for calculations
    rws <- md$nValidObs
    # Predicted values
    rxPredict(md,"CE4_temp.xdf",type="response",predVarNames="pred",overwrite=T,reportProgress=0)
    prd <- rxDataStep("CE4_temp.xdf",varsToKeep=c(dep_var,"pred",weight),reportProgress=0)
    
    full_stats <- data.frame(file = file,
                             variable = "Full File",
                             Estimate = NA,
                             StdEst = NA,
                             StdErr = NA,
                             zValue = NA,
                             zProb = NA,
                             VIF = NA,
                             RelImp = NA,
                             AIC = md$aic,
                             SC = md$deviance + (length(names(md$coefficients))*log(rws)),
                             LogL2 = md$deviance,
                             Rsquare = 1 - exp((-2/nrow(prd))*(md$deviance/-2 - null.deviance/-2)),
                             f.stats_log(prd,weight),
                             sord = 1,
                             stringsAsFactors=F)
    
    coef <- as.data.frame(coefficients(summary(md)))
    se_var <- coef[,"Estimate"] / ((pi/sqrt(3))/c(NA,apply(rxDataStep("CE4_temp.xdf",varsToKeep=varlist,reportProgress=0),2,sd)))
    # Variance inflation factors
    vif <- vif(lm(frm,
                  rxDataStep("CE4_temp.xdf",
                             transforms=list(vif_wgt=pred*(1-pred)),
                             reportProgress=0),
                  weights=vif_wgt))
    
    var_stats <- data.frame(file = file,
                            variable = names(md$coefficients),
                            Estimate = coef[,"Estimate"],
                            StdEst = se_var,
                            StdErr = coef[,"Std. Error"],
                            zValue = coef[,3],
                            zProb = coef[,4],
                            VIF = c(NA,vif),
                            RelImp = abs(se_var) / sum(abs(se_var),na.rm=T),
                            stringsAsFactors=F)
    
    # No intercept model
    cat("No intercept model \n")
    frm <- as.formula(paste(dep_var, "~", paste(c(-1,varlist), collapse = "+")))
    md <- rxLogit(formula = frm, data = "CE4_temp.xdf", pweights = weight, 
                  initialValues = NA, coeffTolerance = .1, reportProgress=0)
    # Predicted values
    rxPredict(md,"CE4_temp.xdf",type="response",predVarNames="pred",overwrite=T,reportProgress=0)
    prd <- rxDataStep("CE4_temp.xdf",varsToKeep=c(dep_var,"pred",weight),reportProgress=0)
    tmp_stats <- data.frame(variable = "(Intercept)",
                            AIC = md$aic,
                            SC = md$deviance + (length(names(md$coefficients))*log(rws)),
                            LogL2 = md$deviance,
                            Rsquare = 1 - exp((-2/nrow(prd))*(md$deviance/-2 - null.deviance/-2)),
                            f.stats_log(prd,weight),
                            sord = 2,
                            stringsAsFactors=F)
    
    # Models removing one variable at a time
    cat("Subset models - remove one variable at a time \n")
    for (i in 1:length(varlist)) {
      frm <- as.formula(paste(dep_var, "~", paste(varlist[-i], collapse = "+")))
      md <- rxLogit(formula = frm, data = "CE4_temp.xdf", pweights = weight,
                    initialValues = NA, coeffTolerance = .1, reportProgress=0)
      # Predicted values
      rxPredict(md,"CE4_temp.xdf",type="response",predVarNames="pred",overwrite=T,reportProgress=0)
      prd <- rxDataStep("CE4_temp.xdf",varsToKeep=c(dep_var,"pred",weight),reportProgress=0)
      tmp_stats <- rbind(tmp_stats,
                         data.frame(variable = varlist[i],
                                    AIC = md$aic,
                                    SC = md$deviance + (length(names(md$coefficients))*log(rws)),
                                    LogL2 = md$deviance,
                                    Rsquare = 1 - exp((-2/nrow(prd))*(md$deviance/-2 - null.deviance/-2)),
                                    f.stats_log(prd,weight),
                                    sord = 3,
                                    stringsAsFactors=F))
    }
    var_stats <- merge(var_stats,tmp_stats,by="variable")
    
    # Subtract full model metrics from sub-model metrics
    # This shows the impact of each variable
    fin_stats <- rbind(full_stats,
                       data.frame(file = var_stats$file,
                                  variable = var_stats$variable,
                                  Estimate = var_stats$Estimate,
                                  StdEst = var_stats$StdEst,
                                  StdErr = var_stats$StdErr,
                                  zValue = var_stats$zValue,
                                  zProb = var_stats$zProb,
                                  VIF = var_stats$VIF,
                                  RelImp = var_stats$RelImp,
                                  AIC = var_stats$AIC - full_stats$AIC,
                                  SC = var_stats$SC - full_stats$SC,
                                  LogL2 = var_stats$LogL2 - full_stats$LogL2,
                                  Rsquare = var_stats$Rsquare - full_stats$Rsquare,
                                  SomersD = var_stats$SomersD - full_stats$SomersD,
                                  Gamma = var_stats$Gamma - full_stats$Gamma,
                                  TauA = var_stats$TauA - full_stats$TauA,
                                  c = var_stats$c - full_stats$c,
                                  Concord = var_stats$Concord - full_stats$Concord,
                                  Discon = var_stats$Discon - full_stats$Discon,
                                  LackFit = var_stats$LackFit - full_stats$LackFit,
                                  Lift_Index = var_stats$Lift_Index - full_stats$Lift_Index,
                                  infv = var_stats$infv - full_stats$infv,
                                  ks = var_stats$ks - full_stats$ks,
                                  sord = var_stats$sord,
                                  stringsAsFactors=F))
    fin_stats <- fin_stats[order(fin_stats$sord, -fin_stats$RelImp),]
    file.remove('CE4_temp.xdf')
    
    return(fin_stats)
  }
  
  # Binary model statistics
  f.stats_log <- function(y2,weight,g = 10){
    f.ks <- function(p,a) {
      a0 <- p[a==0]
      a1 <- p[a==1]
      n.a1 <- as.double(length(a1))
      n.a0 <- as.double(length(a0))
      n <- n.a1 * n.a0 / (n.a1 + n.a0)
      w <- c(a1, a0)
      d <- max(abs(cumsum(ifelse(order(w) <= n.a1, 1 / n.a1, - 1 / n.a0))))
      return(d)
    }
    
    if (is.null(weight)){
      names(y2) <- c("dv","PScore")
      y2$rank <- floor((rank(y2[,"PScore"])*g)/(nrow(y2)+1)+1)
      y3 <- ddply(y2,"rank",summarize,
                  mactual=mean(dv),
                  badcnt=sum(ifelse(dv==1,1,0)),
                  goodcnt=sum(ifelse(dv==0,1,0)))
      norm <- mean(y2$dv)
      nbad <- nrow(y2[y2$dv==1,])
      ngood <- nrow(y2[y2$dv==0,])
    } else {
      names(y2) <- c("dv","PScore","wgt")
      y2 <- y2[order(y2$PScore),]
      y2$rank <- floor((cumsum(y2$wgt)*g)/(sum(y2$wgt)+1)+1)
      y3 <- ddply(y2,"rank",summarize,
                  mactual=weighted.mean(x=dv,w=wgt),
                  badcnt=sum(ifelse(dv==1,wgt,0)),
                  goodcnt=sum(ifelse(dv==0,wgt,0)))
      norm <- weighted.mean(x=y2$dv,w=y2$wgt)
      nbad <- sum(y2[y2$dv==1,"wgt"])
      ngood <- sum(y2[y2$dv==0,"wgt"])
    }
    fix <- data.frame(rank=seq(1:g))
    y3 <- merge(y3,fix,by="rank",all=T,sort=T)
    y3$mactual <- ifelse(is.na(y3$mactual),0,y3$mactual)
    y3$badcnt <- ifelse(is.na(y3$badcnt),0,y3$badcnt)
    y3$goodcnt <- ifelse(is.na(y3$goodcnt),0,y3$goodcnt)
    
    y3$goodratio <- ifelse(y3$goodcnt * y3$badcnt == 0,(y3$goodcnt+.5)/(ngood+.5),
                           y3$goodcnt / ngood)
    y3$badratio <- ifelse(y3$goodcnt * y3$badcnt == 0,(y3$badcnt+.5)/(nbad+.5),
                          y3$badcnt / nbad)
    y3$infv <- cumsum((y3$badratio-y3$goodratio)*log(y3$badratio/y3$goodratio))
    
    obs <- xtabs(cbind(1 - y2$dv, y2$dv) ~ y2$rank)
    expect <- xtabs(cbind(1 - y2$PScore, y2$PScore) ~ y2$rank)
    chisq <- sum((obs - expect)^2/expect)
    
    f1 <- floor(y2$PScore*500)
    event <- f1[y2$dv==1]
    nonevent <- f1[y2$dv==0]
    
    conc <- sum(sapply(event, function(x) sum(x > nonevent)))
    disc <- sum(sapply(event, function(x) sum(x < nonevent)))
    cnt <- length(event) * length(nonevent)
    tied <- cnt - conc - disc
    
    keep <- data.frame(SomersD = (conc-disc)/cnt,
                       Gamma = (conc-disc)/(conc+disc),
                       TauA = (conc-disc)/(.5*length(f1)*(length(f1)-1)),
                       c = (conc+.5*tied)/cnt,
                       Concord = (conc/cnt)*100,
                       Discon = (disc/cnt)*100,
                       LackFit = 1 - pchisq(chisq, g - 2),
                       Lift_Index = (y3[y3$rank==g,"mactual"] / norm)*100,
                       infv =  y3[y3$rank==g,"infv"],
                       ks =f.ks(y2$PScore,y2$dv),
                       stringsAsFactors=F)
    
    return(keep)
  }
  
  # Master function for continuous dependent variable
  f.Model_Val_Reg <- function(insdn,
                              regvlist, 
                              includelist = NULL,
                              startlist = NULL,
                              excludelist = NULL,
                              sle = 0.05,
                              sls = 0.05){
    
    
    cat('Get variable list from forward & backward select \n')
    origlst <- regvlist
    
    # Remove excluded variables from lists
    if (length(excludelist) > 0) {
      regvlist <- setdiff(regvlist,excludelist)
      if (length(startlist) > 0) {startlist <- setdiff(startlist,excludelist)}
      if (length(includelist) > 0) {includelist <- setdiff(includelist,excludelist)}
    }
    
    # Build base variable list
    if (length(startlist) > 0 & length(includelist) > 0) {basevars <- unique(startlist,includelist)}
    if (length(startlist) > 0 & length(includelist) == 0) {basevars <- unique(startlist)}
    if (length(startlist) == 0 & length(includelist) > 0) {basevars <- unique(includelist)}
    if (length(startlist) == 0 & length(includelist) == 0) {basevars <- 1}
    
    # Forward selection
    form <- as.formula(paste(dep_var, "~", paste(basevars, collapse = "+")))
    if (length(regvlist) > 1){
      scp <- as.formula(paste("~", paste(regvlist, collapse = "+")))
    } else {
      scp <- as.formula(paste("~", regvlist))
    }
    stp <- rxStepControl(method="forward",scope = scp,stepCriterion = "SigLevel",
                         maxSigLevelToAdd = sle, test="F", refitEachStep = F)
    # Model
    cat("Forward - Model \n")
    md <- rxLinMod(form,
                   data = insdn,
                   rowSelection = mod_val_test == 1,
                   variableSelection = stp,
                   reportProgress=0)
    Tforward <- row.names(md$coefficients)
    Tforward <- Tforward[Tforward != "(Intercept)"]
    rm(md)
    
    # Validation
    cat("Forward - Validation \n")
    md <- rxLinMod(form,
                   data = insdn,
                   rowSelection = mod_val_test == 2,
                   variableSelection = stp,
                   reportProgress=0)
    Vforward <- row.names(md$coefficients)
    Vforward <- Vforward[Vforward != "(Intercept)"]
    rm(md)
    
    # Backward selection
    if (length(regvlist) > 1){
      form <- as.formula(paste(dep_var, "~", paste(regvlist, collapse = "+")))
    } else {
      form <- as.formula(paste(dep_var, "~", regvlist))
    }
    scp <- as.formula(paste("~", paste(basevars, collapse = "+")))
    stp <- rxStepControl(method="backward",scope = list(lower=scp),stepCriterion = "SigLevel",
                         minSigLevelToDrop = sls,test="F", refitEachStep = F)
    # Model
    cat("Backward - Model \n")
    md <- rxLinMod(form,
                   data = insdn,
                   rowSelection = mod_val_test == 1,
                   variableSelection = stp,
                   reportProgress=0)
    Tbackward <- row.names(md$coefficients)
    Tbackward <- Tbackward[Tbackward != "(Intercept)"]
    rm(md)
    
    # Validation
    cat("Backward - Validation \n")
    md <- rxLinMod(form,
                   data = insdn,
                   rowSelection = mod_val_test == 2,
                   variableSelection = stp,
                   reportProgress=0)
    Vbackward <- row.names(md$coefficients)
    Vbackward <- Vbackward[Vbackward != "(Intercept)"]
    rm(md) 
    
    # Get statistics - model
    
    if (max(length(Tforward),length(Tbackward),length(Vforward),length(Vbackward)) == 0) {
      stop("No variables selected")
    } else {
      CE4_Model_Metric_Mod <- rbind(f.ctl_Stats_Reg(insdn, 1, dep_var, union(Tforward,Tbackward),"Dev"),
                                    f.ctl_Stats_Reg(insdn, 2, dep_var, union(Tforward,Tbackward),"Val"))
      save(CE4_Model_Metric_Mod, file = "CE4_Model_Metric_Mod.RData",compress=T)
      
      # Get statistics - validation
      CE4_Model_Metric_Val <- rbind(f.ctl_Stats_Reg(insdn, 2, dep_var, union(Vforward,Vbackward),"Dev"),
                                    f.ctl_Stats_Reg(insdn, 1, dep_var, union(Vforward,Vbackward),"Val"))
      save(CE4_Model_Metric_Val, file = "CE4_Model_Metric_Val.RData",compress=T)
      
      CE4_Variables <- data.frame(variable = origlst,
                                  Tforward = ifelse(origlst %in% Tforward,"Y",NA),
                                  Tbackward = ifelse(origlst %in% Tbackward,"Y",NA),
                                  Vforward = ifelse(origlst %in% Vforward,"Y",NA),
                                  Vbackward = ifelse(origlst %in% Vbackward,"Y",NA),
                                  stringsAsFactors=F)
      save(CE4_Variables, file = "CE4_Variables.RData",compress=T)
    }
  
  }
  # Continuous model control
  f.ctl_Stats_Reg <- function(insdn, part, dep_var, varlist, file) {
    
    cat('Get statistics for mod_val_test = ',part,'\n')
    
    # Separate data needed into dataframe
    f.tmp <- function(lst) {
      lst$sel <- lst$mod_val_test == ft
      return(lst)
    }
    rxDataStep(inData = insdn,
               outFile = "CE4_temp.xdf",
               overwrite = T,
               rowSelection = sel,
               varsToKeep  = c(dep_var,weight,varlist),
               transformFunc = f.tmp,
               transformObjects = list(ft = part),
               transformVars = "mod_val_test",
               reportProgress=0)
    
    # Full model
    cat("Full model \n")
    frm <- as.formula(paste(dep_var, "~", paste(varlist, collapse = "+")))
    md <- rxLinMod(formula = frm, data = "CE4_temp.xdf", pweights = weight, reportProgress=0)
    # Predicted values
    rxPredict(md,"CE4_temp.xdf",predVarNames="pred",overwrite=T,reportProgress=0)
    prd <- rxDataStep("CE4_temp.xdf",varsToKeep=c(dep_var,"pred",weight),reportProgress=0)
    # Elements needed for metric calculations
    sse <- md$residual.squares
    sst <- md$y.var
    #n <- md$nValidObs
    n <- nrow(prd)
    p <- nrow(md$coefficients)
    if (is.null(weight)) {mdv <- mean(prd[,dep_var])} else 
    {mdv <- weighted.mean(x=prd[,dep_var],w=prd[,weight])}
    
    full_stats <- data.frame(file = file,
                             variable = "Full File",
                             Estimate = NA,
                             StdEst = NA,
                             StdErr = NA,
                             tValue = NA,
                             Probt = NA,
                             VIF = NA,
                             RelImp = NA,
                             AIC = n*log(sse/n) + (2*p),
                             SBC = n * log(sse/n) + (p*log(n)),
                             JP = ((n + p) / n) * (sse / (n - p)),
                             Rsquare = md$r.squared,
                             AdjRsq = md$adj.r.squared,
                             RMSE = sqrt(sse / (n - p)),
                             CoeffVar = (sqrt(sse / (n - p)) * 100) / mdv,
                             PC = ((n+p)/(n-p))*(1-md$r.squared),
                             f.stats_reg(prd,weight),
                             sord = 1,
                             stringsAsFactors=F)
    
    se_var <- (md$coefficients * c(NA,apply(rxDataStep("CE4_temp.xdf",varsToKeep=varlist,reportProgress=0),2,sd)))/
      mean(prd[,dep_var])
    # Variance inflation factors
    f.tmp <- function(lst) {
      lst$vif_wgt <- lst[[1]]
      return(lst)
    }
    vif <- vif(lm(frm,
                  rxDataStep("CE4_temp.xdf",
                             transformFunc = f.tmp,
                             transformVars = weight,
                             reportProgress=0),
                  weights=vif_wgt))
    
    var_stats <- data.frame(file = file,
                            variable = row.names(md$coefficients),
                            Estimate = md$coefficients,
                            StdEst = se_var,
                            StdErr = md$coef.std.error,
                            tValue = md$coef.t.value,
                            Probt = md$coef.p.value,
                            VIF = c(NA,vif),
                            RelImp = abs(se_var) / sum(abs(se_var),na.rm=T),
                            stringsAsFactors=F)
    names(var_stats) <- c("file","variable","Estimate","StdEst","StdErr","tValue","Probt","VIF","RelImp")
    
    # No intercept model
    cat("No intercept model \n")
    frm <- as.formula(paste(dep_var, "~", paste(c(-1,varlist), collapse = "+")))
    md <- rxLinMod(formula = frm, data = "CE4_temp.xdf", pweights = weight, reportProgress=0)
    # Predicted values
    rxPredict(md,"CE4_temp.xdf",predVarNames="pred",overwrite=T,reportProgress=0)
    prd <- rxDataStep("CE4_temp.xdf",varsToKeep=c(dep_var,"pred",weight),reportProgress=0)
    # Elements needed for metric calculations
    sse <- md$residual.squares
    sst <- md$y.var
    n <- md$nValidObs
    p <- nrow(md$coefficients)
    if (is.null(weight)) {mdv <- mean(prd[,dep_var])} else 
    {mdv <- weighted.mean(x=prd[,dep_var],w=prd[,weight])}
    
    tmp_stats <- data.frame(variable = "(Intercept)",
                            AIC = n*log(sse/n) + (2*p),
                            SBC = n * log(sse/n) + (p*log(n)),
                            JP = ((n + p) / n) * (sse / (n - p)),
                            Rsquare = md$r.squared,
                            AdjRsq = md$adj.r.squared,
                            RMSE = sqrt(sse / (n - p)),
                            CoeffVar = (sqrt(sse / (n - p)) * 100) / mdv,
                            PC = ((n+p)/(n-p))*(1-md$r.squared),
                            f.stats_reg(prd,weight),
                            sord = 2,
                            stringsAsFactors=F)
    
    # Models removing one variable at a time
    cat("Subset models - remove one variable at a time \n")
    for (i in 1:length(varlist)) {
      frm <- as.formula(paste(dep_var, "~", paste(varlist[-i], collapse = "+")))
      md <- rxLinMod(formula = frm, data = "CE4_temp.xdf", pweights = weight, reportProgress=0)
      # Predicted values
      rxPredict(md,"CE4_temp.xdf",predVarNames="pred",overwrite=T,reportProgress=0)
      prd <- rxDataStep("CE4_temp.xdf",varsToKeep=c(dep_var,"pred",weight),reportProgress=0)
      # Elements needed for metric calculations
      sse <- md$residual.squares
      sst <- md$y.var
      n <- md$nValidObs
      p <- nrow(md$coefficients)
      if (is.null(weight)) {mdv <- mean(prd[,dep_var])} else 
      {mdv <- weighted.mean(x=prd[,dep_var],w=prd[,weight])}
      
      tmp_stats <- rbind(tmp_stats,
                         data.frame(variable = varlist[i],
                                    AIC = n*log(sse/n) + (2*p),
                                    SBC = n * log(sse/n) + (p*log(n)),
                                    JP = ((n + p) / n) * (sse / (n - p)),
                                    Rsquare = md$r.squared,
                                    AdjRsq = md$adj.r.squared,
                                    RMSE = sqrt(sse / (n - p)),
                                    CoeffVar = (sqrt(sse / (n - p)) * 100) / mdv,
                                    PC = ((n+p)/(n-p))*(1-md$r.squared),
                                    f.stats_reg(prd,weight),
                                    sord = 3,
                                    stringsAsFactors=F))
    }
    var_stats <- merge(var_stats,tmp_stats,by="variable")
    
    # Subtract full model metrics from sub-model metrics
    # This shows the impact of each variable
    fin_stats <- rbind(full_stats,
                       data.frame(file = var_stats$file,
                                  variable = var_stats$variable,
                                  Estimate = var_stats$Estimate,
                                  StdEst = var_stats$StdEst,
                                  StdErr = var_stats$StdErr,
                                  tValue = var_stats$tValue,
                                  Probt = var_stats$Probt,
                                  VIF = var_stats$VIF,
                                  RelImp = var_stats$RelImp,
                                  AIC = var_stats$AIC - full_stats$AIC,
                                  SBC = var_stats$SBC - full_stats$SBC,
                                  JP = var_stats$JP - full_stats$JP,
                                  Rsquare = var_stats$Rsquare - full_stats$Rsquare,
                                  AdjRsq = var_stats$AdjRsq - full_stats$AdjRsq,
                                  RMSE = var_stats$RMSE - full_stats$RMSE,
                                  CoeffVar = var_stats$CoeffVar - full_stats$CoeffVar,
                                  PC = var_stats$PC - full_stats$PC,
                                  Lift_Index = var_stats$Lift_Index - full_stats$Lift_Index,
                                  gini = var_stats$gini - full_stats$gini,
                                  infv = var_stats$infv - full_stats$infv,
                                  ks = var_stats$ks - full_stats$ks,
                                  sord = var_stats$sord,
                                  stringsAsFactors=F))
    fin_stats <- fin_stats[order(fin_stats$sord, -fin_stats$RelImp),]
    file.remove('CE4_temp.xdf')
    
    return(fin_stats)
  }
  # Continuous model statistics
  f.stats_reg <- function(y2,weight,g = 10){
    
    if (is.null(weight)){
      names(y2) <- c("dv","PScore")
      y2$rank <- floor((rank(y2[,"PScore"])*g)/(nrow(y2)+1)+1)
      y3 <- ddply(y2,"rank",summarize,
                  count=length(dv),
                  mactual=mean(dv),
                  mpred=mean(PScore))
      norm <- mean(y2$dv)
    } else {
      names(y2) <- c("dv","PScore","wgt")
      y2 <- y2[order(y2$PScore),]
      y2$rank <- floor((cumsum(y2$wgt)*g)/(sum(y2$wgt)+1)+1)
      y3 <- ddply(y2,"rank",summarize,
                  count=sum(wgt),
                  mactual=weighted.mean(x=dv,w=wgt),
                  mpred=weighted.mean(x=PScore,w=wgt))
      norm <- weighted.mean(x=y2$dv,w=y2$wgt)
    }
    fix <- data.frame(rank=seq(1:g))
    y3 <- merge(y3,fix,by="rank",all=T,sort=T)
    y3$count <- ifelse(is.na(y3$count),0,y3$count)
    y3$mactual <- ifelse(is.na(y3$mactual),0,y3$mactual)
    y3$mpred <- ifelse(is.na(y3$mpred),0,y3$mpred)
    
    totalbad <- sum(y3$count*(y3$count*y3$mactual/sum(y3$count*y3$mactual)))
    
    cumresp <- cumsum(y3$mactual*y3$count)
    cumpct_resp <- cumresp / sum(y3$count*y3$mactual)
    cumtotal <- cumsum(y3$count)
    cumpct_freq <- cumtotal / sum(y3$count)
    cumavg <- cumresp / cumtotal
    badrate <- y3$count * y3$mactual / sum(y3$count*y3$mactual)
    goodrate <- (1 - badrate) * y3$count / (sum(y3$count) - totalbad)
    prev_dv <- c(NA,cumpct_resp[1:(g-1)])
    prev_pct <- c(NA,cumpct_freq[1:(g-1)])
    infv <- cumsum((badrate-goodrate)*log(badrate/goodrate))
    cumbad <- cumsum(badrate)
    cumgood <- cumsum(goodrate)
    ks <- abs(cumbad-cumgood)
    gini <- 2*((cumpct_freq+prev_pct)/2-(cumpct_resp+prev_dv)/2)*(cumpct_freq-prev_pct)
    gini[1] <- 0
    
    keep <- data.frame(Lift_Index = (y3[y3$rank==g,"mactual"] / norm)*100,
                       gini = sum(gini),
                       infv = infv[g],
                       ks = max(ks),
                       stringsAsFactors=F)
    
    return(keep)
  }
  
  # Variable Tuning and Selection Based on model metric output
  f.vars_tune <- function(dt1,dt2,criteria,threshold,minimp,join){
    cat("Variable Tuning and Selection \n")
    # Change variable names to upper to avoid case issues
    names(dt1) <- tolower(names(dt1))
    names(dt2) <- tolower(names(dt2))
    criteria <- tolower(criteria)
   
    # Check criteria and reset if necessary
    if (toupper(binary_dv) == "Y") {
      vl <- c("aic","sc","logl2","rsquare","somersd","gamma","taua","c","concord","discon",
              "lackfit","lift_index","infv","ks")
      if (!(criteria %in% vl)) {
        criteria <- "c"
        cat("Criteria changed to c \n")
      }
    } else {
      vl <- c("aic","sbc","jp","rsquare","adjrsq","rmse","coeffvar","pc","lift_index","infv","ks")
      if (!(criteria %in% vl)) {
        criteria <- "adjrsq"
        cat("Criteria changed to AdjRsq \n")
      }
    }
    
    if (criteria %in%  c("rsquare","somersd","gamma","taua","c","concord","infv","ks",
                         "adjrsq","lift_index","gini")) {threshold <- -threshold}
    if (toupper(join) == "UNION") {
      vars <- union(dt1[dt1$sord==3 & dt1[,criteria] <= threshold & dt1$relimp >= minimp & dt1$file =="Val","variable"],
                    dt2[dt2$sord==3 & dt2[,criteria] <= threshold & dt2$relimp >= minimp & dt2$file =="Val","variable"])
    } else {
      vars <- intersect(dt1[dt1$sord==3 & dt1[,criteria] <= threshold & dt1$relimp >= minimp & dt1$file =="Val","variable"],
                        dt2[dt2$sord==3 & dt2[,criteria] <= threshold & dt2$relimp >= minimp & dt2$file =="Val","variable"])
    }
    
    # Write out list of variables to keep
    f.catwrite<-function(...,f='CE4_Varlist_Final.txt',app=F) cat(...,file=f,append=app)
    
    if (length(vars)> 1){
      temp <- data.frame(name=vars,len=nchar(vars),splt=NA,outvar=NA,stringsAsFactors=F)
      temp[1,"splt"] <- 1
      temp[1,"outvar"] <- paste("varlist_final <- c('",temp[1,"name"],"'",sep="")
      cnt <- temp[1,"len"]
      ov <- temp[1,"outvar"]
      splt <- 1
      
      for (i in 2:nrow(temp)) {
        if (is.na(temp[i,"splt"])) {
          if (cnt + temp[i,"len"] > 250) {
            splt <- splt + 1
            cnt <- temp[i,"len"]
            ov <- paste("'",temp[i,"name"],"'",sep="")
          } else
          { ov <- paste(ov,",'",temp[i,"name"],"'",sep="")
            cnt <- cnt + temp[i,"len"] + 3
          }
          temp[i,"splt"] <- splt
          temp[i,"outvar"] <- ov
        }
      }
      # Mark first and last record for strings
      temp$slast <- !duplicated( temp[, "splt" ], fromLast=T)
      # Keep lest record for each string
      temp2 <- temp[temp$slast==T,]
      temp2$outvar <- ifelse(temp2$splt==splt,
                             paste(temp2$outvar,")",sep=""),
                             paste(temp2$outvar,",",sep=""))
      
      f.catwrite(paste(temp2$outvar,'\n',sep=""))
      rm(f.catwrite,temp,cnt,ov,splt,temp2,i)
      cat(length(vars)," variables selected \n")
    } else {
      if (length(vars)==1){
        f.catwrite(paste("varlist_final <- c('",vars,')\n',sep=""))  
        cat("1 variable selected \n")
      } else {
        cat("No variables selected \n")
      }
    }
    load("CE4_Variables.RData")
    CE4_Variables$final <- ifelse(CE4_Variables$variable %in% vars,"Y",NA)
    save(CE4_Variables, file = "CE4_Variables.RData",compress=T)
  }
  
  # Graphing function
  f.graph <- function(x,dep_var,vars,cats=20){
    # Set up workbook
    cat("Create graph workbook \n")
    if (file.exists('CE4_Graphs.xlsx')) file.remove('CE4_Graphs.xlsx')
    gr <- loadWorkbook('CE4_Graphs.xlsx', create = TRUE)
    
    csComma <- createCellStyle(gr, name = "comma")
    setDataFormat(csComma, format = "#,##0")
    csDec4 <- createCellStyle(gr, name = "dec4")
    setDataFormat(csDec4, format = "0.0000")
    
    #Loop through variables and write out
    for (i in 1:length(vars)){
      # Create sheet name
      sn <- paste0("var",i)
      # Create temporary dataset with one independent variable
      tmp <- data.frame(dv=x[,dep_var],iv=x[,vars[i]])
      # Summarize data
      if (nrow(table(tmp$iv,useNA="ifany")) <= cats){
        f1 <- ddply(tmp,"iv",summarize,
                    count=length(iv),
                    mean=mean(dv))
        f1$iv <- round(f1$iv, 4)
        names(f1) <- c(vars[i],"count","mean")
      } else {
        tmp$bin <- cut(tmp$iv,b=cats,labels=F,include.lowest=T)
        f1 <- ddply(tmp,"bin",summarize,
                    count=length(iv),
                    lo=min(iv),
                    hi=max(iv),
                    mean=mean(dv))  
      }
      allRows <- seq(length = nrow(f1)) + 3
      createSheet(gr, name = sn)
      writeWorksheet(gr, paste0("Variable: ",vars[i]), sheet = sn, startRow = 1, startCol = 1, header = F)
      
      writeWorksheet(gr, f1, sheet = sn, startRow = 3, startCol = 1)
      if (ncol(f1)==3){
        setColumnWidth(gr, sheet = sn, column = 1, width = -1)
        setCellStyle(gr, sheet = sn, row = allRows, col = 2, cellstyle = csComma)
        setCellStyle(gr, sheet = sn, row = allRows, col = 3, cellstyle = csDec4)
      } else {
        setCellStyle(gr, sheet = sn, row = allRows, col = 2, cellstyle = csComma)
        setCellStyle(gr, sheet = sn, row = allRows, col = 3, cellstyle = csDec4)
        setCellStyle(gr, sheet = sn, row = allRows, col = 4, cellstyle = csDec4)
        setCellStyle(gr, sheet = sn, row = allRows, col = 5, cellstyle = csDec4)
      }
      
      # Create chart
      jpeg(filename = "var.jpeg", width = 600, height = 400)
      layout(rbind(1,2), heights=c(1,7))  # put legend on top 1/8th of the chart
      # setup for no margins on the legend
      par(mar=c(0, 0, 1, 7))
      plot.new()
      legend('right','groups',c(paste0("Mean ",dep_var),"# of Customers"),
             col=c('lightblue','dodgerblue4'),
             lwd=c(10,3),ncol=2,bty ="n")
      # Build plots
      par(mar=c(7, 5, 0, 7) + 0.2)
      if (ncol(f1)==3){
        xx <- barplot(f1$mean, ylab = paste0("Mean ",dep_var), font.lab=2,
                      names.arg = round(f1[,vars[i]],4) , 
                      ylim = c(0, max(f1$mean)+.1), las=2,
                      col="lightblue")
      } else {
        xx <- barplot(f1$mean, ylab = paste0("Mean ",dep_var), font.lab=2,
                      names.arg = round(f1$lo,4) , 
                      ylim = c(0, max(f1$mean)+.1), las=2,
                      col="lightblue")
      }
      
      par(new = T)
      plot(x=xx+.5, f1$count, "l", lwd = 3, col="dodgerblue4", lty=1, lend=2,
           axes=F, ylim=c(0, max(f1$count)), xlim=c(min(xx), max(xx)+1),
           xlab = "", ylab="")
      axis(4, ylim=c(0, max(f1$count)), las=1)
      mtext(vars[i],side=1,line=5,font=2)
      mtext("# of Customers",side=4, line=4,font=2)
      box(which="plot", lty="solid")
      box(which="outer", lty="solid", lwd=3)
      dev.off()
      
      createName(gr, name = "vgr", formula = paste0(sn,"!$G$3"), overwrite=T)
      addImage(gr, filename = "var.jpeg", name = "vgr", originalSize = TRUE)
      removeName(gr, name = "vgr")
    }
    file.remove('var.jpeg')
    saveWorkbook(gr)
  }
  
  if (toupper(binary_dv) == "Y") {
    f.Model_Val_Logistic(insdn,
                         varlist, 
                         includelist = includelist,
                         startlist = startlist,
                         excludelist = excludelist,
                         sle = sel_alpha,
                         sls = sel_alpha)
    
    load("CE4_Model_Metric_Mod.RData")
    load("CE4_Model_Metric_Val.RData")
    
    f.vars_tune(CE4_Model_Metric_Mod,CE4_Model_Metric_Val,criteria,threshold,minimp,SQL_join)
    
    load("CE4_Variables.RData")
    
    CE4_Model_Metric_Mod$sord <- NULL
    CE4_Model_Metric_Val$sord <- NULL
    
    # Write out Excel report
    cat("Create Excel report \n")
    if (file.exists('CE4_Model_Report.xlsx')) file.remove('CE4_Model_Report.xlsx')
    mr <- loadWorkbook('CE4_Model_Report.xlsx', create = TRUE)
    
    csDec2 = createCellStyle(mr, name = "dec2")
    setDataFormat(csDec2, format = "0.00")
    csDec4 = createCellStyle(mr, name = "dec4")
    setDataFormat(csDec4, format = "0.0000")
    csDec6 = createCellStyle(mr, name = "dec6")
    setDataFormat(csDec6, format = "0.000000")
    
    ## Write out report
    createSheet(mr, name = "Variables")
    writeWorksheet(mr, CE4_Variables, sheet = "Variables", startRow = 1, startCol = 1)
    setColumnWidth(mr, sheet = "Variables", column = c(1), width = 8192)
    setColumnWidth(mr, sheet = "Variables", column = c(2,3,4,5,6), width = 2560)
    
    allRows <- seq(length = nrow(CE4_Model_Metric_Mod)) + 3
    cls <- seq(3:23)+2
    createSheet(mr, name = "Model")
    writeWorksheet(mr, "Diagnostic statistics: Model Portion", sheet = "Model", startRow = 1, startCol = 1, header = F)
    writeWorksheet(mr, CE4_Model_Metric_Mod, sheet = "Model", startRow = 3, startCol = 1)
    setColumnWidth(mr, sheet = "Model", column = 2, width = 8192)
    setColumnWidth(mr, sheet = "Model", column = cls, width = 2560)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 3, cellstyle = csDec6)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 4, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 5, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 6, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 7, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 8, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 9, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 10, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 11, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 12, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 13, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 14, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 15, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 16, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 17, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 18, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 19, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 20, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 21, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 22, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 23, cellstyle = csDec4)
    
    allRows <- seq(length = nrow(CE4_Model_Metric_Val)) + 3
    createSheet(mr, name = "Validation")
    writeWorksheet(mr, "Diagnostic statistics: Validation Portion", sheet = "Validation", startRow = 1, startCol = 1, header = F)
    writeWorksheet(mr, CE4_Model_Metric_Val, sheet = "Validation", startRow = 3, startCol = 1)
    setColumnWidth(mr, sheet = "Validation", column = 2, width = 8192)
    setColumnWidth(mr, sheet = "Validation", column = cls, width = 2560)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 3, cellstyle = csDec6)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 4, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 5, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 6, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 7, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 8, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 9, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 10, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 11, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 12, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 13, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 14, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 15, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 16, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 17, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 18, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 19, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 20, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 21, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 22, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 23, cellstyle = csDec4)
    
    saveWorkbook(mr)
    
  } else {
    f.Model_Val_Reg(insdn,
                    varlist, 
                    includelist = includelist,
                    startlist = startlist,
                    excludelist = excludelist,
                    sle = sel_alpha,
                    sls = sel_alpha)
    
    
    load("CE4_Model_Metric_Mod.RData")
    load("CE4_Model_Metric_Val.RData")
    
    f.vars_tune(CE4_Model_Metric_Mod,CE4_Model_Metric_Val,criteria,threshold,minimp,SQL_join)
    
    load("CE4_Variables.RData")
    
    CE4_Model_Metric_Mod$sord <- NULL
    CE4_Model_Metric_Val$sord <- NULL
    
    # Write out Excel report
    cat("Create Excel report \n")
    if (file.exists('CE4_Model_Report.xlsx')) file.remove('CE4_Model_Report.xlsx')
    mr <- loadWorkbook('CE4_Model_Report.xlsx', create = TRUE)
    
    csDec2 = createCellStyle(mr, name = "dec2")
    setDataFormat(csDec2, format = "0.00")
    csDec4 = createCellStyle(mr, name = "dec4")
    setDataFormat(csDec4, format = "0.0000")
    csDec6 = createCellStyle(mr, name = "dec6")
    setDataFormat(csDec6, format = "0.000000")
    
    ## Write out report
    createSheet(mr, name = "Variables")
    writeWorksheet(mr, CE4_Variables, sheet = "Variables", startRow = 1, startCol = 1)
    setColumnWidth(mr, sheet = "Variables", column = c(1), width = 8192)
    setColumnWidth(mr, sheet = "Variables", column = c(2,3,4,5,6), width = 2560)
    
    allRows <- seq(length = nrow(CE4_Model_Metric_Mod)) + 3
    cls <- seq(3:21)+2
    createSheet(mr, name = "Model")
    writeWorksheet(mr, "Diagnostic statistics: Model Portion", sheet = "Model", startRow = 1, startCol = 1, header = F)
    writeWorksheet(mr, CE4_Model_Metric_Mod, sheet = "Model", startRow = 3, startCol = 1)
    setColumnWidth(mr, sheet = "Model", column = 2, width = 8192)
    setColumnWidth(mr, sheet = "Model", column = cls, width = 2560)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 3, cellstyle = csDec6)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 4, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 5, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 6, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 7, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 8, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 9, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 10, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 11, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 12, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 13, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 14, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 15, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 16, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 17, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 18, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 19, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 20, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Model", row = allRows, col = 21, cellstyle = csDec4)
    
    allRows <- seq(length = nrow(CE4_Model_Metric_Val)) + 3
    createSheet(mr, name = "Validation")
    writeWorksheet(mr, "Diagnostic statistics: Validation Portion", sheet = "Validation", startRow = 1, startCol = 1, header = F)
    writeWorksheet(mr, CE4_Model_Metric_Val, sheet = "Validation", startRow = 3, startCol = 1)
    setColumnWidth(mr, sheet = "Validation", column = 2, width = 8192)
    setColumnWidth(mr, sheet = "Validation", column = cls, width = 2560)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 3, cellstyle = csDec6)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 4, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 5, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 6, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 7, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 8, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 9, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 10, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 11, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 12, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 13, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 14, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 15, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 16, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 17, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 18, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 19, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 20, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Validation", row = allRows, col = 21, cellstyle = csDec4)
    
    saveWorkbook(mr)
    
  }
  
  source('CE4_Varlist_Final.txt')
  tst <- rxGetInfo(data=insdn)
  mx <- tst$numRows*length(varlist_final)
  if (toupper(graph_plot) == "Y") {
    f.graph(rxDataStep(inData=insdn,
                       rowSelection = mod_val_test != 3,
                       varsToKeep = c(dep_var,varlist_final),
                       maxRowsByCols = mx,
                       reportProgress=0),
            dep_var,varlist_final,cats=20)
  }
  # Clean up
  rm(list=ls())
}

# 5. FINAL MODEL BUILD ##########################

## *************************************************************
## ***** 5.Final Model build and validation on test sample *****
## *************************************************************

## *** Macro 5: Final Model selection ***;
## %let fin_alpha = .05;                 ** Alpha level when selecting/removing variable into the model;
## %let method = stepwise;               ** Model build method;
## %let fin_num_category = 10;           ** Maximum number of categories for profiling variables;
## %let fin_equal_dist = N;              ** Use equal distance for dividing variables into groups for profiling (Y/N);

f.model_lift <- function(insdn,
                         varlist,
                         dep_var,
                         binary_dv,
                         weight = NULL,
                         fin_alpha = .05,
                         refit = F,
                         method = "stepwise",
                         fin_num_category = 10,
                         fin_equal_dist = "N") {
  
  varlist <- varlist[!is.na(varlist) & varlist!=""] # Make sure there are no extraneous items
  
  # Separate data needed
  rxDataStep(inData = insdn,
             outFile = "CE5_mod.xdf",
             overwrite = T,
             rowSelection = mod_val_test != 3,
             varsToKeep  = c(dep_var,weight,"mod_val_test",varlist),
             reportProgress=0)
  rxDataStep(inData = insdn,
             outFile = "CE5_val.xdf",
             overwrite = T,
             rowSelection = mod_val_test == 3,
             varsToKeep  = c(dep_var,weight,"mod_val_test",varlist),
             reportProgress=0)
  
  # Master function for binary dependent variable
  f.Model_Fin_Logistic <- function(insdn,
                                   regvlist,
                                   dep_var,
                                   weight = NULL,
                                   method = "stepwise",
                                   sle = 0.05,
                                   sls = 0.05){
    
    cat('Build final model \n')
    
    # Model
    if (tolower(method) %in% c("stepwise","forward")){
      form <- as.formula(paste(dep_var,"~1"))
      if (length(regvlist) > 1){
        scp <- as.formula(paste("~", paste(regvlist, collapse = "+")))  
      } else {
        scp <- as.formula(paste("~", regvlist))
      }
      stp <- rxStepControl(method=method, scope = scp, stepCriterion = "SigLevel",
                           maxSigLevelToAdd = sle, minSigLevelToDrop = sls,
                           test="Chisq", refitEachStep = refit)
    } else {
      if (length(regvlist) > 1){
        form <- as.formula(paste(dep_var, "~", paste(regvlist, collapse = "+")))
      } else {
        form <- as.formula(paste(dep_var, "~", regvlist))
      }
      scp <- as.formula("~1")
      stp <- rxStepControl(method=method,scope = list(lower=scp),stepCriterion = "SigLevel",
                           minSigLevelToDrop = sls,test="Chisq", refitEachStep = refit)
    }
    
    md0 <- rxLogit(form,
                   data = "CE5_mod.xdf",
                   pweights = weight,
                   variableSelection = stp,
                   initialValues = NA,
                   coeffTolerance = .1,
                   reportProgress=0)
    nvarlist <- names(md0$coefficients)
    nvarlist <- nvarlist[nvarlist != "(Intercept)"]
    
    # Create files with predicted score
    rxPredict(md0,"CE5_mod.xdf",type="response",predVarNames="pred",overwrite=T,reportProgress=0)
    rxPredict(md0,"CE5_val.xdf",type="response",predVarNames="pred",overwrite=T,reportProgress=0)
    rxMerge(inData1="CE5_mod.xdf",inData2="CE5_val.xdf",outFile="CE5_scored.xdf",
            type="union",overwrite=T,xdfCompressionLevel = 1,reportProgress=0)
    
    # Get statistics
    CE5_Model_Metric_All <- rbind(f.ctl_Stats_Log("CE5_mod.xdf", dep_var, nvarlist,"Dev"),
                                  f.ctl_Stats_Log("CE5_val.xdf", dep_var, nvarlist,"Val"))
    save(CE5_Model_Metric_All, file = "CE5_Model_Metric_All.RData",compress=T)
    
    # Parameters
    CE5_ParameterEstimates <- CE5_Model_Metric_All[CE5_Model_Metric_All$file=="Dev" & 
                                                     CE5_Model_Metric_All$variable != "Full File" ,
                                                   c("variable","Estimate","StdErr",
                                                     "zValue","zProb","StdEst")]
    save(CE5_ParameterEstimates, file = "CE5_ParameterEstimates.RData",compress=T)
    
    # Write out scoring equation
    f.rwrite<-function(...,f='CE5_Scoring_equation.txt',app=T) cat(...,'\n',file=f,append=app)
    f.rwrite('f.score <- function(x){',app=F)
    f.rwrite('#####  Model  #####')
    
    var <- CE5_ParameterEstimates$variable
    est <- round(CE5_ParameterEstimates$Estimate,6)
    state <- ifelse(var=="(Intercept)",paste0("x$logit <- (1 * (",est,") + \n"),
                    paste0("x$",var," * (",est,") + \n"))
    state[length(var)] <- paste0("x$",var[length(var)]," * (",est[length(var)],")) \n")
    
    f.rwrite(state)
    f.rwrite("x$score <- 1/(1+exp(-x$logit)) \n")
    f.rwrite("return(x)} \n")
    
    # Performance and variable validation
    tmp <- f.gains(rxDataStep("CE5_mod.xdf",varsToKeep=c(dep_var,"pred",weight,nvarlist),
                              reportProgress=0),
                   weight=weight, binary_dv=binary_dv, vars=nvarlist)
    CE5_Gainstable_Train <- tmp[[1]]
    CE5_Profiles_Train <- tmp[[2]]
    tmp <- f.gains(rxDataStep("CE5_val.xdf",varsToKeep=c(dep_var,"pred",weight,nvarlist)
                              ,reportProgress=0),
                   weight=weight, binary_dv=binary_dv, vars=nvarlist)
    CE5_Gainstable_Test <- tmp[[1]]
    CE5_Profiles_Test <- tmp[[2]]
    rm(tmp)
    save(CE5_Gainstable_Train, file = "CE5_Gainstable_Train.RData",compress=T)
    save(CE5_Profiles_Train, file = "CE5_Profiles_Train.RData",compress=T)
    save(CE5_Gainstable_Test, file = "CE5_Gainstable_Test.RData",compress=T)
    save(CE5_Profiles_Test, file = "CE5_Profiles_Test.RData",compress=T)
    
    # Correlations
    CE5_Corr <- data.frame(rxCor(as.formula(paste("~", paste(c(dep_var,nvarlist), collapse = "+"))),
                                 "CE5_scored.xdf",reportProgress=0))
    Variable <- row.names(CE5_Corr)
    CE5_Corr <- cbind(Variable,CE5_Corr)
    save(CE5_Corr, file = "CE5_Corr.RData",compress=T)
    
    # Profiling
    CE5_Profile <- f.profile(rxDataStep("CE5_mod.xdf",
                                        varsToKeep=c(dep_var,nvarlist),
                                        reportProgress=0),
                             dep_var = dep_var,
                             vars=nvarlist,
                             num_category = fin_num_category,
                             equal_dist = fin_equal_dist)
    save(CE5_Profile, file = "CE5_Profile.RData",compress=T)
    
    # Write out variable list
    f.mod_write(nvarlist)
  }
  # Binary model control
  f.ctl_Stats_Log <- function(insdn, dep_var, varlist, file) {
    
    cat('Get statistics for ',file,'\n')
    
    # Full model
    cat("Full model \n")
    null.deviance <- rxLogit(formula = as.formula(paste(dep_var,"~1")),
                             data = insdn,
                             pweights = weight,
                             reportProgress=0)$deviance
    frm <- as.formula(paste(dep_var, "~", paste(varlist, collapse = "+")))
    md <- rxLogit(formula = frm, data = insdn, pweights = weight, initialValues = NA,
                  coeffTolerance = .1, reportProgress=0)
    # Number of valid observations for calculations
    rws <- md$nValidObs
    # Predicted values
    rxPredict(md,insdn,type="response",predVarNames="pred",overwrite=T,reportProgress=0)
    prd <- rxDataStep(insdn,varsToKeep=c(dep_var,"pred",weight),reportProgress=0)
    
    full_stats <- data.frame(file = file,
                             variable = "Full File",
                             Estimate = NA,
                             StdEst = NA,
                             StdErr = NA,
                             zValue = NA,
                             zProb = NA,
                             VIF = NA,
                             RelImp = NA,
                             AIC = md$aic,
                             SC = md$deviance + (length(names(md$coefficients))*log(rws)),
                             LogL2 = md$deviance,
                             Rsquare = 1 - exp((-2/nrow(prd))*(md$deviance/-2 - null.deviance/-2)),
                             f.stats_log(prd,weight),
                             sord = 1,
                             stringsAsFactors=F)
    
    coef <- as.data.frame(coefficients(summary(md)))
    se_var <- coef[,"Estimate"] / ((pi/sqrt(3))/c(NA,apply(rxDataStep(insdn,varsToKeep=varlist,reportProgress=0),2,sd)))
    # Variance inflation factors
    vif <- vif(lm(frm,
                  rxDataStep(insdn,
                             transforms=list(vif_wgt=pred*(1-pred)),
                             reportProgress=0),
                  weights=vif_wgt))
    
    var_stats <- data.frame(file = file,
                            variable = names(md$coefficients),
                            Estimate = coef[,"Estimate"],
                            StdEst = se_var,
                            StdErr = coef[,"Std. Error"],
                            zValue = coef[,3],
                            zProb = coef[,4],
                            VIF = c(NA,vif),
                            RelImp = abs(se_var) / sum(abs(se_var),na.rm=T),
                            stringsAsFactors=F)
    
    # No intercept model
    cat("No intercept model \n")
    frm <- as.formula(paste(dep_var, "~", paste(c(-1,varlist), collapse = "+")))
    md <- rxLogit(formula = frm, data = insdn, pweights = weight, initialValues = NA,
                  coeffTolerance = .1, reportProgress=0)
    # Predicted values
    rxPredict(md,insdn,type="response",predVarNames="pred",overwrite=T,reportProgress=0)
    prd <- rxDataStep(insdn,varsToKeep=c(dep_var,"pred",weight),reportProgress=0)
    tmp_stats <- data.frame(variable = "(Intercept)",
                            AIC = md$aic,
                            SC = md$deviance + (length(names(md$coefficients))*log(rws)),
                            LogL2 = md$deviance,
                            Rsquare = 1 - exp((-2/nrow(prd))*(md$deviance/-2 - null.deviance/-2)),
                            f.stats_log(prd,weight),
                            sord = 2,
                            stringsAsFactors=F)
    
    # Models removing one variable at a time
    cat("Subset models - remove one variable at a time \n")
    for (i in 1:length(varlist)) {
      frm <- as.formula(paste(dep_var, "~", paste(varlist[-i], collapse = "+")))
      md <- rxLogit(formula = frm, data = insdn, pweights = weight, initialValues = NA,
                    coeffTolerance = .1, reportProgress=0)
      # Predicted values
      rxPredict(md,insdn,type="response",predVarNames="pred",overwrite=T,reportProgress=0)
      prd <- rxDataStep(insdn,varsToKeep=c(dep_var,"pred",weight),reportProgress=0)
      tmp_stats <- rbind(tmp_stats,
                         data.frame(variable = varlist[i],
                                    AIC = md$aic,
                                    SC = md$deviance + (length(names(md$coefficients))*log(rws)),
                                    LogL2 = md$deviance,
                                    Rsquare = 1 - exp((-2/nrow(prd))*(md$deviance/-2 - null.deviance/-2)),
                                    f.stats_log(prd,weight),
                                    sord = 3,
                                    stringsAsFactors=F))
    }
    var_stats <- merge(var_stats,tmp_stats,by="variable")
    
    # Subtract full model metrics from sub-model metrics
    # This shows the impact of each variable
    fin_stats <- rbind(full_stats,
                       data.frame(file = var_stats$file,
                                  variable = var_stats$variable,
                                  Estimate = var_stats$Estimate,
                                  StdEst = var_stats$StdEst,
                                  StdErr = var_stats$StdErr,
                                  zValue = var_stats$zValue,
                                  zProb = var_stats$zProb,
                                  VIF = var_stats$VIF,
                                  RelImp = var_stats$RelImp,
                                  AIC = var_stats$AIC - full_stats$AIC,
                                  SC = var_stats$SC - full_stats$SC,
                                  LogL2 = var_stats$LogL2 - full_stats$LogL2,
                                  Rsquare = var_stats$Rsquare - full_stats$Rsquare,
                                  SomersD = var_stats$SomersD - full_stats$SomersD,
                                  Gamma = var_stats$Gamma - full_stats$Gamma,
                                  TauA = var_stats$TauA - full_stats$TauA,
                                  c = var_stats$c - full_stats$c,
                                  Concord = var_stats$Concord - full_stats$Concord,
                                  Discon = var_stats$Discon - full_stats$Discon,
                                  LackFit = var_stats$LackFit - full_stats$LackFit,
                                  Lift_Index = var_stats$Lift_Index - full_stats$Lift_Index,
                                  infv = var_stats$infv - full_stats$infv,
                                  ks = var_stats$ks - full_stats$ks,
                                  sord = var_stats$sord,
                                  stringsAsFactors=F))
    fin_stats <- fin_stats[order(fin_stats$sord, -fin_stats$RelImp),]
    
    return(fin_stats)
  }
  
  # Binary model statistics
  f.stats_log <- function(y2,weight,g = 10){
    f.ks <- function(p,a) {
      a0 <- p[a==0]
      a1 <- p[a==1]
      n.a1 <- as.double(length(a1))
      n.a0 <- as.double(length(a0))
      n <- n.a1 * n.a0 / (n.a1 + n.a0)
      w <- c(a1, a0)
      d <- max(abs(cumsum(ifelse(order(w) <= n.a1, 1 / n.a1, - 1 / n.a0))))
      return(d)
    }
    
    if (is.null(weight)){
      names(y2) <- c("dv","PScore")
      y2$rank <- floor((rank(y2[,"PScore"])*g)/(nrow(y2)+1)+1)
      y3 <- ddply(y2,"rank",summarize,
                  mactual=mean(dv),
                  badcnt=sum(ifelse(dv==1,1,0)),
                  goodcnt=sum(ifelse(dv==0,1,0)))
      norm <- mean(y2$dv)
      nbad <- nrow(y2[y2$dv==1,])
      ngood <- nrow(y2[y2$dv==0,])
    } else {
      names(y2) <- c("dv","PScore","wgt")
      y2 <- y2[order(y2$PScore),]
      y2$rank <- floor((cumsum(y2$wgt)*g)/(sum(y2$wgt)+1)+1)
      y3 <- ddply(y2,"rank",summarize,
                  mactual=weighted.mean(x=dv,w=wgt),
                  badcnt=sum(ifelse(dv==1,wgt,0)),
                  goodcnt=sum(ifelse(dv==0,wgt,0)))
      norm <- weighted.mean(x=y2$dv,w=y2$wgt)
      nbad <- sum(y2[y2$dv==1,"wgt"])
      ngood <- sum(y2[y2$dv==0,"wgt"])
    }
    fix <- data.frame(rank=seq(1:g))
    y3 <- merge(y3,fix,by="rank",all=T,sort=T)
    y3$mactual <- ifelse(is.na(y3$mactual),0,y3$mactual)
    y3$badcnt <- ifelse(is.na(y3$badcnt),0,y3$badcnt)
    y3$goodcnt <- ifelse(is.na(y3$goodcnt),0,y3$goodcnt)
    
    y3$goodratio <- ifelse(y3$goodcnt * y3$badcnt == 0,(y3$goodcnt+.5)/(ngood+.5),
                           y3$goodcnt / ngood)
    y3$badratio <- ifelse(y3$goodcnt * y3$badcnt == 0,(y3$badcnt+.5)/(nbad+.5),
                          y3$badcnt / nbad)
    y3$infv <- cumsum((y3$badratio-y3$goodratio)*log(y3$badratio/y3$goodratio))
    
    obs <- xtabs(cbind(1 - y2$dv, y2$dv) ~ y2$rank)
    expect <- xtabs(cbind(1 - y2$PScore, y2$PScore) ~ y2$rank)
    chisq <- sum((obs - expect)^2/expect)
    
    f1 <- floor(y2$PScore*500)
    event <- f1[y2$dv==1]
    nonevent <- f1[y2$dv==0]
    
    conc <- sum(sapply(event, function(x) sum(x > nonevent)))
    disc <- sum(sapply(event, function(x) sum(x < nonevent)))
    cnt <- length(event) * length(nonevent)
    tied <- cnt - conc - disc
    
    keep <- data.frame(SomersD = (conc-disc)/cnt,
                       Gamma = (conc-disc)/(conc+disc),
                       TauA = (conc-disc)/(.5*length(f1)*(length(f1)-1)),
                       c = (conc+.5*tied)/cnt,
                       Concord = (conc/cnt)*100,
                       Discon = (disc/cnt)*100,
                       LackFit = 1 - pchisq(chisq, g - 2),
                       Lift_Index = (y3[y3$rank==g,"mactual"] / norm)*100,
                       infv =  y3[y3$rank==g,"infv"],
                       ks =f.ks(y2$PScore,y2$dv),
                       stringsAsFactors=F)
    
    return(keep)
  }
  
  # Master function for continuous dependent variable
  f.Model_Fin_Reg <- function(insdn,
                              regvlist, 
                              dep_var,
                              weight = NULL,
                              method = "stepwise",
                              sle = 0.05,
                              sls = 0.05){
    
    cat('Build final model \n')
    
    # Model
    if (method %in% c("stepwise","forward")){
      form <- as.formula(paste(dep_var,"~1"))
      if (length(regvlist) > 1){
        scp <- as.formula(paste("~", paste(regvlist, collapse = "+")))  
      } else {
        scp <- as.formula(paste("~", regvlist))
      }
      stp <- rxStepControl(method=method, scope = scp, stepCriterion = "SigLevel",
                           maxSigLevelToAdd = sle, minSigLevelToDrop = sls,
                           test="F", refitEachStep = F)
    } else {
      if (length(regvlist) > 1){
        form <- as.formula(paste(dep_var, "~", paste(regvlist, collapse = "+")))
      } else {
        form <- as.formula(paste(dep_var, "~", regvlist))
      }
      scp <- as.formula("~1")
      stp <- rxStepControl(method=method,scope = list(lower=scp),stepCriterion = "SigLevel",
                           minSigLevelToDrop = sls,test="F", refitEachStep = F)
    }
    
    md0 <- rxLinMod(form,
                    data = "CE5_mod.xdf",
                    pweights = weight,
                    variableSelection = stp,
                    reportProgress=0)
    nvarlist <- row.names(md0$coefficients)
    nvarlist <- nvarlist[nvarlist != "(Intercept)"]
    
    # Create file with predicted score
    rxPredict(md0,"CE5_mod.xdf",type="response",predVarNames="pred",overwrite=T,reportProgress=0)
    rxPredict(md0,"CE5_val.xdf",type="response",predVarNames="pred",overwrite=T,reportProgress=0)
    rxMerge(inData1="CE5_mod.xdf",inData2="CE5_val.xdf",outFile="CE5_scored.xdf",
            type="union",overwrite=T,xdfCompressionLevel = 1,reportProgress=0)
    
    # Get statistics
    CE5_Model_Metric_All <- rbind(f.ctl_Stats_Reg("CE5_mod.xdf", dep_var, nvarlist,"Dev"),
                                  f.ctl_Stats_Reg("CE5_val.xdf", dep_var, nvarlist,"Val"))
    save(CE5_Model_Metric_All, file = "CE5_Model_Metric_All.RData",compress=T)
    
    # Parameters
    CE5_ParameterEstimates <- CE5_Model_Metric_All[CE5_Model_Metric_All$file=="Dev" & 
                                                     CE5_Model_Metric_All$variable != "Full File" ,
                                                   c("variable","Estimate","StdErr",
                                                     "tValue","Probt","StdEst")]
    save(CE5_ParameterEstimates, file = "CE5_ParameterEstimates.RData",compress=T)
    
    # Write out scoring equation
    f.rwrite<-function(...,f='CE5_Scoring_equation.txt',app=T) cat(...,'\n',file=f,append=app)
    f.rwrite('f.score <- function(x){',app=F)
    f.rwrite('#####  Model  #####')
    
    var <- CE5_ParameterEstimates$variable
    est <- round(CE5_ParameterEstimates$Estimate,6)
    state <- ifelse(var=="(Intercept)",paste0("x$score <- (1 * (",est,") + \n"),
                    paste0("x$",var," * (",est,") + \n"))
    state[length(var)] <- paste0("x$",var[length(var)]," * (",est[length(var)],")) \n")
    
    f.rwrite(state)
    f.rwrite("return(x)}")
    
    # Performance and variable validation
    tmp <- f.gains(rxDataStep("CE5_mod.xdf",varsToKeep=c(dep_var,"pred",weight,nvarlist),
                              reportProgress=0),
                   weight=weight, binary_dv=binary_dv, vars=nvarlist)
    CE5_Gainstable_Train <- tmp[[1]]
    CE5_Profiles_Train <- tmp[[2]]
    tmp <- f.gains(rxDataStep("CE5_val.xdf",varsToKeep=c(dep_var,"pred",weight,nvarlist)
                              ,reportProgress=0),
                   weight=weight, binary_dv=binary_dv, vars=nvarlist)
    CE5_Gainstable_Test <- tmp[[1]]
    CE5_Profiles_Test <- tmp[[2]]
    save(CE5_Gainstable_Train, file = "CE5_Gainstable_Train.RData",compress=T)
    save(CE5_Profiles_Train, file = "CE5_Profiles_Train.RData",compress=T)
    save(CE5_Gainstable_Test, file = "CE5_Gainstable_Test.RData",compress=T)
    save(CE5_Profiles_Test, file = "CE5_Profiles_Test.RData",compress=T)
    rm(tmp)
    
    # Correlations
    CE5_Corr <- data.frame(rxCor(as.formula(paste("~", paste(c(dep_var,nvarlist), collapse = "+"))),
                                 "CE5_scored.xdf",reportProgress=0))
    Variable <- row.names(CE5_Corr)
    CE5_Corr <- cbind(Variable,CE5_Corr)
    save(CE5_Corr, file = "CE5_Corr.RData",compress=T)
    
    # Profiling
    CE5_Profile <- f.profile(rxDataStep("CE5_mod.xdf",
                                        varsToKeep=c(dep_var,nvarlist),
                                        reportProgress=0),
                             dep_var = dep_var,
                             vars=nvarlist,
                             num_category = fin_num_category,
                             equal_dist = fin_equal_dist)
    save(CE5_Profile, file = "CE5_Profile.RData",compress=T)
    
    # Write out variable list
    f.mod_write(nvarlist)
  }
  # Continuous model control
  f.ctl_Stats_Reg <- function(insdn, dep_var, varlist, file) {
    
    cat('Get statistics for ',file,'\n')
    
    # Full model
    cat("Full model \n")
    frm <- as.formula(paste(dep_var, "~", paste(varlist, collapse = "+")))
    md <- rxLinMod(formula = frm, data = insdn, pweights = weight, reportProgress=0)
    # Predicted values
    rxPredict(md,insdn,predVarNames="pred",overwrite=T,reportProgress=0)
    prd <- rxDataStep(insdn,varsToKeep=c(dep_var,"pred",weight),reportProgress=0)
    # Elements needed for metric calculations
    sse <- md$residual.squares
    sst <- md$y.var
    n <- nrow(prd)
    p <- nrow(md$coefficients)
    if (is.null(weight)) {mdv <- mean(prd[,dep_var])} else 
    {mdv <- weighted.mean(x=prd[,dep_var],w=prd[,weight])}
    
    full_stats <- data.frame(file = file,
                             variable = "Full File",
                             Estimate = NA,
                             StdEst = NA,
                             StdErr = NA,
                             tValue = NA,
                             Probt = NA,
                             VIF = NA,
                             RelImp = NA,
                             AIC = n*log(sse/n) + (2*p),
                             SBC = n * log(sse/n) + (p*log(n)),
                             JP = ((n + p) / n) * (sse / (n - p)),
                             Rsquare = md$r.squared,
                             AdjRsq = md$adj.r.squared,
                             RMSE = sqrt(sse / (n - p)),
                             CoeffVar = (sqrt(sse / (n - p)) * 100) / mdv,
                             PC = ((n+p)/(n-p))*(1-md$r.squared),
                             f.stats_reg(prd,weight),
                             sord = 1,
                             stringsAsFactors=F)
    
    se_var <- (md$coefficients * c(NA,apply(rxDataStep(insdn,varsToKeep=varlist,reportProgress=0),2,sd)))/
      mean(prd[,dep_var])
    # Variance inflation factors
    f.tmp <- function(lst) {
      lst$vif_wgt <- lst[[1]]
      return(lst)
    }
    vif <- vif(lm(frm,
                  rxDataStep(insdn,
                             transformFunc = f.tmp,
                             transformVars = weight,
                             reportProgress=0),
                  weights=vif_wgt))
    
    var_stats <- data.frame(file = file,
                            variable = row.names(md$coefficients),
                            Estimate = md$coefficients,
                            StdEst = se_var,
                            StdErr = md$coef.std.error,
                            tValue = md$coef.t.value,
                            Probt = md$coef.p.value,
                            VIF = c(NA,vif),
                            RelImp = abs(se_var) / sum(abs(se_var),na.rm=T),
                            stringsAsFactors=F)
    names(var_stats) <- c("file","variable","Estimate","StdEst","StdErr","tValue","Probt","VIF","RelImp")
    
    # No intercept model
    cat("No intercept model \n")
    frm <- as.formula(paste(dep_var, "~", paste(c(-1,varlist), collapse = "+")))
    md <- rxLinMod(formula = frm, data = insdn, pweights = weight, reportProgress=0)
    # Predicted values
    rxPredict(md,insdn,predVarNames="pred",overwrite=T,reportProgress=0)
    prd <- rxDataStep(insdn,varsToKeep=c(dep_var,"pred",weight),reportProgress=0)
    # Elements needed for metric calculations
    sse <- md$residual.squares
    sst <- md$y.var
    n <- md$nValidObs
    p <- nrow(md$coefficients)
    if (is.null(weight)) {mdv <- mean(prd[,dep_var])} else 
    {mdv <- weighted.mean(x=prd[,dep_var],w=prd[,weight])}
    
    tmp_stats <- data.frame(variable = "(Intercept)",
                            AIC = n*log(sse/n) + (2*p),
                            SBC = n * log(sse/n) + (p*log(n)),
                            JP = ((n + p) / n) * (sse / (n - p)),
                            Rsquare = md$r.squared,
                            AdjRsq = md$adj.r.squared,
                            RMSE = sqrt(sse / (n - p)),
                            CoeffVar = (sqrt(sse / (n - p)) * 100) / mdv,
                            PC = ((n+p)/(n-p))*(1-md$r.squared),
                            f.stats_reg(prd,weight),
                            sord = 2,
                            stringsAsFactors=F)
    
    # Models removing one variable at a time
    cat("Subset models - remove one variable at a time \n")
    for (i in 1:length(varlist)) {
      frm <- as.formula(paste(dep_var, "~", paste(varlist[-i], collapse = "+")))
      md <- rxLinMod(formula = frm, data = insdn, pweights = weight, reportProgress=0)
      # Predicted values
      rxPredict(md,insdn,predVarNames="pred",overwrite=T,reportProgress=0)
      prd <- rxDataStep(insdn,varsToKeep=c(dep_var,"pred",weight),reportProgress=0)
      # Elements needed for metric calculations
      sse <- md$residual.squares
      sst <- md$y.var
      n <- md$nValidObs
      p <- nrow(md$coefficients)
      if (is.null(weight)) {mdv <- mean(prd[,dep_var])} else 
      {mdv <- weighted.mean(x=prd[,dep_var],w=prd[,weight])}
      
      tmp_stats <- rbind(tmp_stats,
                         data.frame(variable = varlist[i],
                                    AIC = n*log(sse/n) + (2*p),
                                    SBC = n * log(sse/n) + (p*log(n)),
                                    JP = ((n + p) / n) * (sse / (n - p)),
                                    Rsquare = md$r.squared,
                                    AdjRsq = md$adj.r.squared,
                                    RMSE = sqrt(sse / (n - p)),
                                    CoeffVar = (sqrt(sse / (n - p)) * 100) / mdv,
                                    PC = ((n+p)/(n-p))*(1-md$r.squared),
                                    f.stats_reg(prd,weight),
                                    sord = 3,
                                    stringsAsFactors=F))
    }
    var_stats <- merge(var_stats,tmp_stats,by="variable")
    
    # Subtract full model metrics from sub-model metrics
    # This shows the impact of each variable
    fin_stats <- rbind(full_stats,
                       data.frame(file = var_stats$file,
                                  variable = var_stats$variable,
                                  Estimate = var_stats$Estimate,
                                  StdEst = var_stats$StdEst,
                                  StdErr = var_stats$StdErr,
                                  tValue = var_stats$tValue,
                                  Probt = var_stats$Probt,
                                  VIF = var_stats$VIF,
                                  RelImp = var_stats$RelImp,
                                  AIC = var_stats$AIC - full_stats$AIC,
                                  SBC = var_stats$SBC - full_stats$SBC,
                                  JP = var_stats$JP - full_stats$JP,
                                  Rsquare = var_stats$Rsquare - full_stats$Rsquare,
                                  AdjRsq = var_stats$AdjRsq - full_stats$AdjRsq,
                                  RMSE = var_stats$RMSE - full_stats$RMSE,
                                  CoeffVar = var_stats$CoeffVar - full_stats$CoeffVar,
                                  PC = var_stats$PC - full_stats$PC,
                                  Lift_Index = var_stats$Lift_Index - full_stats$Lift_Index,
                                  gini = var_stats$gini - full_stats$gini,
                                  infv = var_stats$infv - full_stats$infv,
                                  ks = var_stats$ks - full_stats$ks,
                                  sord = var_stats$sord,
                                  stringsAsFactors=F))
    fin_stats <- fin_stats[order(fin_stats$sord, -fin_stats$RelImp),]
    
    return(fin_stats)
  }
  # Continuous model statistics
  f.stats_reg <- function(y2,weight,g = 10){
    
    if (is.null(weight)){
      names(y2) <- c("dv","PScore")
      y2$rank <- floor((rank(y2[,"PScore"])*g)/(nrow(y2)+1)+1)
      y3 <- ddply(y2,"rank",summarize,
                  count=sum(1),
                  mactual=mean(dv),
                  mpred=mean(PScore))
      norm <- mean(y2$dv)
    } else {
      names(y2) <- c("dv","PScore","wgt")
      y2 <- y2[order(y2$PScore),]
      y2$rank <- floor((cumsum(y2$wgt)*g)/(sum(y2$wgt)+1)+1)
      y3 <- ddply(y2,"rank",summarize,
                  count=sum(wgt),
                  mactual=weighted.mean(x=dv,w=wgt),
                  mpred=weighted.mean(x=PScore,w=wgt))
      norm <- weighted.mean(x=y2$dv,w=y2$wgt)
    }
    fix <- data.frame(rank=seq(1:g))
    y3 <- merge(y3,fix,by="rank",all=T,sort=T)
    y3$count <- ifelse(is.na(y3$count),0,y3$count)
    y3$mactual <- ifelse(is.na(y3$mactual),0,y3$mactual)
    y3$mpred <- ifelse(is.na(y3$mpred),0,y3$mpred)
    
    totalbad <- sum(y3$count*(y3$count*y3$mactual/sum(y3$count*y3$mactual)))
    
    cumresp <- cumsum(y3$mactual*y3$count)
    cumpct_resp <- cumresp / sum(y3$count*y3$mactual)
    cumtotal <- cumsum(y3$count)
    cumpct_freq <- cumtotal / sum(y3$count)
    cumavg <- cumresp / cumtotal
    badrate <- y3$count * y3$mactual / sum(y3$count*y3$mactual)
    goodrate <- (1 - badrate) * y3$count / (sum(y3$count) - totalbad)
    prev_dv <- c(NA,cumpct_resp[1:(g-1)])
    prev_pct <- c(NA,cumpct_freq[1:(g-1)])
    infv <- cumsum((badrate-goodrate)*log(badrate/goodrate))
    cumbad <- cumsum(badrate)
    cumgood <- cumsum(goodrate)
    ks <- abs(cumbad-cumgood)
    gini <- 2*((cumpct_freq+prev_pct)/2-(cumpct_resp+prev_dv)/2)*(cumpct_freq-prev_pct)
    gini[1] <- 0
    
    keep <- data.frame(Lift_Index = (y3[y3$rank==g,"mactual"] / norm)*100,
                       gini = sum(gini),
                       infv = infv[g],
                       ks = max(ks),
                       stringsAsFactors=F)
    
    return(keep)
  }
  # Function to write out final variable list
  f.mod_write <- function(vars){
    # Write out list of variables to keep
    f.catwrite<-function(...,f='CE5_Varlist_Model.txt',app=F) cat(...,file=f,append=app)
    
    if (length(vars)> 1){
      temp <- data.frame(name=vars,len=nchar(vars),splt=NA,outvar=NA,stringsAsFactors=F)
      temp[1,"splt"] <- 1
      temp[1,"outvar"] <- paste("varlist_final <- c('",temp[1,"name"],"'",sep="")
      cnt <- temp[1,"len"]
      ov <- temp[1,"outvar"]
      splt <- 1
      
      for (i in 2:nrow(temp)) {
        if (is.na(temp[i,"splt"])) {
          if (cnt + temp[i,"len"] > 250) {
            splt <- splt + 1
            cnt <- temp[i,"len"]
            ov <- paste("'",temp[i,"name"],"'",sep="")
          } else
          { ov <- paste(ov,",'",temp[i,"name"],"'",sep="")
            cnt <- cnt + temp[i,"len"] + 3
          }
          temp[i,"splt"] <- splt
          temp[i,"outvar"] <- ov
        }
      }
      # Mark first and last record for strings
      temp$slast <- !duplicated( temp[, "splt" ], fromLast=T)
      # Keep lest record for each string
      temp2 <- temp[temp$slast==T,]
      temp2$outvar <- ifelse(temp2$splt==splt,
                             paste(temp2$outvar,")",sep=""),
                             paste(temp2$outvar,",",sep=""))
      
      f.catwrite(paste(temp2$outvar,'\n',sep=""))
      rm(f.catwrite,temp,cnt,ov,splt,temp2,i)
      cat(length(vars)," variables selected \n")
    } else {
      if (length(vars)==1){
        f.catwrite(paste("varlist_final <- c('",vars,')\n',sep=""))  
        cat("1 variable selected \n")
      } else {
        cat("No variables selected \n")
      }
    }
  }
  # Performance and variable validation function
  f.gains <- function(y2, weight, binary_dv, vars){
    if (is.null(weight)){
      names(y2) <- c("dv","PScore", vars)
      y2$wgt <- 1
    } else {
      names(y2) <- c("dv","PScore","wgt", vars)
    }
    y2 <- y2[order(y2$PScore,decreasing=T),]
    y2$decile <- as.factor(floor((cumsum(y2$wgt)*10)/(sum(y2$wgt)+1))+1)
    y3 <- ddply(y2,"decile",summarize,
                count=sum(wgt),
                numresp=sum(dv*wgt),  # only need if binary_dv
                mean_dv=weighted.mean(x=dv,w=wgt),
                mean_pred=weighted.mean(x=PScore,w=wgt),
                minscore=min(PScore),
                maxscore=max(PScore))
    if (binary_dv == "Y") {
      y3$avg_rate <- y3$numresp  / y3$count
      y3$lift <- (y3$avg_rate / weighted.mean(x=y2$dv, w=y2$wgt))*100
      y3$cum_index <- ((cumsum(y3$numresp)/cumsum(y3$count)) / weighted.mean(x=y2$dv, w=y2$wgt)) * 100
      y3$resp_pct <- y3$numresp / sum(y3$numresp)
      y3$cum_resp_pct <- cumsum(y3$numresp) / sum(y3$numresp)
    } else {
      y3$numresp <- NULL
      y3$lift = (y3$mean_dv / weighted.mean(x=y2$dv, w=y2$wgt))*100
    } 
    y4 <- ddply(y2, "decile", function(x) colwise(weighted.mean, w = x$wgt)(x[vars]))
    return(list(y3,y4))
  }
  # Profiling function
  f.profile <- function(x,dep_var,vars,num_category = 10,equal_dist = "N") {
    dv <- x[,dep_var]
    # Get overall average and count
    cnt <- nrow(x)
    mn <- mean(dv)
    
    # Get variables and unique counts
    f.cnt <- function(z) {nrow(table(z,useNA="ifany"))}
    cnts <- apply(x[vars],2,f.cnt)
    v1 <- vars[cnts <= num_category]
    v2 <- vars[cnts > num_category]
    
    # If number of unique values <= number of categories
    f.f2 <- function(y,num_category=10){
      f1 <- as.data.frame(table(y))
      f1$Percent <- f1$Freq / cnt
      f2 <- as.data.frame(as.table(by(dv,y,mean)))
      f3 <- cbind(f1,f2)
      f3 <- f3[,-4]
      names(f3) <- c("Category","Count","Percent","MeanDV")
      f3$Index <- (f3[,4] / mn)*100
      ov <- data.frame(Category="Overall",Count=cnt,Percent=1,MeanDV=mn,Index=100,stringsAsFactors=F)
      f4 <- rbind(ov,f3)
      return(f4)
    }
    # If number of unique values > number of categories
    f.f3 <- function(y,num_category=10){
      if (equal_dist == 'Y') {
        f1 <- as.data.frame(table(cut(y,b=num_category,include.lowest=T)))
        f2 <- as.data.frame(as.table(by(dv,cut(y,b=num_category,include.lowest=T),mean)))
      } else {
        f1 <- as.data.frame(table(cut(y,b=unique(quantile(y, seq(0, 1, 1/num_category),na.rm=T)),include.lowest=T)))
        f2 <- as.data.frame(as.table(by(dv,cut(y,b=unique(quantile(y, seq(0, 1, 1/num_category),na.rm=T)),include.lowest=T),mean)))
      }
      f1$Percent <- f1$Freq / cnt
      f3 <- cbind(f1,f2)
      f3 <- f3[,-4]
      names(f3) <- c("Category","Count","Percent","MeanDV")
      f3$Index <- (f3[,4] / mn)*100
      ov <- data.frame(Category="Overall",Count=cnt,Percent=1,MeanDV=mn,Index=100,stringsAsFactors=F)
      f4 <- rbind(ov,f3)
      return(f4)
    }
    
    # Summarize data
    if (length(v2) > 1) {
      t1<-data.frame(do.call(rbind,lapply(x[,v2],f.f3,num_category)))
      rn <- row.names(t1)
      t1$Variable <- ifelse(substr(rn,nchar(rn)-1,nchar(rn)-1)==".",
                            substr(rn,1,nchar(rn)-2),
                            substr(rn,1,nchar(rn)-3))
    } else {
      if (length(v2) == 1) {
        t1<-f.f3(x[,v2],num_category)
        rn <- row.names(t1)
        t1$Variable <- ifelse(substr(rn,nchar(rn)-1,nchar(rn)-1)==".",
                              substr(rn,1,nchar(rn)-2),
                              substr(rn,1,nchar(rn)-3))
      }
    }
    if (length(v1) > 1) {
      t2<-data.frame(do.call(rbind,lapply(x[,v1],f.f2,num_category)))
      rn <- row.names(t2)
      t2$Variable <- ifelse(substr(rn,nchar(rn)-1,nchar(rn)-1)==".",
                            substr(rn,1,nchar(rn)-2),
                            substr(rn,1,nchar(rn)-3))
      
    } else {
      if (length(v1) == 1) {
        t2<-f.f2(x[,v1],num_category)
        rn <- row.names(t2)
        t2$Variable <- ifelse(substr(rn,nchar(rn)-1,nchar(rn)-1)==".",
                              substr(rn,1,nchar(rn)-2),
                              substr(rn,1,nchar(rn)-3))
        
      }
    }
    
    if (length(v1) > 0 & length(v2) > 0) {
      tmp<-rbind(t2,t1)  
    } else {
      if (length(v1) > 0) {
        tmp<-t2
      } else {
        tmp<-t1
      }
    }
    
    tmp2 <- tmp[,c("Variable","Category","Count","Percent","MeanDV","Index")]
    tmp2$Star <- ifelse(tmp2$Index >= 110,'* (+)',
                        ifelse(tmp2$Index > 100,'  (+)',
                               ifelse(tmp2$Index <= 90,'* (-)',
                                      ifelse(tmp2$Index < 100,'  (-)','  (0)'))))
    
    return(tmp2)
  }
  
  if (toupper(binary_dv) == "Y") {
    f.Model_Fin_Logistic(insdn,
                         varlist,
                         dep_var,
                         weight,
                         method,
                         sle = fin_alpha,
                         sls = fin_alpha)
  } else {
    f.Model_Fin_Reg(insdn,
                    varlist,
                    dep_var,
                    weight,
                    method,
                    sle = fin_alpha,
                    sls = fin_alpha)
  }
  
  # Load data needed for report
  load("CE5_Model_Metric_All.RData")
  CE5_Model_Metric_All$sord <- NULL
  load("CE5_ParameterEstimates.RData")
  load("CE5_Gainstable_Train.RData")
  if (toupper(binary_dv) == "Y") {CE5_Gainstable_Train$numresp <- NULL}
  load("CE5_Profiles_Train.RData")
  load("CE5_Gainstable_Test.RData")
  if (toupper(binary_dv) == "Y") {CE5_Gainstable_Test$numresp <- NULL}
  load("CE5_Profiles_Test.RData")
  load("CE5_Corr.RData")
  load("CE5_Profile.RData")
  
  # Write out Excel report
  cat("Create Excel report \n")
  if (file.exists('CE5_Model_report.xlsx')) file.remove('CE5_Model_report.xlsx')
  mr <- loadWorkbook('CE5_Model_report.xlsx', create = TRUE)
  
  # Create formats
  csPercent = createCellStyle(mr, name = "pct")
  setDataFormat(csPercent, format = "0.00%")
  
  csPercent1 = createCellStyle(mr, name = "pct1")
  setDataFormat(csPercent1, format = "0.0%")
  
  csComma = createCellStyle(mr, name = "comma")
  setDataFormat(csComma, format = "#,##0")
  
  csDec0 = createCellStyle(mr, name = "dec0")
  setDataFormat(csDec0, format = "0")
  
  csDec2 = createCellStyle(mr, name = "dec2")
  setDataFormat(csDec2, format = "0.00")
  
  csDec4 = createCellStyle(mr, name = "dec4")
  setDataFormat(csDec4, format = "0.0000")
  
  csDec6 = createCellStyle(mr, name = "dec6")
  setDataFormat(csDec6, format = "0.000000")
  
  csLine = createCellStyle(mr, name = "line")
  setBorder(csLine, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")
  
  cslComma = createCellStyle(mr, name = "lcomma")
  setDataFormat(cslComma, format = "#,##0")
  setBorder(cslComma, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")
  
  cslPercent1 = createCellStyle(mr, name = "lpct1")
  setDataFormat(cslPercent1, format = "0.0%")
  setBorder(cslPercent1, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")
  
  cslDec2 = createCellStyle(mr, name = "ldec2")
  setDataFormat(cslDec2, format = "#,##0.00")
  setBorder(cslDec2, side = "top", type = XLC$"BORDER.MEDIUM", color = XLC$"COLOR.BLACK")
  
  ## Write out report
  #createSheet(mr, name = "Selection")
  
  allRows <- seq(length = nrow(CE5_ParameterEstimates)) + 1
  createSheet(mr, name = "Parameters")
  writeWorksheet(mr, CE5_ParameterEstimates, sheet = "Parameters", startRow = 1, startCol = 1)
  setColumnWidth(mr, sheet = "Parameters", column = c(1), width = 8192)
  setCellStyle(mr, sheet = "Parameters", row = allRows, col = 2, cellstyle = csDec4)
  setCellStyle(mr, sheet = "Parameters", row = allRows, col = 3, cellstyle = csDec4)
  setCellStyle(mr, sheet = "Parameters", row = allRows, col = 4, cellstyle = csDec4)
  setCellStyle(mr, sheet = "Parameters", row = allRows, col = 5, cellstyle = csDec4)
  setCellStyle(mr, sheet = "Parameters", row = allRows, col = 6, cellstyle = csDec4)
  
  allRows <- seq(length = nrow(CE5_Model_Metric_All)) + 1
  cls <- seq(3:21)+2
  createSheet(mr, name = "Statistics")
  writeWorksheet(mr, CE5_Model_Metric_All, sheet = "Statistics", startRow = 1, startCol = 1)
  setColumnWidth(mr, sheet = "Statistics", column = 2, width = 8192)
  setColumnWidth(mr, sheet = "Statistics", column = cls, width = 2560)
  if (toupper(binary_dv) == "Y") {
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 3, cellstyle = csDec6)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 4, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 5, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 6, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 7, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 8, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 9, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 10, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 11, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 12, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 13, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 14, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 15, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 16, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 17, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 18, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 19, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 20, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 21, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 22, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 23, cellstyle = csDec4)
  } else {
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 3, cellstyle = csDec6)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 4, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 5, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 6, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 7, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 8, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 9, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 10, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 11, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 12, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 13, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 14, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 15, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 16, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 17, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 18, cellstyle = csDec2)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 19, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 20, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Statistics", row = allRows, col = 21, cellstyle = csDec4)
  }
  
  allRows <- seq(length = nrow(CE5_Gainstable_Train)) + 3
  createSheet(mr, name = "Performance")
  writeWorksheet(mr, "Train", sheet = "Performance", startRow = 1, startCol = 1, header = F)
  writeWorksheet(mr, CE5_Gainstable_Train, sheet = "Performance", startRow = 3, startCol = 1)
  if (toupper(binary_dv) == "Y") {
    cls <- seq(2:10)+1
    setColumnWidth(mr, sheet = "Performance", column = cls, width = 2560)
    setColumnWidth(mr, sheet = "Performance", column = 11, width = 3072)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 1, cellstyle = csDec0)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 2, cellstyle = csComma)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 3, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 4, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 5, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 6, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 7, cellstyle = csPercent)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 8, cellstyle = csDec0)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 9, cellstyle = csDec0)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 10, cellstyle = csPercent)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 11, cellstyle = csPercent)
  } else {
    cls <- seq(2:7)+1
    setColumnWidth(mr, sheet = "Performance", column = cls, width = 2560)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 1, cellstyle = csDec0)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 2, cellstyle = csComma)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 3, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 4, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 5, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 6, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 7, cellstyle = csDec0)
  }
  
  allRows <- seq(length = nrow(CE5_Gainstable_Test)) + 17
  writeWorksheet(mr, "Test", sheet = "Performance", startRow = 15, startCol = 1, header = F)
  writeWorksheet(mr, CE5_Gainstable_Test, sheet = "Performance", startRow = 17, startCol = 1)
  if (toupper(binary_dv) == "Y") {
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 1, cellstyle = csDec0)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 2, cellstyle = csComma)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 3, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 4, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 5, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 6, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 7, cellstyle = csPercent)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 8, cellstyle = csDec0)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 9, cellstyle = csDec0)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 10, cellstyle = csPercent)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 11, cellstyle = csPercent)
  } else {
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 1, cellstyle = csDec0)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 2, cellstyle = csComma)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 3, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 4, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 5, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 6, cellstyle = csDec4)
    setCellStyle(mr, sheet = "Performance", row = allRows, col = 7, cellstyle = csDec0)
  }
  
  allRows <- seq(length = nrow(CE5_Profiles_Train)) + 3
  allCols <- nchar(names(CE5_Profiles_Train))
  createSheet(mr, name = "Variable Validation")
  writeWorksheet(mr, "Train", sheet = "Variable Validation", startRow = 1, startCol = 1, header = F)
  writeWorksheet(mr, CE5_Profiles_Train, sheet = "Variable Validation", startRow = 3, startCol = 1)
  setCellStyle(mr, sheet = "Variable Validation", row = allRows, col = 1, cellstyle = csDec0)
  for (i in 2:ncol(CE5_Profiles_Train)){
    setColumnWidth(mr, sheet = "Variable Validation", column = i, width = -1)
    setCellStyle(mr, sheet = "Variable Validation", row = allRows, col = i, cellstyle = csDec4)  
  }
  
  allRows <- seq(length = nrow(CE5_Profiles_Test)) + 17
  writeWorksheet(mr, "Test", sheet = "Variable Validation", startRow = 15, startCol = 1, header = F)
  writeWorksheet(mr, CE5_Profiles_Test, sheet = "Variable Validation", startRow = 17, startCol = 1)
  setCellStyle(mr, sheet = "Variable Validation", row = allRows, col = 1, cellstyle = csDec0)
  for (i in 2:ncol(CE5_Profiles_Test)){
    setCellStyle(mr, sheet = "Variable Validation", row = allRows, col = i, cellstyle = csDec4)  
  }
  
  allRows <- seq(length = nrow(CE5_Corr)) + 1
  createSheet(mr, name = "Correlations")
  writeWorksheet(mr, CE5_Corr, sheet = "Correlations", startRow = 1, startCol = 1,rownames=T)
  setColumnWidth(mr, sheet = "Correlations", column = 1, width = 8192)
  for (i in 2:ncol(CE5_Corr)){
    setColumnWidth(mr, sheet = "Correlations", column = i, width = -1)
    setCellStyle(mr, sheet = "Correlations", row = allRows, col = i, cellstyle = csDec4)  
  }
  
  allRows <- seq(length = nrow(CE5_Profile)) + 1
  rows <- as.character(CE5_Profile[,1])
  lrows <- c(NA,rows[1:(length(rows)-1)])
  chg <- which(rows != lrows)+1
  createSheet(mr, name = "Profiles")
  writeWorksheet(mr, CE5_Profile, sheet = "Profiles", startRow = 1, startCol = 1)
  setColumnWidth(mr, sheet = "Profiles", column = c(1,2), width = 8192)
  setColumnWidth(mr, sheet = "Profiles", column = c(3,4,5,6,7), width = 3072)
  if (length(chg) > 0) {
    setCellStyle(mr, sheet = "Profiles", row = chg, col = 1, cellstyle = csLine)
    setCellStyle(mr, sheet = "Profiles", row = chg, col = 2, cellstyle = csLine)  
  }
  setCellStyle(mr, sheet = "Profiles", row = allRows, col = 3, cellstyle = csComma)
  setCellStyle(mr, sheet = "Profiles", row = allRows, col = 4, cellstyle = csPercent1)
  setCellStyle(mr, sheet = "Profiles", row = allRows, col = 5, cellstyle = csDec2)
  setCellStyle(mr, sheet = "Profiles", row = allRows, col = 6, cellstyle = csDec2)
  if (length(chg) > 0) {
    setCellStyle(mr, sheet = "Profiles", row = chg, col = 3, cellstyle = cslComma)
    setCellStyle(mr, sheet = "Profiles", row = chg, col = 4, cellstyle = cslPercent1)
    setCellStyle(mr, sheet = "Profiles", row = chg, col = 5, cellstyle = cslDec2)
    setCellStyle(mr, sheet = "Profiles", row = chg, col = 6, cellstyle = cslDec2)
    setCellStyle(mr, sheet = "Profiles", row = chg, col = 7, cellstyle = csLine)  
  }
  
  saveWorkbook(mr)
  
  # Clean up
  rm(list=ls())
  file.remove("CE5_mod.xdf")
  file.remove("CE5_val.xdf")
}