#################################################
#####             Sampling Code             #####
#################################################

f.sample <- function(x,dv=lab.targ,samp.size=30000,samp.strat='Y',samp.mult=5,samp.lv=.7,
                     samp.cut=0,seed1=123456789,seed2=987654321,out="N",out.max=99) {
  ## First throw out missing target values
  f1 <- x[!is.na(x[dv]),]
  ## Determine if binary or continuous
  if (nrow(table(f1[dv]))==2) mod.typ <- 'Binary' else mod.typ <- 'Continuous'
  if (mod.typ == 'Continuous' & dv.type == 'O') mod.typ <- 'Ordinal'
  cat("Model type is", mod.typ, "\n")
  
  ## Handle outliers for continuous
  if (mod.typ=='Continuous' & out != 'N'){
    if (samp.strat == 'N') {limit <- quantile(f1[,dv],probs=out.max/100)} else
    {limit <- quantile(f1[f1[,dv]>samp.cut,dv],probs=out.max/100)}
    
    if (out == 'Y') {cat(length(f1[f1[,dv]>limit,dv]),"records over",limit,"removed\n")
                     f1 <- f1[f1[,dv]<=limit,]} else
                     {cat(length(f1[f1[,dv]>limit,dv]),"records over",limit,"recoded\n")
                      f1[,dv] <- ifelse(f1[,dv]<=limit,f1[,dv],limit)}
  }
  
  ## Create sample
  if (nrow(f1) < samp.size) cat("Cannot create sample. File has only",nrow(f1),"rows.","\n") else {
    good = 0
    set.seed(seed1)
    if (samp.strat == 'N') {                                      ## Not stratified
      f2 <- f1[sample(1:nrow(f1),samp.size,replace=F),]
      good = 1
    } else if (mod.typ == 'Binary' & samp.strat == 'Y') {         ## Binary, stratified
      zeros <- f1[f1[dv]==0,] 
      ones <- f1[f1[dv]==1,]
      if (nrow(zeros) < floor(samp.size*samp.mult/(samp.mult+1))) {
        cat("Cannot create sample. File has only",nrow(zeros),"zeros.","\n")} else
          if (nrow(ones) < floor(samp.size/(samp.mult+1))) {
            cat("Cannot create sample. File has only",nrow(ones),"ones.","\n")} else {
              good = 1
              zeros.sel <- zeros[sample(1:nrow(zeros),floor(samp.size*samp.mult/(samp.mult+1)),replace=F),]
              ones.sel <- ones[sample(1:nrow(ones),floor(samp.size/(samp.mult+1)),replace=F),]
              f2 <- rbind(ones.sel,zeros.sel)
            }
    } else if (mod.typ == 'Continuous' & samp.strat == 'Y') {     ## Continuous, stratified
      f1$stratind <- ifelse(f1[dv] <= samp.cut,1,0)
      zeros <- f1[f1$stratind==0,] 
      ones <- f1[f1$stratind==1,]
      if (nrow(zeros) < floor(samp.size*samp.mult/(samp.mult+1))) {
        cat("Cannot create sample. File has only",nrow(zeros),"above cutpoint.","\n")} else
          if (nrow(ones) < floor(samp.size/(samp.mult+1))) {
            cat("Cannot create sample. File has only",nrow(ones),"below cutpoint.","\n")} else {
              good = 1
              zeros.sel <- zeros[sample(1:nrow(zeros),floor(samp.size*samp.mult/(samp.mult+1)),replace=F),]
              ones.sel <- ones[sample(1:nrow(ones),floor(samp.size/(samp.mult+1)),replace=F),]
              f2 <- rbind(ones.sel,zeros.sel)
            }
    } else if (mod.typ == 'Ordinal' & samp.strat == 'Y') {     ## Ordinal, stratified
      if (length(dv.vals) != length(dv.cnts)) {
        cat("Cannot create sample. dv.vals and dv.cnts do not have the same number of elements.\n")
      } else {
        good = 1
        for (i in 1:length(dv.vals)) {
          f1a <- f1[f1[dv]==dv.vals[i],]
          f1a <- f1a[sample(1:nrow(f1a),dv.cnts[i],replace=F),]
          if (i==1) f2 <- f1a else f2 <- rbind(f2,f1a)
          rm(f1a)
        }
      }
    }
    if (good == 1) {
      ## Split file between learning and validation
      set.seed(seed2)
      f2$lrnval <- 'V'
      f2[sample(nrow(f2),floor(nrow(f2)*samp.lv),replace=F),]$lrnval <- 'L'
      
      ## Set weights
      if (samp.strat == 'N') {f2$wgt <- nrow(f1)/nrow(f2) 
      } else if (mod.typ == 'Binary'){
        f2$wgt <- ifelse (f2[dv]==0, nrow(f1[f1[dv]==0,])/nrow(f2[f2[dv]==0 ,]),
                          ifelse (f2[dv]==1, nrow(f1[f1[dv]==1,])/nrow(f2[f2[dv]==1 ,]),NA))
      } else if (mod.typ == 'Ordinal'){
        wgt <- merge(as.data.frame(table(f2[dv],useNA="ifany"),stringsAsFactors=F),
                     as.data.frame(table(f1[dv],useNA="ifany"),stringsAsFactors=F),
                     by="Var1",all=T)
        wgt$wgt <- wgt[,3] / wgt[,2]
        names(wgt) <- c(dv,"smp","full",lab.wgt)
        f2[,lab.wgt] <- wgt[,lab.wgt][match(f2[,dv],wgt[,dv])]
      } else {
        f2$wgt <- ifelse (f2$stratind==0, nrow(f1[f1$stratind==0,])/nrow(f2[f2$stratind==0 ,]),
                          ifelse (f2$stratind==1, nrow(f1[f1$stratind==1,])/nrow(f2[f2$stratind==1 ,]),NA))
      }
      
      ## Output summary data
      if (mod.typ == 'Binary'){
        out <- rbind(data.frame(class = '0',
                                FullCnt = nrow(f1[f1[lab.targ]==0,]),
                                FullPct = nrow(f1[f1[lab.targ]==0,]) / nrow(f1),
                                SampleCnt = nrow(f2[f2[lab.targ]==0,]),
                                SamplePct = nrow(f2[f2[lab.targ]==0,]) / nrow(f2),
                                LearnCnt = nrow(f2[f2[lab.targ]==0 & f2$lrnval=='L',]),
                                LearnPct = nrow(f2[f2[lab.targ]==0 & f2$lrnval=='L',]) / nrow(f2[f2$lrnval=='L',]),
                                ValidCnt = nrow(f2[f2[lab.targ]==0 & f2$lrnval=='V',]),
                                ValidPct = nrow(f2[f2[lab.targ]==0 & f2$lrnval=='V',]) / nrow(f2[f2$lrnval=='V',])),
                     data.frame(class = '1',
                                FullCnt = nrow(f1[f1[lab.targ]==1,]),
                                FullPct = nrow(f1[f1[lab.targ]==1,]) / nrow(f1),
                                SampleCnt = nrow(f2[f2[lab.targ]==1,]),
                                SamplePct = nrow(f2[f2[lab.targ]==1,]) / nrow(f2),
                                LearnCnt = nrow(f2[f2[lab.targ]==1 & f2$lrnval=='L',]),
                                LearnPct = nrow(f2[f2[lab.targ]==1 & f2$lrnval=='L',]) / nrow(f2[f2$lrnval=='L',]),
                                ValidCnt = nrow(f2[f2[lab.targ]==1 & f2$lrnval=='V',]),
                                ValidPct = nrow(f2[f2[lab.targ]==1 & f2$lrnval=='V',]) / nrow(f2[f2$lrnval=='V',])))
        
        out2 <- rbind(data.frame(class = 'Learning',
                                 Count = nrow(f2[f2$lrnval=='L',]),
                                 Percent = nrow(f2[f2$lrnval=='L',]) / nrow(f2),
                                 AvgDV = mean(f2[f2$lrnval=='L',dv])),
                      data.frame(class = 'Validation',
                                 Count = nrow(f2[f2$lrnval=='V',]),
                                 Percent = nrow(f2[f2$lrnval=='V',]) / nrow(f2),
                                 AvgDV = mean(f2[f2$lrnval=='V',dv])))
        
        createSheet(ff, name = "File Statistics")
        writeWorksheet(ff, out, sheet = "File Statistics", startRow = 1, startCol = 1)
        setColumnWidth(ff, sheet = "File Statistics", column = 1, width = 5120)
        setColumnWidth(ff, sheet = "File Statistics", column = c(2,3,4,5,6,7,8,9), width = 3072)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 2, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 3, cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 4, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 5, cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 6, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 7, cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 8, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 9, cellstyle = csPercent)
        writeWorksheet(ff, out2, sheet = "File Statistics", startRow = 5, startCol = 1)
        setCellStyle(ff, sheet = "File Statistics", row = c(6,7), col = c(2), cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = c(6,7), col = c(3), cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = c(6,7), col = c(4), cellstyle = csDec4)
        
      } else if (mod.typ == 'Ordinal'){
        out <- data.frame(FullCnt = as.data.frame(table(f1[,dv])),
                          FullPct = as.data.frame(prop.table(table(f1[,dv])))[,2],
                          SampleCnt = as.data.frame(table(f2[,dv]))[,2],
                          SamplePct = as.data.frame(prop.table(table(f2[,dv])))[,2],
                          LearnCnt = as.data.frame(table(f2[f2[lab.split]=='L',dv]))[,2],
                          LearnPct = as.data.frame(prop.table(table(f2[f2[lab.split]=='L',dv])))[,2],
                          ValidCnt = as.data.frame(table(f2[f2[lab.split]=='V',dv]))[,2],
                          ValidPct = as.data.frame(prop.table(table(f2[f2[lab.split]=='V',dv])))[,2],
                          stringsAsFactors=F)
        names(out) <- c(dv,"FullCnt","FullPct","SampleCnt","SamplePct","LearnCnt","LearnPct",
                        "ValidCnt","ValidPct")
        
        out2 <- rbind(data.frame(class = 'Learning',
                                 Count = nrow(f2[f2$lrnval=='L',]),
                                 Percent = nrow(f2[f2$lrnval=='L',]) / nrow(f2),
                                 AvgDV = mean(f2[f2$lrnval=='L',dv])),
                      data.frame(class = 'Validation',
                                 Count = nrow(f2[f2$lrnval=='V',]),
                                 Percent = nrow(f2[f2$lrnval=='V',]) / nrow(f2),
                                 AvgDV = mean(f2[f2$lrnval=='V',dv])))
        
        allRows <- seq(length = nrow(out)) + 1
        createSheet(ff, name = "File Statistics")
        writeWorksheet(ff, out, sheet = "File Statistics", startRow = 1, startCol = 1)
        setColumnWidth(ff, sheet = "File Statistics", column = 1, width = 5120)
        setColumnWidth(ff, sheet = "File Statistics", column = c(2,3,4,5,6,7,8,9), width = 3072)
        setCellStyle(ff, sheet = "File Statistics", row = allRows, col = 2, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = allRows, col = 3, cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = allRows, col = 4, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = allRows, col = 5, cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = allRows, col = 6, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = allRows, col = 7, cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = allRows, col = 8, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = allRows, col = 9, cellstyle = csPercent)
        writeWorksheet(ff, out2, sheet = "File Statistics", startRow = nrow(out)+4, startCol = 1)
        setCellStyle(ff, sheet = "File Statistics", row = c(nrow(out)+5,nrow(out)+6), col = c(2), cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = c(nrow(out)+5,nrow(out)+6), col = c(3), cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = c(nrow(out)+5,nrow(out)+6), col = c(4), cellstyle = csDec4)
        
      } else if (samp.strat == 'Y') {
        out <- rbind(data.frame(class = 'Below Cutpoint',
                                FullCnt = nrow(f1[f1$stratind==1,]),
                                FullPct = nrow(f1[f1$stratind==1,]) / nrow(f1),
                                SampleCnt = nrow(f2[f2$stratind==1,]),
                                SamplePct = nrow(f2[f2$stratind==1,]) / nrow(f2),
                                LearnCnt = nrow(f2[f2$stratind==1 & f2$lrnval=='L',]),
                                LearnPct = nrow(f2[f2$stratind==1 & f2$lrnval=='L',]) / nrow(f2[f2$lrnval=='L',]),
                                ValidCnt = nrow(f2[f2$stratind==1 & f2$lrnval=='V',]),
                                ValidPct = nrow(f2[f2$stratind==1 & f2$lrnval=='V',]) / nrow(f2[f2$lrnval=='V',])),
                     data.frame(class = 'Above Cutpoint',
                                FullCnt = nrow(f1[f1$stratind==0,]),
                                FullPct = nrow(f1[f1$stratind==0,]) / nrow(f1),
                                SampleCnt = nrow(f2[f2$stratind==0,]),
                                SamplePct = nrow(f2[f2$stratind==0,]) / nrow(f2),
                                LearnCnt = nrow(f2[f2$stratind==0 & f2$lrnval=='L',]),
                                LearnPct = nrow(f2[f2$stratind==0 & f2$lrnval=='L',]) / nrow(f2[f2$lrnval=='L',]),
                                ValidCnt = nrow(f2[f2$stratind==0 & f2$lrnval=='V',]),
                                ValidPct = nrow(f2[f2$stratind==0 & f2$lrnval=='V',]) / nrow(f2[f2$lrnval=='V',])))
        
        out2 <- rbind(data.frame(class = 'Learning',
                                 Count = nrow(f2[f2$lrnval=='L',]),
                                 Percent = nrow(f2[f2$lrnval=='L',]) / nrow(f2),
                                 AvgDV = mean(f2[f2$lrnval=='L',dv])),
                      data.frame(class = 'Validation',
                                 Count = nrow(f2[f2$lrnval=='V',]),
                                 Percent = nrow(f2[f2$lrnval=='V',]) / nrow(f2),
                                 AvgDV = mean(f2[f2$lrnval=='V',dv])))
        
        createSheet(ff, name = "File Statistics")
        writeWorksheet(ff, out, sheet = "File Statistics", startRow = 1, startCol = 1)
        setColumnWidth(ff, sheet = "File Statistics", column = 1, width = 5120)
        setColumnWidth(ff, sheet = "File Statistics", column = c(2,3,4,5,6,7,8,9), width = 3072)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 2, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 3, cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 4, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 5, cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 6, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 7, cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 8, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = c(2,3), col = 9, cellstyle = csPercent)
        writeWorksheet(ff, out2, sheet = "File Statistics", startRow = 5, startCol = 1)
        setCellStyle(ff, sheet = "File Statistics", row = c(6,7), col = c(2), cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = c(6,7), col = c(3), cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = c(6,7), col = c(4), cellstyle = csDec2)
        f2$stratind <- NULL
      } else  {
        out <- data.frame(class = 'File',
                          FullCnt = nrow(f1),
                          FullAvg = mean(f1[,dv]),
                          SampleCnt = nrow(f2),
                          SampleAvg = mean(f2[,dv]))
        
        out2 <- rbind(data.frame(class = 'Learning',
                                 Count = nrow(f2[f2$lrnval=='L',]),
                                 Percent = nrow(f2[f2$lrnval=='L',]) / nrow(f2),
                                 AvgDV = mean(f2[f2$lrnval=='L',dv])),
                      data.frame(class = 'Validation',
                                 Count = nrow(f2[f2$lrnval=='V',]),
                                 Percent = nrow(f2[f2$lrnval=='V',]) / nrow(f2),
                                 AvgDV = mean(f2[f2$lrnval=='V',dv])))
        
        createSheet(ff, name = "File Statistics")
        writeWorksheet(ff, out, sheet = "File Statistics", startRow = 1, startCol = 1)
        setColumnWidth(ff, sheet = "File Statistics", column = 1, width = 5120)
        setColumnWidth(ff, sheet = "File Statistics", column = c(2,3,4,5), width = 3072)
        setCellStyle(ff, sheet = "File Statistics", row = 2, col = 2, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = 2, col = 3, cellstyle = csDec2)
        setCellStyle(ff, sheet = "File Statistics", row = 2, col = 4, cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = 2, col = 5, cellstyle = csDec2)
        writeWorksheet(ff, out2, sheet = "File Statistics", startRow = 5, startCol = 1)
        setCellStyle(ff, sheet = "File Statistics", row = c(6,7), col = c(2), cellstyle = csComma)
        setCellStyle(ff, sheet = "File Statistics", row = c(6,7), col = c(3), cellstyle = csPercent)
        setCellStyle(ff, sheet = "File Statistics", row = c(6,7), col = c(4), cellstyle = csDec2)
      }
      
      if (mod.typ == 'Continuous'){
        createSheet(ff, name = "Dependent Variable")
        
        jpeg(filename = "dv.jpeg", width = 600, height = 400, bg="grey90")
        hist(f2[,dv],main=paste("Histogram of",dv),xlab=dv)
        box(which="outer", lty="solid")
        dev.off()
        createName(ff, name = "dgraph", formula = "'Dependent Variable'!$B$2")
        addImage(ff, filename = "dv.jpeg", name = "dgraph", originalSize = TRUE)
      }
      
      return(f2)
    }
  }
}

#################################################
#####          No Sampling Code             #####
#################################################
f.nosample <- function(x,dv=lab.targ,samp.lv=.7,seed2=987654321) {
  f2 <- x
  if (is.null(lab.split)) {
    set.seed(seed2)
    f2$lrnval <- 'V'
    f2[sample(nrow(f2),floor(nrow(f2)*samp.lv),replace=F),]$lrnval <- 'L'
    lab.split <<- 'lrnval'
  }
  if (is.null(lab.wgt)) {
    f2$wgt <- 1
    lab.wgt <<- 'wgt'
  }
  out2 <- rbind(data.frame(class = 'Learning',
                            Count = nrow(f2[f2[,lab.split]=='L',]),
                            Percent = nrow(f2[f2[,lab.split]=='L',]) / nrow(f2),
                            AvgDV = mean(f2[f2[,lab.split]=='L',dv])),
                 data.frame(class = 'Validation',
                            Count = nrow(f2[f2[,lab.split]=='V',]),
                            Percent = nrow(f2[f2[,lab.split]=='V',]) / nrow(f2),
                            AvgDV = mean(f2[f2[,lab.split]=='V',dv])))
  
  createSheet(ff, name = "File Statistics")
  writeWorksheet(ff, out2, sheet = "File Statistics", startRow = 5, startCol = 1)
  setColumnWidth(ff, sheet = "File Statistics", column = 1, width = 5120)
  setColumnWidth(ff, sheet = "File Statistics", column = c(2,3,4,5,6,7,8,9), width = 3072)
  setCellStyle(ff, sheet = "File Statistics", row = c(6,7), col = c(2), cellstyle = csComma)
  setCellStyle(ff, sheet = "File Statistics", row = c(6,7), col = c(3), cellstyle = csPercent)
  setCellStyle(ff, sheet = "File Statistics", row = c(6,7), col = c(4), cellstyle = csDec4)
  
  ## Determine if binary or continuous
  if (nrow(table(f2[dv]))==2) mod.typ <- 'Binary' else mod.typ <- 'Continuous'
  if (mod.typ == 'Continuous' & dv.type == 'O') mod.typ <- 'Ordinal'
  if (mod.typ == 'Continuous'){
    createSheet(ff, name = "Dependent Variable")
    
    jpeg(filename = "dv.jpeg", width = 600, height = 400, bg="grey90")
    hist(f2[,dv],main=paste("Histogram of",dv),xlab=dv)
    box(which="outer", lty="solid")
    dev.off()
    createName(ff, name = "dgraph", formula = "'Dependent Variable'!$B$2")
    addImage(ff, filename = "dv.jpeg", name = "dgraph", originalSize = TRUE)
  }
  return(f2)
}

#################################################
#####  DataSource variable recode function  #####
#################################################
f.DSVars <- function(x) {
  vars <- names(x)
  DSvars1 <- c('AU003','DS946','EQ034','EQ036','MS924','MS925','MS957','XB044','XB107','XB108',
               'XB005','XB875','EQ018','PD016','EQ014','XM915','FD027','EC001','EC002','EC003',
               'EC004','EC005','EC006','EC007','EQ015','EQ025','MS456','XM008','AU012','XB111',
               'XB112','XB113','XB114','XB115','XB116','XB117','XB118','XB119','XB120','XB121',
               'XB123','XB124','XB125','XB126','XB127','XB128')
  
  keep <- DSvars1[DSvars1 %in% vars]
  if (length(keep) > 0) {
    keepN <- paste(keep,"N",sep="")
    f.DSchg <- function(y) {as.numeric(y)}
    new <- as.data.frame(apply(x[keep],2,f.DSchg))
    colnames(new) <- keepN
    vars <- c(setdiff(vars,keep),keepN)
    x <- data.frame(x,new)[,vars]
  }
  
  DSvars2 <- c('XB040','XB043','XB101','XB102','XB103','XB105')
  keep2 <- DSvars2[DSvars2 %in% vars]
  if (length(keep2) > 0) {
    keep2N <- paste(keep2,"N",sep="")
    f.DSchg2 <- function(y) {as.numeric(ifelse(y=='Y',10,y))}
    new <- as.data.frame(apply(x[keep2],2,f.DSchg2))
    colnames(new) <- keep2N
    vars <- c(setdiff(vars,keep2),keep2N)
    x <- data.frame(x,new)[,vars]
  }
  
  if ("DS915" %in% vars) {
    tmp <- data.frame(DS915=c('01','02','03','04','05','06','07','08','09','10','11','12'),
                      DS915N=c(8000,20000,30000,42500,62500,87500,112500,137500,162500,187500,225000,275000))
    x$DS915N <- tmp$DS915N[match(x$DS915,tmp$DS915)]
    x$DS915 <- NULL
  }
  
  if ("DS922" %in% vars) {
    tmp <- data.frame(DS922=c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V'),
                      DS922N=c(12500,37500,62500,87500,112500,137500,162500,187500,212500,237500,262500,287500,325000,
                               375000,425000,475000,550000,650000,750000,850000,950000,1050000))
    x$DS922N <- tmp$DS922N[match(x$DS922,tmp$DS922)]
    x$DS922 <- NULL
  }
  
  if ("PD017" %in% vars) {
    tmp <- data.frame(PD017=c('A','B','C','D','E','F','G','H'),
                      PD017N=c(25,75,175,375,750,1750,6500,15000))
    x$PD017N <- tmp$PD017N[match(x$PD017,tmp$PD017)]
    x$PD017 <- NULL
  }
  
  if ("EQ013" %in% vars) {
    tmp <- data.frame(EQ013=c(0,1,2,3,4,5,6,7),
                      EQ013N=c(0,1,2,2,3,3,4,4))
    x$EQ013N <- tmp$EQ013N[match(x$EQ013,tmp$EQ013)]
    x$EQ013 <- NULL
  }
  
  
  if ("EQ019" %in% vars) {x$EQ019C <- substr(x$EQ019,1,1); x$EQ019 <- NULL;}
  if ("MS459" %in% vars) {x$MS459C <- substr(x$MS459,1,1); x$MS459 <- NULL;}
  
  if ("AU009" %in% vars & "MS074" %in% vars) {
    x$MS074_INF <- ifelse(is.na(x$MS074) | x$MS074==0,x$AU009,x$MS074)
    x$MS074_INF_FLAG <- ifelse(is.na(x$MS074) | x$MS074==0,1,0)
  }
  
  if ("AU016" %in% vars & "MS919" %in% vars) {
    x$MS919_INF <- ifelse(!is.na(x$MS919),x$MS919,x$AU016)
    x$MS919_INF_FLAG <- ifelse(!is.na(x$MS919),0,1)
  }
  
  if ("ET012" %in% vars & "ET022" %in% vars & "ET032" %in% vars & "ET042" %in% vars) {
    x$DV001 <- ifelse(x$ET012=='20' | x$ET022=='20' | x$ET032=='20' | x$ET042=='20','1','0')
    x$ET022 <- NULL
    x$ET032 <- NULL
    x$ET042 <- NULL
  }
  
  if ("ET014" %in% vars & "ET024" %in% vars & "ET034" %in% vars & "ET044" %in% vars) {
    x$DV002 <- ifelse(x$ET014=='A' | x$ET024=='A' | x$ET034=='A' | x$ET044=='A','1','0')
    x$DV003 <- ifelse(x$ET014 %in% c('B','C','D') | x$ET024 %in% c('B','C','D') | 
                        x$ET034 %in% c('B','C','D') | x$ET044 %in% c('B','C','D'),'1','0')
    
    x$ET024 <- NULL
    x$ET034 <- NULL
    x$ET044 <- NULL
  }
  
  if ("FT050" %in% vars) {
    x$FT050_AGE <- ifelse(!is.na(x$FT050) & x$FT050 != "0000",year(Sys.Date())-as.numeric(x$FT050),NA)
    x$FT050 <- NULL
  }
  
  if ("FT024" %in% vars) {
    x$FT024_AGE <- ifelse(!is.na(x$FT024) & x$FT024 != "0000",year(Sys.Date())-as.numeric(x$FT024),NA)
    x$FT024 <- NULL
  }
  
  if ("XB037" %in% vars) {
    x$XB037_AGE <- ifelse(!is.na(x$XB037) & x$XB037 != 0,year(Sys.Date())-as.numeric(x$XB037),NA)
    x$XB037 <- NULL
  }
  
  return(x)
}

########################################
#####  Numeric variables function  #####
########################################
f.numvars <- function(x,id,targ,split,wgt,misspct,ql=.05,qu=.95) {
  vars <- names(x)
  nvars<-setdiff(vars[sapply(x,is.numeric)],c(id,targ,split,wgt))
  cat("Starting number of numeric variables is",length(nvars),"\n")
  dv0<-x[,targ]
  
  ## Basic statistics needed for transformations
  f.test <- function(z) {sum(dv0[!is.na(z)])}
  f.stats <- function(y){
    tmp <- data.frame(
      var = names(y),
      Obs = apply(!is.na(y), 2, sum),
      MissingObs = apply(is.na(y), 2, sum),
      SumDV =  apply(y,2,f.test),
      Minimum = round(apply(y, 2, min, na.rm=TRUE), 4),
      LowerPct = round(apply(y, 2,quantile, probs=ql, na.rm=TRUE),4),
      Median = round(apply(y, 2,quantile, probs=0.50, na.rm=TRUE),4),
      UpperPct = round(apply(y, 2,quantile, probs=qu, na.rm=TRUE),4),
      Maximum = round(apply(y, 2, max, na.rm=TRUE),4)
    )
  }
  
  stats <- f.stats(x[nvars])
  ustats <- subset(stats, subset = (MissingObs / (Obs+MissingObs) >= misspct | LowerPct == UpperPct | SumDV == 0 | SumDV == Obs))
  stats <- subset(stats, subset = (MissingObs / (Obs+MissingObs) < misspct & LowerPct != UpperPct & SumDV != 0 & SumDV != Obs))
  nvars2  <- as.character(stats$var)
  mstats <- subset(stats, subset = MissingObs > 0)
  mvars  <- as.character(mstats$var)
  cat("Unusable numeric variable count is",length(nvars)-length(nvars2),"\n")
  
  ## Output statistics
  if (nrow(stats) > 0) {
    stats$SumDV <- NULL
    allRows <- seq(length = nrow(stats)) + 1
    createSheet(ff, name = "Numeric Variables")
    writeWorksheet(ff, stats, sheet = "Numeric Variables", startRow = 1, startCol = 1)
    setColumnWidth(ff, sheet = "Numeric Variables", column = 1, width = 8192)
    setColumnWidth(ff, sheet = "Numeric Variables", column = c(2,3,4,5,6,7,8), width = 3072)
    setCellStyle(ff, sheet = "Numeric Variables", row = allRows, col = 2, cellstyle = csComma)
    setCellStyle(ff, sheet = "Numeric Variables", row = allRows, col = 3, cellstyle = csComma)
    setCellStyle(ff, sheet = "Numeric Variables", row = allRows, col = 4, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Numeric Variables", row = allRows, col = 5, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Numeric Variables", row = allRows, col = 6, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Numeric Variables", row = allRows, col = 7, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Numeric Variables", row = allRows, col = 8, cellstyle = csDec2)
  }
  
  if (nrow(ustats) > 0) {
    ustats$SumDV <- NULL
    allRows <- seq(length = nrow(ustats)) + 1
    createSheet(ff, name = "Unusable Numeric Variables")
    writeWorksheet(ff, ustats, sheet = "Unusable Numeric Variables", startRow = 1, startCol = 1)
    setColumnWidth(ff, sheet = "Unusable Numeric Variables", column = 1, width = 8192)
    setColumnWidth(ff, sheet = "Unusable Numeric Variables", column = c(2,3,4,5,6,7,8), width = 3072)
    setCellStyle(ff, sheet = "Unusable Numeric Variables", row = allRows, col = 2, cellstyle = csComma)
    setCellStyle(ff, sheet = "Unusable Numeric Variables", row = allRows, col = 3, cellstyle = csComma)
    setCellStyle(ff, sheet = "Unusable Numeric Variables", row = allRows, col = 4, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Unusable Numeric Variables", row = allRows, col = 5, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Unusable Numeric Variables", row = allRows, col = 6, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Unusable Numeric Variables", row = allRows, col = 7, cellstyle = csDec2)
    setCellStyle(ff, sheet = "Unusable Numeric Variables", row = allRows, col = 8, cellstyle = csDec2)
  }
  
  ## Transforms
  f.bounds<-function(y){
    quants<-quantile(y,c(0,.01,.25,.5,.75,.99,1),na.rm=T)
    iqr <- pmax((quants[5]-quants[3]),(quants[6]-quants[5]),(quants[3]-quants[2]))
    lb <- pmin(pmax((quants[3] - (1.5*iqr)), quants[1]), quants[2])
    ub <- pmin(pmax((quants[5] + (1.5*iqr)), quants[7]), quants[6])
    lb <- ifelse(lb == ub, quants[1], lb)
    ub <- ifelse(lb == ub, quants[7], ub)
    RW = pmin(pmax(y,lb),ub)
  }
  f.transform<-function(y,sqrt_lo=0,log_lo=0.0001){
    cbind(
      RW = y,
      SR = sqrt(pmax(y,sqrt_lo,na.rm=T)),
      SQ = y^2,
      LN = log(pmax(y,log_lo,na.rm=T)),
      EP = exp(pmin(y,100,na.rm=T)),
      IV = ifelse(y==0, 0, 1/y))
  }
  f.miss<- function(y){MS <- as.numeric(is.na(y))}
  
  x2<-data.frame(do.call(cbind,lapply(x[nvars2],f.bounds)))                               ## Apply upper and lower bounds
  x2[is.na(x2)] <- matrix(stats$Median,nrow(x2),length(nvars2),byrow=T)[is.na(x2)]        ## Replace NAs
  x3<-do.call(cbind,lapply(x2,f.transform))                                               ## Transforms
  lab.trans<-apply(expand.grid(colnames(f.transform(NA)),nvars2)[,2:1],1,paste,collapse='.')
  colnames(x3)<-lab.trans                                                                 ## Label Vars
  x4<-scale(x3,center=T,scale=T)                                                          ## Standardize
  ## Eliminate variables with 0 standard deviation
  ck <- round(apply(x3,2,sd),2)
  nvars4 <- names(ck[ck>0])
  nvars4 <- nvars4[!is.na(nvars4)]
  x4 <- x4[,nvars4]
  
  ## Choose best form for each variable
  ck <- as.data.frame(cor(x4,dv0))
  ck$var <- row.names(ck)
  ck$basevar <- substr(ck$var,1,nchar(ck$var)-3)
  ck$abscor <- abs(ck$V1)
  ck <- ck[order(ck$basevar,ck$abscor),]
  ck <- ck[!duplicated( ck[, "basevar" ], fromLast=T),]
  nvars3 <- ck$var
  x4 <- x4[,nvars3]
  
  ## Create missing flags
  x5<-apply(x[mvars],2,f.miss)
  lab.miss<-paste(mvars,".MS",sep="")
  colnames(x5)<-lab.miss
  
  x6 <- cbind(x4,x5)
  x7 <- cbind(x3,x5)
  
  return(list(x7,x6))
}

#########################################
#####  Character variable function  #####
#########################################
f.charvars <- function(x,id,targ,split,wgt,misspct,maxbinpct) {
  
  vars <- names(x)
  nvars<-setdiff(vars[sapply(x,is.numeric)],c(id,targ,split,wgt))
  cvars<-setdiff(vars,c(targ,nvars,id,split,wgt))
  cat("Starting number of character variables is",length(cvars),"\n")
  dv0<-x[,targ]
  
  if (length(cvars)<=0 ) {
    f.catwrite<-function(...,f='Categorical Recodes.txt',app=T) cat(...,file=f,append=app)
    
    f.catwrite('######################################################################\n',app=F)
    f.catwrite('#####               Categorical Variable Recodes                 #####\n')
    f.catwrite('######################################################################\n')
    f.catwrite('\n')
    dat <- NA
  } else {
    ## Statistics to clean character variables
    f.max <- function(z) {max(table(z,useNA="ifany"))}
    f.cnt <- function(z) {nrow(table(z,useNA="ifany"))}
    f.test <- function(z) {sum(dv0[!is.na(z)])}
    f.cstats <- function(y){
      tmp <- data.frame(
        var = names(y),
        Obs = apply(!is.na(y), 2, sum),
        MissingObs = apply(is.na(y), 2, sum),
        SumDV =  apply(y,2,f.test),
        MaxCnt = apply(y,2,f.max),
        VarCnt = apply(y,2,f.cnt)
      )
    }
    
    cstats <- f.cstats(x[cvars])
    ucstats <- subset(cstats, subset = (MissingObs / (Obs+MissingObs) >= misspct | 
                                          MaxCnt / (Obs+MissingObs) >= maxbinpct | 
                                          SumDV == 0 | 
                                          SumDV == Obs |
                                          VarCnt >= 50))
    cstats <- subset(cstats, subset = (MissingObs / (Obs+MissingObs) < misspct & 
                                         MaxCnt / (Obs+MissingObs) < maxbinpct & 
                                         SumDV != 0 & 
                                         SumDV != Obs &
                                         VarCnt < 50))
    cvars2  <- as.character(cstats$var)
    cat("Unusable character variable count is",length(cvars)-length(cvars2),"\n")
    
    ## Output statistics
    if (nrow(cstats) > 0) {
      cstats$SumDV <- NULL
      allRows <- seq(length = nrow(cstats)) + 1
      createSheet(ff, name = "Character Variables")
      writeWorksheet(ff, cstats, sheet = "Character Variables", startRow = 1, startCol = 1)
      setColumnWidth(ff, sheet = "Character Variables", column = 1, width = 8192)
      setColumnWidth(ff, sheet = "Character Variables", column = c(2,3,4,5), width = 3072)
      setCellStyle(ff, sheet = "Character Variables", row = allRows, col = 2, cellstyle = csComma)
      setCellStyle(ff, sheet = "Character Variables", row = allRows, col = 3, cellstyle = csComma)
      setCellStyle(ff, sheet = "Character Variables", row = allRows, col = 4, cellstyle = csComma)
      setCellStyle(ff, sheet = "Character Variables", row = allRows, col = 5, cellstyle = csComma)
    }
    
    if (nrow(ucstats) > 0) {
      ucstats$SumDV <- NULL
      allRows <- seq(length = nrow(ucstats)) + 1
      createSheet(ff, name = "Unusable Character Variables")
      writeWorksheet(ff, ucstats, sheet = "Unusable Character Variables", startRow = 1, startCol = 1)
      setColumnWidth(ff, sheet = "Unusable Character Variables", column = 1, width = 8192)
      setColumnWidth(ff, sheet = "Unusable Character Variables", column = c(2,3,4,5), width = 3072)
      setCellStyle(ff, sheet = "Unusable Character Variables", row = allRows, col = 2, cellstyle = csComma)
      setCellStyle(ff, sheet = "Unusable Character Variables", row = allRows, col = 3, cellstyle = csComma)
      setCellStyle(ff, sheet = "Unusable Character Variables", row = allRows, col = 4, cellstyle = csComma)
      setCellStyle(ff, sheet = "Unusable Character Variables", row = allRows, col = 5, cellstyle = csComma)
    }
    
    ##  Categorical recodes
    f.catwrite<-function(...,f='Categorical Recodes.txt',app=T) cat(...,file=f,append=app)
    
    f.catwrite('######################################################################\n',app=F)
    f.catwrite('#####               Categorical Variable Recodes                 #####\n')
    f.catwrite('######################################################################\n')
    f.catwrite('\n')
    
    if (nrow(cstats) > 0) {
      ### Imputation
      impute_char<-function(pred,value="Unk"){
        if (is.character(pred))  {pred[is.na(pred)]<- value}
        if (is.factor(pred)) {
          pred <- factor(pred, levels = c(levels(pred), value))
          pred[is.na(pred)]<- value
        }
        return(pred)
      }
      dat<-data.frame(apply(x[,cvars2],2,impute_char))
      
      mincount <- nrow(x) * .05
      
      for (j in 1:length(cvars2)) {
        varname <- dat[,cvars2[j]]
        temp1 <- data.frame(
          freq = as.data.frame(table(dat[,cvars2[j]])),
          mean = as.data.frame(aggregate(dv0, by=data.frame(varname), mean)),
          var = as.data.frame(aggregate(dv0, by=data.frame(varname), var)))
        temp1 <- temp1[c(1,2,4,6)]
        colnames(temp1)<-c("value","freq","mean","var") 
        temp1$var[is.na(temp1$var)] <- 0
        temp1 <- temp1[order(temp1$mean, decreasing=T),]
        temp1$group <- NA
        temp1[1,5] <- 1
        temp1$gfreq <- NA
        temp1[1,6] <- temp1[1,2]
        
        grpfreq <- temp1[1,2]
        grpmean <- temp1[1,3]
        grpvar <- temp1[1,4]
        grp <- 1
        
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
        temp1 <- temp1[temp1$group < grp,]
        temp1$first <- !duplicated( temp1[, "group" ]) 
        temp1$last <- !duplicated( temp1[, "group" ], fromLast=T)
        f.catwrite("### Recode: ",cvars2[j],"\n")
        temp1$OutVar <- ifelse(temp1$first & temp1$last,paste("dat$",cvars2[j],".",temp1$group," <- as.numeric(dat$",cvars2[j]," %in% c('",temp1$value,"'))\n",sep=""),
                               ifelse(temp1$first,paste("dat$",cvars2[j],".",temp1$group," <- as.numeric(dat$",cvars2[j]," %in% c('",temp1$value,"'",sep=""),
                                      ifelse(temp1$last,paste(",'",temp1$value,"'))\n",sep=""), paste(",'",temp1$value,"'",sep=""))))
        f.catwrite(temp1$OutVar)
        f.catwrite('\n')
      }
      source("Categorical Recodes.txt", local=T,echo = F)
      cvars3 <- setdiff(names(dat),cvars2)
      dat <- dat[cvars3]
    } else {dat <- NA}
  }
  
  
  
  return(dat)
}

#####################################################
#####   Remove univariate correlation problems  #####
#####################################################
f.cleancorr <- function(x,dv,maxcorr=.75) {
  cor1 <- cor(x)                  ## Correlations between independent and w/ dependent
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
    keep <- c(dv, colnames(cor3[rsums==0,csums==0]),stored.names[is.na(drop)])  ## final list of variables
  } else {keep <- colnames(cor2)}
  return(keep)
}

#################################################
#####              Model Prep               #####
#################################################
f.modelprep<-function(m){
  f.sym<-function(x) x*lower.tri(x,diag=T)+t(x*lower.tri(x))
  b <- coef(m)[-1]
  Sb <- cov.shrink(vcov(m)[-1,-1])
  h <- solve(Sb,b)
  H <- h %o% h
  FI <- solve(Sb)
  FI <- f.sym(FI)
  attr(FI,'FisherI')<-TRUE
  list(mat=FI,H=H,r=1,call=match.call())
}

#################################################
#####          Variable Reduction           #####
#################################################
f.var_redu <- function(ds,dv,mc=.75,tol=.000001) {
  cat("Starting number of variables is",ncol(ds)-1,"\n")
  ##  Remove univariate correlation problems and multicolliniarity
  k1 <- f.cleancorr(x=ds, dv=dv, maxcorr=mc)
  cat("Univariate correlation cleanup removed",ncol(ds)-length(k1) ,"variables","\n")
  trim <- trim.matrix(cor(ds[,setdiff(k1,dv)]),tol)
  k2 <- row.names(as.data.frame(trim["trimmedmat"]))
  ds2 <- ds[,c(dv,k2)]
  cat("Multicolliniarity cleanup removed",length(k1)-ncol(ds2) ,"variables","\n")
  
  if (length(k2) > 150) {
    ##  Build preliminary model
    f.formula<-function(y,x) as.formula(paste(y, "~", paste(x, collapse = "+")))
    form <- f.formula(dv, k2)
    mod1 <- biglm(form, ds2)
    
    ##  Extract top 400 variables  
    t1 <- summary(mod1)
    t2 <- t1[[2]][-1,5]
    t3 <- t2[order(abs(t2))]
    k3 <- names(t3[1:min(400,length(t3))])
    cat("High level variable reduction eliminated",length(t3)-length(k3),"variables","\n")
    
    ##  Build model on top 400 variables  
    form2 <- f.formula(dv, k3)
    mod2 <- biglm(form2, ds2)
    x2 <- f.modelprep(mod2)
    
    ##  Run subsetting options
    z2 <- anneal(x2$mat,H=x2$H,r=x2$r,kmin=150,kmax=150,criterion="Wald")
    z3 <- improve(x2$mat,H=x2$H,r=x2$r,kmin=150,kmax=150,criterion="Wald")
    z4 <- genetic(x2$mat,H=x2$H,r=x2$r,kmin=150,kmax=150,criterion="Wald")
    
    ##  Select subset of variables
    best1 <- data.frame(var=k3[1:150],lm=rep(1,150))
    best2 <- data.frame(var=k3[t(z2$bestsets)],anneal=rep(1,150))
    best3 <- data.frame(var=k3[t(z3$bestsets)],improve=rep(1,150))
    best4 <- data.frame(var=k3[t(z4$bestsets)],genetic=rep(1,150))
    t1 <- merge(best1, best2, by = "var", all=TRUE)
    t2 <- merge(best3, best4, by = "var", all=TRUE)
    t3 <- merge(t1, t2, by = "var", all=TRUE)
    f.cnt <- function(x) {sum(x,na.rm=T)}
    t3$cnt <- apply(t3[2:5],1,f.cnt)
    k4 <- as.character(t3[t3$cnt>=3,1])
  } else {k4 <- k2}
  cat("Final number of variables is",length(k4),"\n")
  return( k4)
}

#################################################
#####          Backward Selection           #####
#################################################
f.backward <- function(dat,dv,id,split='lrnval',wgt='wgt',sls=.05) {
  f.formula<-function(y,x) as.formula(paste(y, "~", paste(x, collapse = "+")))
  newvars <- setdiff(names(dat),c(dv,id,split,wgt))
  cnt <- 1
  stp <- 0
  while (cnt == 1){
    form <- f.formula(dv, newvars)
    mod1 <- biglm(form, dat)
    t1 <- summary(mod1)
    t2 <- t1[[2]][-1,5]
    t3 <- t2[order(abs(t2),decreasing=T)]
    if (t3[1] > sls) {
      d <- names(t3[1])
      newvars <- setdiff(newvars,d)
      stp <- stp + 1
      cat("Step # ",stp," dropped ",d," with p-value of ",t3[1],"\n")
    } else cnt = 0
  }
  newvars <- newvars[!is.na(newvars)]
  return(newvars)
}

#################################################
#####           Forward Selection           #####
#################################################
f.forward <- function(dat,dv,id,split='lrnval',wgt='wgt',sls=.05,maxsteps=20) {
  f.formula<-function(y,x) as.formula(paste(y, "~", paste(x, collapse = "+")))
  newvars <- setdiff(names(dat),c(dv,id,split,wgt))
  
  ## Add first variable
  cont <- 1
  minp <- sls
  minvar <- NA
  for (i in 1:length(newvars)) {
    form <- f.formula(dv,newvars[i])
    mod1 <- biglm(form, dat)
    t1 <- summary(mod1)
    t2 <- t1[[2]][-1,5]
    if (t2 < minp) {
      minp <- t2
      minvar <- newvars[i]
    } 
  }
  if (!is.na(minvar)){
    modvars <- minvar
    newvars <- setdiff(newvars,minvar)
    cat("Step # 1 added ",minvar," with p-value of ",minp,"\n")
  } else cont <- 0
  
  ## Add additional variables to maxsteps
  for (j in 2:maxsteps) {
    if (cont == 1){
      minp <- sls
      minvar <- NA
      for (i in 1:length(newvars)) {
        form <- f.formula(dv,c(modvars,newvars[i]))
        mod1 <- biglm(form, dat)
        t1 <- summary(mod1)
        t2 <- t1[[2]][-1,5]
        if (!is.na(t2[j]) & t2[j] < minp) {
          minp <- t2[j]
          minvar <- newvars[i]
        } 
      }
      if (!is.na(minvar)){
        modvars <- c(modvars,minvar)
        newvars <- setdiff(newvars,minvar)
        cat("Step #",j," added ",minvar," with p-value of ",minp,"\n")
      } else cont <- 0
    }
  }
  modvars <- modvars[!is.na(modvars)]
  return(modvars)
}

#######################################################
#####  Code to create continuous and categorical  #####
#####  recodes code to score new file. Will run   #####
#####  after final model creation.                #####
#######################################################

f.Recodes <- function(fvars,df1=dat3, df2=dat4, cvars) {
  
  final.numb <- setdiff(fvars,cvars)
  final.char <- setdiff(fvars,final.numb)
  base.numb <- substr(final.numb,1,nchar(final.numb)-3)
  base.char <- substr(final.char,1,nchar(final.char)-2)
  
  f.rwrite<-function(...,f='Final Recodes.txt',app=T) cat(...,'\n',file=f,append=app)
  f.rwrite('f.recode <- function(x){',app=F)
  f.rwrite('######################################################################')
  f.rwrite('#####                DataSource Variable Recodes                 #####')
  f.rwrite('######################################################################')
  
  DSvars1 <- c('AU003N','DS946N','EQ034N','EQ036N','MS924N','MS925N','MS957N','XB044N','XB107N','XB108N',
               'XB005N','XB875N','EQ018N','PD016N','EQ014N','XM915N','FD027N','EC001N','EC002N','EC003N',
               'EC004N','EC005N','EC006N','EC007N','EQ015N','EQ025N','MS456N','XM008N','AU012N','XB111N',
               'XB112N','XB113N','XB114N','XB115N','XB116N','XB117N','XB118N','XB119N','XB120N','XB121N',
               'XB123N','XB124N','XB125N','XB126N','XB127N','XB128N')
  keep <- DSvars1[DSvars1 %in% base.numb]
  if (length(keep) > 0) {
    for (i in 1:length(keep)) {
      f.rwrite(paste("x$",keep[i]," <- as.numeric(x$",substr(keep[i],1,nchar(keep[i])-1),")\n",sep=""))
    }
  }
  
  DSvars2 <- c('XB040N','XB043N','XB101N','XB102N','XB103N','XB105N')
  keep2 <- DSvars2[DSvars2 %in% base.numb]
  if (length(keep2) > 0) {
    for (i in 1:length(keep2)) {
      f.rwrite(paste("x$",keep2[i]," <- as.numeric(ifelse(x$",substr(keep2[i],1,nchar(keep[i])-1),"=='Y',10,x$",substr(keep2[i],1,nchar(keep[i])-1),"))\n",sep=""))
    }
  }
  
  if ("DS915N" %in% base.numb) {
    f.rwrite("tmp <- data.frame(DS915=c('01','02','03','04','05','06','07','08','09','10','11','12'),")
    f.rwrite("                  DS915N=c(8000,20000,30000,42500,62500,87500,112500,137500,162500,187500,225000,275000))")
    f.rwrite("x$DS915N <- tmp$DS915N[match(x$DS915,tmp$DS915)]")
    f.rwrite("rm(tmp)\n")
  }
  
  if ("DS922N" %in% base.numb) {
    f.rwrite("tmp <- data.frame(DS922=c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V'),")
    f.rwrite("                  DS922N=c(12500,37500,62500,87500,112500,137500,162500,187500,212500,237500,262500,287500,325000,")
    f.rwrite("                           375000,425000,475000,550000,650000,750000,850000,950000,1050000))")
    f.rwrite("x$DS922N <- tmp$DS922N[match(x$DS922,tmp$DS922)]")
    f.rwrite("rm(tmp)\n")
  }
  
  if ("PD017N" %in% base.numb) {
    f.rwrite("tmp <- data.frame(PD017=c('A','B','C','D','E','F','G','H'),")
    f.rwrite("                  PD017N=c(25,75,175,375,750,1750,6500,15000))")
    f.rwrite("x$PD017N <- tmp$PD017N[match(x$PD017,tmp$PD017)]")
    f.rwrite("rm(tmp)\n")
  }
  
  if ("EQ013N" %in% base.numb) {
    f.rwrite("tmp <- data.frame(EQ013=c(0,1,2,3,4,5,6,7),")
    f.rwrite("                  EQ013N=c(0,1,2,2,3,3,4,4))")
    f.rwrite("x$EQ013N <- tmp$EQ013N[match(x$EQ013,tmp$EQ013)]")
    f.rwrite("rm(tmp)\n")
  }
  
  if ("EQ019C" %in% base.char) f.rwrite("x$EQ019C <- substr(x$EQ019,1,1)\n")
  if ("MS459C" %in% base.char) f.rwrite("x$MS459C <- substr(x$MS459,1,1)\n")
  if ("MS074_INF" %in% base.numb) f.rwrite("x$MS074_INF <- ifelse(is.na(x$MS074) | x$MS074==0,x$AU009,x$MS074)\n")
  if ("MS074_INF_FLAG" %in% base.numb) f.rwrite("x$MS074_INF_FLAG <- ifelse(is.na(x$MS074) | x$MS074==0,1,0)\n")
  if ("DV001" %in% base.numb) f.rwrite("x$DV001 <- ifelse(x$ET012=='20' | x$ET022=='20' | x$ET032=='20' | x$ET042=='20','1','0')\n")
  if ("DV002" %in% base.numb) f.rwrite("x$DV002 <- ifelse(x$ET014=='A' | x$ET024=='A' | x$ET034=='A' | x$ET044=='A','1','0')\n")
  if ("DV003" %in% base.numb) {
    f.rwrite("x$DV003 <- ifelse(x$ET014 %in% c('B','C','D') | x$ET024 %in% c('B','C','D') | ")
    f.rwrite("                    x$ET034 %in% c('B','C','D') | x$ET044 %in% c('B','C','D'),'1','0')\n")
  }
  if ("FT050_AGE" %in% base.numb) f.rwrite("x$FT050_AGE <- ifelse(!is.na(x$FT050) & x$FT050 != '0000',year(Sys.Date())-as.numeric(x$FT050),NA)\n")
  if ("FT024_AGE" %in% base.numb) f.rwrite("x$FT024_AGE <- ifelse(!is.na(x$FT024) & x$FT024 != '0000',year(Sys.Date())-as.numeric(x$FT024),NA)\n")
  if ("XB037_AGE" %in% base.numb) f.rwrite("x$XB037_AGE <- ifelse(!is.na(x$XB037) & x$XB037 != 0,year(Sys.Date())-as.numeric(x$XB037),NA)\n")
  
  ###############################
  #####  Numeric Variables  #####
  ###############################
  
  ## Basic statistics needed for transformations
  f.fstats <- function(y){
    tmp <- data.frame(
      var = names(y),
      Minimum = round(apply(y, 2, min, na.rm=TRUE), 4),
      P1 = round(apply(y, 2,quantile, probs=0.01, na.rm=TRUE),4),
      P25 = round(apply(y, 2,quantile, probs=0.25, na.rm=TRUE),4),
      Median = round(apply(y, 2,quantile, probs=0.50, na.rm=TRUE),4),
      P75 = round(apply(y, 2,quantile, probs=0.75, na.rm=TRUE),4),
      P99 = round(apply(y, 2,quantile, probs=0.99, na.rm=TRUE),4),
      Maximum = round(apply(y, 2, max, na.rm=TRUE),4)
    )
  }
  f.fstats2 <- function(y){
    tmp <- data.frame(
      var2 = names(y),
      Mean = round(apply(y, 2, mean, na.rm=TRUE), 4),
      StDev = round(apply(y, 2, sd, na.rm=TRUE), 4)
    )
  }
  
  fstats <- f.fstats(df1[unique(base.numb)])
  fstats2 <- f.fstats2(as.data.frame(df2[,final.numb]))
  fstats2$var <- substr(fstats2$var2,1,nchar(as.vector(fstats2$var2))-3)
  fstats <- merge(fstats,fstats2,by="var",all=T)
  fstats$iqr <- pmax((fstats$P75-fstats$P25),(fstats$P99-fstats$P75),(fstats$P25-fstats$P1))
  fstats$lb <- pmin(pmax((fstats$P25 - (1.5*fstats$iqr)), fstats$Minimum), fstats$P1)
  fstats$ub <- pmin(pmax((fstats$P75 + (1.5*fstats$iqr)), fstats$Maximum), fstats$P99)
  fstats$lb <- ifelse(fstats$lb == fstats$ub, fstats$Minimum, fstats$lb)
  fstats$ub <- ifelse(fstats$lb == fstats$ub, fstats$Maximum, fstats$ub)
  fstats$type <- substr(fstats$var2,nchar(as.vector(fstats$var2))-1,nchar(as.vector(fstats$var2)))
  
  f.trans <- function(y) {
    tmp <- ifelse (y["type"] == 'SR', paste("x$",y["var2"]," <- sqrt(pmax(x$",y["var"],".RW,0,na.rm=T))",sep=""),
             ifelse (y["type"] == 'SQ', paste("x$",y["var2"]," <- x$",y["var"],".RW^2",sep=""),
               ifelse (y["type"] == 'LN', paste("x$",y["var2"]," <- log(pmax(x$",y["var"],".RW,.00001,na.rm=T))",sep=""),
                 ifelse(y["type"] == 'EP', paste("x$",y["var2"]," <- exp(pmin(x$",y["var"],".RW,100,na.rm=T))",sep=""),
                   ifelse(y["type"] == 'IV', paste("x$",y["var2"]," <- ifelse(x$",y["var"],".RW==0,0,1/x$",y["var"],".RW)",sep=""),
                     ifelse (y["type"] == 'MS', paste("x$",y["var2"]," <- as.numeric(is.na(x$",y["var"],"))",sep=""),
                                   ""))))))

  }
  fstats$trans <- apply(fstats,1,f.trans)
  
  ##  Create code to do transformations
  
  f.rwrite('######################################################################')
  f.rwrite('#####                Continuous Variable Recodes                 #####')
  f.rwrite('######################################################################')
  f.rwrite('')
  
  f.cont <- function(y){
    f.rwrite(paste("### Recode: ",y["var"],sep="")) 
    f.rwrite(paste("x$",y["var"],".RW <- x$",y["var"],sep=""))                                           # Create new variable
    f.rwrite(paste("x$",y["var"],".RW[is.na(x$",y["var"],".RW)] <- ",y["Median"],sep=""))                # Set missing values to median
    f.rwrite(paste("x$",y["var"],".RW <- pmin(pmax(x$",y["var"],".RW,",y["lb"],"),",y["ub"],")",sep="")) # Set variable bounds
    f.rwrite(y["trans"]) 
    f.rwrite(paste("x$",y["var2"]," <- (x$",y["var2"]," - ",y["Mean"],") / ",y["StDev"],sep=""))         # transform
    f.rwrite('')
  }   
  f.cont2 <- function(y){
    f.rwrite(paste("### Missing Flag: ",y["var"],sep="")) 
    f.rwrite(y["trans"]) 
    f.rwrite('')
  }   
  
  sub1 <- subset(fstats, subset = type != 'MS')
  sub2 <- subset(fstats, subset = type == 'MS')
  if (nrow(sub1) > 0) apply(sub1, 1, f.cont)
  if (nrow(sub2) > 0) apply(sub2, 1, f.cont2)
  
  for (i in 1:length(final.numb)) {
    if (i == 1) {nlist <- paste("final.numb <- c('",final.numb[i],"'",sep="")} else
    {nlist <- paste(nlist,",'",final.numb[i],"'",sep="")}
    if (i == length(final.numb)) {nlist <- paste(nlist,")",sep="")}
  }
  f.rwrite(nlist)  
  
  #################################
  #####  Character Variables  #####
  #################################
  if (length(final.char) > 0){
    basechar <- paste("### Recode:  ",substr(final.char,1,regexpr(".$",final.char)-2)," ",sep="")
    recodes <- read.fwf(file="Categorical Recodes.txt",width=1000,comment.char="",stringsAsFactors=F)
    state <- na.omit(recodes$V1)
    var <- substr(state,regexpr("dat",state)+4,regexpr(" <- ",state)-1)
    f1 <- ifelse(var %in% final.char,2, ifelse(state %in% basechar,1,0))
    state2 <- na.omit(ifelse(f1 == 1,paste("\n",state,sep=""),
                             ifelse(f1==2,gsub("dat","x",state),NA)))
    state2 <- gsub("'Unk'","NA",state2)
    
    ##  Categorical recodes
    f.rwrite('######################################################################')
    f.rwrite('#####               Categorical Variable Recodes                 #####')
    f.rwrite('######################################################################')
    
    f.rwrite(paste(state2,'\n',sep=""))
    
    for (i in 1:length(final.char)) {
      if (i == 1) {clist <- paste("final.char <- c('",final.char[i],"'",sep="")} else
      {clist <- paste(clist,",'",final.char[i],"'",sep="")}
      if (i == length(final.char)) {clist <- paste(clist,")",sep="")}
    }
    f.rwrite(clist)
  } else f.rwrite("final.char <- NA")
  
  f.rwrite("return(list(x,final.numb=final.numb,final.char=final.char))}")
}

#################################################
#####            Model Creation             #####
#################################################
f.model <- function(x,y,z,dv,id,split='lrnval',wgt='wgt',modtyp='forward',sls=.05,maxsteps=20,cvars) {
  if (nrow(table(x[dv]))==2) mod.typ <- 'Binary' else mod.typ <- 'Continuous'
  lrn <- x[x[split]=='L',]
  val <- x[x[split]=='V',]
  if (modtyp == 'forward') finvars <- f.forward(lrn, dv, id, split, wgt, sls, maxsteps) else
    finvars <- f.backward(lrn, dv, id, split, wgt, sls)

  wgts <- lrn[,wgt]
  f.formula<-function(y,x) as.formula(paste(y, "~", paste(x, collapse = "+")))
  form <- f.formula(dv,finvars)
  mod1 <- lm(form, lrn, w=wgts)
  
  ## Create code to score new files
  f.Recodes(finvars,y,z,cvars)
  
  f.rwrite<-function(...,f='Scoring Code.txt',app=T) cat(...,'\n',file=f,append=app)
  f.rwrite('f.score <- function(x){',app=F)
  f.rwrite('######################################################################')
  f.rwrite('#####                Scoring Code                                #####')
  f.rwrite('######################################################################')
  cf <- data.frame(coefficients(mod1))
  cf$state <- ifelse (rownames(cf)=="(Intercept)",paste("x$PScore = ((1 * (",cf$coefficients.mod1,"))\n",sep=""),
                      paste("+ ( x$",row.names(cf),"* (",cf$coefficients.mod1,"))\n",sep=""))
  f.rwrite(cf$state)
  f.rwrite(")\n")
  f.rwrite("return(x)}")
  
  ## Summary data
  ck <- summary(mod1)
  coef <- as.data.frame(coefficients(ck))
  coef$variable <- row.names(coef)
  coef <- coef[,c(5,1,2,3,4)]
  stats <- data.frame(
    RSquare = ck[9],
    AdjRSquare = ck[10],
    Sigma = ck[7],
    F = ck[11],
    LogLik =logLik(mod1),
    AIC = AIC(mod1)
  )
  
  allRows <- seq(length = nrow(coefficients(ck))) + 1
  newrow <- nrow(coefficients(ck)) + 3
  createSheet(fm, name = "Summary")
  writeWorksheet(fm, coef, sheet = "Summary", startRow = 1, startCol = 1)
  writeWorksheet(fm, stats[1,], sheet = "Summary", startRow = newrow, startCol = 2)
  setColumnWidth(fm, sheet = "Summary", column = 1, width = 7680)
  setColumnWidth(fm, sheet = "Summary", column = c(2,3,4,5,6,7), width = 3072)
  setCellStyle(fm, sheet = "Summary", row = allRows, col = 2, cellstyle = csDec6)
  setCellStyle(fm, sheet = "Summary", row = allRows, col = 3, cellstyle = csDec4)
  setCellStyle(fm, sheet = "Summary", row = allRows, col = 4, cellstyle = csDec4)
  setCellStyle(fm, sheet = "Summary", row = allRows, col = 5, cellstyle = csDec6)
  setCellStyle(fm, sheet = "Summary", row = newrow+1, col = 2, cellstyle = csDec4)
  setCellStyle(fm, sheet = "Summary", row = newrow+1, col = 3, cellstyle = csDec4)
  setCellStyle(fm, sheet = "Summary", row = newrow+1, col = 4, cellstyle = csDec4)
  setCellStyle(fm, sheet = "Summary", row = newrow+1, col = 5, cellstyle = csDec4)
  setCellStyle(fm, sheet = "Summary", row = newrow+1, col = 6, cellstyle = csDec2)
  setCellStyle(fm, sheet = "Summary", row = newrow+1, col = 7, cellstyle = csDec2)
  
  ## Base reporting files
  source("Scoring Code.txt", echo = F)
  lrn2 <- f.score(lrn)
  totcnt <- sum(lrn2[,wgt])
  lrn3 <- lrn2[order(lrn2$PScore,decreasing=T),c(dv,"PScore",wgt,finvars)]
  lrn3$decile <- as.factor(floor((cumsum(lrn3[,wgt])*10)/(totcnt+1))+1)
  
  val2 <- f.score(val)
  totcnt <- sum(val2[,wgt])
  val3 <- val2[order(val2$PScore,decreasing=T),c(dv,"PScore",wgt)]
  val3$decile <- as.factor(floor((cumsum(val3[,wgt])*10)/(totcnt+1))+1)
  
  if (mod.typ == 'Binary') {
    ## Binary performance
    lrn4 <- data.frame(
      count = as.table(by(lrn3[,wgt],lrn3$decile,sum)),
      responders = as.table(by(lrn3[lrn3[dv]==1,wgt],lrn3[lrn3[dv]==1,"decile"],sum)))
    names(lrn4)<-c("Decile","Count","responders.deciles","Responders")
    lrn4$responders.deciles <- NULL
    lrn4$ResponseRate <- lrn4$Responders / lrn4$Count
    lrn4$Lift <- (lrn4$ResponseRate / (sum(lrn4$Responders)/sum(lrn4$Count))) * 100
    lrn4$CumLift <- ((cumsum(lrn4$Responders)/cumsum(lrn4$Count)) / (sum(lrn4$Responders)/sum(lrn4$Count))) * 100
    lrn4$ResponderPct <- lrn4$Responders / sum(lrn4$Responders)
    lrn4$CumRsponderPct <- cumsum(lrn4$ResponderPct)
  
    val4 <- data.frame(
      count = as.table(by(val3[,wgt],val3$decile,sum)),
      responders = as.table(by(val3[val3[dv]==1,wgt],val3[val3[dv]==1,"decile"],sum)))
    names(val4)<-c("Decile","Count","responders.deciles","Responders")
    val4$responders.deciles <- NULL
    val4$ResponseRate <- val4$Responders / val4$Count
    val4$Lift <- (val4$ResponseRate / (sum(val4$Responders)/sum(val4$Count))) * 100
    val4$CumLift <- ((cumsum(val4$Responders)/cumsum(val4$Count)) / (sum(val4$Responders)/sum(val4$Count))) * 100
    val4$ResponderPct <- val4$Responders / sum(val4$Responders)
    val4$CumRsponderPct <- cumsum(val4$ResponderPct)
  
    allRows <- seq(length = nrow(lrn4)) + 3
    createSheet(fm, name = "Performance")
    writeWorksheet(fm, "Learning", sheet = "Performance", startRow = 1, startCol = 1, header=F)
    writeWorksheet(fm, lrn4, sheet = "Performance", startRow = 3, startCol = 1)
    setColumnWidth(fm, sheet = "Performance", column = c(1,2,3,4,5,6,7,8), width = 3072)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 2, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 3, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 4, cellstyle = csPercent1)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 5, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 6, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 7, cellstyle = csPercent0)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 8, cellstyle = csPercent0)
  
    writeWorksheet(fm, "Validation", sheet = "Performance", startRow = 1, startCol = 10, header=F)
    writeWorksheet(fm, val4, sheet = "Performance", startRow = 3, startCol = 10)
    setColumnWidth(fm, sheet = "Performance", column = c(10,11,12,13,14,15,16,17), width = 3072)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 11, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 12, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 13, cellstyle = csPercent1)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 14, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 15, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 16, cellstyle = csPercent0)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 17, cellstyle = csPercent0)
    
    jpeg(filename = "perf.jpeg", width = 600, height = 400, bg="grey90")
    barplot(rbind(lrn4$Lift,val4$Lift),
            names.arg=lrn4$Decile,
            legend.text=c("Learning","Validation"),
            beside=T,
            col=c("darkgreen","lightgreen"),
            xlab="Decile",
            ylab="Lift",
            ylim=c(0,round(max(max(lrn4$Lift),max(val4$Lift))+50,0)),
            axis.lty=1,
            main="Model Lift"
    )
    abline(100,0)
    box(which="outer", lty="solid")
    dev.off()
    createName(fm, name = "pgraph", formula = "Performance!$B$15")
    addImage(fm, filename = "perf.jpeg", name = "pgraph", originalSize = TRUE)
    
    jpeg(filename = "cperf.jpeg", width = 600, height = 400, bg="grey90")
    barplot(rbind(lrn4$CumLift,val4$CumLift),
            names.arg=lrn4$Decile,
            legend.text=c("Learning","Validation"),
            beside=T,
            col=c("darkgreen","lightgreen"),
            xlab="Decile",
            ylab="Lift",
            ylim=c(0,round(max(max(lrn4$Lift),max(val4$Lift))+50,0)),
            axis.lty=1,
            main="Cumulative Model Lift"
    )
    abline(100,0)
    box(which="outer", lty="solid")
    dev.off()
    createName(fm, name = "pgraph2", formula = "Performance!$J$15")
    addImage(fm, filename = "cperf.jpeg", name = "pgraph2", originalSize = TRUE)
    
  
  } else {
    ## Continuous performance
    lrn4 <- data.frame(
      count = as.table(by(lrn3[,wgt],lrn3$decile,sum)),
      avgact = as.table(by(lrn3[,dv],lrn3$decile,mean)),
      avgprd = as.table(by(lrn3$PScore,lrn3$decile,mean)),
      minprd = as.table(by(lrn3$PScore,lrn3$decile,min)),
      maxprd = as.table(by(lrn3$PScore,lrn3$decile,max)))
    lrn4 <- lrn4[,c(1,2,4,6,8,10)]
    names(lrn4)<-c("Decile","Count","AvgAct","AvgPred","MinPred","MaxPred")
    lrn4$Lift <- (lrn4$AvgAct / mean(lrn3[,dv])*100)
    
    val4 <- data.frame(
      count = as.table(by(val3[,wgt],val3$decile,sum)),
      avgact = as.table(by(val3[,dv],val3$decile,mean)),
      avgprd = as.table(by(val3$PScore,val3$decile,mean)),
      minprd = as.table(by(val3$PScore,val3$decile,min)),
      maxprd = as.table(by(val3$PScore,val3$decile,max)))
    val4 <- val4[,c(1,2,4,6,8,10)]
    names(val4)<-c("Decile","Count","AvgAct","AvgPred","MinPred","MaxPred")
    val4$Lift <- (val4$AvgAct / mean(val3[,dv])*100)
    
    allRows <- seq(length = nrow(lrn4)) + 3
    createSheet(fm, name = "Performance")
    writeWorksheet(fm, "Learning", sheet = "Performance", startRow = 1, startCol = 1, header=F)
    writeWorksheet(fm, lrn4, sheet = "Performance", startRow = 3, startCol = 1)
    setColumnWidth(fm, sheet = "Performance", column = c(1,2,3,4,5,6,7), width = 3072)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 2, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 3, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 4, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 5, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 6, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 7, cellstyle = csComma)
    
    writeWorksheet(fm, "Validation", sheet = "Performance", startRow = 1, startCol = 9, header=F)
    writeWorksheet(fm, val4, sheet = "Performance", startRow = 3, startCol = 9)
    setColumnWidth(fm, sheet = "Performance", column = c(9,10,11,12,13,14,15), width = 3072)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 10, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 11, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 12, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 13, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 14, cellstyle = csComma)
    setCellStyle(fm, sheet = "Performance", row = allRows, col = 15, cellstyle = csComma)
    
    jpeg(filename = "perf.jpeg", width = 600, height = 400, bg="grey90")
    barplot(rbind(lrn4$Lift,val4$Lift),
            names.arg=lrn4$Decile,
            legend.text=c("Learning","Validation"),
            beside=T,
            col=c("darkgreen","lightgreen"),
            xlab="Decile",
            ylab="Lift",
            ylim=c(0,round(max(max(lrn4$Lift),max(val4$Lift))+50,0)),
            axis.lty=1,
            main="Model Lift"
    )
    abline(100,0)
    box(which="outer", lty="solid")
    dev.off()
    createName(fm, name = "pgraph", formula = "Performance!$B$15")
    addImage(fm, filename = "perf.jpeg", name = "pgraph", originalSize = TRUE)

  }
  ## Correlations
  cor <- as.data.frame(cor(dat8[,c(dv,finvars)]))
  variable <- row.names(cor)
  cor2 <- cbind(variable,cor)
  allRows <- seq(length = nrow(cor2)) + 1
  allCols <- seq(length = nrow(cor2)) + 1
  createSheet(fm, name = "Correlation")
  writeWorksheet(fm, cor2, sheet = "Correlation", startRow = 1, startCol = 1)
  setColumnWidth(fm, sheet = "Correlation", column = 1, width = 7680)
  setColumnWidth(fm, sheet = "Correlation", column = allCols, width = 3072)
  for (i in 2:ncol(cor2)){
    setCellStyle(fm, sheet = "Correlation", row = allRows, col = i, cellstyle = csDec6)
  }
  
  ## Variable Validation
  f.val <- function(x,dec) {out <- by(x,dec,mean)}
  vval <- as.data.frame(apply(lrn3[,finvars],2,f.val,lrn3$decile))
  decile <- seq(1:10)
  vval2 <- cbind(decile,vval)
  
  createSheet(fm, name = "Variable Validation")
  for (i in 2:ncol(vval2)){
    p <- vval2[,c(1,i)]
    if (i == 2) nrow <- 1 else nrow <- ((i-2)*12)+1
    writeWorksheet(fm, p, sheet = "Variable Validation", startRow = nrow, startCol = 1)
    allRows <- seq(from=(nrow+1),length = 10)
    setCellStyle(fm, sheet = "Variable Validation", row = allRows, col = 2, cellstyle = csDec2)
    
    jpeg(filename = paste("val",i,".jpeg",sep=""), width = 330, height = 220, bg="grey90")
    par(pin=c(4,2.75))
    barplot(p[,2],
            col="darkblue",
            axes=F,
            space=.6
    )
    box(which="outer", lty="solid")
    dev.off()
    createName(fm, name = paste("val",i,sep=""), formula = paste("'Variable Validation'",idx2cref(c(nrow,5)),sep="!"))
    addImage(fm, filename = paste("val",i,".jpeg",sep=""), name = paste("val",i,sep=""), originalSize = TRUE)
  }
  setColumnWidth(fm, sheet = "Variable Validation", column = 1, width = 3072)
  setColumnWidth(fm, sheet = "Variable Validation", column = 2, width = 7680)
  
  ## Profiling
  profile <- f.profile(lrn3,dv=dv,finvars,b=5,fm)
  
  ## Residuals
  createSheet(fm, name = "Residuals")
  
  jpeg(filename = "qqnorm.jpeg", width = 600, height = 400)
  qqnorm(rstandard(mod1), ylab = "Standardized Residuals")
  qqline(rstandard(mod1))
  dev.off()
  createName(fm, name = "rgraph1", formula = "Residuals!$B$2")
  addImage(fm, filename = "qqnorm.jpeg", name = "rgraph1", originalSize = TRUE)
  
  jpeg(filename = "resid.jpeg", width = 600, height = 400)
  plot(fitted(mod1), residuals(mod1), xlab = "Fitted values",
       ylab = "Residuals", type = "n",
       ylim = max(abs(residuals(mod1))) * c(-1, 1))
  abline(h = 0, lty = 2)
  text(fitted(mod1), residuals(mod1),labels="*")
  dev.off()
  createName(fm, name = "rgraph2", formula = "Residuals!$B$23")
  addImage(fm, filename = "resid.jpeg", name = "rgraph2", originalSize = TRUE)
  
  return(mod1)
  
}

#################################################
#####              Profiling                #####
#################################################
f.profile <- function(x,dv,vars,b=5,eo){
  f.cnt <- function(y) {nrow(table(y,useNA="ifany"))}
  cnts <- data.frame(var = names(x[vars]),
                     VarCnt = apply(x[vars],2,f.cnt))
  
  catvars <- as.character(cnts[cnts$VarCnt <= b,"var"])
  cntvars <- as.character(cnts[cnts$VarCnt > b,"var"])
  
  f.f1 <- function(y,cnt,vars,dv,mn){
    for (i in 1:length(vars)){
      f1 <- data.frame(Variable=vars[i],as.data.frame(table(y[,vars[i]])))
      f1$Percent <- f1$Freq / cnt
      f2 <- as.data.frame(as.table(by(dv,y[,vars[i]],mean)))
      f3 <- cbind(f1,f2)
      if (i==1) out <- f3 else out <- rbind(out,f3)
    }
    out <- out[,-5]
    names(out) <- c("Variable","Values","Count","Percent","MeanDV")
    out$Index <- (out$MeanDV / mn)*100
    return(out)
  }
  test <- f.f1(x,nrow(x),catvars,x[,dv],mean(x[,dv]))
  
  f.f2 <- function(y,cnt,vars,dv,mn,b=5){
    for (i in 1:length(vars)){
      f1 <- data.frame(Variable=vars[i],as.data.frame(table(cut(y[,vars[i]],b=b,include.lowest=T))))
      f1$Percent <- f1$Freq / cnt
      f2 <- as.data.frame(as.table(by(dv,cut(y[,vars[i]],b=b,include.lowest=T),mean)))
      f3 <- cbind(f1,f2)
      if (i==1) out <- f3 else out <- rbind(out,f3)
    }
    out <- out[,-5]
    names(out) <- c("Variable","Values","Count","Percent","MeanDV")
    out$Index <- (out$MeanDV / mn)*100
    return(out)
  }
  test2 <- f.f2(x,nrow(x),cntvars,x[,dv],mean(x[,dv]),b=b)
  profile <- rbind(test,test2)
  profile$Star <- ifelse(profile$Index >= 110,'     * (+)',
                    ifelse(profile$Index > 100,'        (+)',
                      ifelse(profile$Index <= 90,'     * (-)',
                        ifelse(profile$Index < 100,'        (-)','        (0)'))))
  
  ## Write out report
  allRows <- seq(length = nrow(profile)) + 1
  rows <- as.character(profile[,1])
  lrows <- c(NA,rows[1:(length(rows)-1)])
  chg <- which(rows != lrows)+1
  createSheet(eo, name = "Profiling")
  writeWorksheet(eo, profile, sheet = "Profiling", startRow = 1, startCol = 1)
  setColumnWidth(eo, sheet = "Profiling", column = c(1,2), width = 7680)
  setColumnWidth(eo, sheet = "Profiling", column = c(3,4,5,6,7), width = 3072)
  setCellStyle(eo, sheet = "Profiling", row = chg, col = 1, cellstyle = csLine)
  setCellStyle(eo, sheet = "Profiling", row = chg, col = 2, cellstyle = csLine)
  setCellStyle(eo, sheet = "Profiling", row = allRows, col = 3, cellstyle = csComma)
  setCellStyle(eo, sheet = "Profiling", row = allRows, col = 4, cellstyle = csPercent1)
  setCellStyle(eo, sheet = "Profiling", row = allRows, col = 5, cellstyle = csDec2)
  setCellStyle(eo, sheet = "Profiling", row = allRows, col = 6, cellstyle = csDec2)
  setCellStyle(eo, sheet = "Profiling", row = chg, col = 3, cellstyle = cslComma)
  setCellStyle(eo, sheet = "Profiling", row = chg, col = 4, cellstyle = cslPercent1)
  setCellStyle(eo, sheet = "Profiling", row = chg, col = 5, cellstyle = cslDec2)
  setCellStyle(eo, sheet = "Profiling", row = chg, col = 6, cellstyle = cslDec2)
  setCellStyle(eo, sheet = "Profiling", row = chg, col = 7, cellstyle = csLine)
  
  return(profile)
}
