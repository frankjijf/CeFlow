require(XLConnect)
# Set working directory to write out Excel output
setwd("/mnt/projects/locked/pst_qmgisi/Modeling/CE Rebuild/binary/R_work")

# Location of file to analyze
lib <- "/mnt/projects/locked/pst_qmgisi/Modeling/CE Rebuild/sasdata"

# File to analyze
inData <- file.path(lib, "bin_file.xdf")

### Do not make changes below this line ###
# Get variable names
tst <- rxGetVarInfo(data=inData)

# Get variable descriptions
ds <- sapply(tst,getElement,"description")
ds2 <- character(length=length(ds))
for (i in 1:length(ds)) {
  ds2[i] <- if(is.null(ds[[i]])) {NA} else {ds[[i]]}
}

# Build dataframe with variable information
contents <- data.frame(Variable=names(tst),
                       Type=sapply(tst,getElement,"varType"),
                       Label=ds2,
                       stringsAsFactors=F)

# Write out workbook with variable information
l <- max(nchar(contents$Label))*256  # column width for label field

if (file.exists('var_contents.xlsx')) file.remove('var_contents.xlsx')
vc <- loadWorkbook('var_contents.xlsx', create = TRUE)
createSheet(vc, name = "Variables")
writeWorksheet(vc, contents, sheet = "Variables", startRow = 1, startCol = 1)
setColumnWidth(vc, sheet = "Variables", column = 1, width = 8192)
setColumnWidth(vc, sheet = "Variables", column = 2, width = 2304)
setColumnWidth(vc, sheet = "Variables", column = 3, width = l)
saveWorkbook(vc)

# Clean up
rm(tst,ds,ds2,i,l,vc)
