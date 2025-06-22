%let loc=/mnt/projects/locked/pst_qmgisi/Modeling/CE Rebuild/;     ** Project folder;
libname lib   "&loc.sasdata";                                      ** SAS dataset libname;
%let inds=lib.bin_file;               ** your SAS input dataset name;

ods listing close;
ods Tagsets.ExcelxP body="&loc.vars.xls" style=sasweb;

proc contents data=&inds; run;

ods Tagsets.ExcelxP close;
ods listing;
