#################################################
#####  DataSource variable recode function  #####
#################################################
f.DSVars <- function(x) {
vars <- names(x)
if ("AU003" %in% vars) {x$AU003N <- suppressWarnings(as.numeric(x$AU003))}
if ("XB005" %in% vars) {x$XB005N <- suppressWarnings(as.numeric(x$XB005))}
if ("XM915" %in% vars) {x$XM915N <- suppressWarnings(as.numeric(x$XM915))}
if ("FD027" %in% vars) {x$FD027r <- suppressWarnings(as.numeric(x$FD027))}
if ("EC001" %in% vars) {x$EC001r <- suppressWarnings(as.numeric(x$EC001))}
if ("EC002" %in% vars) {x$EC002r <- suppressWarnings(as.numeric(x$EC002))}
if ("EC003" %in% vars) {x$EC003r <- suppressWarnings(as.numeric(x$EC003))}
if ("EC004" %in% vars) {x$EC004r <- suppressWarnings(as.numeric(x$EC004))}
if ("EC005" %in% vars) {x$EC005r <- suppressWarnings(as.numeric(x$EC005))}
if ("EC006" %in% vars) {x$EC006r <- suppressWarnings(as.numeric(x$EC006))}
if ("EC007" %in% vars) {x$EC007r <- suppressWarnings(as.numeric(x$EC007))}
if ("MS456" %in% vars) {x$MS456r <- suppressWarnings(as.numeric(x$MS456))}
if ("XM008" %in% vars) {x$XM008r <- suppressWarnings(as.numeric(x$XM008))}
if ("XB044" %in% vars) {x$XB044r <- suppressWarnings(as.numeric(x$XB044))}
if ("XB107" %in% vars) {x$XB107r <- suppressWarnings(as.numeric(x$XB107))}
if ("XB108" %in% vars) {x$XB108r <- suppressWarnings(as.numeric(x$XB108))}
if ("XB111" %in% vars) {x$XB111r <- suppressWarnings(as.numeric(x$XB111))}
if ("XB112" %in% vars) {x$XB112r <- suppressWarnings(as.numeric(x$XB112))}
if ("XB113" %in% vars) {x$XB113r <- suppressWarnings(as.numeric(x$XB113))}
if ("XB114" %in% vars) {x$XB114r <- suppressWarnings(as.numeric(x$XB114))}
if ("XB115" %in% vars) {x$XB115r <- suppressWarnings(as.numeric(x$XB115))}
if ("XB116" %in% vars) {x$XB116r <- suppressWarnings(as.numeric(x$XB116))}
if ("XB117" %in% vars) {x$XB117r <- suppressWarnings(as.numeric(x$XB117))}
if ("XB118" %in% vars) {x$XB118r <- suppressWarnings(as.numeric(x$XB118))}
if ("XB119" %in% vars) {x$XB119r <- suppressWarnings(as.numeric(x$XB119))}
if ("XB120" %in% vars) {x$XB120r <- suppressWarnings(as.numeric(x$XB120))}
if ("XB121" %in% vars) {x$XB121r <- suppressWarnings(as.numeric(x$XB121))}
if ("XB123" %in% vars) {x$XB123r <- suppressWarnings(as.numeric(x$XB123))}
if ("XB124" %in% vars) {x$XB124r <- suppressWarnings(as.numeric(x$XB124))}
if ("XB125" %in% vars) {x$XB125r <- suppressWarnings(as.numeric(x$XB125))}
if ("XB126" %in% vars) {x$XB126r <- suppressWarnings(as.numeric(x$XB126))}
if ("XB127" %in% vars) {x$XB127r <- suppressWarnings(as.numeric(x$XB127))}
if ("XB128" %in% vars) {x$XB128r <- suppressWarnings(as.numeric(x$XB128))}
if ("VW526" %in% vars) {x$VW526r <- suppressWarnings(as.numeric(x$VW526))}
if ("MS496" %in% vars) {x$MS496r <- suppressWarnings(as.numeric(x$MS496))}
if ("MS112" %in% vars) {x$MS112r <- suppressWarnings(as.numeric(x$MS112))}
if ("MS127" %in% vars) {x$MS127r <- suppressWarnings(as.numeric(x$MS127))}
if ("MV012" %in% vars) {x$MV012r <- suppressWarnings(as.numeric(x$MV012))}
if ("MS052" %in% vars) {x$MS052r <- suppressWarnings(as.numeric(x$MS052))}
if ("XB040" %in% vars) {x$XB040r <- suppressWarnings(as.numeric(ifelse(x$XB040=='Y',10,x$XB040)))}
if ("XB043" %in% vars) {x$XB043r <- suppressWarnings(as.numeric(ifelse(x$XB043=='Y',10,x$XB043)))}
if ("XB101" %in% vars) {x$XB101r <- suppressWarnings(as.numeric(ifelse(x$XB101=='Y',10,x$XB101)))}
if ("XB102" %in% vars) {x$XB102r <- suppressWarnings(as.numeric(ifelse(x$XB102=='Y',10,x$XB102)))}
if ("XB103" %in% vars) {x$XB103r <- suppressWarnings(as.numeric(ifelse(x$XB103=='Y',10,x$XB103)))}
if ("XB105" %in% vars) {x$XB105r <- suppressWarnings(as.numeric(ifelse(x$XB105=='Y',10,x$XB105)))}
if ("DS915" %in% vars) {
  x$DS915N <- recode(x$DS915,"'01'=8000;'02'=20000;'03'=30000;'04'=42500;
                     '05'=62500;'06'=87500;'07'=112500;'08'=137500;'09'=162500;'10'=187500;
                     '11'=225000;'12'=275000;else=NA")
}
if ("DS922" %in% vars) {
  x$DS922N <- recode(x$DS922,"'A'=12500;'B'=37500;'C'=62500;'D'=87500;
                     'E'=112500;'F'=137500;'G'=162500;'H'=187500;'I'=212500;'J'=237500;
                     'K'=262500;'L'=287500;'M'=325000;'N'=375000;'O'=425000;'P'=475000;
                     'Q'=550000;'R'=650000;'S'=750000;'T'=850000;'U'=950000;'V'=1050000;
                     else=NA")
}
if ("PD017" %in% vars) {
  x$PD017N <- recode(x$PD017,"'A'=25;'B'=75;'C'=175;'D'=375;
                     'E'=750;'F'=1750;'G'=6500;'H'=15000;else=NA")
}
if ("MS459" %in% vars) {x$MS459C <- suppressWarnings(substr(x$MS459,1,1))}
if ("MS525" %in% vars) {x$MS525C <- suppressWarnings(substr(x$MS525,1,1))}
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
}
if ("ET014" %in% vars & "ET024" %in% vars & "ET034" %in% vars & "ET044" %in% vars) {
  x$DV002 <- ifelse(x$ET014=='A' | x$ET024=='A' | x$ET034=='A' | x$ET044=='A','1','0')
  x$DV003 <- ifelse(x$ET014 %in% c('B','C','D') | x$ET024 %in% c('B','C','D') | 
                      x$ET034 %in% c('B','C','D') | x$ET044 %in% c('B','C','D'),'1','0')
}
if ("FT050" %in% vars) {
  x$FT050_AGE <- ifelse(!is.na(x$FT050) & x$FT050 != "0000",
                        as.numeric(substr(Sys.Date(),1,4))-as.numeric(x$FT050),NA)
}
if ("FT024" %in% vars) {
  x$FT024_AGE <- ifelse(!is.na(x$FT024) & x$FT024 != "0000",
                        as.numeric(substr(Sys.Date(),1,4))-as.numeric(x$FT024),NA)
}
if ("XB036" %in% vars) {
  x$XB036_AGE <- ifelse(!is.na(x$XB036) & x$XB036 != 0,
                        as.numeric(substr(Sys.Date(),1,4))-as.numeric(x$XB036),NA)
}

if ("XB037" %in% vars) {
  x$XB037_AGE <- ifelse(!is.na(x$XB037) & x$XB037 != 0,
                        as.numeric(substr(Sys.Date(),1,4))-as.numeric(x$XB037),NA)
}
if ("VW052" %in% vars) {
  x$VW052r <- recode(x$VW052,"'A'=0;'B'=2500;'C'=7500;'D'=17500;
                     'E'=37500;'F'=75000;'G'=175000;'H'=375000;'I'=625000;
                     else=NA")
}
if ("VW107" %in% vars) {
  x$VW107r <- recode(x$VW107,"'A'=800;'B'=750;'C'=700;'D'=650;
                     'E'=600;'F'=550;'G'=500;'H'=450;else=NA")
}
if ("PD012" %in% vars) {
  x$PD012r <- recode(x$PD012,"'A'=1;'B'=2;'C'=3;'D'=4;'E'=5;'F'=6;
                     'G'=7;'H'=8;'I'=9;'J'=10;'K'=11;'L'=12;'M'=13;'N'=14;
                     'O'=15;'P'=16;'Q'=17;'R'=18;'S'=19;'T'=20;'U'=21;
                     else=NA")
}
if ("PD013" %in% vars) {
  x$PD013r <- recode(x$PD013,"'A'=1;'B'=2;'C'=3;'D'=4;'E'=5;'F'=6;
                     'G'=7;'H'=8;'I'=9;'J'=10;'K'=11;'L'=12;'M'=13;'N'=14;
                     'O'=15;'P'=16;'Q'=17;'R'=18;'S'=19;'T'=20;'U'=21;
                     else=NA")
}
if ("PD014" %in% vars) {
  x$PD014r <- recode(x$PD014,"'A'=1;'B'=2;'C'=3;'D'=4;'E'=5;'F'=6;
                     'G'=7;'H'=8;'I'=9;'J'=10;'K'=11;'L'=12;'M'=13;'N'=14;
                     'O'=15;'P'=16;'Q'=17;'R'=18;'S'=19;'T'=20;'U'=21;
                     else=NA")
}
if ("PD015" %in% vars) {
  x$PD015r <- recode(x$PD015,"'A'=1;'B'=2;'C'=3;'D'=4;'E'=5;'F'=6;
                     'G'=7;'H'=8;'I'=9;'J'=10;'K'=11;'L'=12;'M'=13;'N'=14;
                     'O'=15;'P'=16;'Q'=17;'R'=18;'S'=19;'T'=20;'U'=21;
                     else=NA")
}
if ("MV002" %in% vars) {
  x$MV002r <- recode(x$MV002,"'01'=8000;'02'=20000;'03'=30000;'04'=42500;
                     '05'=62500;'06'=87500;'07'=112500;'08'=137500;'09'=162500;'10'=187500;
                     '11'=225000;'12'=275000;else=NA")
}
if ("MV010" %in% vars) {
  x$MV010r <- recode(x$MV010,"'A'=12500;'B'=37500;'C'=62500;'D'=87500;
                     'E'=112500;'F'=137500;'G'=162500;'H'=187500;'I'=212500;'J'=237500;
                     'K'=262500;'L'=287500;'M'=325000;'N'=375000;'O'=425000;'P'=475000;
                     'Q'=550000;'R'=650000;'S'=750000;'T'=850000;'U'=950000;'V'=1050000;
                     else=NA")
}
if ("IN057" %in% vars) {
  x$IN057r <- recode(x$IN057,"'A'=5000;'B'=15000;'C'=25000;'D'=35000;
                     'E'=45000;'F'=55000;'G'=65000;'H'=75000;'I'=85000;'J'=95000;
                     'K'=105000;'L'=115000;'M'=125000;'N'=135000;'O'=145000;'P'=155000;
                     'Q'=165000;'R'=175000;'S'=185000;'T'=195000;'U'=205000;'V'=215000;
                     'W'=225000;'X'=235000;'Y'=245000;'Z'=255000;else=NA")
}
if ("IN077" %in% vars) {
  x$IN077r <- recode(x$IN077,"'A'=5000;'B'=15000;'C'=25000;'D'=35000;
                     'E'=45000;'F'=55000;'G'=65000;'H'=75000;'I'=85000;'J'=95000;
                     'K'=105000;'L'=115000;'M'=125000;'N'=135000;'O'=145000;'P'=155000;
                     'Q'=165000;'R'=175000;'S'=185000;'T'=195000;'U'=205000;'V'=215000;
                     'W'=225000;'X'=235000;'Y'=245000;'Z'=255000;else=NA")
}
if ("IN097" %in% vars) {
  x$IN097r <- recode(x$IN097,"'A'=5000;'B'=15000;'C'=25000;'D'=35000;
                     'E'=45000;'F'=55000;'G'=65000;'H'=75000;'I'=85000;'J'=95000;
                     'K'=105000;'L'=115000;'M'=125000;'N'=135000;'O'=145000;'P'=155000;
                     'Q'=165000;'R'=175000;'S'=185000;'T'=195000;'U'=205000;'V'=215000;
                     'W'=225000;'X'=235000;'Y'=245000;'Z'=255000;else=NA")
}
if ("IN117" %in% vars) {
  x$IN117r <- recode(x$IN117,"'A'=5000;'B'=15000;'C'=25000;'D'=35000;
                     'E'=45000;'F'=55000;'G'=65000;'H'=75000;'I'=85000;'J'=95000;
                     'K'=105000;'L'=115000;'M'=125000;'N'=135000;'O'=145000;'P'=155000;
                     'Q'=165000;'R'=175000;'S'=185000;'T'=195000;'U'=205000;'V'=215000;
                     'W'=225000;'X'=235000;'Y'=245000;'Z'=255000;else=NA")
}
if ("IN137" %in% vars) {
  x$IN137r <- recode(x$IN137,"'A'=5000;'B'=15000;'C'=25000;'D'=35000;
                     'E'=45000;'F'=55000;'G'=65000;'H'=75000;'I'=85000;'J'=95000;
                     'K'=105000;'L'=115000;'M'=125000;'N'=135000;'O'=145000;'P'=155000;
                     'Q'=165000;'R'=175000;'S'=185000;'T'=195000;'U'=205000;'V'=215000;
                     'W'=225000;'X'=235000;'Y'=245000;'Z'=255000;else=NA")
}

if ("XB151" %in% vars) {x$XB151r <- suppressWarnings(as.numeric(x$XB151))}
if ("XB152" %in% vars) {x$XB152r <- suppressWarnings(as.numeric(x$XB152))}
if ("XB153" %in% vars) {x$XB153r <- suppressWarnings(as.numeric(x$XB153))}
if ("XB154" %in% vars) {x$XB154r <- suppressWarnings(as.numeric(x$XB154))}
if ("XB155" %in% vars) {x$XB155r <- suppressWarnings(as.numeric(x$XB155))}
if ("XB156" %in% vars) {x$XB156r <- suppressWarnings(as.numeric(x$XB156))}
if ("XB157" %in% vars) {x$XB157r <- suppressWarnings(as.numeric(x$XB157))}
if ("XB158" %in% vars) {x$XB158r <- suppressWarnings(as.numeric(x$XB158))}
if ("XB159" %in% vars) {x$XB159r <- suppressWarnings(as.numeric(x$XB159))}
if ("XB160" %in% vars) {x$XB160r <- suppressWarnings(as.numeric(x$XB160))}
if ("XB161" %in% vars) {x$XB161r <- suppressWarnings(as.numeric(x$XB161))}
if ("XB162" %in% vars) {x$XB162r <- suppressWarnings(as.numeric(x$XB162))}
if ("XB163" %in% vars) {x$XB163r <- suppressWarnings(as.numeric(x$XB163))}
if ("XB164" %in% vars) {x$XB164r <- suppressWarnings(as.numeric(x$XB164))}
if ("XB165" %in% vars) {x$XB165r <- suppressWarnings(as.numeric(x$XB165))}
if ("XB166" %in% vars) {x$XB166r <- suppressWarnings(as.numeric(x$XB166))}
if ("XB167" %in% vars) {x$XB167r <- suppressWarnings(as.numeric(x$XB167))}
if ("XB168" %in% vars) {x$XB168r <- suppressWarnings(as.numeric(x$XB168))}
if ("XB169" %in% vars) {x$XB169r <- suppressWarnings(as.numeric(x$XB169))}
if ("XB170" %in% vars) {x$XB170r <- suppressWarnings(as.numeric(x$XB170))}
if ("XB171" %in% vars) {x$XB171r <- suppressWarnings(as.numeric(x$XB171))}
if ("XB172" %in% vars) {x$XB172r <- suppressWarnings(as.numeric(x$XB172))}
if ("XB173" %in% vars) {x$XB173r <- suppressWarnings(as.numeric(x$XB173))}
if ("XB174" %in% vars) {x$XB174r <- suppressWarnings(as.numeric(x$XB174))}
if ("XB175" %in% vars) {x$XB175r <- suppressWarnings(as.numeric(x$XB175))}
if ("XB176" %in% vars) {x$XB176r <- suppressWarnings(as.numeric(x$XB176))}

return(x)
}
