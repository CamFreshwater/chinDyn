#Merging data into one data frame and attempting to redo analysis from 

### Installing and loading libraries----
#install.packages("rmarkdown")
#install.packages("plyr")
#install.packages("visreg")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("MCMCpack")
library(rmarkdown)
library(plyr)
library(visreg)
library(ggplot2)
library(dplyr)
library(MCMCpack)


### Loading data----

setwd("U:orca analysis")

#If want to see just since 2007
#NRKW_FecuTrimmed <- read.csv("NRKW_fecunidty_trimmed.csv",
 #                            stringsAsFactors = FALSE,
  #                           strip.white = TRUE,
   #                          na.strings =c("NA"))

#SRKW_FecuTrimmed <- read.csv("SRKW_fecunidty_trimmed.csv",
 #                            stringsAsFactors = FALSE,
  #                           strip.white = TRUE,
   #                          na.strings =c("NA"))


#For since-1980 data
NRKW_Live_Females <- read.csv("NRKW_Live_Females_since_1980.csv",
                              stringsAsFactors = FALSE,
                              strip.white = TRUE,
                              na.strings =c("NA"))

NRKW_Females_YOB <- read.csv("NRKW_Females_YOB.csv",
                             stringsAsFactors = FALSE,
                             strip.white = TRUE,
                             na.strings =c("NA"))

NRKW_Calved <- read.csv("NRKW_Calved_since_1980.csv",
                        stringsAsFactors = FALSE,
                        strip.white = TRUE,
                        na.strings =c("NA"))

SRKW_Live_Females <- read.csv("SRKW_Live_Females_since_1980.csv",
                              stringsAsFactors = FALSE,
                              strip.white = TRUE,
                              na.strings =c("NA"))


SRKW_Females_YOB <- read.csv("SRKW_Females_YOB.csv",
                             stringsAsFactors = FALSE,
                             strip.white = TRUE,
                             na.strings =c("NA"))

SRKW_Calved <- read.csv("SRKW_Calved_since_1980.csv",
                        stringsAsFactors = FALSE,
                        strip.white = TRUE,
                        na.strings =c("NA"))

Salmon_since_1980 <- read.csv("Salmon_since_1980.csv",
                              stringsAsFactors = FALSE,
                              strip.white = TRUE,
                              na.strings =c("NA"))

### Checking data----
head(NRKW_Females_YOB)
head(NRKW_Calved)
head(NRKW_Live_Females)

head(SRKW_Females_YOB)
head(SRKW_Calved)
head(SRKW_Live_Females)

head(Salmon_since_1980)


### Create "NRKW_fecundity" dataframe by adding birth, calf data, age, and salmon to live females and trimming out females on years adjacent births and those removed from Ward et al. ----

##Add year of birth to live females
#dfsNRKW1<- list(NRKW_Live_Females, NRKW_Females_YOB)
#join_all(dfsNRKW1, type="left", match = "all")

NRKW_fecundity3 <- merge(NRKW_Live_Females, NRKW_Females_YOB, by=c("ID"))
head(NRKW_fecundity3)

##Add calf data to fecundity3
#change column names of calf data from "mother" to "ID" to allow merging.  We're interested in whether the mother Calved, rather than who mothered each whale
head(NRKW_Calved)
colnames(NRKW_Calved) <- c("ID", "Year", "Calved")

dfsNRKW1980 <- list(NRKW_fecundity3, NRKW_Calved)
NRKW_fecundity2 <- join_all(dfsNRKW1980, type="left", match = "all") 


head(NRKW_fecundity2)

NRKW_fecundity2$Calved[is.na(NRKW_fecundity2$Calved)]<-0 #Changing NAs to 0s to allow for binomial analysis

#adding age

head(NRKW_fecundity2)
NRKW_fecundity2$Age <- (NRKW_fecundity2$Year-NRKW_fecundity2$YoB)
NRKW_fecundity2 <- NRKW_fecundity2[!(NRKW_fecundity2$Age<8),] #removing those less than 8 years old
head(NRKW_fecundity2)


### Remove NRKW whales that gave birth in an adjacent year----
#ie. if calved in 2009, not counted in 2008 or 2010

#copy dataframe into new one for trimming
NRKW_FecuTrimmed<-NRKW_fecundity2


###A14, B07, C06, G03, I11, I20, R05, R07, I17*  had calves in 1979    *I17 was "1979.5" but included here
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1980" & NRKW_FecuTrimmed$ID %in% c("A14", "B07", "C06", "G03", "I11", "I20", "R05", "R07", "I17")),]

### A08, G17, I01, I15, I18, I22, I31 had calves in 1980
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1981" & NRKW_FecuTrimmed$ID %in% c("A08", "G17", "I01", "I15", "I18", "I22", "I31")),]

### A23, A24, A30, G02, G18, H03 had calves in 1981
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1980" & NRKW_FecuTrimmed$ID %in% c("A23", "A24", "A30", "G02", "G18", "H03")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1982" & NRKW_FecuTrimmed$ID %in% c("A23", "A24", "A30", "G02", "G18", "H03")),]

### A36, D08, G03, R18, I17* had calves in 1982    *I17 was "1982.5" but included here
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1981" & NRKW_FecuTrimmed$ID %in% c("A36", "D08", "G03", "R18", "I17")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1983" & NRKW_FecuTrimmed$ID %in% c("A36", "D08", "G03", "R18", "I17")),]

### A10, A11, A24, I07, I11, I16, I19 had calves in 1983
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1982" & NRKW_FecuTrimmed$ID %in% c("A10", "A11", "A24", "I07", "I11", "I16", "I19")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1984" & NRKW_FecuTrimmed$ID %in% c("A10", "A11", "A24", "I07", "I11", "I16", "I19")),]

### A30, B07, D07, G20, R04 had calves in 1984
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1983" & NRKW_FecuTrimmed$ID %in% c("A30", "B07", "D07", "G20", "R04")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1985" & NRKW_FecuTrimmed$ID %in% c("A30", "B07", "D07", "G20", "R04")),]

### A24, C06, C10, G02, G08, I12, I15, I31, I33, R17 had calves in 1985
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1984" & NRKW_FecuTrimmed$ID %in% c("A24", "C06", "C10", "G02", "G08", "I12", "I15", "I31", "I33", "R17")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1986" & NRKW_FecuTrimmed$ID %in% c("A24", "C06", "C10", "G02", "G08", "I12", "I15", "I31", "I33", "R17")),]

### A25, G12, G16, G17, H03, I16, I18, I19, I20 had calves in 1986
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1985" & NRKW_FecuTrimmed$ID %in% c("A25", "G12", "G16", "G17", "H03", "I16", "I18", "I19", "I20")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1987" & NRKW_FecuTrimmed$ID %in% c("A25", "G12", "G16", "G17", "H03", "I16", "I18", "I19", "I20")),]

### A35, B07, D08, D09, D11, G25, G27, I22, I26, R05, R18 had calves in 1987
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1986" & NRKW_FecuTrimmed$ID %in% c("A35", "B07", "D08", "D09", "D11", "G25", "G27", "I22", "I26", "R05", "R18")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1988" & NRKW_FecuTrimmed$ID %in% c("A35", "B07", "D08", "D09", "D11", "G25", "G27", "I22", "I26", "R05", "R18")),]

### A24, G29, H05, I07, I33, I35, R04 had calves in 1988
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1987" & NRKW_FecuTrimmed$ID %in% c("A24", "G29", "H05", "I07", "I33", "I35", "R04")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1989" & NRKW_FecuTrimmed$ID %in% c("A24", "G29", "H05", "I07", "I33", "I35", "R04")),]

### A30, C08, C10, G02, G03, I17, I18 had calves in 1989
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1988" & NRKW_FecuTrimmed$ID %in% c("A30", "C08", "C10", "G02", "G03", "I17", "I18")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1990" & NRKW_FecuTrimmed$ID %in% c("A30", "C08", "C10", "G02", "G03", "I17", "I18")),]

### A11, A34, D11, G08, G20, G34, I11, I15, I27, R17 had calves in 1990
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1989" & NRKW_FecuTrimmed$ID %in% c("A11", "A34", "D11", "G08", "G20", "G34", "I11", "I15", "I27", "R17")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1991" & NRKW_FecuTrimmed$ID %in% c("A11", "A34", "D11", "G08", "G20", "G34", "I11", "I15", "I27", "R17")),]

### A42, B07, C06, C08, G23, G29, I20, I31, I35, I40 had calves in 1991
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1990" & NRKW_FecuTrimmed$ID %in% c("A42", "B07", "C06", "C08", "G23", "G29", "I20", "I31", "I35", "I40")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1992" & NRKW_FecuTrimmed$ID %in% c("A42", "B07", "C06", "C08", "G23", "G29", "I20", "I31", "I35", "I40")),]

### A23, A24, A35, G22, G25, I18, I26, R04, G16* had calves in 1992    *G16 was "1992.5" but included here
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1991" & NRKW_FecuTrimmed$ID %in% c("A23", "A24", "A35", "G22", "G25", "I18", "I26", "R04", "G16")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1993" & NRKW_FecuTrimmed$ID %in% c("A23", "A24", "A35", "G22", "G25", "I18", "I26", "R04", "G16")),]

### C10, I16, I19, I54 had calves in 1993
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1992" & NRKW_FecuTrimmed$ID %in% c("C10", "I16", "I19", "I54")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1994" & NRKW_FecuTrimmed$ID %in% c("C10", "I16", "I19", "I54")),]

### A25, A34, C08, G02, G20, G31, G34, H03, I35, I40, R05, R22 had calves in 1994
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1993" & NRKW_FecuTrimmed$ID %in% c("A25", "A34", "C08", "G02", "G20", "G31", "G34", "H03", "I35", "I40", "R05", "R22")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1995" & NRKW_FecuTrimmed$ID %in% c("A25", "A34", "C08", "G02", "G20", "G31", "G34", "H03", "I35", "I40", "R05", "R22")),]

### A24, A43, B07, D08, G23, I13, R04 had calves in 1995
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1994" & NRKW_FecuTrimmed$ID %in% c("A24", "A43", "B07", "D08", "G23", "I13", "R04")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1996" & NRKW_FecuTrimmed$ID %in% c("A24", "A43", "B07", "D08", "G23", "I13", "R04")),]

### A42, A48, G27, I17, I21, I26, I50, R17 had calves in 1996
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1995" & NRKW_FecuTrimmed$ID %in% c("A42", "A48", "G27", "I17", "I21", "I26", "I50", "R17")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1997" & NRKW_FecuTrimmed$ID %in% c("A42", "A48", "G27", "I17", "I21", "I26", "I50", "R17")),]

### A34, A43, A45, C08, I04, I12, I20, I27, I33, I36, I54, R20, R24 had calves in 1997
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1996" & NRKW_FecuTrimmed$ID %in% c("A34", "A43", "A45", "C08", "I04", "I12", "I20", "I27", "I33", "I36", "I54", "R20", "R24")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1998" & NRKW_FecuTrimmed$ID %in% c("A34", "A43", "A45", "C08", "I04", "I12", "I20", "I27", "I33", "I36", "I54", "R20", "R24")),]

### C10, D13, G20, G25, I16, I35, R04, R22 had calves in 1998
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1997" & NRKW_FecuTrimmed$ID %in% c("C10", "D13", "G20", "G25", "I16", "I35", "R04", "R22")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1999" & NRKW_FecuTrimmed$ID %in% c("C10", "D13", "G20", "G25", "I16", "I35", "R04", "R22")),]

### A24, A35, A50, D12, G08, G31, G37, I19, I40, R18, R23* had calves in 1999    *R23 was "1999.5" but included here
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1998" & NRKW_FecuTrimmed$ID %in% c("A24, A35, A50, D12, G08, G31, G37, I19, I40, R18, R23")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2000" & NRKW_FecuTrimmed$ID %in% c("A24, A35, A50, D12, G08, G31, G37, I19, I40, R18, R23")),]

### A45, C13, G27, H05, I17 had calves in 2000
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="1999" & NRKW_FecuTrimmed$ID %in% c("A45", "C13", "G27", "H05", "I17")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2001" & NRKW_FecuTrimmed$ID %in% c("A45", "C13", "G27", "H05", "I17")),]

### A34, G20, G22, G23, I21, I33, R22, R24 had calves in 2001
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2000" & NRKW_FecuTrimmed$ID %in% c("A34", "G20", "G22", "G23", "I21", "I33", "R22", "R24")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2002" & NRKW_FecuTrimmed$ID %in% c("A34", "G20", "G22", "G23", "I21", "I33", "R22", "R24")),]

### A52, A54, H09, I16, I24, I35, I54, R04, R17, R23 had calves in 2002
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2001" & NRKW_FecuTrimmed$ID %in% c("A52", "A54", "H09", "I16", "I24", "I35", "I54", "R04", "R17", "R23")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2003" & NRKW_FecuTrimmed$ID %in% c("A52", "A54", "H09", "I16", "I24", "I35", "I54", "R04", "R17", "R23")),]

### A24, A35, G34, I04, I17, I20, I50, I65, C08* had calves in 2003    *C08 was "2003.5" but included here
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2002" & NRKW_FecuTrimmed$ID %in% c("A24", "A35", "G34", "I04", "I17", "I20", "I50", "I65", "C08")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2004" & NRKW_FecuTrimmed$ID %in% c("A24", "A35", "G34", "I04", "I17", "I20", "I50", "I65", "C08")),]

### A34, A42, A52, A59, B14, C10, D11, D17, G20, G25, H05, I12, I27, I51, R18, R22, R24, G46* had calves in 2004    *G46 was "2004.5" but included here
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2003" & NRKW_FecuTrimmed$ID %in% c("A34", "A42", "A52", "A59", "B14", "C10", "D11", "D17", "G20", "G25", "H05", "I12", "I27", "I51", "R18", "R22", "R24", "G46")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2005" & NRKW_FecuTrimmed$ID %in% c("A34", "A42", "A52", "A59", "B14", "C10", "D11", "D17", "G20", "G25", "H05", "I12", "I27", "I51", "R18", "R22", "R24", "G46")),]

### A50, A62, D13, G02, G22, G23, G27, G41, G49, H08, I13, I19, I40, R13 had calves in 2005
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2004" & NRKW_FecuTrimmed$ID %in% c("A50", "A62", "D13", "G02", "G22", "G23", "G27", "G41", "G49", "H08", "I13", "I19", "I40", "R13")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2006" & NRKW_FecuTrimmed$ID %in% c("A50", "A62", "D13", "G02", "G22", "G23", "G27", "G41", "G49", "H08", "I13", "I19", "I40", "R13")),]

### A51, A54, C21, G37, G51, I16, I21, I26, I33, I54, I68, I69, I71, I82, R23, R29 had calves in 2006
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2005" & NRKW_FecuTrimmed$ID %in% c("A51", "A54", "C21", "G37", "G51", "I16", "I21", "I26", "I33", "I54", "I68", "I69", "I71", "I82", "R23", "R29")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2007" & NRKW_FecuTrimmed$ID %in% c("A51", "A54", "C21", "G37", "G51", "I16", "I21", "I26", "I33", "I54", "I68", "I69", "I71", "I82", "R23", "R29")),]

### 13: A56, C13, D10, D11, G31, G48, G53, I22, I24, I35, I57, I65, R22 had calves in 2007
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2006" & NRKW_FecuTrimmed$ID %in% c("A56", "C13", "D10", "D11", "G31", "G48", "G53", "I22", "I24", "I35", "I57", "I65", "R22")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2008" & NRKW_FecuTrimmed$ID %in% c("A56", "C13", "D10", "D11", "G31", "G48", "G53", "I22", "I24", "I35", "I57", "I65", "R22")),]

### 17*: A35, A42, A64, B14, G22, G25, G29, G34, G41, G52, G54, H09, I17, I83, I92, R24 in 2008 *B14 calved twice according to the data... so only 16 producing females in this list
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2007" & NRKW_FecuTrimmed$ID %in% c("A35","A42", "A64", "B14", "G22", "G25", "G29", "G34", "G41", "G52", "G54", "H09", "I17", "I83", "I92", "R24")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2009" & NRKW_FecuTrimmed$ID %in% c("A35","A42", "A64", "B14", "G22", "G25", "G29", "G34", "G41", "G52", "G54", "H09", "I17", "I83", "I92", "R24")),]

### 13: A24	A54	A62	A67	A69	C10	C19	C23	G20	I16	I51	R23	R38 in 2009
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2008" & NRKW_FecuTrimmed$ID %in% c("A24", "A54", "A62", "A67", "A69", "C10", "C19", "C23", "G20", "I16", "I51", "R23", "R38")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2010" & NRKW_FecuTrimmed$ID %in% c("A24", "A54", "A62", "A67", "A69", "C10", "C19", "C23", "G20", "I16", "I51", "R23", "R38")),]

### 13: A34	A51	A56	D17	G52	G63	I19	I21	I40	I50	I71	R29	R35 in 2010
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2009" & NRKW_FecuTrimmed$ID %in% c("A34", "A51", "A56", "D17", "G52", "G63", "I19", "I21", "I40", "I50", "I71", "R29", "R35")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2011" & NRKW_FecuTrimmed$ID %in% c("A34", "A51", "A56", "D17", "G52", "G63", "I19", "I21", "I40", "I50", "I71", "R29", "R35")),]

### 10: A50	A70	G31	G49	G50	G58	I13	R13	R22	R24 in 2011
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2010" & NRKW_FecuTrimmed$ID %in% c("A50", "A70", "G31", "G49", "G50", "G58", "I13", "R13", "R22", "R24")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2012" & NRKW_FecuTrimmed$ID %in% c("A50", "A70", "G31", "G49", "G50", "G58", "I13", "R13", "R22", "R24")),]

### 14: A67	A75	D13	G22	G27	G37	G53	G54	G64	H09	I12	I35	I68	I90 in 2012
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2011" & NRKW_FecuTrimmed$ID %in% c("A67", "A75", "D13", "G22", "G27", "G37", "G53", "G54", "G64", "H09", "I12", "I35", "I68", "I90")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2013" & NRKW_FecuTrimmed$ID %in% c("A67", "A75", "D13", "G22", "G27", "G37", "G53", "G54", "G64", "H09", "I12", "I35", "I68", "I90")),]

### 15: A42	A52	A54	A62	A73	C13	D11	G48	G50	H12	I54	I69	I83	R38	R39 in 2013
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2012" & NRKW_FecuTrimmed$ID %in% c("A42", "A52", "A54", "A62", "A73", "C13", "D11", "G48", "G50", "H12", "I54", "I69", "I83", "R38", "R39")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2014" & NRKW_FecuTrimmed$ID %in% c("A42", "A52", "A54", "A62", "A73", "C13", "D11", "G48", "G50", "H12", "I54", "I69", "I83", "R38", "R39")),]

### 18: A69	A72	B14	C08	C26	D20	G51	G69	I104	I21	I26	I50	I51	I63	I65	I99	R40	R41 in 2014
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2013" & NRKW_FecuTrimmed$ID %in% c("A69", "A72", "B14", "C08", "C26", "D20", "G51", "G69", "I104", "I21", "I26", "I50", "I51", "I63", "I65", "I99", "R40", "R41")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2015" & NRKW_FecuTrimmed$ID %in% c("A69", "A72", "B14", "C08", "C26", "D20", "G51", "G69", "I104", "I21", "I26", "I50", "I51", "I63", "I65", "I99", "R40", "R41")),]

### 12: A64	C23	D17	D19	G50	G63	I100	I96	R22	R24	R29	R35 in 2015
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2014" & NRKW_FecuTrimmed$ID %in% c("A64", "C23", "D17", "D19", "G50", "G63", "I100", "I96", "R22", "R24", "R29", "R35")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2016" & NRKW_FecuTrimmed$ID %in% c("A64", "C23", "D17", "D19", "G50", "G63", "I100", "I96", "R22", "R24", "R29", "R35")),]

###  7: A67	A70	A75	G41	I16	I27	R23 in 2016
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2015" & NRKW_FecuTrimmed$ID %in% c("A67", "A70", "A75", "G41", "I16", "I27", "R23")),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2017" & NRKW_FecuTrimmed$ID %in% c("A67", "A70", "A75", "G41", "I16", "I27", "R23")),]

### 10: A35	A42	A50	A73	G22	G58 G77	G62	I40	I92 in 2017
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$Year=="2016" & NRKW_FecuTrimmed$ID %in% c("A35", "A42", "A50", "A73", "G22", "G58 G77", "G62", "I40", "I92")),]




### Removal of NRKW females with assigned ages, to follow Ward et al. 2009 ----
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$ID=="A35"),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$ID=="B17"),] #did not give birth in our data
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$ID=="I17"),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$ID=="I18"),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$ID=="I31"),]
NRKW_FecuTrimmed<-NRKW_FecuTrimmed[!(NRKW_FecuTrimmed$ID=="R5"),]


#just in case, why not have a hard copy?
#write.csv(NRKW_FecuTrimmed, file="NRKW_fecunidty_trimmed_since_1980.csv")





### Create "SRKW_fecundity" dataframe by adding birth, calf data, age, and salmon to live females and trimming out females on years adjacent births and those removed from Ward et al. ----

##Add year of birth to live females
head(SRKW_Females_YOB)
colnames(SRKW_Females_YOB) <- c("ID", "YoB")
SRKW_fecundity3 <- merge(SRKW_Live_Females, SRKW_Females_YOB, by=c("ID"))
head(SRKW_fecundity3)

##Add calf data to fecundity3
#change column names of calf data from "MOTHER" to "ID" to allow merging.  We're interested in whether the mother Calved, rather than who mothered each whale
head(SRKW_Calved)
colnames(SRKW_Calved) <- c("Year", "ID", "Calved")


dfsSRKW <- list(SRKW_fecundity3, SRKW_Calved)
SRKW_fecundity2 <- join_all(dfsSRKW, type="left", match = "all") 


head(SRKW_fecundity2)

SRKW_fecundity2$Calved[is.na(SRKW_fecundity2$Calved)]<-0 #Changing NAs to 0s to allow for binomial analysis

#adding age

SRKW_fecundity2$Age <- (SRKW_fecundity2$Year-SRKW_fecundity2$YoB)
SRKW_fecundity2 <- SRKW_fecundity2[!(SRKW_fecundity2$Age<8),] #removing those less than 8 years old
head(SRKW_fecundity2)


### Remove SRKW whales that gave birth in an adjacent year----
#ie. if Calved in 2009, not counted in 2008 or 2010

#copy dataframe into new one for trimming
SRKW_FecuTrimmed<-SRKW_fecundity2


#J04, L03, L23 had calves in 1979
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1980" & SRKW_FecuTrimmed$ID %in% c("J04", "L03", "L23")),]

#L05, L26, L27 had calves in 1980
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1981" & SRKW_FecuTrimmed$ID %in% c("L05", "L26", "L27")),]

#L10 had a calf in 1981
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1980" & SRKW_FecuTrimmed$ID %in% c("L10")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1982" & SRKW_FecuTrimmed$ID %in% c("L10")),]

#J04 had a calf in 1982
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1981" & SRKW_FecuTrimmed$ID %in% c("J04")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1983" & SRKW_FecuTrimmed$ID %in% c("J04")),]

#no calves in 1983

#L28, L32, L35 had calves in 1984
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1983" & SRKW_FecuTrimmed$ID %in% c("L28", "L32", "L35")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1985" & SRKW_FecuTrimmed$ID %in% c("L28", "L32", "L35")),]

#J10, K03, L02, L11, L27 had calves in 1985
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1984" & SRKW_FecuTrimmed$ID %in% c("J10", "K03", "L02", "L11", "L27")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1986" & SRKW_FecuTrimmed$ID %in% c("J10", "K03", "L02", "L11", "L27")),]

#K13, K18, L03, L05, L22, L26, L43 had calves in 1986
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1985" & SRKW_FecuTrimmed$ID %in% c("K13", "K18", "L03", "L05", "L22", "L26", "L43")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1987" & SRKW_FecuTrimmed$ID %in% c("K13", "K18", "L03", "L05", "L22", "L26", "L43")),]

#J14, K12, L07, L11 had calves in 1987
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1986" & SRKW_FecuTrimmed$ID %in% c("J14", "K12", "L07", "L11")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1988" & SRKW_FecuTrimmed$ID %in% c("J14", "K12", "L07", "L11")),]

#J11, K14 had calves in 1988
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1987" & SRKW_FecuTrimmed$ID %in% c("J11", "K14")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1989" & SRKW_FecuTrimmed$ID %in% c("J11", "K14")),]

#L02, L22 had calves in 1989
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1988" & SRKW_FecuTrimmed$ID %in% c("L02", "L22")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1990" & SRKW_FecuTrimmed$ID %in% c("L02", "L22")),]

#K14, L27, L47, L51, L55, L60 had calves in 1990
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1989" & SRKW_FecuTrimmed$ID %in% c("K14", "L27", "L47", "L51", "L55", "L60")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1991" & SRKW_FecuTrimmed$ID %in% c("K14", "L27", "L47", "L51", "L55", "L60")),]

#J11, J16, K13, L04, L28 had calves in 1991
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1990" & SRKW_FecuTrimmed$ID %in% c("J11", "J16", "K13", "L04", "L28")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1992" & SRKW_FecuTrimmed$ID %in% c("J11", "J16", "K13", "L04", "L28")),]

#L32 had calves in 1992
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1991" & SRKW_FecuTrimmed$ID %in% c("L32")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1993" & SRKW_FecuTrimmed$ID %in% c("L32")),]

#J0A, J17, J19, K14, L02, L22, L26 had calves in 1993
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1992" & SRKW_FecuTrimmed$ID %in% c("J0A", "J17", "J19", "K14", "L02", "L22", "L26")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1994" & SRKW_FecuTrimmed$ID %in% c("J0A", "J17", "J19", "K14", "L02", "L22", "L26")),]

#K12, K13 had calves in 1994
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1993" & SRKW_FecuTrimmed$ID %in% c("K12", "K13")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1995" & SRKW_FecuTrimmed$ID %in% c("K12", "K13")),]

#J11, J14, L11, L27, L47, L60 had calves in 1995
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1994" & SRKW_FecuTrimmed$ID %in% c("J11", "J14", "L11", "L27", "L47", "L60")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1996" & SRKW_FecuTrimmed$ID %in% c("J11", "J14", "L11", "L27", "L47", "L60")),]

#J16, J20, K03, L43, L55 had calves in 1996
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1995" & SRKW_FecuTrimmed$ID %in% c("J16", "J20", "K03", "L43", "L55")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1997" & SRKW_FecuTrimmed$ID %in% c("J16", "J20", "K03", "L43", "L55")),]

# No calves in 1997

#J17, J22 had calves in 1998
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1997" & SRKW_FecuTrimmed$ID %in% c("J17", "J22")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1999" & SRKW_FecuTrimmed$ID %in% c("J17", "J22")),]

#J16, K12, L51, L67 had calves in 1999
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1998" & SRKW_FecuTrimmed$ID %in% c("J16", "K12", "L51", "L67")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2000" & SRKW_FecuTrimmed$ID %in% c("J16", "K12", "L51", "L67")),]

#K16, L47 had calves in 2000
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="1999" & SRKW_FecuTrimmed$ID %in% c("K16", "L47")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2001" & SRKW_FecuTrimmed$ID %in% c("K16", "L47")),]

#J14, K13, K22, L54 had calves in 2001
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2000" & SRKW_FecuTrimmed$ID %in% c("J14", "K13", "K22", "L54")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2002" & SRKW_FecuTrimmed$ID %in% c("J14", "K13", "K22", "L54")),]

#K16, L47, L67 had calves in 2002
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2001" & SRKW_FecuTrimmed$ID %in% c("K16", "L47", "L67")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2003" & SRKW_FecuTrimmed$ID %in% c("K16", "L47", "L67")),]

#J11, J22, K12, K14, L55 had calves in 2003
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2002" & SRKW_FecuTrimmed$ID %in% c("J11", "J22", "K12", "K14", "L55")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2004" & SRKW_FecuTrimmed$ID %in% c("J11", "J22", "K12", "K14", "L55")),]

#J14, K20, L43, L72 had calves in 2004
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2003" & SRKW_FecuTrimmed$ID %in% c("J14", "K20", "L43", "L72")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2005" & SRKW_FecuTrimmed$ID %in% c("J14", "K20", "L43", "L72")),]

#J19, L47, L86 had calves in 2005
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2004" & SRKW_FecuTrimmed$ID %in% c("J19", "L47", "L86")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2006" & SRKW_FecuTrimmed$ID %in% c("J19", "L47", "L86")),]

#K22, K28, L54 had calves in 2006
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2005" & SRKW_FecuTrimmed$ID %in% c("K22", "K28", "L54")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2007" & SRKW_FecuTrimmed$ID %in% c("K22", "K28", "L54")),]

#J14, J16, L55, L83 had calves in 2007
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2006" & SRKW_FecuTrimmed$ID %in% c("J14", "J16", "L55", "L83")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2008" & SRKW_FecuTrimmed$ID %in% c("J14", "J16", "L55", "L83")),]

#K14, L47 had calves in 2008
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2007" & SRKW_FecuTrimmed$ID %in% c("K14", "L47")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2009" & SRKW_FecuTrimmed$ID %in% c("K14", "L47")),]

#J14, J17, J28, L86, L94 had calves in 2009
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2008" & SRKW_FecuTrimmed$ID %in% c("J14", "J17", "J28", "L86", "L94")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2010" & SRKW_FecuTrimmed$ID %in% c("J14", "J17", "J28", "L86", "L94")),]

#J35, K12, L47, L54, L77, L82 had calves in 2010
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2009" & SRKW_FecuTrimmed$ID %in% c("J35", "K12", "L47", "L54", "L77", "L82")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2011" & SRKW_FecuTrimmed$ID %in% c("J35", "K12", "L47", "L54", "L77", "L82")),]

#J16, K27, L55 had calves in 2011
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2010" & SRKW_FecuTrimmed$ID %in% c("J16", "K27", "L55")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2012" & SRKW_FecuTrimmed$ID %in% c("J16", "K27", "L55")),]

#J37, L77 had calves in 2012
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2011" & SRKW_FecuTrimmed$ID %in% c("J37", "L77")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2013" & SRKW_FecuTrimmed$ID %in% c("J37", "L77")),]

# No calves in 2013

#J16, L86 had calves in 2014
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2013" & SRKW_FecuTrimmed$ID %in% c("J16", "L86")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2015" & SRKW_FecuTrimmed$ID %in% c("J16", "L86")),]

#J17, J28, J36, J41, L103, L91, L94 had calves in 2015
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2014" & SRKW_FecuTrimmed$ID %in% c("J17", "J28", "J36", "J41", "L103", "L91", "L94")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2016" & SRKW_FecuTrimmed$ID %in% c("J17", "J28", "J36", "J41", "L103", "L91", "L94")),]

#J40 had calves in 2016
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2015" & SRKW_FecuTrimmed$ID %in% c("J40")),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2017" & SRKW_FecuTrimmed$ID %in% c("J40")),]

# No calves in 2017
#SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$Year=="2016" & SRKW_FecuTrimmed$ID %in% c("")),]

### Removal of females with assigned ages, to follow Ward et al. 2009 ----
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$ID=="J05"),] #No births after 1980
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$ID=="J12"),] #No births after 1980
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$ID=="K03"),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$ID=="L02"),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$ID=="L03"),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$ID=="L11"),]
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$ID=="L21"),] #No births after 1980
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$ID=="L23"),] #No births after 1980
SRKW_FecuTrimmed<-SRKW_FecuTrimmed[!(SRKW_FecuTrimmed$ID=="L45"),] #No births after 1980


#just in case, why not have a hard copy?
#write.csv(SRKW_FecuTrimmed, file="SRKW_fecunidty_trimmed_since_1980.csv")



### Joining NRKW and SRKW ----

## Add region so we can tell them apart
NRKW_FecuTrimmed$Region <- "North"
SRKW_FecuTrimmed$Region <- "South"

## Add in birthrates and change from average birthrate (1980-2017) for each year
NRBirthrate <- merge(aggregate(NRKW_FecuTrimmed$Calved, by=list(Category=NRKW_FecuTrimmed$Year), FUN=sum), aggregate(NRKW_FecuTrimmed$ID, by=list(Category=NRKW_FecuTrimmed$Year), FUN=length), by="Category")
head(NRBirthrate)
colnames(NRBirthrate) <- c("Year", "SumCalved","N")
head(NRBirthrate)
NRBirthrate$Birthrate <- NRBirthrate$SumCalved/NRBirthrate$N
NRBirthrate$BRDiff <- NRBirthrate$Birthrate/mean(NRBirthrate$Birthrate)
NRKW_FecuTrimmed <- merge(NRKW_FecuTrimmed, NRBirthrate, by="Year")

SRBirthrate <- merge(aggregate(SRKW_FecuTrimmed$Calved, by=list(Category=SRKW_FecuTrimmed$Year), FUN=sum), aggregate(SRKW_FecuTrimmed$ID, by=list(Category=SRKW_FecuTrimmed$Year), FUN=length), by="Category")
head(SRBirthrate)
colnames(SRBirthrate) <- c("Year", "SumCalved","N")
head(SRBirthrate)
SRBirthrate$Birthrate <- SRBirthrate$SumCalved/SRBirthrate$N
SRBirthrate$BRDiff <- SRBirthrate$Birthrate/mean(SRBirthrate$Birthrate)
SRKW_FecuTrimmed <- merge(SRKW_FecuTrimmed, SRBirthrate, by="Year")
  
head(NRKW_FecuTrimmed)
head(SRKW_FecuTrimmed)

#Check if "BIRTH_YEAR" is column 4 or 5, if still present (if loaded from merged fecutrimmed files, it'll be 5)
#colnames(SRKW_FecuTrimmed)[4] <- "YoB"


## and bind!
FecuTrimmed_allyears2 <- rbind(SRKW_FecuTrimmed, NRKW_FecuTrimmed)



##adding salmon data to fecundity2  
dfsALL2 <- list(FecuTrimmed_allyears2, Salmon_since_1980)
FecuTrimmed_allyears <- join_all(dfsALL2, type="full", match = "all") 

#FecuTrimmed_allyears <- merge(FecuTrimmed_allyears2, Salmon_since_1980, by="Year", no.dups = FALSE)
head(FecuTrimmed_allyears)


#just in case, why not have a hard copy?
#write.csv(FecuTrimmed_allyears, file="Combined_fecunidty_trimmed_allyears.csv")


### Binomial GLM of fecunidty data with region ----

# First looking at data


plot(FecuTrimmed_allyears$Calved~FecuTrimmed_allyears$Age)
is.numeric(FecuTrimmed_allyears$Age)

plot(FecuTrimmed_allyears$BRDiff~FecuTrimmed_allyears$WCVI_Index)
plot(FecuTrimmed_allyears$BRDiff~FecuTrimmed_allyears$WCVI_Index.lag1)
#summary(lm(FecuTrimmed_allyears$BRDiff~FecuTrimmed_allyears$WCVI_Index.lag1))

plot(Salmon_since_1980$WCVI_Index~Salmon_since_1980$Year)
plot(FecuTrimmed_allyears$WCVI_Index~FecuTrimmed_allyears$Year)


## WCVI index


plot(FecuTrimmed_allyears$Calved~FecuTrimmed_allyears$WCVI_Index)

ALL_WCVI_Index_fit <- glm(Calved~poly(Age,4)+WCVI_Index+Region, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_WCVI_Index_fit)
hist(resid(ALL_WCVI_Index_fit))

#FTexperiment <- FecuTrimmed_allyears
#FTexperiment$WIFP<- fitted(ALL_WCVI_Index_fit)
#FTexperiment$WIFPdiff <- FTexperiment$WIFP/mean(FTexperiment$WIFP)
#boxplot(FTexperiment$WIFPdiff~FTexperiment$WCVI_Index)

FecuTrimmed_allyears$colour[FecuTrimmed_allyears$Region=="North"]<-"blue"
FecuTrimmed_allyears$colour[FecuTrimmed_allyears$Region=="South"]<-"red"
FecuTrimmed_allyears$dot[FecuTrimmed_allyears$Region=="North"]<-1
FecuTrimmed_allyears$dot[FecuTrimmed_allyears$Region=="South"]<-19
FecuTrimmed_allyears$dot<-as.numeric(FecuTrimmed_allyears$dot)


plot(fitted(ALL_WCVI_Index_fit)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="WCVI Index", xlab="Age", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_WCVI_Index_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="WCVI Index") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))


boxplot(fitted(ALL_WCVI_Index_fit)~FecuTrimmed_allyears$WCVI_Index*FecuTrimmed_allyears$Region, ylim=c(0,0.4), main="WCVI Index", xlab="WCVI Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)

plot(fitted(ALL_WCVI_Index_fit)~jitter(FecuTrimmed_allyears$WCVI_Index, 2.5), ylim=c(0,0.4), main="WCVI Index", xlab="WCVI Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


#visreg(ALL_WCVI_Index_fit, scale='response')


#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$WCVIfact<-as.factor(FecuTrimmed_allyears$WCVI_Index)

ggplot(FecuTrimmed_allyears, aes(x=WCVIfact, y=fitted(ALL_WCVI_Index_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.301) +
  scale_x_discrete(breaks=c(0.2240855, 0.623705813, 0.931366345, 1.335899735, 2.379985586), labels=c(0.22, 0.62, 0.93, 1.34, 2.38)) +
  labs(x="WCVI Index", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))

  


length(FecuTrimmed_allyears$Age)  
length(predict(ALL_WCVI_Index_fit))

 #allRgraph+geom_point(aes(size = 2)) + theme_bw()+theme(axis.text.x=element_text(size=12, angle=90, vjust=0.5, hjust=1), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), legend.position="none")+ labs(size=12, x=NULL, y="Risk Score\n") + scale_y_continuous(limits=c(1,9), breaks=c(1:9))+  annotate("text", x = 59.45, y=8.9, label = "All Regions", size=6)
#ggsave("All Regions AM algae with sterr - calculated without natives.pdf", width=13, height= 5)





### WCVI_Index lagged 1 year


ALL_WCVI_Index.lag1_fit <- glm(Calved~poly(Age,4)+WCVI_Index.lag1+Region, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_WCVI_Index.lag1_fit)

plot(fitted(ALL_WCVI_Index.lag1_fit)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="WCVI Index lag 1 year", xlab="Age", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_WCVI_Index.lag1_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="WCVI Index lag 1 year") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))




boxplot(fitted(ALL_WCVI_Index.lag1_fit)~FecuTrimmed_allyears$WCVI_Index.lag1, ylim=c(0,0.4), main="WCVI Index lag 1 year", xlab="WCVI Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)

plot(fitted(ALL_WCVI_Index.lag1_fit)~jitter(FecuTrimmed_allyears$WCVI_Index.lag1, 2.5), ylim=c(0,0.4), main="WCVI Index lag 1 year", xlab="WCVI Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 

#visreg(ALL_WCVI_Index.lag1_fit, scale='response')


#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$WCVILag1fact<-as.factor(FecuTrimmed_allyears$WCVI_Index.lag1)

ggplot(FecuTrimmed_allyears, aes(x=WCVILag1fact, y=fitted(ALL_WCVI_Index.lag1_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.301) +
  scale_x_discrete(breaks=c(0.2240855, 0.623705813, 0.931366345, 1.335899735, 2.379985586), labels=c(0.22, 0.62, 0.93, 1.34, 2.38)) +
  labs(x="WCVI Index lag 1 year", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))



#FTexperiment <- FecuTrimmed_allyears
#FTexperiment$WIFP.1<- fitted(ALL_WCVI_Index.lag1_fit)
#FTexperiment$WIFP.1diff <- FTexperiment$WIFP.1/mean(FTexperiment$WIFP.1)
#boxplot(FTexperiment$WIFPdiff~FTexperiment$WCVI_Index.lag1)





### Coastwide index

plot(FecuTrimmed_allyears$Calved~FecuTrimmed_allyears$Coastwide_Index)

ALL_Coastwide_index_fit <- glm(Calved~poly(Age,4)+Coastwide_Index+Region, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_Coastwide_index_fit)

plot(fitted(ALL_Coastwide_index_fit, scale='Calved')~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.5))
#visreg(ALL_Coastwide_index_fit, scale='response')

plot(fitted(ALL_Coastwide_index_fit)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="Coastwide Index", xlab="Age", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_Coastwide_index_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="Coastwide Index") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))


boxplot(fitted(ALL_Coastwide_index_fit)~FecuTrimmed_allyears$Coastwide_Index, ylim=c(0,0.4), main="Coastwide Index", xlab="Coastwide Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)

plot(fitted(ALL_Coastwide_index_fit)~jitter(FecuTrimmed_allyears$Coastwide_Index, 2.5), ylim=c(0,0.4), main="Coastwide Index", xlab="Coastwide Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 

#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$Coastwidefact<-as.factor(FecuTrimmed_allyears$Coastwide_Index)

ggplot(FecuTrimmed_allyears, aes(x=Coastwidefact, y=fitted(ALL_Coastwide_index_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.34) +
  scale_x_discrete(breaks=c(0.60754998, 0.781029037, 1.049534295, 1.147386131, 1.445263092), labels=c(0.61, 0.78, 1.05, 1.15, 1.15)) +
  labs(x="Coastwide Index", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))




### Coastwide index lag 1 year

plot(FecuTrimmed_allyears$Calved~FecuTrimmed_allyears$Coastwide_Index.lag1)

ALL_Coastwide_Index.lag1_fit <- glm(Calved~poly(Age,4)+Coastwide_Index.lag1+Region, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_Coastwide_Index.lag1_fit)

plot(fitted(ALL_Coastwide_Index.lag1_fit, scale='Calved')~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.5))
#visreg(ALL_Coastwide_Index.lag1_fit, scale='response')

plot(fitted(ALL_Coastwide_Index.lag1_fit)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="Coastwide Index lag 1 year", xlab="Age", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_Coastwide_Index.lag1_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="Coastwide Index lag 1 year") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))



boxplot(fitted(ALL_Coastwide_Index.lag1_fit)~FecuTrimmed_allyears$Coastwide_Index.lag1, ylim=c(0,0.4), main="Coastwide Index lag 1 year", xlab="Coastwide Index lag 1 year", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)

plot(fitted(ALL_Coastwide_Index.lag1_fit)~jitter(FecuTrimmed_allyears$Coastwide_Index.lag1, 2.5), ylim=c(0,0.4), main="Coastwide Index lag 1 year", xlab="Coastwide Index lag 1 year", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 

#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$CoastwideLag1fact<-as.factor(FecuTrimmed_allyears$Coastwide_Index.lag1)

ggplot(FecuTrimmed_allyears, aes(x=CoastwideLag1fact, y=fitted(ALL_Coastwide_Index.lag1_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.34) +
  scale_x_discrete(breaks=c(0.60754998, 0.781029037, 1.049534295, 1.147386131, 1.445263092), labels=c(0.61, 0.78, 1.05, 1.15, 1.15)) +
  labs(x="Coastwide Index lag 1 year", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))


### SRKW index

plot(FecuTrimmed_allyears$Calved~FecuTrimmed_allyears$SRKW_Index)

ALL_SRKW_Index_fit <- glm(Calved~poly(Age,4)+SRKW_Index+Region, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_SRKW_Index_fit)

plot(fitted(ALL_SRKW_Index_fit, scale='Calved')~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.5))
#visreg(ALL_SRKW_Index_fit, scale='response')

plot(fitted(ALL_SRKW_Index_fit)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="SRKW Index", xlab="Age", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)
                                        
ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_SRKW_Index_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="SRKW Index") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))
                                        

boxplot(fitted(ALL_SRKW_Index_fit)~FecuTrimmed_allyears$SRKW_Index, ylim=c(0,0.4), main="SRKW Index lag 1 year", xlab="SRKW Index lag 1 year", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)

plot(fitted(ALL_SRKW_Index_fit)~jitter(FecuTrimmed_allyears$SRKW_Index, 2.5), ylim=c(0,0.4), main="SRKW Index", xlab="SRKW Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$SRKWfact<-as.factor(FecuTrimmed_allyears$SRKW_Index)

ggplot(FecuTrimmed_allyears, aes(x=SRKWfact, y=fitted(ALL_SRKW_Index_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.33) +
  scale_x_discrete(breaks=c(0.466835142, 0.717500422, 1.01140825, 1.345621478, 1.610760791), labels=c(0.47, 0.72, 1.01, 1.34, 1.61)) +
  labs(x="SRKW Index", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))
                                 
                                 

### SRKW Index lag 1 year

ALL_SRKW_Index.lag1_fit <- glm(Calved~poly(Age,4)+SRKW_Index.lag1+Region, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_SRKW_Index.lag1_fit)

plot(fitted(ALL_SRKW_Index.lag1_fit)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="SRKW Index lag 1 year", xlab="Age", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 



#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_SRKW_Index.lag1_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="SRKW Index lag 1 year") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))




boxplot(fitted(ALL_SRKW_Index.lag1_fit)~FecuTrimmed_allyears$SRKW_Index.lag1, ylim=c(0,0.4), main="SRKW Index lag 1 year", xlab="SRKW Index lag 1 year", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)

plot(fitted(ALL_SRKW_Index.lag1_fit)~jitter(FecuTrimmed_allyears$SRKW_Index.lag1, 2.5), ylim=c(0,0.4), main="SRKW Index lag 1 year", xlab="SRKW Index lag 1 year", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$SRKWLag1fact<-as.factor(FecuTrimmed_allyears$SRKW_Index.lag1)

ggplot(FecuTrimmed_allyears, aes(x=SRKWLag1fact, y=fitted(ALL_SRKW_Index.lag1_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.33) +
  scale_x_discrete(breaks=c(0.466835142, 0.717500422, 1.01140825, 1.345621478, 1.610760791), labels=c(0.47, 0.72, 1.01, 1.34, 1.61)) +
  labs(x="SRKW Index lag 1 year", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))




### NRKW index

plot(FecuTrimmed_allyears$Calved~FecuTrimmed_allyears$NRKW_Index)

ALL_NRKW_Index_fit <- glm(Calved~poly(Age,4)+NRKW_Index+Region, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_NRKW_Index_fit)

plot(fitted(ALL_NRKW_Index_fit, scale='Calved')~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.5))
#visreg(ALL_NRKW_Index_fit, scale='response')

plot(fitted(ALL_NRKW_Index_fit)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="NRKW Index", xlab="Age", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_NRKW_Index_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="NRKW Index") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))



boxplot(fitted(ALL_NRKW_Index_fit)~FecuTrimmed_allyears$NRKW_Index, ylim=c(0,0.4), main="NRKW Index", xlab="NRKW Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)

plot(fitted(ALL_NRKW_Index_fit)~jitter(FecuTrimmed_allyears$NRKW_Index, 2.5), ylim=c(0,0.4), main="NRKW Index", xlab="NRKW Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


### Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$NRKWfact<-as.factor(FecuTrimmed_allyears$NRKW_Index)

ggplot(FecuTrimmed_allyears, aes(x=NRKWfact, y=fitted(ALL_NRKW_Index_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.355) +
  scale_x_discrete(breaks=c(0.583906655, 0.747866325, 0.940862357, 1.19058655, 1.817588067), labels=c(0.58, 0.73, 0.94, 1.19, 1.82)) +
  labs(x="NRKW Index", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))




#NRKW Index lag 1 year

ALL_NRKW_Index.lag1_fit <- glm(Calved~poly(Age,4)+NRKW_Index.lag1+Region, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_NRKW_Index.lag1_fit)

plot(fitted(ALL_NRKW_Index.lag1_fit)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="NRKW Index lag 1 year", xlab="Age", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


#Make a agefactor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_NRKW_Index.lag1_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="NRKW Index lag 1 year") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))


boxplot(fitted(ALL_NRKW_Index.lag1_fit)~FecuTrimmed_allyears$NRKW_Index.lag1, ylim=c(0,0.4), main="NRKW Index lag 1 year", xlab="NRKW Index lag 1 year", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)

plot(fitted(ALL_NRKW_Index.lag1_fit)~jitter(FecuTrimmed_allyears$NRKW_Index.lag1, 2.5), ylim=c(0,0.4), main="NRKW Index lag 1 year", xlab="NRKW Index lag 1 year", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$NRKWLag1fact<-as.factor(FecuTrimmed_allyears$NRKW_Index.lag1)

ggplot(FecuTrimmed_allyears, aes(x=NRKWLag1fact, y=fitted(ALL_NRKW_Index.lag1_fit))) + 
  geom_boxplot(aes(fill=Region)) +
  ylim(-0.001, 0.355) +
  scale_x_discrete(breaks=c(0.583906655, 0.747866325, 0.940862357, 1.19058655, 1.817588067), labels=c(0.58, 0.73, 0.94, 1.19, 1.82)) +
  labs(x="NRKW Index lag 1 year", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))



### JC calculated WCVI index with average calculated from 1980-2017 ----

plot(FecuTrimmed_allyears$Calved~FecuTrimmed_allyears$JC_WCVI_Index)

ALL_JC_WCVI_Index_fit <- glm(Calved~poly(Age,4)+JC_WCVI_Index+Region, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_JC_WCVI_Index_fit)

plot(fitted(ALL_JC_WCVI_Index_fit, scale='Calved')~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.5))
#visreg(ALL_JC_WCVI_Index_fit, scale='response')

plot(fitted(ALL_JC_WCVI_Index_fit)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="JC_WCVI_Index", xlab="Age", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 

boxplot(fitted(ALL_JC_WCVI_Index_fit)~FecuTrimmed_allyears$JC_WCVI_Index, ylim=c(0,0.4), main="JC_WCVI_Index", xlab="JC_WCVI_Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)

plot(fitted(ALL_JC_WCVI_Index_fit)~jitter(FecuTrimmed_allyears$JC_WCVI_Index, 2.5), ylim=c(0,0.4), main="JC_WCVI_Index", xlab="JC_WCVI_Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 



### JC calculated WCVI index with one year lag

plot(FecuTrimmed_allyears$Calved~FecuTrimmed_allyears$JC_WCVI_Index.lag1)

ALL_JC_WCVI_Index.lag1_fit <- glm(Calved~poly(Age,4)+JC_WCVI_Index.lag1+Region, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_JC_WCVI_Index.lag1_fit)

plot(fitted(ALL_JC_WCVI_Index.lag1_fit, scale='Calved')~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.5))
#visreg(ALL_JC_WCVI_Index.lag1_fit, scale='response')

plot(fitted(ALL_JC_WCVI_Index.lag1_fit)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="JC_WCVI_Index.lag1", xlab="Age", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 

boxplot(fitted(ALL_JC_WCVI_Index.lag1_fit)~FecuTrimmed_allyears$JC_WCVI_Index.lag1, ylim=c(0,0.4), main="JC_WCVI_Index.lag1", xlab="JC_WCVI_Index.lag1", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)

plot(fitted(ALL_JC_WCVI_Index.lag1_fit)~jitter(FecuTrimmed_allyears$JC_WCVI_Index.lag1, 2.5), ylim=c(0,0.4), main="JC_WCVI_Index.lag1", xlab="JC_WCVI_Index.lag1", ylab="Probability of calf", col=FecuTrimmed_allyears$colour, pch=FecuTrimmed_allyears$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 





### Model comparison ----

AIC(ALL_WCVI_Index_fit, ALL_WCVI_Index.lag1_fit, ALL_Coastwide_index_fit, ALL_Coastwide_Index.lag1_fit, ALL_SRKW_Index_fit, ALL_SRKW_Index.lag1_fit, ALL_NRKW_Index_fit, ALL_NRKW_Index.lag1_fit)

#AIC(ALL_Coastwide_index_fit, ALL_NRKW_Index_fit, ALL_SRKW_Index_fit, ALL_WCVI_Index_fit, ALL_WCVI_Index.lag1_fit)#, ALL_JC_WCVI_Index_fit, ALL_JC_WCVI_Index.lag1_fit)

# AICc using parameter estimation recommended by bolker et al on glmm.wikidot.com/faq
q <- 0 # number of components in each random effect. We didn't use any random effects 
K <- function(x) {length(coef(x)) + 1}
AICc.mem <- function(x) {-2*as.numeric(logLik(x)) + 2*K(x)*(length(FecuTrimmed_allyears$Calved)/(length(FecuTrimmed_allyears$Calved)-K(x)-1))}

AIC.sum<- as.data.frame(cbind(AICc.mem(ALL_Coastwide_index_fit), AICc.mem(ALL_NRKW_Index_fit), AICc.mem(ALL_SRKW_Index_fit), AICc.mem(ALL_WCVI_Index_fit), AICc.mem(ALL_WCVI_Index.lag1_fit)))#, AICc.mem(ALL_JC_WCVI_Index_fit), AICc.mem(ALL_JC_WCVI_Index.lag1_fit)))
names(AIC.sum) <- c('Coastwide Index', 'NRKW Index', 'SRKW Index', 'WCVI Index', 'ALL_WCVI_Index.lag1_fit')#, 'ALL_JC_WCVI_Index_fit','ALL_JC_WCVI_Index_fit.lag1')
AIC.sum



#### Why is Pcalf lower now that data is since 1980? ----
L2008 <- subset(FecuTrimmed_allyears[FecuTrimmed_allyears$Year<2008,])
G2007 <- subset(FecuTrimmed_allyears[FecuTrimmed_allyears$Year>2007,])


L2008_WCVI_Index_fit <- glm(Calved~poly(Age,4)+WCVI_Index+Region, family=binomial(link="logit"), data=L2008)
summary(L2008_WCVI_Index_fit)
#plot(residuals(L2008_WCVI_Index_fit)~L2008$WCVI_Index)

plot(fitted(L2008_WCVI_Index_fit)~L2008$Age, xlim=c(0,60), ylim=c(0,0.4), main="WCVI Index <2008", xlab="Age", ylab="Probability of calf", col=L2008$colour, pch=L2008$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 


G2007_WCVI_Index_fit <- glm(Calved~poly(Age,4)+WCVI_Index+Region, family=binomial(link="logit"), data=G2007)
summary(G2007_WCVI_Index_fit)
#plot(residuals(G2007_WCVI_Index_fit)~G2007$WCVI_Index)

plot(fitted(G2007_WCVI_Index_fit)~G2007$Age, xlim=c(0,60), ylim=c(0,0.4), main="WCVI Index >2007", xlab="Age", ylab="Probability of calf", col=G2007$colour, pch=G2007$dot)
legend(
  x ="topright",
  legend = c("North","South"),
  col = c("blue", 'red'),
  pch = 19, 
  cex = .7) 





#### Models without region ----

### Binomial GLM of fecunidty data----


## WCVI index


ALL_WCVI_Index_fit_noRegion <- glm(Calved~poly(Age,4)+WCVI_Index, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_WCVI_Index_fit_noRegion)
hist(resid(ALL_WCVI_Index_fit_noRegion))


plot(fitted(ALL_WCVI_Index_fit_noRegion)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="WCVI Index", xlab="Age", ylab="Probability of calf")



#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_WCVI_Index_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="WCVI Index") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))


boxplot(fitted(ALL_WCVI_Index_fit_noRegion)~FecuTrimmed_allyears$WCVI_Index*FecuTrimmed_allyears$Region, ylim=c(0,0.4), main="WCVI Index", xlab="WCVI Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)


#visreg(ALL_WCVI_Index_fit_noRegion, scale='response')


#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$WCVIfact<-as.factor(FecuTrimmed_allyears$WCVI_Index)

ggplot(FecuTrimmed_allyears, aes(x=WCVIfact, y=fitted(ALL_WCVI_Index_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.301) +
  scale_x_discrete(breaks=c(0.2240855, 0.623705813, 0.931366345, 1.335899735, 2.379985586), labels=c(0.22, 0.62, 0.93, 1.34, 2.38)) +
  labs(x="WCVI Index", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))



#ggsave("All Regions AM algae with sterr - calculated without natives.pdf", width=13, height= 5)




### WCVI_Index lagged 1 year


ALL_WCVI_Index.lag1_fit_noRegion <- glm(Calved~poly(Age,4)+WCVI_Index.lag1, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_WCVI_Index.lag1_fit_noRegion)

plot(fitted(ALL_WCVI_Index.lag1_fit_noRegion)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="WCVI Index lag 1 year", xlab="Age", ylab="Probability of calf")



#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_WCVI_Index.lag1_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="WCVI Index lag 1 year") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))




boxplot(fitted(ALL_WCVI_Index.lag1_fit_noRegion)~FecuTrimmed_allyears$WCVI_Index.lag1, ylim=c(0,0.4), main="WCVI Index lag 1 year", xlab="WCVI Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)



#visreg(ALL_WCVI_Index.lag1_fit_noRegion, scale='response')


#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$WCVILag1fact<-as.factor(FecuTrimmed_allyears$WCVI_Index.lag1)

ggplot(FecuTrimmed_allyears, aes(x=WCVILag1fact, y=fitted(ALL_WCVI_Index.lag1_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.301) +
  scale_x_discrete(breaks=c(0.2240855, 0.623705813, 0.931366345, 1.335899735, 2.379985586), labels=c(0.22, 0.62, 0.93, 1.34, 2.38)) +
  labs(x="WCVI Index lag 1 year", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))



#FTexperiment <- FecuTrimmed_allyears
#FTexperiment$WIFP.1<- fitted(ALL_WCVI_Index.lag1_fit_noRegion)
#FTexperiment$WIFP.1diff <- FTexperiment$WIFP.1/mean(FTexperiment$WIFP.1)
#boxplot(FTexperiment$WIFPdiff~FTexperiment$WCVI_Index.lag1)





### Coastwide index

ALL_Coastwide_index_fit_noRegion <- glm(Calved~poly(Age,4)+Coastwide_Index, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_Coastwide_index_fit_noRegion)

plot(fitted(ALL_Coastwide_index_fit_noRegion, scale='Calved')~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.5))
#visreg(ALL_Coastwide_index_fit_noRegion, scale='response')

plot(fitted(ALL_Coastwide_index_fit_noRegion)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="Coastwide Index", xlab="Age", ylab="Probability of calf")



#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_Coastwide_index_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="Coastwide Index") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))


boxplot(fitted(ALL_Coastwide_index_fit_noRegion)~FecuTrimmed_allyears$Coastwide_Index, ylim=c(0,0.4), main="Coastwide Index", xlab="Coastwide Index", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)

plot(fitted(ALL_Coastwide_index_fit_noRegion)~jitter(FecuTrimmed_allyears$Coastwide_Index, 2.5), ylim=c(0,0.4), main="Coastwide Index", xlab="Coastwide Index", ylab="Probability of calf")

#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$Coastwidefact<-as.factor(FecuTrimmed_allyears$Coastwide_Index)

ggplot(FecuTrimmed_allyears, aes(x=Coastwidefact, y=fitted(ALL_Coastwide_index_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.34) +
  scale_x_discrete(breaks=c(0.60754998, 0.781029037, 1.049534295, 1.147386131, 1.445263092), labels=c(0.61, 0.78, 1.05, 1.15, 1.15)) +
  labs(x="Coastwide Index", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))




### Coastwide index lag 1 year

ALL_Coastwide_Index.lag1_fit_noRegion <- glm(Calved~poly(Age,4)+Coastwide_Index.lag1, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_Coastwide_Index.lag1_fit_noRegion)

plot(fitted(ALL_Coastwide_Index.lag1_fit_noRegion, scale='Calved')~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.5))
#visreg(ALL_Coastwide_Index.lag1_fit_noRegion, scale='response')

plot(fitted(ALL_Coastwide_Index.lag1_fit_noRegion)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="Coastwide Index lag 1 year", xlab="Age", ylab="Probability of calf")



#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_Coastwide_Index.lag1_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="Coastwide Index lag 1 year") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))



boxplot(fitted(ALL_Coastwide_Index.lag1_fit_noRegion)~FecuTrimmed_allyears$Coastwide_Index.lag1, ylim=c(0,0.4), main="Coastwide Index lag 1 year", xlab="Coastwide Index lag 1 year", ylab="Probability of calf", col=FecuTrimmed_allyears$colour)

plot(fitted(ALL_Coastwide_Index.lag1_fit_noRegion)~jitter(FecuTrimmed_allyears$Coastwide_Index.lag1, 2.5), ylim=c(0,0.4), main="Coastwide Index lag 1 year", xlab="Coastwide Index lag 1 year", ylab="Probability of calf")


#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$CoastwideLag1fact<-as.factor(FecuTrimmed_allyears$Coastwide_Index.lag1)

ggplot(FecuTrimmed_allyears, aes(x=CoastwideLag1fact, y=fitted(ALL_Coastwide_Index.lag1_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.34) +
  scale_x_discrete(breaks=c(0.60754998, 0.781029037, 1.049534295, 1.147386131, 1.445263092), labels=c(0.61, 0.78, 1.05, 1.15, 1.15)) +
  labs(x="Coastwide Index lag 1 year", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))


### SRKW index

ALL_SRKW_Index_fit_noRegion <- glm(Calved~poly(Age,4)+SRKW_Index, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_SRKW_Index_fit_noRegion)

plot(fitted(ALL_SRKW_Index_fit_noRegion, scale='Calved')~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.5))
#visreg(ALL_SRKW_Index_fit_noRegion, scale='response')

plot(fitted(ALL_SRKW_Index_fit_noRegion)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="SRKW Index", xlab="Age", ylab="Probability of calf")



#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_SRKW_Index_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="SRKW Index") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))


boxplot(fitted(ALL_SRKW_Index_fit_noRegion)~FecuTrimmed_allyears$SRKW_Index, ylim=c(0,0.4), main="SRKW Index lag 1 year", xlab="SRKW Index lag 1 year", ylab="Probability of calf")

plot(fitted(ALL_SRKW_Index_fit_noRegion)~jitter(FecuTrimmed_allyears$SRKW_Index, 2.5), ylim=c(0,0.4), main="SRKW Index", xlab="SRKW Index", ylab="Probability of calf")



#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$SRKWfact<-as.factor(FecuTrimmed_allyears$SRKW_Index)

ggplot(FecuTrimmed_allyears, aes(x=SRKWfact, y=fitted(ALL_SRKW_Index_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.33) +
  scale_x_discrete(breaks=c(0.466835142, 0.717500422, 1.01140825, 1.345621478, 1.610760791), labels=c(0.47, 0.72, 1.01, 1.34, 1.61)) +
  labs(x="SRKW Index", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))



### SRKW Index lag 1 year

ALL_SRKW_Index.lag1_fit_noRegion <- glm(Calved~poly(Age,4)+SRKW_Index.lag1, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_SRKW_Index.lag1_fit_noRegion)

plot(fitted(ALL_SRKW_Index.lag1_fit_noRegion)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="SRKW Index lag 1 year", xlab="Age", ylab="Probability of calf")



#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_SRKW_Index.lag1_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="SRKW Index lag 1 year") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))




boxplot(fitted(ALL_SRKW_Index.lag1_fit_noRegion)~FecuTrimmed_allyears$SRKW_Index.lag1, ylim=c(0,0.4), main="SRKW Index lag 1 year", xlab="SRKW Index lag 1 year", ylab="Probability of calf")

plot(fitted(ALL_SRKW_Index.lag1_fit_noRegion)~jitter(FecuTrimmed_allyears$SRKW_Index.lag1, 2.5), ylim=c(0,0.4), main="SRKW Index lag 1 year", xlab="SRKW Index lag 1 year", ylab="Probability of calf")


#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$SRKWLag1fact<-as.factor(FecuTrimmed_allyears$SRKW_Index.lag1)

ggplot(FecuTrimmed_allyears, aes(x=SRKWLag1fact, y=fitted(ALL_SRKW_Index.lag1_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.33) +
  scale_x_discrete(breaks=c(0.466835142, 0.717500422, 1.01140825, 1.345621478, 1.610760791), labels=c(0.47, 0.72, 1.01, 1.34, 1.61)) +
  labs(x="SRKW Index lag 1 year", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))




### NRKW index

ALL_NRKW_Index_fit_noRegion <- glm(Calved~poly(Age,4)+NRKW_Index, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_NRKW_Index_fit_noRegion)

plot(fitted(ALL_NRKW_Index_fit_noRegion, scale='Calved')~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.5))
#visreg(ALL_NRKW_Index_fit_noRegion, scale='response')

plot(fitted(ALL_NRKW_Index_fit_noRegion)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="NRKW Index", xlab="Age", ylab="Probability of calf")



#Make a factor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_NRKW_Index_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="NRKW Index") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))



boxplot(fitted(ALL_NRKW_Index_fit_noRegion)~FecuTrimmed_allyears$NRKW_Index, ylim=c(0,0.4), main="NRKW Index", xlab="NRKW Index", ylab="Probability of calf")

plot(fitted(ALL_NRKW_Index_fit_noRegion)~jitter(FecuTrimmed_allyears$NRKW_Index, 2.5), ylim=c(0,0.4), main="NRKW Index", xlab="NRKW Index", ylab="Probability of calf")


### Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$NRKWfact<-as.factor(FecuTrimmed_allyears$NRKW_Index)

ggplot(FecuTrimmed_allyears, aes(x=NRKWfact, y=fitted(ALL_NRKW_Index_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.355) +
  scale_x_discrete(breaks=c(0.583906655, 0.747866325, 0.940862357, 1.19058655, 1.817588067), labels=c(0.58, 0.73, 0.94, 1.19, 1.82)) +
  labs(x="NRKW Index", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))




#NRKW Index lag 1 year

ALL_NRKW_Index.lag1_fit_noRegion <- glm(Calved~poly(Age,4)+NRKW_Index.lag1, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_NRKW_Index.lag1_fit_noRegion)

plot(fitted(ALL_NRKW_Index.lag1_fit_noRegion)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="NRKW Index lag 1 year", xlab="Age", ylab="Probability of calf")



#Make a agefactor-copy of the index so you don't accidentally run the analysis on a agefactor
FecuTrimmed_allyears$Agefact<-as.factor(FecuTrimmed_allyears$Age)

ggplot(FecuTrimmed_allyears, aes(x=Agefact, y=fitted(ALL_NRKW_Index.lag1_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.4) +
  scale_x_discrete(limits=4:50) +
  labs(x="Age", y="Probability of calving", title="NRKW Index lag 1 year") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))


boxplot(fitted(ALL_NRKW_Index.lag1_fit_noRegion)~FecuTrimmed_allyears$NRKW_Index.lag1, ylim=c(0,0.4), main="NRKW Index lag 1 year", xlab="NRKW Index lag 1 year", ylab="Probability of calf")

plot(fitted(ALL_NRKW_Index.lag1_fit_noRegion)~jitter(FecuTrimmed_allyears$NRKW_Index.lag1, 2.5), ylim=c(0,0.4), main="NRKW Index lag 1 year", xlab="NRKW Index lag 1 year", ylab="Probability of calf")



#Make a factor-copy of the index so you don't accidentally run the analysis on a factor
FecuTrimmed_allyears$NRKWLag1fact<-as.factor(FecuTrimmed_allyears$NRKW_Index.lag1)

ggplot(FecuTrimmed_allyears, aes(x=NRKWLag1fact, y=fitted(ALL_NRKW_Index.lag1_fit_noRegion))) + 
  geom_boxplot(fill="light grey") +
  ylim(-0.001, 0.355) +
  scale_x_discrete(breaks=c(0.583906655, 0.747866325, 0.940862357, 1.19058655, 1.817588067), labels=c(0.58, 0.73, 0.94, 1.19, 1.82)) +
  labs(x="NRKW Index lag 1 year", y="Probability of calving") +
  theme_classic()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))




### Model comparison

AIC(ALL_WCVI_Index_fit_noRegion, ALL_WCVI_Index.lag1_fit_noRegion, ALL_Coastwide_index_fit_noRegion, ALL_Coastwide_Index.lag1_fit_noRegion, ALL_SRKW_Index_fit_noRegion, ALL_SRKW_Index.lag1_fit_noRegion, ALL_NRKW_Index_fit_noRegion, ALL_NRKW_Index.lag1_fit_noRegion)




### JC calculated WCVI index with average calculated from 1980-2017 ----

plot(FecuTrimmed_allyears$Calved~FecuTrimmed_allyears$JC_WCVI_Index)

ALL_JC_WCVI_Index_fit_noRegion <- glm(Calved~poly(Age,4)+JC_WCVI_Index, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_JC_WCVI_Index_fit_noRegion)

plot(fitted(ALL_JC_WCVI_Index_fit_noRegion, scale='Calved')~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.5))
#visreg(ALL_JC_WCVI_Index_fit_noRegion, scale='response')

plot(fitted(ALL_JC_WCVI_Index_fit_noRegion)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="JC_WCVI_Index", xlab="Age", ylab="Probability of calf")

boxplot(fitted(ALL_JC_WCVI_Index_fit_noRegion)~FecuTrimmed_allyears$JC_WCVI_Index, ylim=c(0,0.4), main="JC_WCVI_Index", xlab="JC_WCVI_Index", ylab="Probability of calf")

plot(fitted(ALL_JC_WCVI_Index_fit_noRegion)~jitter(FecuTrimmed_allyears$JC_WCVI_Index, 2.5), ylim=c(0,0.4), main="JC_WCVI_Index", xlab="JC_WCVI_Index", ylab="Probability of calf")



### JC calculated WCVI index with one year lag

plot(FecuTrimmed_allyears$Calved~FecuTrimmed_allyears$JC_WCVI_Index.lag1)

ALL_JC_WCVI_Index.lag1_fit_noRegion <- glm(Calved~poly(Age,4)+JC_WCVI_Index.lag1, family=binomial(link="logit"), data=FecuTrimmed_allyears)
summary(ALL_JC_WCVI_Index.lag1_fit_noRegion)

plot(fitted(ALL_JC_WCVI_Index.lag1_fit_noRegion, scale='Calved')~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.5))
#visreg(ALL_JC_WCVI_Index.lag1_fit_noRegion, scale='response')

plot(fitted(ALL_JC_WCVI_Index.lag1_fit_noRegion)~FecuTrimmed_allyears$Age, xlim=c(0,60), ylim=c(0,0.4), main="JC_WCVI_Index.lag1", xlab="Age", ylab="Probability of calf")

boxplot(fitted(ALL_JC_WCVI_Index.lag1_fit_noRegion)~FecuTrimmed_allyears$JC_WCVI_Index.lag1, ylim=c(0,0.4), main="JC_WCVI_Index.lag1", xlab="JC_WCVI_Index.lag1", ylab="Probability of calf")

plot(fitted(ALL_JC_WCVI_Index.lag1_fit_noRegion)~jitter(FecuTrimmed_allyears$JC_WCVI_Index.lag1, 2.5), ylim=c(0,0.4), main="JC_WCVI_Index.lag1", xlab="JC_WCVI_Index.lag1", ylab="Probability of calf")



