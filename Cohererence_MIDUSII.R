# RStudio Version 1.0.136
# Author: Sasha L. Sommerfeldt 
# University of Wisconsin - Madison 
# October 2016 - January 2018 

# Clear the workspace
rm(list=ls())

# Date 
today = '20180115'
#### Directories ####
dir = '~/Desktop/UWMadison/MIDUS'
ddir = paste(dir, '/data', sep='')
adir = paste(dir, '/analysis', sep='')
# Directory with raw data downloaded from http://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/29282 and /04652
myddir = '/Volumes/sommerfeldt/fyp/data/publicMIDUS'
# Set working directory 
setwd(ddir)

#### Packages ####
library('foreign') # to read spss file 
library('data.table')
library(multilevel)
library(lme4)
library(AICcmodavg)
library(pbkrtest)
library(boot)
source("~/Desktop/Statistics710/homework/indirectMLM.R")

#### Read in Data ####
# Project 4 (biomarker)
P4file = paste(myddir,"/ICPSR_29282_m2p4/DS0001/29282-0001-Data.rda", sep='')
dfP4 = load(P4file)
dfP4 = da29282.0001
names(dfP4)

# Project 1 (surveys)
P1file = paste(myddir, '/ICPSR_04652_m2/DS0001/04652-0001-Data.sav', sep='')
dfP1 = read.spss(P1file, to.data.frame=TRUE)
names(dfP1)
summary(dfP1$SAMPLMAJ)
length(dfP1$M2ID)




#####################################################
#### Clean & organize Project 4 (Biomarker) data ####
#####################################################
## Self reports of stress
stressSRO = c('B4VSRB1', 'B4VSRCS1', 'B4VSRR1', 'B4VSRCS2', 'B4VSRR2', 'B4VSRU1') # old names
stressSR = c('stress1', 'stress2', 'stress3', 'stress4', 'stress5', 'stress6') # new names

## Heart rate
ecg_HRO = c('B4VB1HR', 'B4VM1HR', 'B4VR1HR', 'B4VS1HR', 'B4VR2HR', 'B4VU1HR') # old names
ecg_HR = c('hr1', 'hrM', 'hr3', 'hrS', 'hr5', 'hr6') # new names


## ECG quality 
# Translate ECG quality ratings to number format 
old = c('(1) EXCELLENT', '(2) GOOD', '(3) SCOREABLE', '(4) UNSCOREABLE', '(5) NO DATA')
new = c(1, 2, 3, 4, 5)
dfP4$B4VBEQn = varRecode(dfP4$B4VBEQ, old, new)
dfP4$B4VMEQn = varRecode(dfP4$B4VMEQ, old, new)
dfP4$B4VPEQn = varRecode(dfP4$B4VPEQ, old, new)
dfP4$B4VSEQn = varRecode(dfP4$B4VSEQ, old, new)
dfP4$B4VR1EQn = varRecode(dfP4$B4VR1EQ, old, new)
dfP4$B4VR2EQn = varRecode(dfP4$B4VR2EQ, old, new)
dfP4$B4VUEQn = varRecode(dfP4$B4VUEQ, old, new)

# c('B4VBEQ', 'B4VMEQ', 'B4VPEQ', 'B4VSEQ', 'B4VR1EQ', 'B4VR2EQ', 'B4VUEQ') # original names 
ecg_QO = c('B4VBEQn', 'B4VMEQn', 'B4VR1EQn', 'B4VSEQn', 'B4VR2EQn', 'B4VUEQn') # original names, in number format
ecg_Q = c('ecgQ1', 'ecgQM', 'ecgQ3', 'ecgQS', 'ecgQ5', 'ecgQ6') # new names 

## Miscellaneous variables to rename
miscO = c('B1PGENDER', 'B1PAGE_M2', 'B4QTA_AX', 'B4QCESD', 'B4H1I', 'B4PBMI')
#d[miscO]
misc = c('gender', 'age', 'P4_STAItrait', 'P4_CESD', 'P4_diabetes', 'P4_BMI')

## Remove whitespace from task variable values
dfP4$B4VTASK1str = varRecode(dfP4$B4VTASK1, c("STROOP    ", "MATH      ", "INAPPLIC  ", "PASAT     "), c('STROOP', 'MATH', 'INAPPLIC', 'PASAT'))

## Rename variables with my intuitive names
setnames(dfP4, old=stressSRO, new=stressSR)
setnames(dfP4, old=ecg_HRO, new=ecg_HR)
setnames(dfP4, old=ecg_QO, new=ecg_Q)
setnames(dfP4, old=miscO, new=misc)
sort(names(dfP4))

# Check some quality ratings before removing any bad data
varDescribeBy(dfP4$hr1, dfP4$ecgQ1) # 1 = 742, 2 = 364, 3 = 47
varDescribeBy(dfP4$hrM, dfP4$ecgQM) # 1 = 667, 2 = 363, 3 = 66
varDescribeBy(dfP4$hr3, dfP4$ecgQ3) # 1 = 718, 2 = 365, 3 = 42
varDescribeBy(dfP4$hrS, dfP4$ecgQS) # 1 = 663, 2 = 374, 3 = 76
varDescribeBy(dfP4$hr5, dfP4$ecgQ5) # 1 = 698, 2 = 378, 3 = 43
varDescribeBy(dfP4$hr6, dfP4$ecgQ6) # 1 = 716, 2 = 345, 3 = 43 
count(is.na(dfP4$hr1)) # 102
count(is.na(dfP4$hrM)) # 159
count(is.na(dfP4$hr3)) # 130
count(is.na(dfP4$hrS)) # 142
count(is.na(dfP4$hr5)) # 136
count(is.na(dfP4$hr6)) # 151


#### Clean ECG Data ####

## Replace all ecg data that is not good (2) or excellent (1) with NA 
dfP4['hr1'][ (dfP4['ecgQ1'] != 1) & (dfP4['ecgQ1'] != 2) ] = NA
dfP4['hrM'][ (dfP4['ecgQM'] != 1) & (dfP4['ecgQM'] != 2) ] = NA
dfP4['hr3'][ (dfP4['ecgQ3'] != 1) & (dfP4['ecgQ3'] != 2) ] = NA
dfP4['hrS'][ (dfP4['ecgQS'] != 1) & (dfP4['ecgQS'] != 2) ] = NA
dfP4['hr5'][ (dfP4['ecgQ5'] != 1) & (dfP4['ecgQ5'] != 2) ] = NA
dfP4['hr6'][ (dfP4['ecgQ6'] != 1) & (dfP4['ecgQ6'] != 2) ] = NA

# How many removed? 
count(is.na(dfP4$hr1)) # 144 - 102 = 42
count(is.na(dfP4$hrM)) # 225 - 159 = 66
count(is.na(dfP4$hr3)) # 172 - 130 = 42
count(is.na(dfP4$hrS)) # 218 - 142 = 76
count(is.na(dfP4$hr5)) # 179 - 136 = 43
count(is.na(dfP4$hr6)) # 194 - 151 = 43

# Subset variables I want into a separate dataframe (dfP4ss)
P4cols = c("M2ID", "M2FAMNUM", "SAMPLMAJ", 'B4VTASK1str', misc, stressSR, ecg_HR, ecg_Q,'B4BIL6', 'B4BCRP')
names(dfP4[P4cols])
dfP4ss = dfP4[P4cols]


#### Sort out order of stressor tasks ####
# ECG variables are tied to the task they were collected during/after, but tasks were counterbalanced. 
# Stress self-reports in order they were measured.
# I need heart rate in order measured, task irrelevant, so can look at coherence with stress self-report.  
summary(dfP4ss$B4VTASK1str)

## Subset stroop-firsts and maths-firsts into 2 separate dataframes
# dS = stroop first
dS = subset(dfP4ss, B4VTASK1str == 'STROOP')
length(dS$M2ID) # 603

# dM = math first
dM = subset(dfP4ss, B4VTASK1str == 'MATH')
length(dM$M2ID) # 592

## Stroop first (stroop = 2, math = 4)
names(dS)[names(dS) == 'hrS'] = 'hr2'
names(dS)[names(dS) == 'hrM'] = 'hr4'
names(dS)[names(dS) == 'ecgQS'] = 'ecgQ2'
names(dS)[names(dS) == 'ecgQM'] = 'ecgQ4'

# Check did not edit other data frame like data.table sometimes does
names(dS)
names(dM) # good, still M and S

## Math first (math = 2, stroop = 4)
names(dM)[names(dM) == 'hrM'] = 'hr2'
names(dM)[names(dM) == 'hrS'] = 'hr4'
names(dM)[names(dM) == 'ecgQM'] = 'ecgQ2'
names(dM)[names(dM) == 'ecgQS'] = 'ecgQ4'

# Merge stroop-first and math-first dataframes back together
dfP4ss2 = merge.data.frame(dM, dS, all=TRUE)

# New ECG variable names with numbers instead of 's' or 'm'
ecg_Qn = c('ecgQ1', 'ecgQ2', 'ecgQ3', 'ecgQ4', 'ecgQ5', 'ecgQ6')
ecg_HRn = c('hr1', 'hr2', 'hr3', 'hr4', 'hr5', 'hr6')


#### Within-subj correlations between stress and heart rate ####
## Compute within subject stress x heart-rate correlation magnitudes and store as a variable (initial check of distribution of coherence)
## Also plot each subjects stress and heart rate data individually (for QA) 

# Start a PDF file for the individual subject plots (X=stress, Y=heart rate)
pdf(paste(adir,"/coherence_hr_subjplots_",today,".pdf",sep=''))
# For each subject
for (s in dfP4ss2$M2ID) {
  s = 10231
  print(s)
  # Subset and transpose subject's stress SRs and heart rate
  SUBstress = t(dfP4ss2[stressSR[1:5]][dfP4ss2$M2ID == s,])
  SUBhr = t(dfP4ss2[ecg_HRn[1:5]][dfP4ss2$M2ID == s,])
  # Put subjects stress and heart rate into their own dataframe
  SUBdf = data.frame(SUBstress, SUBhr)
  names(SUBdf)[1] = "stress"
  names(SUBdf)[2] = "hr"
  
  # Count number of NAs in stress self-reports
  cS = count((is.na(SUBdf$stress)))
  stressNotNA = as.numeric(cS[cS$x == 'FALSE',][2])
  
  # Count number of NAs in heart rate
  cH = count((is.na(SUBdf$hr)))
  hrNotNA = as.numeric(cH[cH$x == 'FALSE',][2])
  
  # Title for plot 
  tit = paste("Subject",s)
  
  # If there is at least 1 set of values (otherwise code bombs)
  if ( (!is.na(hrNotNA)) & (!is.na(stressNotNA)) ) {
    # Compute subject's stressSR x heart rate correlation
    SUBcor = cor(SUBdf, method="pearson", use="complete.obs")
    # Plot it 
    myplot = ggplot(SUBdf,aes(x=stress,y=hr))+
      geom_point(stat="identity",size=4.5, shape=1) +
      ggtitle(tit) +
      scale_x_continuous(name = "Self-reported stress", breaks = 1:10, limits = c(1,10)) +
      scale_y_continuous(name = "Heart rate", limits = c(40,130)) +
      theme_bw()+
      theme(panel.grid.minor = element_blank(), axis.text=element_text(size=14),
            axis.title=element_text(size=16), plot.title=element_text(size=16, hjust=.5)) +
      annotate("text", x=9.3, y=120, label = paste("r =",round(SUBcor[1,2],3)) ) # put correlation magnitude in plot
    print(myplot)
    # 
    dfP4ss2$coherence_as_r[dfP4ss2$M2ID==s] = as.numeric(SUBcor[1,2]) # add subject correlation magnitude to dataframe
    
    # Compute correlations only for subjects with all 5 timepoints 
    if ( (hrNotNA == 5) & (stressNotNA == 5) ) {
      SUBcor5 = cor(SUBdf, method="pearson", use="all.obs")
      dfP4ss2$coherence_as_r5[dfP4ss2$M2ID==s] = as.numeric(SUBcor5[1,2])
    }
    if ( (hrNotNA != 5) & (stressNotNA != 5) ) {
      dfP4ss2$coherence_as_r5[dfP4ss2$M2ID==s] = NA
    }
    
  }
  rm(SUBstress, SUBhr, SUBdf, SUBcor, SUBcor5, cS, cH, stressNotNA, hrNotNA, myplot) # clear variables before looping for next subject
}
dev.off() # Finish PDF with subject plots

# Check distribution
varDescribe(dfP4ss2$coherence_as_r) # n = 1036, mean = .49, sd = .47, median = .65, skew = -1.18, kurt = .54
hist(dfP4ss2$coherence_as_r, main='Coherence', xlab='Coherence (stress x HR correlation)', col='#1F618D', cex.axis=1.7, cex.main=1.6)
varDescribe(dfP4ss2$coherence_as_r5) # n = 944, mean = .5, sd = .46, median = .66, skew = -1.16, kurt = .46
hist(dfP4ss2$coherence_as_r5, main='Coherence - complete data only', xlab='Coherence (stress x HR correlation)', col='#1F618D', cex.axis=1.7, cex.main=1.6)


##################################################
#### Clean & organize Project 1 (Survey) data ####
##################################################
## MIDUS2 version PWB (7-items per sub-scale)
PWB2O = c('B1SPWBA2', 'B1SPWBE2', 'B1SPWBG2', 'B1SPWBR2', 'B1SPWBU2', 'B1SPWBS2') # old names
PWB2 = c('autonomy2', 'envMast2', 'persGrow2', 'posRela2', 'purpLife2', 'selfAcce2') # new names

COPEO = c('B1SEMCOP', 'B1SPRCOP', 'B1SDENIA', 'B1SVENT', 'B1SDISEN', 'B1SREINT', 'B1SACTIV', 'B1SPLAN') # old names
COPE = c('COPEem', 'COPEprob', 'COPE_denial', 'COPE_vent', 'COPE_disengage', 'COPE_posReGrow', 'COPE_active', 'COPE_plan') # new names

# Rename columns
setnames(dfP1, old=PWB2O, new=PWB2)
length(dfP1$M2ID)

## Composite PWB
dfP1$pwb2 = varScore(dfP1, Forward = PWB2, MaxMiss = .0)

# Miscellaneous
P1miscO = c('B1PBYEAR','B1PRSEX','pwb2')
P1misc = c('birth_year','P1_sex','pwb2')
setnames(dfP1, old=P1miscO, new=P1misc)

P1cols = c("M2ID", P1misc, PWB2)
# Subset P1 data
dfP1ss = dfP1[P1cols]

########################

#####################
####  Merge data #### 
#####################
dfTemp = merge.data.frame(dfP4ss2, dfP1ss, by='M2ID', all=TRUE)
names(dfTemp)

df = dfTemp

################################
### Convert to long format  #### 
################################
## Long format 
# Make a copy of the dataframe to convert to a data table 
# because data tables link to all versions of table
dt = df
# All variables except stress[1-6], hr[1-6], ecgQ[1-6]
dt$stress1

varsid = c(names(dt[1:10]), names(dt[29:40]))
dfLTemp = melt.data.table(setDT(dt),
                       id.vars = varsid, # ID variables - all the variables to keep but not split apart on
                       measure = patterns("^stress", "^hr", "^ecgQ"),
                       variable.name = "timepoint")
OLD = c("value1", "value2", "value3")
NEW = c("stress", "hr", "ecgQ")
setnames(dfLTemp, old=OLD, new=NEW)
# convert to beloved data.frame
dfLTemp = data.frame(dfLTemp)
names(dfLTemp)
View(dfLTemp)

#### Remove orthostatic stress from long format df ####
dfLnoO = dfLTemp[dfLTemp$timepoint != 6,] 

dfL = dfLnoO


################################
#### Write out my data file #### 
################################

today='20180115'

# Wide format
fnameW = paste("coh_",today,".csv",sep='')
fpathW = paste(ddir,"/",fnameW, sep='')
write.csv(df, file=fpathW)
fpathW = paste(myddir,"/",fnameW, sep='')
write.csv(df, file=fpathW)

# Long format
fnameL = paste("cohLong_",today,".csv",sep='')
fpathL = paste(ddir,"/",fnameL, sep='')
write.csv(dfL, file=fpathL)
fpathL = paste(myddir,"/",fnameL, sep='')
write.csv(dfL, file=fpathL)


#### Read in data in future without all that trouble #### 
df = read.csv(fpathW)
dfL = read.csv(fpathL)

#################
#### ANALYZE #### 
#################

#######################################
####   LMER/LONG format analysis   #### 
#######################################

#### Cluster Mean Center ####
#varDescribeBy(dfL$stress, dfL$timepoint)
dfL$stress_CMC = dfL$stress - ave(dfL$stress, dfL$M2ID, na.rm=T)

#varDescribeBy(dfL$hr, dfL$timepoint)
dfL$hr_CMC = dfL$hr - ave(dfL$hr, dfL$M2ID, na.rm=T)

#varDescribeBy(dfL$stress_CMC, dfL$timepoint)
#varDescribeBy(dfL$hr_CMC, dfL$timepoint)
#varDescribe(dfL$timepoint)

dfL$pwb2_C = dfL$pwb2 - mean(dfL$pwb2, na.rm=T)
dfL$P4_CESD_C = dfL$P4_CESD- mean(dfL$P4_CESD, na.rm=T)
dfL$P4_STAItrait_C = dfL$P4_STAItrait - mean(dfL$P4_STAItrait, na.rm=T)
dfL$P5_ERQ_s_C = dfL$P5_ERQ_s - mean(dfL$P5_ERQ_s, na.rm=T)
dfL$B4BIL6_C = dfL$B4BIL6 - mean(dfL$B4BIL6, na.rm=T)
dfL$B4BCRP_C = dfL$B4BCRP - mean(dfL$B4BCRP, na.rm=T)
dfL$B1SPRCOP_C = dfL$B1SPRCOP - mean(dfL$B1SPRCOP, na.rm=T)
dfL$B1SEMCOP_C = dfL$B1SEMCOP - mean(dfL$B1SEMCOP, na.rm=T)
dfL$coherence_slope_C = dfL$coherence_slope - mean(dfL$coherence_slope, na.rm=T)
dfL$coherence_slopeNoAge_C = dfL$coherence_slopeNoAge - mean(dfL$coherence_slopeNoAge, na.rm=T)

# Log transformed inflammatory markers
dfL$B4BIL6_T = log2(dfL$B4BIL6)
dfL$B4BCRP_T = log(dfL$B4BCRP, base=10)
dfL$B4BIL6_T_C = dfL$B4BIL6_T - mean(dfL$B4BIL6_T, na.rm=T)
dfL$B4BCRP_T_C = dfL$B4BCRP_T - mean(dfL$B4BCRP_T, na.rm=T)

dfL$age_C = dfL$age - mean(dfL$age, na.rm=T)
dfL$gender_C = varRecode(dfL$gender, c('(1) MALE', '(2) FEMALE'), c(-.5,.5))

names(dfL)



