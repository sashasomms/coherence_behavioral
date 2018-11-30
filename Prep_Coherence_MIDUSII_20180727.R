# RStudio version 1.1.453
# R version 3.5 
# Author: Sasha L. Sommerfeldt 
# University of Wisconsin - Madison 
# October 2016 - July 2018 

# Clear the workspace
rm(list=ls())

# Date 
today = '20180710'
#### Directories ####
dir = '~/Desktop/UWMadison/MIDUS'
ddir = paste(dir, '/data', sep='')
adir = paste(dir, '/analysis', sep='')
# Directory with raw data downloaded from http://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/29282 and /04652
# And from http://midus.colectica.org/ for MIDUS 2 Milwaukee subsample 
myddir = '/Volumes/sommerfeldt/fyp/data/publicMIDUS'
# Set working directory 
setwd(ddir)

#### Packages ####
library(data.table)
library(plyr)
library(ggplot2)
library(multilevel)
library(lme4)
library(lmSupport)
library(AICcmodavg)
library(pbkrtest)
library(boot)


#### Read in Data ####
# Project 4 (biomarker)
P4file = paste(myddir,"/ICPSR_29282_m2p4/DS0001/29282-0001-Data.rda", sep='')
dfP4 = load(P4file)
dfP4 = da29282.0001
names(dfP4)

# Project 1 (surveys)
P1file = paste(myddir, '/ICPSR_04652_m2/DS0001/MIDUS2_P1_colectica.csv', sep='')
dfP1 = read.csv(P1file)
names(dfP1)
summary(dfP1$SAMPLMAJ)
length(dfP1$M2ID)

# Project 1 (surveys) for Milwaukee subsample (also zygotic category from M1)
P1Mfile = paste(myddir, '/mke.sav', sep='')
#dfP1M = read.csv(P1Mfile, header=TRUE)
dfP1M = read.spss(P1Mfile, to.data.frame=TRUE)
names(dfP1M)
summary(dfP1M$SAMPLMAJ)
length(dfP1M$M2ID)
varDescribe(dfP1M)

# Zygosity data from Midus I 
# Downloaded from colectica with variables: M2ID, ZYGCAT, TOT_SIBS
PZfile = paste(myddir, '/zygoSibs.sav', sep='')
dfPZ = read.spss(PZfile, to.data.frame=TRUE)
names(dfPZ)
length(dfPZ$M2ID)




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
miscO = c('B1PGENDER', 'B1PAGE_M2', 'B4ZAGE', 'B4ZB1PLG', 'B4ZB1SLG', 'B4ZB1CLG', 'B4QTA_AX', 'B4QCESD', 'B4H1I', 'B4PBMI', 'B4BIL6', 'B4BCRP')
#d[miscO]
misc = c('gender', 'P1_PIage', 'P4_age', 'months_P1PI_to_P4', 'months_P1SAQ_to_P4', 'months_P1cog_to_P4', 'P4_STAItrait', 'P4_CESD', 'P4_diabetes', 'P4_BMI', 'IL6', 'CRP')


## Remove whitespace from task variable values
dfP4$B4VTASK1str = varRecode(dfP4$B4VTASK1, c("STROOP    ", "MATH      ", "INAPPLIC  ", "PASAT     "), c('STROOP', 'MATH', 'INAPPLIC', 'PASAT'))

## Rename variables with my intuitive names
setnames(dfP4, old=stressSRO, new=stressSR)
setnames(dfP4, old=ecg_HRO, new=ecg_HR)
setnames(dfP4, old=ecg_QO, new=ecg_Q)
setnames(dfP4, old=miscO, new=misc)
sort(names(dfP4))

sort(names(dfP4))[2400:2600]




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
count(is.na(dfP4$hr1)) # 149 - 102 = 47
count(is.na(dfP4$hrM)) # 225 - 159 = 66
count(is.na(dfP4$hr3)) # 172 - 130 = 42
count(is.na(dfP4$hrS)) # 218 - 142 = 76
count(is.na(dfP4$hr5)) # 179 - 136 = 43
count(is.na(dfP4$hr6)) # 194 - 151 = 43

#### Clean diabetes data ####
dfP4$P4_diabetes = varRecode(dfP4$P4_diabetes, c("(1) YES", "(2) NO", "(3) BORDERLINE"), c(3, 1, 2)) # 1 = not diabetic, 2 = borderline diabetic, 3 = diabetic
dfP4$P4_diabetes = as.numeric(dfP4$P4_diabetes)
varDescribeBy(dfP4$M2ID, dfP4$P4_diabetes)

#### Score BMI data ####
# Underweight: BMI is less than 18.5
# Normal weight: BMI is 18.5 to 24.9
# Overweight: BMI is 25 to 29.9
# Obese: BMI is 30 or more
# via: https://www.cancer.org/cancer/cancer-causes/diet-physical-activity/body-weight-and-cancer-risk/adult-bmi.html
varDescribe(dfP4$P4_BMI)
# Exclude underweight
dfP4['P4_BMI'][ (dfP4['P4_BMI'] < 18.5 ) ] = NA

# Subset variables I want into a separate dataframe (dfP4ss)
P4cols = c("M2ID", "M2FAMNUM", "SAMPLMAJ", 'B4VTASK1str', misc, stressSR, ecg_HR, ecg_Q)
names(dfP4[P4cols])
dfP4ss = dfP4[P4cols]


#### Sort out order of stressor tasks ####
# ECG variables are tied to the task they were collected during/after, but tasks were counterbalanced. 
# Stress self-reports are tied to the order they were measured.
# I need heart rate in order measured, task irrelevant, so can look at relation to stress self-report.  
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
## Compute within subject stress x heart-rate correlation magnitudes and store as a variable
## Also plot each subjects stress and heart rate data individually (for QA) 

# Start a PDF file for the individual subject plots (X=stress, Y=heart rate)
pdf(paste(adir,"/coherence_hr_subjplots_",today,".pdf",sep=''))
# For each subject
for (s in dfP4ss2$M2ID) {
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

## change 98s to NA
dfP1[PWB2O][dfP1[PWB2O] == 98] = NA
dfP1[COPEO][dfP1[COPEO] == 98] = NA

# Look-see
varDescribe(dfP1[PWB2O])
varDescribe(dfP1[COPEO])

# Rename columns
setnames(dfP1, old=PWB2O, new=PWB2)
setnames(dfP1, old=COPEO, new=COPE)
length(dfP1$M2ID)

# Look-see
varDescribe(dfP1[PWB2])
varDescribe(dfP1[COPE])
sapply(dfP1[PWB2], class) # factors
sapply(dfP1[COPE], class) # factors
# Convert to numeric 
dfP1[PWB2] = sapply(dfP1[PWB2], as.numeric)
dfP1[COPE] = sapply(dfP1[COPE], as.numeric)

## Composite PWB
dfP1$pwb2 = varScore(dfP1, Forward = PWB2, MaxMiss = .0)
varDescribe(dfP1$pwb2)

# Miscellaneous
P1miscO = c('B1PBYEAR','B1PRSEX','B1PF7A', 'B1PF2A', 'pwb2')
P1misc = c('birth_year','P1_sex','P1_race', 'P1_ethnicity', 'pwb2')
setnames(dfP1, old=P1miscO, new=P1misc)

P1cols = c("M2ID", P1misc, PWB2, COPE)
# Subset P1 data
dfP1ss = dfP1[P1cols]

########################

############################################################
#### Clean & organize Project 1 (Survey) MILWAUKEE data ####
############################################################

## MIDUS2 version PWB (7-items per sub-scale)
PWB2O = c('BASPWBA2', 'BASPWBE2', 'BASPWBG2', 'BASPWBR2', 'BASPWBU2', 'BASPWBS2') # old names
PWB2 = c('autonomy2', 'envMast2', 'persGrow2', 'posRela2', 'purpLife2', 'selfAcce2') # new names

COPEO = c('BASEMCOP', 'BASPRCOP', 'BASDENIA', 'BASVENT', 'BASDISEN', 'BASREINT', 'BASACTIV', 'BASPLAN') # old names
COPE = c('COPEem', 'COPEprob', 'COPE_denial', 'COPE_vent', 'COPE_disengage', 'COPE_posReGrow', 'COPE_active', 'COPE_plan') # new names

# Rename columns
setnames(dfP1M, old=PWB2O, new=PWB2)
setnames(dfP1M, old=COPEO, new=COPE)
length(dfP1M$M2ID)

# Look-see
varDescribe(dfP1M[PWB2])
varDescribe(dfP1M[COPE])
sapply(dfP1M[PWB2], class) # numeric
sapply(dfP1M[COPE], class) # numeric

## Calculate composite PWB
dfP1M$pwb2 = varScore(dfP1M, Forward = PWB2, MaxMiss = .0)

# Miscellaneous
P1MmiscO = c( 'BACBYR', 'BACRSEX', 'BACF7A', 'BACF2A', 'pwb2')
P1Mmisc = c('birth_year', 'P1_sex', 'P1_race', 'P1_ethnicity', 'pwb2')
setnames(dfP1M, old=P1MmiscO, new=P1Mmisc)

class(dfP1$birth_year)
dfP1M$birth_year = sapply(dfP1M$birth_year, as.factor)


P1Mcols = c("M2ID", P1Mmisc, PWB2, COPE)
# Subset P1 data
dfP1Mss = dfP1M[P1Mcols]


varDescribe(dfP1Mss)
varDescribe(dfP1ss)
# make sure no IDs are duplicated
summary(duplicated(c(dfP1ss$M2ID, dfP1Mss$M2ID)))

#### Merge P1 and P1 Milwaukee ####
dfP1_P1M = rbind(dfP1ss, dfP1Mss)
varDescribe(dfP1_P1M)
########################

#####################
####  Merge data #### 
#####################
# Merge zyogsity and twins data to P1 data
dfTemp2 = merge.data.frame(dfP1_P1M, dfPZ,  by='M2ID', all=TRUE)
# Merge that zygosity/twins & P1 data to P4 data
dfTemp = merge.data.frame(dfTemp2, dfP4ss2,  by='M2ID', all=TRUE)


df = dfTemp

############################################
####  Create a M2FAMNUM for MKE subjects #### 
############################################
# MKE subjects don't have a FAMNUM assigned, but none were related. 
# Give each MKE subject a unique M2FAMNUM for analyses.
# Start numbering at 130000
famnum = 130000
for ( mke in df$M2ID[df$SAMPLMAJ == '(13) MILWAUKEE'] ) {
  if ( !is.na(mke) ) {
    print(mke)
    df$M2FAMNUM[df$M2ID == mke] = famnum
    famnum = famnum + 1
  }
} 


################################
### Convert to long format  #### 
################################
## Long format 
# Make a copy of the dataframe to convert to a data table 
# because data tables link to all versions of table, which causes issues when trying to change column names in only one copy
dt = df
# All variables except stress[1-6], hr[1-6], ecgQ[1-6]
names(dt)
varsid = c(names(dt[1:37]), names(dt[56:57]))
dfLTemp = melt.data.table(setDT(dt),
                       id.vars = varsid, # ID variables - all the variables to keep but not split apart on
                       measure = patterns("^stress", "^hr", "^ecgQ"),
                       variable.name = "timepoint")
OLD = c("value1", "value2", "value3")
NEW = c("stress", "hr", "ecgQ")
setnames(dfLTemp, old=OLD, new=NEW)
# Convert back to beloved data frame
dfLTemp = data.frame(dfLTemp)
names(dfLTemp)

#### Remove orthostatic stress timepoint from long format df ####
dfLnoO = dfLTemp[dfLTemp$timepoint != 6,] 

dfL = dfLnoO

#### Create a 'Complete' variable: Find subjects with only all 5 timepoints of stress SRs and HR data ####
completeSubjs = NA
for ( s in unique(dfL$M2ID) ) {
  #print(s)
  # Subset and transpose subject's stress SRs and heart rate
  SUBstress = dfL$stress[dfL$M2ID == s]
  SUBhr = dfL$hr[dfL$M2ID == s]
  # Put subjects stress and heart rate into their own dataframe
  SUBdf = data.frame(SUBstress, SUBhr)
  names(SUBdf)[1] = "stress"
  names(SUBdf)[2] = "hr"
  # Count number of NAs in stress self-reports
  cS = count((is.na(SUBdf$stress)))
  stressNotNA = as.numeric(cS[cS$x == 'FALSE',][2])
  dfL$stressNotNA[dfL$M2ID == s] = stressNotNA
  
  # Count number of NAs in heart rate
  cH = count((is.na(SUBdf$hr)))
  hrNotNA = as.numeric(cH[cH$x == 'FALSE',][2])
  dfL$hrNotNA[dfL$M2ID == s] = hrNotNA
  
  if ( (!is.na(hrNotNA)) & (!is.na(stressNotNA)) ) {
    # Compute correlations only for subjects with all 5 timepoints 
    if ( (hrNotNA == 5) & (stressNotNA == 5) ) {
      print(s)
      completeSubjs = c(completeSubjs, s)
      dfL$complete[dfL$M2ID == s] = 1
    }
    if ( (hrNotNA != 5) & (stressNotNA != 5) ) {
      dfL$complete[dfL$M2ID == s] = 0
    }
    
  }
  
  rm(cS, cH, stressNotNA, hrNotNA) # clear variables before looping for next subject
}

length(completeSubjs) # 968
count(dfL, 'hrNotNA')
count(dfL, 'stressNotNA')

##################################
#### Extract slopes from LMEM ####
##################################
# cluster mean-center stress
dfL$stress_CMC = dfL$stress - ave(dfL$stress, dfL$M2ID, na.rm=T)
lmerS = lmer(hr ~ stress_CMC + (1+ stress_CMC|M2ID), data=dfL)
Anova(lmerS, type=3, test="F")
modelSummary(lmerS)
# The slopes -- actually eBLUPs
slopes = coef(lmerS)$M2ID
#names(slopes)
#slopes
dfSlope = data.frame(slopes)
#names(dfSlope)
#dfSlope$stress_CMC
#rownames(dfSlope)

dfSlope2 = data.frame(M2ID = row.names(dfSlope), dfSlope$stress_CMC)
#dfSlope2
names(dfSlope2)[names(dfSlope2) == 'dfSlope.stress_CMC'] = 'coherence_slope'
# Merge slopes into long format dataframe
dfL = merge.data.frame(dfL, dfSlope2, by='M2ID', all=TRUE)
# Merge slopes into wide format dataframe
df = merge.data.frame(df, dfSlope2, by='M2ID', all=TRUE)
varDescribe(df$coherence_slope)

################################
#### Write out my data file #### 
################################

today='20181124'

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

