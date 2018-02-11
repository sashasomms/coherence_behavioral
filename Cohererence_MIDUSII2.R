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
# And from http://midus.colectica.org/ for MIDUS 2 Milwaukee subsample 
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
setnames(dfP1, old=COPEO, new=COPE)
length(dfP1$M2ID)

## Composite PWB
dfP1$pwb2 = varScore(dfP1, Forward = PWB2, MaxMiss = .0)

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



## Calculate composite PWB
dfP1M$pwb2 = varScore(dfP1M, Forward = PWB2, MaxMiss = .0)

# Miscellaneous
P1MmiscO = c( 'BACBYR', 'BACRSEX', 'BACF7A', 'BACF2A', 'pwb2')
P1Mmisc = c('birth_year', 'P1_sex', 'P1_race', 'P1_ethnicity', 'pwb2')
setnames(dfP1M, old=P1MmiscO, new=P1Mmisc)

P1Mcols = c("M2ID", P1Mmisc, PWB2, COPE)
# Subset P1 data
dfP1Mss = dfP1M[P1Mcols]


varDescribe(dfP1Mss)
varDescribe(dfP1ss)
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
# because data tables link to all versions of table
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
# convert to beloved data.frame
dfLTemp = data.frame(dfLTemp)
names(dfLTemp)

#### Remove orthostatic stress timepoint from long format df ####
dfLnoO = dfLTemp[dfLTemp$timepoint != 6,] 

dfL = dfLnoO

#### 'Complete' variable: Find subjects with only all 5 timepoints of stress SRs and HR data ####
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
  
  # Count number of NAs in heart rate
  cH = count((is.na(SUBdf$hr)))
  hrNotNA = as.numeric(cH[cH$x == 'FALSE',][2])
  
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

##################################
#### Extract slopes from LMEM ####
##################################
# cluster mean-center stress
dfL$stress_CMC = dfL$stress - ave(dfL$stress, dfL$M2ID, na.rm=T)
lmerS = lmer(hr ~ stress_CMC + (1+ stress_CMC|M2ID), data=dfL)
Anova(lmerS, type=3, test="F")
modelSummary(lmerS)
# The slopes
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

today='20180207'

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

#### Plot individual subject slopes ####
dfL$stressMC = dfL$stress - ave(dfL$stress, dfL$M2ID)
dfL$hrM = ave(dfL$stress, dfL$M2ID)

ggplot(dfL, aes(stress, hr, color=as.factor(M2ID)))+
  geom_smooth(aes(group=as.factor(M2ID)),method="lm",se=F,size=.1, alpha=.6, position="jitter")+
  xlim(c(0,11))+
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text=element_text(size=14), axis.title=element_text(size=24)) +
  labs(x="Self-Reported Stress", y="Heart Rate")+
  theme(legend.position="none")

## Make with heat map where color is the magnitude of the slope 
mycol = c("#0710C4", "gray",# negative & zero
           "#FFEC00", "#FFC300", "#FF5733", "#C70039", "#900C3F", "#581845")   # positive
mybreaks = c(-.4, 0, 
             .5, 1, 1.5, 2, 3, 4)

ggplot(dfL, aes(stress, hr, color=coherence_slope))+
  geom_smooth(aes(group=as.factor(M2ID)),method="lm",se=F,size=.1, alpha=.6, position="jitter")+
  xlim(c(0,11))+
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text=element_text(size=14), axis.title=element_text(size=24)) +
  labs(x="Self-Reported Stress", y="Heart Rate")+
  scale_colour_gradientn("Slope",colours=mycol, limits=c(-.4, 4.5), values = scales::rescale(c(-0.5, -0.05, 0, 0.05, 0.5,1,2,3,4)), breaks = mybreaks, guide="colourbar")

#### Histogram of coherence slopes ####
ggplot(dfLsW, aes(coherence_slope)) +
  geom_histogram(aes(fill=as.factor(coherence_slope)), binwidth=.2, col="black", fill="#FFC300") +
  #scale_fill_gradientn("Slope", colours=mycol, limits=c(-.4, 4.5), values = scales::rescale(c(-0.5, -0.05, 0, 0.05, 0.5,1,2,3,4)), breaks = mybreaks, guide="colourbar")+
  labs(x="Slope", y="Number of Subjects") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=24)) +
  theme(legend.position="none")


###########################
####   LMER analysis   #### 
###########################

#### Create a separate dataframe for analysis - excluding the many survey/P1 people without biomarker/P4/coherence data
dfLs = dfL[!is.na(dfL$coherence_slope),]
length(unique(dfLs$M2ID)) # 1065

#### Demographics ####

# Move long data back to wide format
dfLsW = reshape(dfLs, idvar = "M2ID", v.names=c('hr', 'stress', 'stress_CMC', 'ecgQ'), drop=c('X', 'stressMC'), timevar = "timepoint", direction = "wide")

names(dfLsW)

summary(dfLsW$gender) # Male = 455, Female = 610
varDescribe(dfLsW$months_P1SAQ_to_P4) # M = 25.89, SD = 14.19, range = 0 - 62

varDescribe(dfLsW$P4_age) # M = 56.4, SD = 11.21, range = 35-86
varDescribe(dfLsW$P1_PIage) # M = 53.55, SD = 11.4, range = 34-83
varDescribe(dfLsW$months_P1SAQ_to_P4) # M = 25.89, SD = 14.19, range = 0-62
varDescribe(dfLsW$months_P1PI_to_P4) # M = 28.4, SD = 13.93, range = 5-63
varDescribe(dfLsW$months_P1cog_to_P4) # M = 23.62, SD = 12.64, range = 1-61 (N = 973)
varDescribe(dfLsW$pwb2) # N = 1061, M = 232.81 (35.25)
varDescribe(dfLsW$P4_CESD) # N = 1057, M = 8.61 (8.1)
varDescribe(dfLsW$P4_STAItrait) # N = 1057, M = 34.2 (8.98)
varDescribe(dfLsW$IL6) # N = 1058, M = 2.96 (2.89)
varDescribe(dfLsW$CRP) # N = 1052, M = 2.85 (4.26)
varDescribeBy(dfLsW$M2ID, dfLsW$P4_diabetes) # N = 1063, M = 1.25 (.66). 1 = 929, 2 = 2, 3 = 132
varDescribe(dfLsW$P4_BMI) # N = 1062, M = 29.83 (6.63)
varDescribe(dfLsW$COPEem) # N = 1060, M = 22.28 (5.51)
varDescribe(dfLsW$COPEprob) # N = 1059, M = 38.05 (5.94)

summary(dfLsW$P1_race) # Asian = 3, black = 193, Native american or alaska native aleutian islander/eskimo = 14, other = 27, white = 825, DK = 1, REFUSED =1, NA=1
summary(dfLsW$P1_ethnicity)

summary(dfLsW$SAMPLMAJ) # main = 521, Sibling = 6, twin = 337, city oversample = 19, milwaukee=182
## Twins
# Create a data frame with a list of family IDs and the number of times they occurred.
n_occur = data.frame(table(dfLsW$M2FAMNUM))
n_occur[n_occur$Freq > 1,]
# Which M2FAMNUMs occurred more than once.
dfLsW$M2ID[dfLsW$M2FAMNUM %in% n_occur$Var1[n_occur$Freq > 1]]

# Count up the number of subjects with each repeated family id, each subject's value for sibs column is how many people have same family id (including themself)
dfLsW$sibs = 1
for (family in n_occur$Var1[n_occur$Freq > 1] ) {
  print(family)
  numsibs = n_occur$Freq[n_occur$Var1 == family ]
  dfLsW$sibs[dfLsW$M2FAMNUM == family] = numsibs
}
length(dfLsW$M2ID[dfLsW$SAMPLMAJ == '(03) TWIN' & dfLsW$sibs == 2 & dfLsW$ZYGCAT == 'MONOZYGOTIC']) # 128
length(dfLsW$M2ID[dfLsW$SAMPLMAJ == '(03) TWIN' & dfLsW$sibs == 2 & dfLsW$ZYGCAT == 'DIZYGOTIC - SAME SEX']) # 56
length(dfLsW$M2ID[dfLsW$SAMPLMAJ == '(03) TWIN' & dfLsW$sibs == 2 & dfLsW$ZYGCAT == 'DIZYGOTIC - DIFFERENT SEX']) # 46
length(dfLsW$M2ID[dfLsW$SAMPLMAJ == '(03) TWIN' & dfLsW$sibs == 2 & dfLsW$ZYGCAT == 'UNABLE TO DETERMING ZYGOSITY']) # 2
View(dfLsW[dfLsW$SAMPLMAJ == '(03) TWIN' & dfLsW$sibs == 4,]) # fam of 4 is DZ same sex 

# 1 family with 2 twin pairs, 1 family with 3 siblings
View(dfLsW[dfLsW$SAMPLMAJ != '(03) TWIN',])

# 3 siblings
length(dfLsW$M2ID[dfLsW$SAMPLMAJ != '(03) TWIN' & dfLsW$sibs == 3]) # N = 3 (1 family)
length(dfLsW$M2ID[dfLsW$SAMPLMAJ != '(03) TWIN' & dfLsW$sibs == 2]) # N = 8 (4 families)

varDescribe(unique(dfLsW$M2FAMNUM[dfLsW$SAMPLMAJ == '(02) SIBLING'])) # 5 unique family ids that are siblings, 1 pair of siblings 


#### Prep variables ####
# Have age for everyone (so don't need to recenter well-being variable)
# Stress is centered within cluster (centered around each subject's mean)
# Thus: for each analysis, just need to re-center age based on who has that well-being variable
# For mediation, will need to re-center based on who has complete COPE data AND the outcome var

## Cluster Mean Center Stress ##
dfLs$stress_CMC = dfLs$stress - ave(dfLs$stress, dfLs$M2ID, na.rm=T)

## Mean Center ##
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age, na.rm=T)

dfLs$pwb2_C = dfLs$pwb2 - mean(dfLs$pwb2, na.rm=T)
dfLs$P4_CESD_C = dfLs$P4_CESD- mean(dfLs$P4_CESD, na.rm=T)
dfLs$P4_STAItrait_C = dfLs$P4_STAItrait - mean(dfLs$P4_STAItrait, na.rm=T)
dfLs$IL6_C = dfLs$IL6 - mean(dfLs$IL6, na.rm=T)
dfLs$CRP_C = dfLs$CRP - mean(dfLs$CRP, na.rm=T)
dfLs$COPEprob_C = dfLs$COPEprob - mean(dfLs$COPEprob, na.rm=T)
dfLs$COPEem_C = dfLs$COPEem - mean(dfLs$COPEem, na.rm=T)

## Re-code ##
dfLs$gender_C = varRecode(dfLs$gender, c('(1) MALE', '(2) FEMALE'), c(-.5,.5))

## Log transform inflammatory markers for normal distribution 
dfLs$IL6_T = log2(dfLs$IL6)
hist(dfLs$IL6_T)
dfLs$CRP_T = log(dfLs$CRP, base=10)
hist(dfLs$CRP_T)
# then center the transformed versions
dfLs$IL6_T_C = dfLs$IL6_T - mean(dfLs$IL6_T, na.rm=T)
dfLs$CRP_T_C = dfLs$CRP_T - mean(dfLs$CRP_T, na.rm=T)

names(dfLs)




#### LMER TESTS ####

#### Stress & heart rate ####
lmerM = lmer(hr ~ stress_CMC + P4_age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = .877, F(1, 858.7) = 670.10, p < .0001 (< 2.2e-16)

#### P4_age ####
lmerM = lmer(hr ~ stress_CMC * P4_age_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = -0.008, F(1, 867.0) = 7.757, p = .005


#### Gender ####
lmerM = lmer(hr ~ stress_CMC * gender_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
Anova(lmerM, type=3, test="F") # NS
modelSummary(lmerM)
# b = .051, F(1, 876.8) = .560, p = .455

#### PWB ####
# Center age for subjects in this analysis 
varDescribe(dfLs$pwb2_C)
length(dfLs$P4_age[!is.na(dfLs$pwb2_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$pwb2_C)], na.rm=T)
# Run the test 
lmerM = lmer(hr ~ stress_CMC * pwb2_C + P4_age_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = 4.217e-03, F(1, 840.3) = 1.947e+01, p < .0001 (1.16e-05)
# b = .004217, F(1, 840.3) = 19.47, p < .0001 (p = .0000116)
# SE = 9.555e-04 = .00096

#### Depression ####
# Center age for subjects in this analysis 
varDescribe(dfLs$P4_CESD_C)
length(dfLs$P4_age[!is.na(dfLs$P4_CESD_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$P4_CESD_C)], na.rm=T)
# Run the test
lmerM = lmer(hr ~ stress_CMC * P4_CESD_C + P4_age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = -0.022, F(1, 807.4) = 28.628, p < .0001 (1.14e-07)
# b = -0.022, F(1, 807.4) = 28.628, p < .0001 (p = .000000114)
# SE = 0.004061

#### Anxiety ####
# Center age for subjects in this analysis 
varDescribe(dfLs$P4_STAItrait_C)
length(dfLs$P4_age[!is.na(dfLs$P4_STAItrait_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$P4_STAItrait_C)], na.rm=T)
# Run the test
lmerM = lmer(hr ~ stress_CMC * P4_STAItrait_C + P4_age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = -0.018, F(1, 793.5) = 2.523e+01, p < .0001 (6.29e-07)
# b = -0.018, F(1, 793.5) = 25.23, p < .0001 (p = .000000629)
# SE = 0.003656

#### IL6 ####
# Center age for subjects in this analysis 
varDescribe(dfLs$IL6_T_C)
length(dfLs$P4_age[!is.na(dfLs$IL6_T_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$IL6_T_C)], na.rm=T)
# Run the test
lmerM = lmer(hr ~ stress_CMC * IL6_T_C + P4_age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = -0.155, F(1, 802.8) = 25.49, p < .0001 (5.5e-07)
# b = -0.155, F(1, 802.8) = 25.49, p < .0001 (p = .00000055)
# SE = 0.03066

#### CRP ####
# Center age for subjects in this analysis 
varDescribe(dfLs$CRP_T_C)
length(dfLs$P4_age[!is.na(dfLs$CRP_T_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$CRP_T_C)], na.rm=T)
# Run the test
lmerM = lmer(hr ~ stress_CMC * CRP_T_C + P4_age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = -0.166, F(1, 857.0) = 6.394, p = .0116
# SE = 0.06568


#### Diabetes ####
# Center age for subjects in this analysis 
varDescribe(dfLs$P4_diabetes)
length(dfLs$P4_age[!is.na(dfLs$P4_diabetes)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$P4_diabetes)], na.rm=T)
# Run the test
lmerM = lmer(hr ~ stress_CMC * P4_diabetes + P4_age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = -.120, F(1, 812.9) = 5.673, p = .0175
# SE = 0.05035

#### BMI ####
# Center age for subjects in this analysis 
varDescribe(dfLs$P4_BMI)
length(dfLs$P4_age[!is.na(dfLs$P4_BMI)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$P4_BMI)], na.rm=T)
# Run the test
lmerM = lmer(hr ~ stress_CMC * P4_BMI + P4_age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# With <18.5: b = -0.011, F(1, 824.4) = 5.112, p = .024
# Without <18.5: b = -0.012, F(1, 822.1) = 5.2776, p = .0219
# SE = 0.005041


#### Emotion-focused coping ####
# Center age for subjects in this analysis 
varDescribe(dfLs$COPEem_C)
length(dfLs$P4_age[!is.na(dfLs$COPEem_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$COPEem_C)], na.rm=T)
# Run the test
lmerM = lmer(hr ~ stress_CMC * COPEem_C + P4_age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = -0.032, F(1, 864.7) = 2.728e+01, p < .0001 (2.20e-07)
# b = -0.032, F(1, 864.7) = 27.28, p < .0001 (p = .00000022)

#### Problem-focused coping ####
# Center age for subjects in this analysis 
varDescribe(dfLs$COPEprob_C)
length(dfLs$P4_age[!is.na(dfLs$COPEprob_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$COPEprob_C)], na.rm=T)
# Run the test
lmerM = lmer(hr ~ stress_CMC * COPEprob_C + P4_age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = .009, F(1, 833.6) = 2.881, p = .090


#### All-Together Now ####
dfA = dfLsW
## Mean Center ##
dfA$P4_age_C = dfA$P4_age - mean(dfA$P4_age, na.rm=T)

dfA$pwb2_C = dfA$pwb2 - mean(dfA$pwb2, na.rm=T)
dfA$P4_CESD_C = dfA$P4_CESD- mean(dfA$P4_CESD, na.rm=T)
dfA$P4_STAItrait_C = dfA$P4_STAItrait - mean(dfA$P4_STAItrait, na.rm=T)
dfA$IL6_C = dfA$IL6 - mean(dfA$IL6, na.rm=T)
dfA$CRP_C = dfA$CRP - mean(dfA$CRP, na.rm=T)
dfA$COPEprob_C = dfA$COPEprob - mean(dfA$COPEprob, na.rm=T)
dfA$COPEem_C = dfA$COPEem - mean(dfA$COPEem, na.rm=T)

## Re-code ##
dfA$gender_C = varRecode(dfA$gender, c('(1) MALE', '(2) FEMALE'), c(-.5,.5))

## Log transform inflammatory markers for normal distribution 
dfA$IL6_T = log2(dfA$IL6)
hist(dfA$IL6_T)
dfA$CRP_T = log(dfA$CRP, base=10)
hist(dfA$CRP_T)
# then center the transformed versions
dfA$IL6_T_C = dfA$IL6_T - mean(dfA$IL6_T, na.rm=T)
dfA$CRP_T_C = dfA$CRP_T - mean(dfA$CRP_T, na.rm=T)

# Psych
lmerM = lmer(coherence_slope ~ P4_age_C + pwb2_C + P4_CESD_C + P4_STAItrait_C + (1|M2FAMNUM), data=dfA)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)

# Physical
lmerM = lmer(coherence_slope ~ P4_age_C + IL6_T_C + CRP_T_C + P4_BMI + P4_diabetes + (1|M2FAMNUM), data=dfA)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)

# Literally everything
lmerM = lmer(coherence_slope ~ P4_age_C + COPEem_C + pwb2_C + P4_CESD_C + P4_STAItrait_C + IL6_T_C + CRP_T_C + P4_BMI + P4_diabetes + (1|M2FAMNUM), data=dfA)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)


#### MEDIATION -LMEM ####
#### mLMEM. PWB ####
# New dataframe for analysis subsample (complete data on all variables involved in each mediation)
dfM = dfLsW[!is.na(dfLsW$pwb2) & !is.na(dfLsW$COPEem),]
length(dfM$M2ID) # 1057
### CENTER FOR ANALYSIS SUBSAMPLE
# Age
dfM$P4_age_C = dfM$P4_age - mean(dfM$P4_age, na.rm=T)
# COPE
dfM$COPEem_C = dfM$COPEem - mean(dfM$COPEem, na.rm=T)
# Coherence slope
dfM$coherence_slope_C = dfM$coherence_slope - mean(dfM$coherence_slope, na.rm=T)

### (total effect/c) 1. Coherence -> WB
lmerM1 = lmer(pwb2 ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM1, type=3, test="F")
modelSummary(lmerM1)
# b = 11.52, F(1, 1051.3) = 28.17, p = 1.36e-07 = .000000136
# b = 11.52, F(1, 1051.3) = 28.17, p < .0001
# SE = 2.16610

### (path a) 2. Coherence -> COPE (mediator)
lmerM2 = lmer(COPEem ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM2, type=3, test="F")
modelSummary(lmerM2)
# b = -1.79052, F(1, 1051.1) = 26.809, p = 2.69e-07 = .000000269
# b = -1.79, F(1, 1051.1) = 26.81, p < .0001 
# SE = 0.01541 

### 3. Coherence -> WB controlling for COPE (mediator)
lmerM3 = lmer(pwb2 ~ P4_age_C + coherence_slope_C + COPEem_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM3, type=3, test="F")
modelSummary(lmerM3)
# (direct effect/c') Coherence: b = 6.88616, F(1, 1052.8) = 11.84, p = .000601 # SE = 1.99664
# Coherence: b = 6.89, F(1, 1052.8) = 11.84, p < .001

# (path b) COPEem : b = -2.61036, F(1, 1052.6) = 219.36, p <2.2e-16 # SE = 0.17594
# COPEem : b = -2.61, F(1, 1052.6) = 219.36, p < .0001

### 4. Indirect significant? 
med_pwb2 = mediate(lmerM2, lmerM3, treat = "coherence_slope_C", mediator = "COPEem_C")
summary(med_pwb2)
# ACME/indirect: 4.728, CI = 2.943 to 6.680
# ADE/direct:  6.823, CI = 2.869 to 10.963
# % mediated: 41.1%



#### mLMEM. Depression ####
# New dataframe for analysis subsample (complete data on all variables involved in each mediation)
dfM = dfLsW[!is.na(dfLsW$P4_CESD) & !is.na(dfLsW$COPEem),]
length(dfM$M2ID) # 1052
### CENTER FOR ANALYSIS SUBSAMPLE
# Age
dfM$P4_age_C = dfM$P4_age - mean(dfM$P4_age, na.rm=T)
# COPE
dfM$COPEem_C = dfM$COPEem - mean(dfM$COPEem, na.rm=T)
# Coherence slope
dfM$coherence_slope_C = dfM$coherence_slope - mean(dfM$coherence_slope, na.rm=T)

### (total effect/c) 1. Coherence -> WB 
lmerM1 = lmer(P4_CESD ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM1, type=3, test="F")
modelSummary(lmerM1)
# b = -2.92835, F(1, 1042.8) = 34.18, p = 6.71e-09 = .00000000671
# b = -2.93, F(1, 1042.8) = 34.18, p < .0001
# SE = .49972

### (path a) 2. Coherence -> COPE (mediator)
lmerM2 = lmer(COPEem ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM2, type=3, test="F")
modelSummary(lmerM2)
# b = -1.80120, F(1, 1046.5) = 27.525, p = 1.88e-07 = .000000188
# b = -1.80, F(1, 1046.5) = 27.53, p < .0001
# SE = .34256

### 3. Coherence -> WB controlling for COPE (mediator)
lmerM3 = lmer(P4_CESD ~ P4_age_C + coherence_slope_C + COPEem_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM3, type=3, test="F")
modelSummary(lmerM3)
# (direct effect/c') Coherence: b = -2.12014, F(1, 1047.8) = 19.28, p = 1.25e-05 = .0000125
# Coherence: b = -2.12, F(1, 1047.8) = 19.28, p < .0001

# (path b) COPEem : b = .45536, F(1, 1047.7) = 112.49, p < 2e-16
# COPEem : b = .46, F(1, 1047.7) = 112.49, p < .0001

### 4. Indirect significant? 
med_P4_CESD = mediate(lmerM2, lmerM3, treat = "coherence_slope_C", mediator = "COPEem_C")
summary(med_P4_CESD)
# ACME/indirect: -0.817, CI = -1.176 to -0.498
# ADE/direct: -2.119, CI = -3.023 to -1.089
# % mediated: 27.8%



#### mLMEM. Anxiety ####
# New dataframe for analysis subsample (complete data on all variables involved in each mediation)
dfM = dfLsW[!is.na(dfLsW$P4_STAItrait) & !is.na(dfLsW$COPEem),]
length(dfM$M2ID) # 1052
### CENTER FOR ANALYSIS SUBSAMPLE
# Age
dfM$P4_age_C = dfM$P4_age - mean(dfM$P4_age, na.rm=T)
# COPE
dfM$COPEem_C = dfM$COPEem - mean(dfM$COPEem, na.rm=T)
# Coherence slope
dfM$coherence_slope_C = dfM$coherence_slope - mean(dfM$coherence_slope, na.rm=T)

### (total effect/c) 1. Coherence -> WB 
lmerM1 = lmer(P4_STAItrait ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM1, type=3, test="F")
modelSummary(lmerM1)
# b = -3.14613, F(1, 1039.7) = 32.52, p = 1.53e-08 = .0000000153
# b = -3.15, F(1, 1039.7) = 32.52, p < .0001
# SE = .55038

### (path a) 2. Coherence -> COPE (mediator)
lmerM2 = lmer(COPEem ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM2, type=3, test="F")
modelSummary(lmerM2)
# b = -1.80304, F(1, 1046.6) = 27.578, p = 1.83e-07 = .000000183
# b = -1.80, F(1, 1046.6) = 27.58, p < .0001
# SE = .34258

### 3. Coherence -> WB controlling for COPE (mediator)
lmerM3 = lmer(P4_STAItrait ~ P4_age_C + coherence_slope_C + COPEem_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM3, type=3, test="F")
modelSummary(lmerM3)
# (direct effect/c') Coherence: b = -1.94136, F(1, 1047.1) = 14.66, p = .000136
# Coherence: b = -1.94, F(1, 1047.1) = 14.66, p < .001

# (path b) COPEem : b = .68401, F(1, 1048.0) = 230.23, p < 2e-16
# COPEem : b = .68, F(1, 1048.0) = 230.23, p < .0001

### 4. Indirect significant? 
med_P4_STAItrait = mediate(lmerM2, lmerM3, treat = "coherence_slope_C", mediator = "COPEem_C")
summary(med_P4_STAItrait)
# ACME/indirect: -1.236, -1.725 to -0.752
# ADE/direct: -1.956, -2.923 to -1.058
# % mediated: 38.7%


#### mLMEM. IL-6 ####
# New dataframe for analysis subsample (complete data on all variables involved in each mediation)
dfM = dfLsW[!is.na(dfLsW$IL6) & !is.na(dfLsW$COPEem),]
length(dfM$M2ID) # N = 1053
### CENTER FOR ANALYSIS SUBSAMPLE
# IL6
dfM$IL6_T = log2(dfM$IL6)
hist(dfM$IL6_T)
# Age
dfM$P4_age_C = dfM$P4_age - mean(dfM$P4_age, na.rm=T)
# COPE
dfM$COPEem_C = dfM$COPEem - mean(dfM$COPEem, na.rm=T)
# Coherence slope
dfM$coherence_slope_C = dfM$coherence_slope - mean(dfM$coherence_slope, na.rm=T)

### (total effect/c) 1. Coherence -> WB 
lmerM1 = lmer(IL6_T ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM1, type=3, test="F")
modelSummary(lmerM1)
# b = -0.27484, F(1, 1048.5) = 17.25, p = .0000355
# b = = -0.27, F(1, 1048.5) = 17.25, p < .0001
# SE = .06603

### (path a) 2. Coherence -> COPE (mediator)
lmerM2 = lmer(COPEem ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM2, type=3, test="F")
modelSummary(lmerM2)
# b = -1.80737, F(1, 1046.7) = 26.710, p = 2.83e-07
# b = -1.81, F(1, 1046.7) = 26.71, p < .0001
# SE = .34892

### 3. Coherence -> WB controlling for COPE (mediator)
lmerM3 = lmer(IL6_T ~ P4_age_C + coherence_slope_C + COPEem_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM3, type=3, test="F")
modelSummary(lmerM3)
# (direct effect/c') Coherence: b = -0.238168, F(1, 1049.0) = 12.77, p = .000368
# Coherence: b = -0.24, F(1, 1049.0) = 12.77, p < .001

# (path b) COPEem : b = .020752, F(1, 1049.0) = 12.73, p = .000375
# COPEem : b = .02, F(1, 1049.0) = 12.73, p < .001

### 4. Indirect significant? 
med_IL6 = mediate(lmerM2, lmerM3, treat = "coherence_slope_C", mediator = "COPEem_C")
summary(med_IL6)
# ACME/indirect: -0.0371, CI = -0.0639 to -0.0141
# ADE/direct: -0.2366, CI = -0.3641 to -0.1064
# % mediated: 13.62%



#### mLMEM. CRP ####
# New dataframe for analysis subsample (complete data on all variables involved in each mediation)
dfM = dfLsW[!is.na(dfLsW$CRP) & !is.na(dfLsW$COPEem),]
length(dfM$M2ID) # N = 1047
### CENTER FOR ANALYSIS SUBSAMPLE
# CRP
dfM$CRP_T = log(dfM$CRP, base=10)
hist(dfM$CRP_T)
# Age
dfM$P4_age_C = dfM$P4_age - mean(dfM$P4_age, na.rm=T)
# COPE
dfM$COPEem_C = dfM$COPEem - mean(dfM$COPEem, na.rm=T)
# Coherence slope
dfM$coherence_slope_C = dfM$coherence_slope - mean(dfM$coherence_slope, na.rm=T)

### (total effect/c) 1. Coherence -> WB 
lmerM1 = lmer(CRP_T ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM1, type=3, test="F")
modelSummary(lmerM1)
# b = -0.0636751, F(1, 1029.9) = 3.9268, p = .0478
# b = -0.06, F(1, 1029.9) = 3.93, p < .05
# SE = .0014381

### (path a) 2. Coherence -> COPE (mediator)
lmerM2 = lmer(COPEem ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM2, type=3, test="F")
modelSummary(lmerM2)
# b = -1.76927, F(1, 1040.7) = 25.600, p = 4.96e-07
# b = -1.77, F(1, 1040.7) = 25.60, p < .0001
# SE = .34889

### 3. Coherence -> WB controlling for COPE (mediator)
lmerM3 = lmer(CRP_T ~ P4_age_C + coherence_slope_C + COPEem_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM3, type=3, test="F")
modelSummary(lmerM3)
# (direct effect/c') Coherence: b = -0.0567437, F(1, 1032.1) = 3.0442, p = .0813
# Coherence: b = -0.06, F(1, 1032.1) = 3.04, p = .08

# (path b) COPEem : b = 0.0038934, F(1, 1041.3) = 1.8635, p  = .1725
# COPEem : b = .004, F(1, 1041.3) = 1.86, p = .1725

### 4. Indirect significant? 
med_CRP = mediate(lmerM2, lmerM3, treat = "coherence_slope_C", mediator = "COPEem_C")
summary(med_CRP)
# ACME/indirect: b = -0.007073, CI = -0.017934 to 0.003359, p = .16
# ADE/direct: -0.054810, CI = -0.116582 to 0.004111
# % mediated: 10.83



#### mLMEM. BMI ####
# New dataframe for analysis subsample (complete data on all variables involved in each mediation)
dfM = dfLsW[!is.na(dfLsW$P4_BMI) & !is.na(dfLsW$COPEem),]
length(dfM$M2ID) # N = 1057
### CENTER FOR ANALYSIS SUBSAMPLE
# Age
dfM$P4_age_C = dfM$P4_age - mean(dfM$P4_age, na.rm=T)
# COPE
dfM$COPEem_C = dfM$COPEem - mean(dfM$COPEem, na.rm=T)
# WB
dfM$P4_BMI_C = dfM$P4_BMI - mean(dfM$P4_BMI, na.rm=T)
# Coherence slope
dfM$coherence_slope_C = dfM$coherence_slope - mean(dfM$coherence_slope, na.rm=T)

### (total effect/c) 1. Coherence -> WB 
lmerM1 = lmer(P4_BMI ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM1, type=3, test="F")
modelSummary(lmerM1)
# b = -0.87913, F(1, 781.3) = 781.3, p = .0244
# b = -0.88, F(1, 781.3) = 781.3, p < .05
# SE = .38843

### (path a) 2. Coherence -> COPE (mediator)
lmerM2 = lmer(COPEem ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM2, type=3, test="F")
modelSummary(lmerM2)
# b = -1.77760, F(1, 1051.1) = 26.314, p = 3.45e-07
# b = -1.78, F(1, 1051.1) = 26.31, p < .0001
# SE = .34576

### 3. Coherence -> WB controlling for COPE (mediator)
lmerM3 = lmer(P4_BMI ~ P4_age_C + coherence_slope_C + COPEem_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM3, type=3, test="F")
modelSummary(lmerM3)
# (direct effect/c') Coherence: b = -0.870530, F(1, 791.3) = 4.837, p = .0281
# Coherence: b = -0.87, F(1, 791.3) = 4.84, p < .05

# (path b) COPEem : b = .004719, F(1, 903.8) = .01774, p = .8941
# COPEem : b = .0047, F(1, 903.8) = .018, p = .89

### 4. Indirect significant? 
med_P4_BMI = mediate(lmerM2, lmerM3, treat = "coherence_slope_C", mediator = "COPEem_C")
summary(med_P4_BMI)
# ACME/indirect: -0.012, -0.140 to 0.119, p = .86
# ADE/direct: -0.871, -1.638 to -0.129, p = .02
# % mediated: 1.3%



#### mLMEM. Diabetes ####
# New dataframe for analysis subsample (complete data on all variables involved in each mediation)
dfM = dfLsW[!is.na(dfLsW$P4_diabetes) & !is.na(dfLsW$COPEem),]
length(dfM$M2ID) # N = 1058
### CENTER FOR ANALYSIS SUBSAMPLE
# Age
dfM$P4_age_C = dfM$P4_age - mean(dfM$P4_age, na.rm=T)
# COPE
dfM$COPEem_C = dfM$COPEem - mean(dfM$COPEem, na.rm=T)
# Coherence slope
dfM$coherence_slope_C = dfM$coherence_slope - mean(dfM$coherence_slope, na.rm=T)

### (total effect/c) 1. Coherence -> WB 
lmerM1 = lmer(P4_diabetes ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM1, type=3, test="F")
modelSummary(lmerM1)
# b = -0.061364, F(1, 1043.8) = 2.179, p = .14022
# b = -0.06, F(1, 1043.8) = 2.18, p = .14
# SE = .041473

### (path a) 2. Coherence -> COPE (mediator)
lmerM2 = lmer(COPEem ~ coherence_slope_C + P4_age_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM2, type=3, test="F")
modelSummary(lmerM2)
# b = -1.80400, F(1, 1050.9) = 27.113, p = 2.31e-07
# b = -1.80, F(1, 1050.9) = 27.11, p <.0001
# SE = .34567

### 3. Coherence -> WB controlling for COPE (mediator)
lmerM3 = lmer(P4_diabetes ~ P4_age_C + coherence_slope_C + COPEem_C + (1|M2FAMNUM), data=dfM)
Anova(lmerM3, type=3, test="F")
modelSummary(lmerM3)
# (direct effect/c') Coherence: b = -0.074786, F(1, 1045.2) = 3.163, p = .07562
# Coherence: b = -0.07,  F(1, 1045.2) = 3.16, p = .08

# (path b) COPEem : b = -0.007400, F(1, 1053.2) = 4.001, p = .04573
# COPEem : b = -0.01, F(1, 1053.2) = 4.00, p < .05

### 4. Indirect significant? 
med_P4_diabetes = mediate(lmerM2, lmerM3, treat = "coherence_slope_C", mediator = "COPEem_C")
summary(med_P4_diabetes)
# ACME/indirect: .01357, 0.00109 to 0.02778, p = .03
# ADE/direct: -0.07509, -0.15395 to 0.00459, p = .07

# % mediated: -0.18968







