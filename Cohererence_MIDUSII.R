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

#### Functions ####
# APA writeUp function by Adrienne R. Wood
writeUp = function(myModel) {
  mod = myModel
  coeffs = modelSummary(mod,t=F, Print = F)
  effects = modelEffectSizes(mod, Print = F)
  myVars = c() # Makes a vector of your predictor variables' names
  for (i in attributes(mod$terms)["term.labels"]) {
    myVars = append(myVars,i)
  }
  fstart = (length(myVars)+1)*2+1 # specifies the location in the modelSummary and modelEffect sizes output where our Fs and ps will be found
  pstart = (length(myVars)+1)*3+1
  v = 1
  for (variable in myVars) {
    fstart = fstart+1
    pstart = pstart+1
    v = v+1 ## SLS edit (also below for b)
    thisOutput = paste0(variable,": b = ",round(mod$coefficients[v],3),", F(1, ",mod$df.residual,") = ",
                        round(coeffs$coefficients[fstart],3),", p = ", round(coeffs$coefficients[pstart],3),
                        ", partial η² = ", round(effects$Effects[fstart],3))
    thisOutput = gsub("p = 0,", "p < .001,", thisOutput)
    print(thisOutput)
  }
}

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
miscO = c('B1PGENDER', 'B1PAGE_M2', 'B4QTA_AX', 'B4QCESD', 'B4H1I', 'B4PBMI', 'B4BIL6', 'B4BCRP')
#d[miscO]
misc = c('gender', 'age', 'P4_STAItrait', 'P4_CESD', 'P4_diabetes', 'P4_BMI', 'IL6', 'CRP')

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
count(is.na(dfP4$hr1)) # 149 - 102 = 47
count(is.na(dfP4$hrM)) # 225 - 159 = 66
count(is.na(dfP4$hr3)) # 172 - 130 = 42
count(is.na(dfP4$hrS)) # 218 - 142 = 76
count(is.na(dfP4$hr5)) # 179 - 136 = 43
count(is.na(dfP4$hr6)) # 194 - 151 = 43

#### Clean diabetes data ####
dfP4$P4_diabetes = varRecode(dfP4$P4_diabetes, c("(1) YES", "(2) NO", "(3) BORDERLINE"), c(1, 3, 2))
dfP4$P4_diabetes = as.numeric(dfP4$P4_diabetes)
varDescribe(dfP4$P4_diabetes)

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
P1miscO = c('B1PBYEAR','B1PRSEX','B1PF7A', 'pwb2')
P1misc = c('birth_year','P1_sex','P1_race', 'pwb2')
setnames(dfP1, old=P1miscO, new=P1misc)

P1cols = c("M2ID", P1misc, PWB2, COPE)
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

varsid = c(names(dt[1:12]), names(dt[31:50]))
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

length(completeSubjs)

##################################
#### Extract slopes from LMEM ####
##################################
# cluster mean-center
dfL$stress_CMC = dfL$stress - ave(dfL$stress, dfL$M2ID, na.rm=T)
lmerS = lmer(hr ~ stress_CMC + (1+ stress_CMC|M2ID), data=dfL)
Anova(lmerS, type=3, test="F")
modelSummary(lmerS)
# The slopes
slopes = coef(lmerS)$M2ID
names(slopes)
slopes
dfSlope = data.frame(slopes)
names(dfSlope)
dfSlope$stress_CMC
rownames(dfSlope)

dfSlope2 = data.frame(M2ID = row.names(dfSlope), dfSlope$stress_CMC)
dfSlope2
names(dfSlope2)[names(dfSlope2) == 'dfSlope.stress_CMC'] = 'coherence_slope'
dfL = merge.data.frame(dfL, dfSlope2, by='M2ID', all=TRUE)


################################
#### Write out my data file #### 
################################

today='20180123'

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





###########################
####   LMER analysis   #### 
###########################

#### Save a separate file for analysis - excluding the many survey people without coherence data
dfLs = dfL[!is.na(dfL$coherence_slope),]
length(unique(dfLs$M2ID))
varDescribeBy(dfLs$coherence_slope, dfLs$SAMPLMAJ)
View(dfLs)

dfLc = dfL[dfL$complete == 1,]
length(unique(dfLc$M2ID))
varDescribeBy(dfL$M2ID, dfL$stress)


today = '20180116'
fnameL = paste("cohLong_",today,".csv",sep='')
fpathL = paste(ddir,"/",fnameL, sep='')
write.csv(dfL, file=fpathL)

#### Demographics table ####
varDescribeBy(dfLs$M2ID)
# Move long data back to wide format



#### Prep variables ####
## Cluster Mean Center ##
dfL$stress_CMC = dfL$stress - ave(dfL$stress, dfL$M2ID, na.rm=T)
dfL$hr_CMC = dfL$hr - ave(dfL$hr, dfL$M2ID, na.rm=T)

## Mean Center ##
dfL$age_C = dfL$age - mean(dfL$age, na.rm=T)
dfL$pwb2_C = dfL$pwb2 - mean(dfL$pwb2, na.rm=T)
dfL$P4_CESD_C = dfL$P4_CESD- mean(dfL$P4_CESD, na.rm=T)
dfL$P4_STAItrait_C = dfL$P4_STAItrait - mean(dfL$P4_STAItrait, na.rm=T)
dfL$IL6_C = dfL$IL6 - mean(dfL$IL6, na.rm=T)
dfL$CRP_C = dfL$CRP - mean(dfL$CRP, na.rm=T)
dfL$COPEprob_C = dfL$COPEprob - mean(dfL$COPEprob, na.rm=T)
dfL$COPEem_C = dfL$COPEem - mean(dfL$COPEem, na.rm=T)

## Re-code ##
dfL$gender_C = varRecode(dfL$gender, c('(1) MALE', '(2) FEMALE'), c(-.5,.5))

## Log transform inflammatory markers
dfL$IL6_T = log2(dfL$IL6)
hist(dfL$IL6_T)
dfL$CRP_T = log(dfL$CRP, base=10)
hist(dfL$CRP_T)
dfL$IL6_T_C = dfL$IL6_T - mean(dfL$IL6_T, na.rm=T)
dfL$CRP_T_C = dfL$CRP_T - mean(dfL$CRP_T, na.rm=T)

names(dfL)




#### LMER TESTS ####

#### Stress & heart rate ####
lmerM = lmer(hr ~ stress_CMC + age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = .928, F(1, 731.3) = 610.2, p = .0001

#### Age ####
lmerM = lmer(hr ~ stress_CMC * age_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = -0.010, F(1, 745.2) = 9.458, p = .002


#### Gender ####
lmerM = lmer(hr ~ stress_CMC * gender_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F") # NS
modelSummary(lmerM)


#### PWB ####
lmerM = lmer(hr ~ stress_CMC * pwb2_C + age_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = .003, F(1, 729) = 9.969e+00, p = .001657

#### Depression ####
lmerM = lmer(hr ~ stress_CMC * P4_CESD_C + age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)


#### Anxiety ####
lmerM = lmer(hr ~ stress_CMC * P4_STAItrait_C + age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)


#### IL6 ####
lmerM = lmer(hr ~ stress_CMC * IL6_T_C + age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)


#### CRP ####
varDescribe(dfL$CRP_T_C)
hist(dfL$CRP_T_C)
lmerM = lmer(hr ~ stress_CMC * CRP_T_C + age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = -0.130, F(1, 735.9) = 2.954, p = .086

#### Diabetes ####
dfL$P4_diabetes = as.numeric(dfL$P4_diabetes)
varDescribe(dfL$P4_diabetes)
hist(dfL$P4_diabetes)
lmerM = lmer(hr ~ stress_CMC * P4_diabetes + age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)


#### BMI ####
varDescribe(dfL$P4_BMI)
hist(dfL$P4_BMI)
lmerM = lmer(hr ~ stress_CMC * P4_BMI + age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)



#### Emotion-focused coping ####
lmerM = lmer(hr ~ COPEem_C + age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)


#### Problem-focused coping ####
lmerM = lmer(hr ~ COPEprob_C + age_C + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)



#### MEDIATION -LMEM ####
### 2. Coherence -> COPE 
## Problem-focused coping
lmerM = lmer(hr ~ stress_CMC * COPEprob_C + age_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = .009, F(1, 721.56) = 1.826, p = .177

# Slope extracted
lmerM = lmer(COPEprob ~ coherence_slopeNoAge_C + age_C + (1|M2FAMNUM), data=df)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = .72, F(1,870.4) = 3.24, p = .072

## Emotion-focused coping
lmerM = lmer(hr ~ stress_CMC * COPEem_C + age_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfL)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = -0.033, F(1, 735.51) = 22.204, p < .0001

## Slope extracted
lmerM = lmer(COPEem ~ coherence_slopeNoAge_C + age_C + (1|M2FAMNUM), data=df)
Anova(lmerM, type=3, test="F")
modelSummary(lmerM)
# b = -1.61, F(1, 878.38) = 20.84, p < .0001


#### mLMEM. Depression ####
lmerM1 = lmer(COPEem ~ coherence_slopeNoAge_C + age_C + (1|M2FAMNUM), data=df[df$P4_CESD != 'NA',])
Anova(lmerM1, type=3, test="F")
modelSummary(lmerM1)

lmerM2 = lmer(P4_CESD ~ coherence_slopeNoAge_C + age_C + (1|M2FAMNUM), data=df[df$COPEem_C != 'NA',])
Anova(lmerM2, type=3, test="F")
modelSummary(lmerM2)

lmerM3 = lmer(P4_CESD ~ coherence_slopeNoAge_C + COPEem_C + age_C + (1|M2FAMNUM), data=df)
Anova(lmerM3, type=3, test="F")
modelSummary(lmerM3)
# coherence: b = -1.75, F(1,874.98) = 12.278, p = .0005
# emotion coping: b = .39, F(1, 874.76) = 68.27, p < .0001

med = mediate(lmerM2, lmerM3, treat = "coherence_slopeNoAge_C", mediator = "COPEem_C")
summary(med)




#### mLMEM. Anxiety ####
lmerM1 = lmer(COPEem ~ coherence_slopeNoAge_C + age_C + (1|M2FAMNUM), data=df[df$P4_STAItrait != 'NA',])
Anova(lmerM1, type=3, test="F")
modelSummary(lmerM1)

lmerM2 = lmer(P4_STAItrait ~ coherence_slopeNoAge_C + age_C + (1|M2FAMNUM), data=df[df$COPEem != 'NA',])
Anova(lmerM2, type=3, test="F")
modelSummary(lmerM2)

lmerM3 = lmer(P4_STAItrait ~ coherence_slopeNoAge_C + COPEem_C + age_C + (1|M2FAMNUM), data=df)
Anova(lmerM3, type=3, test="F")
modelSummary(lmerM3)
# coherence: b = -1.66, F(1, 874.60) = 9.51, p = .002
# emotion coping: b = .69, F(1, 874.86) = 175.65, p < .0001

# Mediation
med = mediate(lmerM2, lmerM3, treat = "coherence_slopeNoAge_C", mediator = "COPEem_C")
summary(med)

#### mLMEM. PWB #### 
lmerM1 = lmer(COPEem ~ coherence_slopeNoAge_C + age_C + (1|M2FAMNUM), data=df[df$pwb2 != 'NA',])
Anova(lmerM1, type=3, test="F")
modelSummary(lmerM1)

# 1. Coherence -> PWB
lmerM2 = lmer(pwb2 ~ coherence_slopeNoAge_C + age_C + (1|M2FAMNUM), data=df[df$COPEem != 'NA',])
Anova(lmerM2, type=3, test="F")
modelSummary(lmerM2)
hist(dfL$coherence_slope)
# 3. COPE + Coherence -> PWB
lmerM3 = lmer(pwb2 ~ coherence_slopeNoAge_C + COPEem_C + age_C + (1|M2FAMNUM), data=df)
Anova(lmerM3, type=3, test="F")
modelSummary(lmerM3)
# coherence: b = 5.31, F(1, 877.98) =  6.31, p = .012
# emotion coping: b = -2.51, F(1, 877.95) = 158.09, p < .0001

# Mediation
med = mediate(lmerM2, lmerM3, treat = "coherence_slopeNoAge_C", mediator = "COPEem_C")
summary(med)

#### mLMEM. IL 6 ####
lmerM1 = lmer(COPEem ~ coherence_slopeNoAge_C + age_C + (1|M2FAMNUM), data=df[df$B4BIL6 != 'NA',])
Anova(lmerM1, type=3, test="F")
modelSummary(lmerM1)

lmerM2 = lmer(IL6 ~ coherence_slopeNoAge_C + age_C + (1|M2FAMNUM), data=df[df$COPEem_C != 'NA',])
Anova(lmerM2, type=3, test="F")
modelSummary(lmerM2)

lmerM3 = lmer(IL6 ~ coherence_slopeNoAge_C + COPEem_C + age_C + (1|M2FAMNUM), data=df)
Anova(lmerM3, type=3, test="F")
modelSummary(lmerM3)
# coherence: b = -0.45, F(1, 857.34) =  6.02, p = .014
# emotion coping: b = .04, F(1, 853.81) = 5.82, p = .016

# Mediation
med = mediate(lmerM2, lmerM3, treat = "coherence_slope_C", mediator = "COPEem_C")
summary(med)










