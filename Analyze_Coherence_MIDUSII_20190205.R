# RStudio version 1.1.453
# R version 3.5 
# Author: Sasha L. Sommerfeldt 
# University of Wisconsin - Madison 
# October 2016 - November 2018 

# Clear the workspace
rm(list=ls())

### Directories
# Raw data files downloaded from http://www.icpsr.umich.edu/icpsrweb/ICPSR/studies  
# /29282 (biomarker/project 4) and /04652 (survey/project 1)  
# And from http://midus.colectica.org/ for MIDUS 2 Milwaukee subsample   
  
# Then processed through Prep_Coherence_MIDUSII.R script.

dir = '~/Desktop/UWMadison/MIDUS'
# Data directory
ddir = paste(dir, '/data', sep='')
# Analysis directory (to output plots)
adir = paste(dir, '/analysis', sep='')

setwd(ddir)


### Packages
library(data.table)
library(plyr)
library(stats)
library(car)
library(ggplot2)
library(multilevel)
library(lme4)
library(lmSupport)
library(AICcmodavg)
library(pbkrtest)
library(boot)
library(rmarkdown)
library(broom)
library(pander)
library(broom.mixed)

### Read in processed data files 
#Files generated in Prep_Coherence_MIDUSII.R script


today='20181124'

# Wide format
fnameW = paste("coh_",today,".csv",sep='')
fpathW = paste(ddir,"/",fnameW, sep='')

# Long format
fnameL = paste("cohLong_",today,".csv",sep='')
fpathL = paste(ddir,"/",fnameL, sep='')

# Read in processed data
df = read.csv(fpathW)
dfL = read.csv(fpathL)


# PREP
### Subset dataframe

# A condensed/subsetted dataframe for analysis - excluding the many survey/P1 people without biomarker/P4/coherence data

dfLs = dfL[!is.na(dfL$coherence_slope),]
length(unique(dfLs$M2ID)) # 1065

# Transform that subsetted version to wide format

dfLsW = reshape(dfLs, idvar = "M2ID", v.names=c('hr', 'stress', 'stress_CMC', 'ecgQ'), drop=c('X', 'stressMC'), timevar = "timepoint", direction = "wide")
names(dfLsW)
length(dfLsW$M2ID)

### Summary statistics and demographics

summary(dfLsW$gender) 
varDescribe(dfLsW$months_P1SAQ_to_P4) 

varDescribe(dfLsW$P4_age) 
varDescribe(dfLsW$P1_PIage) 
varDescribe(dfLsW$months_P1SAQ_to_P4) 
varDescribe(dfLsW$months_P1PI_to_P4)
varDescribe(dfLsW$months_P1cog_to_P4) 
varDescribe(dfLsW$pwb2)
varDescribe(dfLsW$P4_CESD)
varDescribe(dfLsW$P4_STAItrait)
varDescribe(dfLsW$IL6) 
varDescribe(dfLsW$CRP) 
varDescribe(dfLsW$COPE_denial)

summary(dfLsW$P1_race) # Asian = 3, black = 193, Native american or alaska native aleutian islander/eskimo = 14, other = 27, white = 825, DK = 1, REFUSED =1, NA=1

summary(dfLsW$SAMPLMAJ) 

### Siblings

# Twins
# Create a data frame with a list of family IDs and the number of times they occurred.
n_occur = data.frame(table(dfLsW$M2FAMNUM))
#n_occur[n_occur$Freq > 1,]
# Which M2FAMNUMs occurred more than once.
#dfLsW$M2ID[dfLsW$M2FAMNUM %in% n_occur$Var1[n_occur$Freq > 1]]

# Count up the number of subjects with each repeated family id, each subject's value for sibs column is how many people have same family id (including themself)
dfLsW$sibs = 1
for (family in n_occur$Var1[n_occur$Freq > 1] ) {
  #print(family)
  numsibs = n_occur$Freq[n_occur$Var1 == family ]
  dfLsW$sibs[dfLsW$M2FAMNUM == family] = numsibs
}
length(dfLsW$M2ID[dfLsW$SAMPLMAJ == '(03) TWIN' & dfLsW$sibs == 2 & dfLsW$ZYGCAT == 'MONOZYGOTIC']) # 128
length(dfLsW$M2ID[dfLsW$SAMPLMAJ == '(03) TWIN' & dfLsW$sibs == 2 & dfLsW$ZYGCAT == 'DIZYGOTIC - SAME SEX']) # 56
length(dfLsW$M2ID[dfLsW$SAMPLMAJ == '(03) TWIN' & dfLsW$sibs == 2 & dfLsW$ZYGCAT == 'DIZYGOTIC - DIFFERENT SEX']) # 46
length(dfLsW$M2ID[dfLsW$SAMPLMAJ == '(03) TWIN' & dfLsW$sibs == 2 & dfLsW$ZYGCAT == 'UNABLE TO DETERMING ZYGOSITY']) # 2
#View(dfLsW[dfLsW$SAMPLMAJ == '(03) TWIN' & dfLsW$sibs == 4,]) # fam of 4 is DZ same sex 

# 1 family with 2 twin pairs, 1 family with 3 siblings
#View(dfLsW[dfLsW$SAMPLMAJ != '(03) TWIN',])

# 3 siblings
length(dfLsW$M2ID[dfLsW$SAMPLMAJ != '(03) TWIN' & dfLsW$sibs == 3]) # N = 3 (1 family)
length(dfLsW$M2ID[dfLsW$SAMPLMAJ != '(03) TWIN' & dfLsW$sibs == 2]) # N = 8 (4 families)

varDescribe(unique(dfLsW$M2FAMNUM[dfLsW$SAMPLMAJ == '(02) SIBLING'])) # 5 unique family ids that are siblings, 1 pair of siblings 



## Prep variables in long format df
# - Have age for everyone (so don't need to recenter well-being variable based on who has age)
# - Stress is centered within cluster (centered around each subject's mean)
# - Thus: for each analysis, just need to re-center age based on who has that well-being variable (this is probably overkill, the mean changes very little, but it's done)

### Cluster mean center
dfLs$stress_CMC = dfLs$stress - ave(dfLs$stress, dfLs$M2ID, na.rm=T)
dfLs$hr_CMC = dfLs$hr - ave(dfLs$hr, dfLs$M2ID, na.rm=T)


### Mean Center
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age, na.rm=T)
# Self reports
dfLs$pwb2_C = dfLs$pwb2 - mean(dfLs$pwb2, na.rm=T)
dfLs$P4_CESD_C = dfLs$P4_CESD- mean(dfLs$P4_CESD, na.rm=T)
dfLs$P4_STAItrait_C = dfLs$P4_STAItrait - mean(dfLs$P4_STAItrait, na.rm=T)
dfLs$COPE_denial_C = dfLs$COPE_denial - mean(dfLs$COPE_denial, na.rm=T)
# Divide pwb, cesd, stai by 10 so SEs larger, interpretable
dfLs$pwb2_C_d10 = dfLs$pwb2_C/10.000000
dfLs$P4_CESD_C_d10 = dfLs$P4_CESD_C/10.000000
dfLs$P4_STAItrait_C_d10 = dfLs$P4_STAItrait_C/10.000000
  
# Inflammatory
dfLs$IL6_C = dfLs$IL6 - mean(dfLs$IL6, na.rm=T)
dfLs$CRP_C = dfLs$CRP - mean(dfLs$CRP, na.rm=T)
# PWB subscales
dfLs$autonomy2_C = dfLs$autonomy2 - mean(dfLs$autonomy2, na.rm=T)
dfLs$envMast2_C = dfLs$envMast2 - mean(dfLs$envMast2, na.rm=T)
dfLs$persGrow2_C = dfLs$persGrow2 - mean(dfLs$persGrow2, na.rm=T)
dfLs$posRela2_C = dfLs$posRela2 - mean(dfLs$posRela2, na.rm=T)
dfLs$purpLife2_C = dfLs$purpLife2 - mean(dfLs$purpLife2, na.rm=T)
dfLs$selfAcce2_C = dfLs$selfAcce2 - mean(dfLs$selfAcce2, na.rm=T)

# Wide data frame
dfLsW$P4_age_C = dfLsW$P4_age - mean(dfLsW$P4_age, na.rm=T)
# Self reports
dfLsW$pwb2_C = dfLsW$pwb2 - mean(dfLsW$pwb2, na.rm=T)
dfLsW$P4_CESD_C = dfLsW$P4_CESD- mean(dfLsW$P4_CESD, na.rm=T)
dfLsW$P4_STAItrait_C = dfLsW$P4_STAItrait - mean(dfLsW$P4_STAItrait, na.rm=T)
dfLsW$COPE_denial_C = dfLsW$COPE_denial - mean(dfLsW$COPE_denial, na.rm=T)
# Inflammatory
dfLsW$IL6_C = dfLsW$IL6 - mean(dfLsW$IL6, na.rm=T)
dfLsW$CRP_C = dfLsW$CRP - mean(dfLsW$CRP, na.rm=T)
# PWB subscales
dfLsW$autonomy2_C = dfLsW$autonomy2 - mean(dfLsW$autonomy2, na.rm=T)
dfLsW$envMast2_C = dfLsW$envMast2 - mean(dfLsW$envMast2, na.rm=T)
dfLsW$persGrow2_C = dfLsW$persGrow2 - mean(dfLsW$persGrow2, na.rm=T)
dfLsW$posRela2_C = dfLsW$posRela2 - mean(dfLsW$posRela2, na.rm=T)
dfLsW$purpLife2_C = dfLsW$purpLife2 - mean(dfLsW$purpLife2, na.rm=T)
dfLsW$selfAcce2_C = dfLsW$selfAcce2 - mean(dfLsW$selfAcce2, na.rm=T)


### Recode dichotomous
dfLs$gender_C = varRecode(dfLs$gender, c('(1) MALE', '(2) FEMALE'), c(-.5,.5))


### Log transform inflammatory markers for normal distribution 
dfLs$IL6_T = log2(dfLs$IL6)
dfLsW$IL6_T = log2(dfLsW$IL6)
hist(dfLs$IL6_T)
dfLs$CRP_T = log(dfLs$CRP, base=10)
dfLsW$CRP_T = log(dfLsW$CRP, base=10)
hist(dfLs$CRP_T)

dfLs$IL6_T_C = dfLs$IL6_T - mean(dfLs$IL6_T, na.rm=T)
dfLs$CRP_T_C = dfLs$CRP_T - mean(dfLs$CRP_T, na.rm=T)


# TESTS

## Stress-heart rate coherence associations

### Age
lmerM = lmer(hr ~ stress_CMC * P4_age_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Gender 
lmerM = lmer(hr ~ stress_CMC * gender_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### PWB
# Center age for subjects in this analysis 
varDescribe(dfLs$pwb2_C)
length(dfLs$P4_age[!is.na(dfLs$pwb2_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$pwb2_C)], na.rm=T)

# Run the test 
lmerM = lmer(hr ~ stress_CMC * pwb2_C_d10 + P4_age_C*stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Depression
# Center age for subjects in this analysis 
varDescribe(dfLs$P4_CESD_C)
length(dfLs$P4_age[!is.na(dfLs$P4_CESD_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$P4_CESD_C)], na.rm=T)

# Run the test
lmerM = lmer(hr ~ stress_CMC * P4_CESD_C_d10 + P4_age_C*stress_CMC + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Anxiety
# Center age for subjects in this analysis 
varDescribe(dfLs$P4_STAItrait_C)
length(dfLs$P4_age[!is.na(dfLs$P4_STAItrait_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$P4_STAItrait_C)], na.rm=T)

# Run the test
lmerM = lmer(hr ~ stress_CMC * P4_STAItrait_C_d10 + P4_age_C*stress_CMC + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### IL6
# Center age for subjects in this analysis 
varDescribe(dfLs$IL6_T_C)
length(dfLs$P4_age[!is.na(dfLs$IL6_T_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$IL6_T_C)], na.rm=T)

# Run the test
lmerM = lmer(hr ~ stress_CMC * IL6_T_C + P4_age_C*stress_CMC + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### CRP
# Center age for subjects in this analysis 
varDescribe(dfLs$CRP_T_C)
length(dfLs$P4_age[!is.na(dfLs$CRP_T_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$CRP_T_C)], na.rm=T)

# Run the test
lmerM = lmer(hr ~ stress_CMC * CRP_T_C + P4_age_C*stress_CMC + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Denial
# Center age for subjects in this analysis 
length(dfLs$P4_age[!is.na(dfLs$COPE_denial_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$COPE_denial_C)], na.rm=T)

# Run the test 
lmerM = lmer(hr ~ stress_CMC * COPE_denial_C + P4_age_C*stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = F)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   # Using pander() to view the created table, with 3 sig figs 
pander(glance_obj, digits = 3) 


## Multiple Comparisons Correction
# Holm-Bonferonni
## p value for each test of a well-being marker
p = c(2.99E-07, 2.06E-09, 1.70E-08, 2.91E-06, 0.00762, 6.18E-06)
## Holm-bonferonni
p.adjust(p, method= 'holm')
# 1.196e-06 1.236e-08 8.500e-08 8.730e-06 7.620e-03 1.236e-05



# Reactivity and Recovery
### Compute reactivity measures
# Stress reactivity
dfLsW$stressChange2to1 = dfLsW$stress.2 - dfLsW$stress.1
varDescribe(dfLsW$stressChange2to1)
dfLsW$stressChange4to1 = dfLsW$stress.4 - dfLsW$stress.1
varDescribe(dfLsW$stressChange4to1)
dfLsW$stressChangeStresstoBase = rowMeans(dfLsW[c('stressChange2to1', 'stressChange4to1')], na.rm=TRUE)
varDescribe(dfLsW$stressChangeStresstoBase) # mean = 2.6, sd = 1.75, min = -7.5, max = 8 

# Heart rate reactivity
dfLsW$hrChange2to1 = dfLsW$hr.2 - dfLsW$hr.1
varDescribe(dfLsW$hrChange2to1)
dfLsW$hrChange4to1 = dfLsW$hr.4 - dfLsW$hr.1
varDescribe(dfLsW$hrChange4to1)
dfLsW$hrChangeStresstoBase = rowMeans(dfLsW[c('hrChange2to1', 'hrChange4to1')], na.rm=TRUE)
varDescribe(dfLsW$hrChangeStresstoBase) # mean = 3.42, sd = 3.81, min = -7.1, max = 30.95



# Center reactivity
dfLsW$stressChangeStresstoBase_C = dfLsW$stressChangeStresstoBase - mean(dfLsW$stressChangeStresstoBase, na.rm=T)
dfLsW$hrChangeStresstoBase_C = dfLsW$hrChangeStresstoBase - mean(dfLsW$hrChangeStresstoBase, na.rm=T)



# Self-reported stress
dfLsW$stressChange3to2 = dfLsW$stress.3 - dfLsW$stress.2
varDescribe(dfLsW$stressChange3to2)
dfLsW$stressChange5to4 = dfLsW$stress.5 - dfLsW$stress.4
varDescribe(dfLsW$stressChange5to4)
dfLsW$stressChangeRecovtoStress = rowMeans(dfLsW[c('stressChange3to2', 'stressChange5to4')], na.rm=TRUE)
varDescribe(dfLsW$stressChangeRecovtoStress)
# center
dfLsW$stressChangeRecovtoStress_C = dfLsW$stressChangeRecovtoStress - mean(dfLsW$stressChangeRecovtoStress, na.rm=T)

# Heart rate
dfLsW$hrChange3to2 = dfLsW$hr.3 - dfLsW$hr.2
varDescribe(dfLsW$hrChange3to2)
dfLsW$hrChange5to4 = dfLsW$hr.5 - dfLsW$hr.4
varDescribe(dfLsW$hrChange5to4)
dfLsW$hrChangeRecovtoStress = rowMeans(dfLsW[c('hrChange3to2', 'hrChange5to4')], na.rm=TRUE)
varDescribe(dfLsW$hrChangeRecovtoStress)
# center
dfLsW$hrChangeRecovtoStress_C = dfLsW$hrChangeRecovtoStress - mean(dfLsW$hrChangeRecovtoStress, na.rm=T)


## Merge reactivity and recovery measures into dfLs
varsToMerge = c('M2ID', 'hrChangeStresstoBase', 'hrChangeStresstoBase_C', 'stressChangeStresstoBase', 'stressChangeStresstoBase_C', 'hrChangeRecovtoStress', 'hrChangeRecovtoStress_C', 'stressChangeRecovtoStress', 'stressChangeRecovtoStress_C')
#dfLsW[varsToMerge]
# dfLs with Reactivity and Recovery data = dfLsRR
dfLsRR = merge.data.frame(dfLs, dfLsW[varsToMerge],  by='M2ID', all=TRUE)
#varDescribe(dfLsRR)


# Center age
dfLsRR$P4_age_C = dfLsRR$P4_age - mean(dfLsRR$P4_age, na.rm=T)


## Is reactivity or recovery associated with coherence?
### Heart rate reactivity 
lmerM = lmer(hr ~ stress_CMC * hrChangeStresstoBase_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Stress reactivity
lmerM = lmer(hr ~ stress_CMC * stressChangeStresstoBase_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Heart rate recovery
lmerM = lmer(hr ~ stress_CMC * hrChangeRecovtoStress_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Stress recovery
lmerM = lmer(hr ~ stress_CMC * stressChangeRecovtoStress_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


## Is stress reactivity associated with heart rate reactivity?
lmerM = lmer(hrChangeStresstoBase_C ~ stressChangeStresstoBase_C + P4_age_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 



## Does coherence predict well-being outcomes when adjusting for reactivity?

### PWB + reactivity
lmerM = lmer(hr ~ stress_CMC * pwb2_C_d10 + P4_age_C * stress_CMC +stressChangeStresstoBase_C + hrChangeStresstoBase_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Depression + reactivity 
lmerM = lmer(hr ~ stress_CMC * P4_CESD_C_d10 + P4_age_C * stress_CMC + stressChangeStresstoBase_C + hrChangeStresstoBase_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Anxiety + reactivity
# P4_STAItrait
lmerM = lmer(hr ~ stress_CMC * P4_STAItrait_C_d10 + P4_age_C * stress_CMC + stressChangeStresstoBase_C + hrChangeStresstoBase_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### IL6 + reactivity
lmerM = lmer(hr ~ stress_CMC * IL6_T_C + P4_age_C * stress_CMC + stressChangeStresstoBase_C + hrChangeStresstoBase_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### CRP + reactivity
lmerM = lmer(hr ~ stress_CMC * CRP_T_C + P4_age_C * stress_CMC + stressChangeStresstoBase_C + hrChangeStresstoBase_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Denial + reactivity
lmerM = lmer(hr ~ stress_CMC * COPE_denial_C + P4_age_C * stress_CMC + stressChangeStresstoBase_C + hrChangeStresstoBase_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


## Does coherence predict well-being outcomes when adjusting for recovery?
### PWB + recovery
lmerM = lmer(hr ~ stress_CMC * pwb2_C + P4_age_C * stress_CMC + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Depression + recovery 
# CESD
lmerM = lmer(hr ~ stress_CMC * P4_CESD_C + P4_age_C * stress_CMC + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Anxiety + recovery
# P4_STAItrait
lmerM = lmer(hr ~ stress_CMC * P4_STAItrait_C + P4_age_C * stress_CMC + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### IL6 + recovery
lmerM = lmer(hr ~ stress_CMC * IL6_T_C + P4_age_C * stress_CMC + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### CRP + recovery
lmerM = lmer(hr ~ stress_CMC * CRP_T_C + P4_age_C * stress_CMC + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### denial + recovery
lmerM = lmer(hr ~ stress_CMC * COPE_denial_C + P4_age_C * stress_CMC + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 



## Does coherence predict well-being outcomes when adjusting for reactivity and recovery?
### PWB + reactivity + recovery
lmerM = lmer(hr ~ stress_CMC * pwb2_C_d10 + P4_age_C * stress_CMC + stressChangeStresstoBase_C + hrChangeStresstoBase_C + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Depression + reactivity + recovery
# CESD
lmerM = lmer(hr ~ stress_CMC * P4_CESD_C_d10 + P4_age_C * stress_CMC + stressChangeStresstoBase_C + hrChangeStresstoBase_C + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Anxiety + reactivity + recovery
# P4_STAItrait
lmerM = lmer(hr ~ stress_CMC * P4_STAItrait_C_d10 + P4_age_C * stress_CMC + stressChangeStresstoBase_C + hrChangeStresstoBase_C + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### IL6 + reactivity + recovery
lmerM = lmer(hr ~ stress_CMC * IL6_T_C + P4_age_C +  P4_age_C * stress_CMC + stressChangeStresstoBase_C + hrChangeStresstoBase_C + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 

### CRP + reactivity + recovery
lmerM = lmer(hr ~ stress_CMC * CRP_T_C + P4_age_C + P4_age_C * stress_CMC + stressChangeStresstoBase_C + hrChangeStresstoBase_C + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Denial + reactivity + recovery
lmerM = lmer(hr ~ stress_CMC * COPE_denial_C + P4_age_C + P4_age_C * stress_CMC + stressChangeStresstoBase_C + hrChangeStresstoBase_C + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLsRR)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsRR) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 




## Does reactivity and/or recovery predict well-being outcomes?
### PWB ~ reactivity + recovery
lmerM = lmer(pwb2 ~ P4_age_C + stressChangeStresstoBase_C + hrChangeStresstoBase_C + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsW) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Depression ~ reactivity + recovery
# CESD
lmerM = lmer(P4_CESD ~ P4_age_C + stressChangeStresstoBase_C + hrChangeStresstoBase_C + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsW) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 

### Anxiety ~ reactivity + recovery
# P4_STAItrait
lmerM = lmer(P4_STAItrait ~ P4_age_C + stressChangeStresstoBase_C + hrChangeStresstoBase_C + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsW) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### IL6 ~ reactivity + recovery
lmerM = lmer(IL6_T ~ P4_age_C + stressChangeStresstoBase_C + hrChangeStresstoBase_C + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsW) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### CRP ~ reactivity + recovery
lmerM = lmer(CRP_T ~ P4_age_C + stressChangeStresstoBase_C + hrChangeStresstoBase_C + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsW) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 

### Denial ~ reactivity + recovery
lmerM = lmer(COPE_denial ~ P4_age_C + stressChangeStresstoBase_C + hrChangeStresstoBase_C + stressChangeRecovtoStress_C + hrChangeRecovtoStress_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLsW) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 



# PLOT

## FIGURE 1: Stress and heart rate by phase histograms

# Facet-wrapped histograms of stress and heart rate at each phase of stress induction

ylimits = c(0, 610)
colcode = "#4ECDC1"
stressHist=ggplot()+
geom_histogram(data=dfL, aes(stress), fill=colcode, binwidth=1, color="black") +
facet_wrap(~timepoint, ncol=5) +
labs(x="Self-Reported Stress", y="Number of Subjects") +
ylim(ylimits)+
scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10) )+
theme_bw(base_size=10)+
theme(axis.text=element_text(size=9, family=""), text=element_text(family="Helvetica", size=10),
panel.grid.minor=element_blank(), 
panel.background=element_rect(fill="transparent"), 
plot.background=element_rect(fill="transparent") )
stressHist
#ggsave(stressHist, filename=paste(adir,"/stressHist.png", sep=''), bg="transparent", height=2.8, width=10.45, units="in")

# Heart rate
xlimits = c(30,150)
ylimits = c(0, 240)
colcode = "#900C3F"
hrHist=ggplot()+
geom_histogram(data=dfL, aes(hr), fill=colcode, binwidth=6, color="black") +
facet_wrap(~timepoint, ncol=5) +
labs(x="Heart Rate", y="Number of Subjects") +
ylim(ylimits)+
xlim(xlimits)+
geom_vline(xintercept=75, size=.5, color="blue")+
#scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10) )+
  theme_bw(base_size=10)+
  theme(axis.text=element_text(size=9, family=""), text=element_text(family="Helvetica", size=10),
panel.grid.minor= element_blank(), 
panel.background=element_rect(fill="transparent"), 
plot.background=element_rect(fill="transparent") )
hrHist
#ggsave(hrHist, filename=paste(adir,"/hrHist.png", sep=''), bg="transparent", height=2.8, width=10.45, units="in")


## FIGURE 2: Interaction plots

# PWB
mod = lmer(hr ~ stress*pwb2 + (1 + stress| M2ID), data=dfLs)
# Prepare independent variables for ggplot
XToPredict = seq(min(dfLs$stress), max(dfLs$stress), length = 100)
pwb2_lo = mean(dfLsW$pwb2, na.rm=T) - sd(dfLsW$pwb2, na.rm=T) 
pwb2_hi = mean(dfLsW$pwb2, na.rm=T) + sd(dfLsW$pwb2, na.rm=T) 

# Use modelPredictions() to generate Y-hats
yHats = expand.grid(stress = XToPredict, pwb2=c(pwb2_lo, pwb2_hi))  # all IVs
yHats = modelPredictions(mod, yHats)

modelplot = ggplot() + 
            geom_smooth(aes(ymin = CILo, ymax = CIHi, x = stress, y = Predicted, 
                  colour=as.factor(pwb2), group=as.factor(pwb2)), 
                  data = yHats, stat = "identity")
#modelplot

pwb2plot = modelplot + scale_x_continuous("Self-Reported Stress", breaks = seq(0, 10, by=1)) + 
             scale_y_continuous("Heart Rate", breaks = seq(60, 100, by=1)) +
             scale_color_manual(name ="PWB",             
                     labels=c("Low ( - 1SD)", "High ( + 1SD)"), values=c("#0679A4","#FDA603")) + 
             theme_bw(base_size = 10) +
             theme(legend.position = c(0.75, 0.15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=9), text=element_text(family="Helvetica", size=10))

pwb2plot

## CESD
mod = lmer(hr ~ stress*P4_CESD + (1 + stress| M2ID), data=dfLs)
# Prepare independent variables for ggplot
XToPredict = seq(min(dfLs$stress), max(dfLs$stress), length = 100)
P4_CESD_lo = mean(dfLsW$P4_CESD, na.rm=T) - sd(dfLsW$P4_CESD, na.rm=T) 
P4_CESD_hi = mean(dfLsW$P4_CESD, na.rm=T) + sd(dfLsW$P4_CESD, na.rm=T) 

# Use modelPredictions() to generate Y-hats
yHats = expand.grid(stress = XToPredict, P4_CESD=c(P4_CESD_lo, P4_CESD_hi))  # all IVs
yHats = modelPredictions(mod, yHats)

modelplot = ggplot() + 
  geom_smooth(aes(ymin = CILo, ymax = CIHi, x = stress, y = Predicted, 
                  colour=as.factor(P4_CESD), group=as.factor(P4_CESD)), 
                  data = yHats, stat = "identity")
#modelplot

P4_CESDplot = modelplot + scale_x_continuous("Self-Reported Stress", breaks = seq(0, 10, by=1)) + 
             scale_y_continuous("Heart Rate", breaks = seq(60, 100, by=1)) +
             scale_color_manual(name ="Depression",             
                     labels=c("Low ( - 1SD)", "High ( + 1SD)"), values=c("#FDA603","#0679A4")) + 
             theme_bw(base_size = 10) +
             theme(legend.position = c(0.75, 0.15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=9), text=element_text(family="Helvetica", size=10))

P4_CESDplot



## P4_STAItrait
mod = lmer(hr ~ stress*P4_STAItrait + (1 + stress| M2ID), data=dfLs)
# Prepare independent variables for ggplot
XToPredict = seq(min(dfLs$stress), max(dfLs$stress), length = 100)
P4_STAItrait_lo = mean(dfLsW$P4_STAItrait, na.rm=T) - sd(dfLsW$P4_STAItrait, na.rm=T) 
P4_STAItrait_hi = mean(dfLsW$P4_STAItrait, na.rm=T) + sd(dfLsW$P4_STAItrait, na.rm=T) 

# Use modelPredictions() to generate Y-hats
yHats = expand.grid(stress = XToPredict, P4_STAItrait=c(P4_STAItrait_lo, P4_STAItrait_hi))  # all IVs
yHats = modelPredictions(mod, yHats)

modelplot = ggplot() + 
  geom_smooth(aes(ymin = CILo, ymax = CIHi, x = stress, y = Predicted, 
                  colour=as.factor(P4_STAItrait), group=as.factor(P4_STAItrait)), 
                  data = yHats, stat = "identity")
#modelplot

P4_STAItraitplot = modelplot + scale_x_continuous("Self-Reported Stress", breaks = seq(0, 10, by=1)) + 
                   scale_y_continuous("Heart Rate", breaks = seq(60, 100, by=1)) +
                   scale_color_manual(name ="Trait Anxiety",             
                     labels=c("Low ( - 1SD)", "High ( + 1SD)"), values=c("#FDA603","#0679A4")) + 
                   theme_bw(base_size = 10) +
                   theme(legend.position = c(0.75, 0.15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=9), text=element_text(family="Helvetica", size=10))

P4_STAItraitplot



## IL6
mod = lmer(hr ~ stress*IL6_T + (1 + stress| M2ID), data=dfLs)
# Prepare independent variables for ggplot
XToPredict = seq(min(dfLs$stress), max(dfLs$stress), length = 100)
IL6_T_lo = mean(dfLsW$IL6_T, na.rm=T) - sd(dfLsW$IL6_T, na.rm=T) 
IL6_T_hi = mean(dfLsW$IL6_T, na.rm=T) + sd(dfLsW$IL6_T, na.rm=T) 

# Use modelPredictions() to generate Y-hats
yHats = expand.grid(stress = XToPredict, IL6_T=c(IL6_T_lo, IL6_T_hi))  # all IVs
yHats = modelPredictions(mod, yHats)

# Starting plot in which we graph regression lines
modelplot = ggplot() + 
  geom_smooth(aes(ymin = CILo, ymax = CIHi, x = stress, y = Predicted, 
                  colour=as.factor(IL6_T), group=as.factor(IL6_T), lineIL6_T=as.factor(IL6_T)), 
                  data = yHats, stat = "identity")
#modelplot

IL6_Tplot = modelplot + scale_x_continuous("Self-Reported Stress", breaks = seq(0, 10, by=1)) + 
             scale_y_continuous("Heart Rate", breaks = seq(60, 100, by=1)) +
             scale_color_manual(name ="IL-6 (log2)",             
                     labels=c("Low ( - 1SD)", "High ( + 1SD)"), values=c("#FDA603","#0679A4")) + 
             theme_bw(base_size = 10) +
             theme(legend.position = c(0.75, 0.15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=9), text=element_text(family="Helvetica", size=10))

IL6_Tplot



## CRP
mod = lmer(hr ~ stress*CRP_T + (1 + stress| M2ID), data=dfLs)
# Prepare independent variables for ggplot
XToPredict = seq(min(dfLs$stress), max(dfLs$stress), length = 100)
CRP_T_lo = mean(dfLsW$CRP_T, na.rm=T) - sd(dfLsW$CRP_T, na.rm=T) 
CRP_T_hi = mean(dfLsW$CRP_T, na.rm=T) + sd(dfLsW$CRP_T, na.rm=T) 

# Use modelPredictions() to generate Y-hats
yHats = expand.grid(stress = XToPredict, CRP_T=c(CRP_T_lo, CRP_T_hi))  # all IVs
yHats = modelPredictions(mod, yHats)

modelplot = ggplot() + 
  geom_smooth(aes(ymin = CILo, ymax = CIHi, x = stress, y = Predicted, 
                  colour=as.factor(CRP_T), group=as.factor(CRP_T), lineCRP_T=as.factor(CRP_T)), 
              data = yHats, stat = "identity")
#modelplot

CRP_Tplot = modelplot + scale_x_continuous("Self-Reported Stress", breaks = seq(0, 10, by=1)) + 
            scale_y_continuous("Heart Rate", breaks = seq(60, 100, by=1)) +
            scale_color_manual(name ="CRP (log10)",             
                     labels=c("Low ( - 1SD)", "High ( + 1SD)"), values=c("#FDA603","#0679A4")) + 
            theme_bw(base_size = 10) +
            theme(legend.position = c(0.75, 0.15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=9), text=element_text(family="Helvetica", size=10))

CRP_Tplot

## Denial
mod = lmer(hr ~ stress*COPE_denial + (1 + stress| M2ID), data=dfLs)
# Prepare independent variables for ggplot
XToPredict = seq(min(dfLs$stress), max(dfLs$stress), length = 100)
COPE_denial_lo = mean(dfLsW$COPE_denial, na.rm=T) - sd(dfLsW$COPE_denial, na.rm=T) 
COPE_denial_hi = mean(dfLsW$COPE_denial, na.rm=T) + sd(dfLsW$COPE_denial, na.rm=T) 

# Use modelPredictions() to generate Y-hats
yHats = expand.grid(stress = XToPredict, COPE_denial=c(COPE_denial_lo, COPE_denial_hi))  # all IVs
yHats = modelPredictions(mod, yHats)

modelplot = ggplot() + 
  geom_smooth(aes(ymin = CILo, ymax = CIHi, x = stress, y = Predicted, 
                  colour=as.factor(COPE_denial), group=as.factor(COPE_denial), lineCOPE_denial=as.factor(COPE_denial)), 
                  data = yHats, stat = "identity")
#modelplot

COPE_denialplot = modelplot + scale_x_continuous("Self-Reported Stress", breaks = seq(0, 10, by=1)) + 
             scale_y_continuous("Heart Rate", breaks = seq(60, 100, by=1)) +
             scale_color_manual(name ="Denial Coping",             
                     labels=c("Low ( - 1SD)", "High ( + 1SD)"), values=c("#FDA603","#0679A4")) + 
             theme_bw(base_size = 10) +
             theme(legend.position = c(0.75, 0.15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=9), text=element_text(family="Helvetica", size=10))

COPE_denialplot



## FIGURE 3: Plot individual subject slopes
dfL$stressMC = dfL$stress - ave(dfL$stress, dfL$M2ID)
dfL$hrM = ave(dfL$stress, dfL$M2ID)

ggplot(dfL, aes(stress, hr, color=as.factor(M2ID)))+
geom_smooth(aes(group=as.factor(M2ID)),method="lm",se=F,size=.2, alpha=.6, position="jitter")+
xlim(c(0,11))+
theme_bw() +
theme(panel.grid.minor = element_blank(), axis.text=element_text(size=14), axis.title=element_text(size=24)) +
labs(x="Self-Reported Stress", y="Heart Rate")+
theme(legend.position="none")


# With heat map where color is the magnitude of the slope 

mycol = c("#0710C4", "gray",# negative & zero
          "#FFEC00", "#FFC300", "#FF5733", "#C70039", "#900C3F", "#581845")   # positive
mybreaks = c(-.4, 0, 
             .5, 1, 1.5, 2, 3, 4)

ggplot(dfL, aes(stress, hr, color=coherence_slope))+
  geom_smooth(aes(group=as.factor(M2ID)),method="lm",se=F,size=.2, alpha=.6, position="jitter")+
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=9,family="Helvetica"), panel.border = element_blank(), text=element_text(family="Helvetica", size=10)) +
  labs(y="Heart Rate")+
  scale_colour_gradientn("",colours=mycol, limits=c(-.4, 4.5), values = scales::rescale(c(-0.5, -0.05, 0, 0.05, 0.5,1,2,3,4)), breaks = mybreaks, guide="colourbar")+
  scale_x_continuous("Self-Reported Stress", breaks = seq(0, 10, by=1))


## FIGURE 3: Histogram of BLUPS
ggplot(dfLsW, aes(coherence_slope)) +
  geom_histogram(aes(fill=as.factor(coherence_slope)), binwidth=.2, col="black",         fill="#FFC300") +
  #scale_fill_gradientn("Slope", colours=mycol, limits=c(-.4, 4.5), values = scales::rescale(c(-0.5, -0.05, 0, 0.05, 0.5,1,2,3,4)), breaks = mybreaks, guide="colourbar")+
  labs(x="EBLUP", y="Number of Subjects") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=9, family="Helvetica"), legend.position="none", panel.border = element_blank(), text=element_text(family="Helvetica", size=10))



# SUPPLEMENTAL 

# I. Correlation (r) as coherence
# See Prep_Coherence_MIDUSII.R for correlation computation. There, each subject's set of heart rate and stress measures are subset to their own data frame and a correlation is computed. The resulting within-subject (i.e., single-subject) r's compose a new variable in the main dataframe. 
### Center correlations variable
varDescribe(dfLsW$coherence_as_r) # .49(.47) median.66 skew = -1.18, kurtosis = .55

# Center age for subjects in this analysis 
dfLsW$P4_age_C = dfLsW$P4_age - mean(dfLsW$P4_age[!is.na(dfLsW$coherence_as_r)], na.rm=T)
# Center correlations
dfLsW$coherence_as_r_C = dfLsW$coherence_as_r - mean(dfLsW$coherence_as_r, na.rm=T)


### PWB ~ coherence as r
# Run the test 
lmerM = lmer(pwb2 ~ coherence_as_r_C + P4_age_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = F)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   # Using pander() to view the created table, with 3 sig figs 
pander(glance_obj, digits = 3) 


### Depression ~ coherence as r 
# Run the test 
lmerM = lmer(P4_CESD ~ coherence_as_r_C + P4_age_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = F)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   # Using pander() to view the created table, with 3 sig figs 
pander(glance_obj, digits = 3) 


### Anxiety ~ coherence as r 
# Run the test 
lmerM = lmer(P4_STAItrait ~ coherence_as_r_C + P4_age_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = F)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   # Using pander() to view the created table, with 3 sig figs 
pander(glance_obj, digits = 3) 


### IL6 ~ coherence as r 
# Run the test 
lmerM = lmer(IL6_T ~ coherence_as_r_C + P4_age_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = F)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   # Using pander() to view the created table, with 3 sig figs 
pander(glance_obj, digits = 3) 


### CRP ~ coherence as r 
# Run the test 
lmerM = lmer(CRP_T ~ coherence_as_r_C + P4_age_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = F)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   # Using pander() to view the created table, with 3 sig figs 
pander(glance_obj, digits = 3) 


### Denial ~ coherence as r 
# Run the test 
lmerM = lmer(COPE_denial ~ coherence_as_r_C + P4_age_C + (1|M2FAMNUM), data=dfLsW)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = F)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   # Using pander() to view the created table, with 3 sig figs 
pander(glance_obj, digits = 3) 


## Multiple Comparisons Correction
# Holm-Bonferonni
## p value for each test of a well-being marker/denial
p = c(2.71E-06, 2.78E-07, 4.05E-07, 3.02E-04, 0.172, 3.19E-05)
## Holm-bonferonni
p.adjust(p, method= 'holm')
# 1.084e-05 1.668e-06 2.025e-06 6.040e-04 1.720e-01 9.570e-05



## FIGURE S1: Correlations histogram
ggplot(dfLsW, aes(coherence_as_r)) +
geom_histogram(aes(fill=as.factor(coherence_as_r)), binwidth=.2, col="black", fill="#FFC300") +
#scale_fill_gradientn("Slope", colours=mycol, limits=c(-.4, 4.5), values = scales::rescale(c(-0.5, -0.05, 0, 0.05, 0.5,1,2,3,4)), breaks = mybreaks, guide="colourbar")+
labs(x="Correlation", y="Number of Subjects") +
theme_bw() +
theme(panel.grid.minor = element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=24)) +
theme(legend.position="none")



# II. Lag from Survey to Biomarker substudies
There was a lag of 0-60 months from the survey to the stress-induction (biomarker) substudies. The COPE and PWB were completed as part of the Survey substudy. All other measures were collected as part of the stress-induction substudy. 

### PWB + lag
# Center age for subjects in this analysis 
length(dfLs$P4_age[!is.na(dfLs$pwb2_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$pwb2_C)], na.rm=T)
# Center lag for subjects in this analysis 
length(dfLs$months_P1SAQ_to_P4[!is.na(dfLs$pwb2_C)])
dfLs$months_P1SAQ_to_P4_C = dfLs$months_P1SAQ_to_P4 - mean(dfLs$months_P1SAQ_to_P4[!is.na(dfLs$pwb2_C)], na.rm=T)

# Lag moderate?
lmerM = lmer(hr ~ stress_CMC * pwb2_C_d10 * months_P1SAQ_to_P4_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 

# Adjust for lag
lmerM = lmer(hr ~ stress_CMC * pwb2_C + months_P1SAQ_to_P4_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Denial + lag
# Center age for subjects in this analysis 
length(dfLs$P4_age[!is.na(dfLs$COPE_denial_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$COPE_denial_C)], na.rm=T)
# Center lag for subjects in this analysis 
length(dfLs$months_P1SAQ_to_P4[!is.na(dfLs$COPE_denial_C)])
dfLs$months_P1SAQ_to_P4_C = dfLs$months_P1SAQ_to_P4 - mean(dfLs$months_P1SAQ_to_P4[!is.na(dfLs$COPE_denial_C)], na.rm=T)

# Lag moderate?
lmerM = lmer(hr ~ stress_CMC * COPE_denial_C * months_P1SAQ_to_P4_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 

# Adjust for lag
lmerM = lmer(hr ~ stress_CMC * COPE_denial_C + months_P1SAQ_to_P4_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 



# III. PWB subscales
# Exploratory analyses investigating individual subscales of the Psychological Well-Being Scales

# Center age for subjects in this analysis 
varDescribe(dfLs$pwb2_C)
length(dfLs$P4_age[!is.na(dfLs$pwb2_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$pwb2_C)], na.rm=T)


### Autonomy 
lmerM = lmer(hr ~ stress_CMC * autonomy2_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Environmental Mastery
lmerM = lmer(hr ~ stress_CMC * envMast2_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Personal Growth
lmerM = lmer(hr ~ stress_CMC * persGrow2_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Positive Relations with Others
lmerM = lmer(hr ~ stress_CMC * posRela2_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Purpose in Life
lmerM = lmer(hr ~ stress_CMC * purpLife2_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Self Acceptance
lmerM = lmer(hr ~ stress_CMC * selfAcce2_C + P4_age_C * stress_CMC + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 




# Non-linear Age  
#Including age^2 in our model did not impact results. 



# Center age (been centered for subsets of participants on different analyses where participants are missing data on well-being indicators)
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age, na.rm=T)

dfLs$P4_age_C2 = dfLs$P4_age_C^2


### HR ~ age^2 
lmerM = lmer(hr ~ P4_age_C + P4_age_C2 + (1|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### PWB + age^2
# Center age for subjects in this analysis 
varDescribe(dfLs$pwb2_C)
length(dfLs$P4_age[!is.na(dfLs$pwb2_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$pwb2_C)], na.rm=T)
dfLs$P4_age_C2 = dfLs$P4_age_C^2
# Run the test 
lmerM = lmer(hr ~ stress_CMC * pwb2_C + P4_age_C + P4_age_C2 + (1+ stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### Depression + age^2
# Center age for subjects in this analysis 
varDescribe(dfLs$P4_CESD_C)
length(dfLs$P4_age[!is.na(dfLs$P4_CESD_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$P4_CESD_C)], na.rm=T)
dfLs$P4_age_C2 = dfLs$P4_age_C^2
# Run the test
lmerM = lmer(hr ~ stress_CMC * P4_CESD_C + P4_age_C + P4_age_C2 + (1 + hr_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 



### Anxiety + age^2
# Center age for subjects in this analysis 
varDescribe(dfLs$P4_STAItrait_C)
length(dfLs$P4_age[!is.na(dfLs$P4_STAItrait_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$P4_STAItrait_C)], na.rm=T)
dfLs$P4_age_C2 = dfLs$P4_age_C^2
# Run the test
lmerM = lmer(hr ~ stress_CMC * P4_STAItrait_C + P4_age_C + P4_age_C2 + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### IL6 + age^2
# Center age for subjects in this analysis 
varDescribe(dfLs$IL6_T_C)
length(dfLs$P4_age[!is.na(dfLs$IL6_T_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$IL6_T_C)], na.rm=T)
dfLs$P4_age_C2 = dfLs$P4_age_C^2
# Run the test
lmerM = lmer(hr ~ stress_CMC * IL6_T_C + P4_age_C + P4_age_C2 + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


### CRP + age^2
# Center age for subjects in this analysis 
varDescribe(dfLs$CRP_T_C)
length(dfLs$P4_age[!is.na(dfLs$CRP_T_C)])
dfLs$P4_age_C = dfLs$P4_age - mean(dfLs$P4_age[!is.na(dfLs$CRP_T_C)], na.rm=T)
dfLs$P4_age_C2 = dfLs$P4_age_C^2
# Run the test
lmerM = lmer(hr ~ stress_CMC * CRP_T_C + P4_age_C + P4_age_C2 + (1 + stress_CMC|M2ID) + (1|M2FAMNUM), data=dfLs)
#Anova(lmerM, type=3, test="F")
modelSummary(lmerM, t = FALSE)
table_obj = broom.mixed::tidy(lmerM, conf.int=TRUE, conf.level=.95, conf.method="Wald", effects = c("ran_pars", "fixed"), data=dfLs) 
glance_obj = broom.mixed::glance(lmerM)
pander(table_obj, digits = 3)   
pander(glance_obj, digits = 3) 


