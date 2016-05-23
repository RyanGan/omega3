#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Title: Comparison of IA at base, Incident IA, and No IA in 9hlth CCP3 pos
# Author: Ryan Gan
# Date: 1/29/15
# R Version: 3.2.2
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Note 5/13/16: I'm rerunning these just to make sure numbers are correct

# loading package --------------------------------------------------------------
# library(foreign)
library(ggplot2) 
# library(MASS) # stats package
library(dplyr)
library(tidyr)
library(nnet) # multinomial logistic regression
library(contrast) # linear contrast statements
library(aod) # tests type III fixed effects
library(haven) # haven package for SAS survival file
library(survival) # survival model
 
# import csv dataset of baseline dataframe -------------------------------------
read_path <- paste0('C:/Users/RGan/Google Drive/9_Health_PUFA/CSV Data Sets/',
                     'ia_baseline_analysis_df.csv')

ia_groups_9hlth_first_vis <- read.csv(read_path)

# first visit analyses ---------------------------------------------------------
# use ia_group_9hlth_first_vis dataframe created in data manage script 

# what predicts those who are not IA_RA at base or who do not convert?
ia_analysis_df <- mutate(ia_groups_9hlth_first_vis, ia_never = ifelse(ia_cat == 0, 1, 0),
                  ia_ever = ifelse(ia_cat != 0, 1, 0), # ia outcome ever
                  select_study = ifelse(ia_cat != 1, 1, 0),
                  bmi_cat = ifelse(BMI_Impute <= 25, 0, 
                            ifelse(BMI_Impute >25 & BMI_Impute <= 30, 1, 2)),
                  overweight_obese = ifelse(BMI_Impute <= 25, 0, 1),
                  # id 09-011-00 is missing n3 supplement use, setting to 0 since
                  # other visits are 0
                  n3_supp_impute = ifelse(is.na(Omega3_bi), 0, Omega3_bi),
                  # id 09-008-00 to incident ia group for sen. analysis
                  ia_cat2 = ifelse(ID=='09-008-00', 2, ia_cat),
                  py10 = ifelse(Packyears_new > 10, 1, 0),
                  age50 = ifelse(Age >= 50, 1, 0))


# check groups
xtabs(~ia_cat + ia_ever, ia_analysis_df)
xtabs(~ia_cat + ia_cat2, ia_analysis_df)
# set as factor
ia_analysis_df$ia_cat <- as.factor(ia_analysis_df$ia_cat)

str(ia_analysis_df$ia_cat)

# table 1 descriptives by ia status at baseline --------------------------------
# incident IA who are not positive at baseline are included in the no ia group

# age 
mod <- lm(Age ~ as.factor(IA_base), data = ia_analysis_df)
summary(mod)
anova(mod)
# check distributions of age
mu <- group_by(ia_analysis_df, ia_cat) %>% 
  summarise(mean_age = mean(Age), std_age = sd(Age), med_age = median(Age))
mu

ggplot(ia_analysis_df, aes(x=Age, color=ia_cat)) + 
  geom_density() + 
  geom_vline(data=mu, aes(xintercept=mean_age, color=ia_cat), linetype='dashed')

# parametric tests with age probably okay

# sex
tab <- xtabs(~ IA_base + Sex, data=ia_analysis_df)
tab
prop.table(tab) # cell %
prop.table(tab, 1) # row %
prop.table(tab, 2) # column %
summary(tab)
fisher.test(tab)

# race binary (nhw vs other)
tab <- xtabs(~IA_base + Race_bi, data=ia_analysis_df)
tab
#prop.table(tab) # cell %
prop.table(tab, 1) # row %
#prop.table(tab, 2) # column %
summary(tab)
fisher.test(tab)

# Education >HS (this is probably okay based on looking at education below)
tab <- xtabs(~IA_base + Educ, data=ia_analysis_df)
tab
# i want row percents the way the table is set up
prop.table(tab, 1) # row %
summary(tab)
fisher.test(tab)

# look at education in more detail
tab <- xtabs(~IA_base + Education, data=ia_analysis_df)
tab
prop.table(tab) # cell %
prop.table(tab, 1) # row %
prop.table(tab, 2) # column %
fisher.test(tab)

# Income (looks okay to treat as binary)
tab <- xtabs(~IA_base + Inc, data=ia_analysis_df)
tab
# i want row percents the way the table is set up
prop.table(tab, 1) # row %
summary(tab)
fisher.test(tab)

# BMI (with some imputed)
mod <- lm(BMI_Impute ~ as.factor(IA_base), data=ia_analysis_df)
summary(mod)
anova(mod)

# checking distribution
mu <- group_by(ia_analysis_df, IA_base) %>% 
  summarise(mean_bmi = mean(BMI_Impute), std_bmi = sd(BMI_Impute), med_bmi = median(BMI_Impute))
mu

ggplot(ia_analysis_df, aes(x=BMI_Impute, color=IA_base)) + 
  geom_density() + 
  geom_vline(data=mu, aes(xintercept=mean_bmi, color=IA_base), linetype='dashed')

# check BMI category
tab <- xtabs(~IA_base + Eversmoke, data=ia_analysis_df)
tab
prop.table(tab, 1) # row %
summary(tab)
fisher.test(tab)

# Ever Smoke (impute)
tab <- xtabs(~IA_base + EverSmoke_Impute, data=ia_analysis_df)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)

# Current smoke (impute)
tab <- xtabs(~IA_base + Cursmoke_Impute, data=ia_analysis_df)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)

# packeyears continous
mod <- lm(PackYears~ as.factor(IA_base), data=ia_analysis_df)
summary(mod)
anova(mod)

xtabs(~ PackYears + Packyears_new, ia_analysis_df)

# checking distribution
mu <- group_by(ia_groups_9hlth_first_vis, IA_base) %>% 
  summarise(mean_py = mean(PackYears), med_py = median(PackYears))
mu

ggplot(ia_analysis_df, aes(x=PackYears, color=IA_base)) + 
  geom_density()

# py10
tab <- xtabs(~IA_base + py10, data=ia_analysis_df)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)

# SE pos
tab <- xtabs(~IA_base + SE_num, data=ia_analysis_df)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)
# SE count
tab <- xtabs(~IA_base + SEcount, data=ia_analysis_df)
tab
prop.table(tab)
prop.table(tab,1)
summary(tab)
fisher.test(tab)

# CCP2 positivity
tab <- xtabs(~IA_base + CCP2, data=ia_analysis_df)
tab
prop.table(tab,1)
fisher.test(tab)

# RF IgM Positivity
tab <- xtabs(~IA_base + RFIgM, data=ia_analysis_df)
tab
prop.table(tab,1)
fisher.test(tab)

# RF IgA Pos
tab <- xtabs(~IA_base + RFIgA, data=ia_analysis_df)
tab
prop.table(tab,1)
fisher.test(tab)

# RF IgG Pos
tab <- xtabs(~IA_base + RFIgG, data=ia_analysis_df)
tab
prop.table(tab,1)
fisher.test(tab)
# could look at antibody distributions too

# swollen joint count
tab <- xtabs(~ IA_base + SwollenWrstMCP, data=ia_analysis_df)
tab
prop.table(tab,1)
fisher.test(tab)

# omega3 supplement use
tab <- xtabs(~ IA_base + Omega3_bi, data=ia_analysis_df)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)
summary(ia_analysis_df)

xtabs(~Omega3_bi + n3_supp_impute, ia_analysis_df)
# omega-3 imput for 09-011-00
tab <- xtabs(~ IA_base + n3_supp_impute, data=ia_analysis_df)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)
summary(ia_analysis_df)

# omega6 supplement use
tab <- xtabs(~ IA_base + Omega6_bi, data=ia_analysis_df)
tab
prop.table(tab,1)
fisher.test(tab)

# vit d supplement use
tab <- xtabs(~ IA_base + VitaminD_bi, data=ia_analysis_df)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)

# anti-oxidant(Antioxidant_bi)
tab <- xtabs(~ IA_base + Antioxidant_bi, data=ia_analysis_df)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)

# Multivits_bi
tab <- xtabs(~ IA_base + Multivits_bi, data=ia_analysis_df)
tab
prop.table(tab,1)
fisher.test(tab)
summary(tab)

# swollen joints
tab <- xtabs(~ IA_base + RASpecficTenderJtCnt_num, data=ia_analysis_df)
tab
prop.table(tab,1)
fisher.test(tab)
summary(tab)

# association with covariates and baseline omegas ------------------------------
# age
n3_mod <- lm(total_n3 ~ Age, data = ia_analysis_df)
summary(n3_mod) # Age is associated

# sex
n3_mod <- lm(total_n3 ~ Sex, data = ia_analysis_df)
summary(n3_mod) # sex is not associated (although femals have lower levels)

# race binary
n3_mod <- lm(total_n3 ~ Race_bi, data = ia_analysis_df)
summary(n3_mod) # race not associated

# education (>HS)
n3_mod <- lm(total_n3 ~ Educ, data = ia_analysis_df)
summary(n3_mod) # ed not associated

# income 
n3_mod <- lm(total_n3 ~ Inc, data = ia_analysis_df)
summary(n3_mod) # income might be associated

# bmi impute 
n3_mod <- lm(total_n3 ~ BMI_Impute, data = ia_analysis_df)
summary(n3_mod) # not really an association with bmi 

# bmi category 
n3_mod <- lm(total_n3 ~ as.factor(bmi_cat), data = ia_analysis_df)
summary(n3_mod) # bmi category associated where both overweight and obese
# subjects have lower total omega fatty acids
# bmi overweight and obese
n3_mod <- lm(total_n3 ~ overweight_obese, data = ia_analysis_df)
summary(n3_mod) # significantly associated

# ever smoke imputed
n3_mod <- lm(total_n3 ~ EverSmoke_Impute, data = ia_analysis_df)
summary(n3_mod) 

# current smoker
n3_mod <- lm(total_n3 ~ Cursmoke_Impute, data = ia_analysis_df)
summary(n3_mod) # significantly associated where smokers have lower levels

# continous packyears
n3_mod <- lm(total_n3 ~ PackYears, data = ia_analysis_df)
summary(n3_mod) # not associated

# py10
n3_mod <- lm(total_n3 ~ py10, data = ia_analysis_df)
summary(n3_mod)

# shared epitope pos
n3_mod <- lm(total_n3 ~ SE_num, data = ia_analysis_df)
summary(n3_mod)
# SE count
n3_mod <- lm(total_n3 ~ as.factor(SEcount), data = ia_analysis_df)
summary(n3_mod) # 2 SE copies have lower omega levels

# CCP2, this might be good to look at excluding base pos folks
n3_mod <- lm(total_n3 ~ CCP2 + as.factor(ia_cat), data = ia_analysis_df)
summary(n3_mod)

# RFIgM
n3_mod <- lm(total_n3 ~ RFIgM, data = ia_analysis_df)
summary(n3_mod)

# omega 3 distirbuionts --------------------------------------------------------
# distribution of omegas by ia cat
mod <- lm(total_n3 ~ as.factor(ia_cat), data=ia_analysis_df)
summary(mod)
anova(mod)
# checking distribution
mu <- group_by(ia_analysis_df, ia_cat) %>% 
  summarise(mean_n3 = mean(total_n3), std_n3 = sd(total_n3), med_n3 = median(total_n3))
mu

ggplot(ia_analysis_df, aes(x=total_n3, color=ia_cat)) + 
  geom_density() + 
  geom_vline(data=mu, aes(xintercept=mean_n3, color=ia_cat), linetype='dashed')

summary(ia_analysis_df$total_n3)

# see if there is a trend in omegas, split by tertile first
tert <- quantile(ia_analysis_df$total_n3, c(0, 1/3, 2/3, 3/3))
tert

# baseline tertile cut points
ia_analysis_df <- mutate(ia_analysis_df, 
                         tot_n3_cat = ifelse(total_n3 > 4.5 & total_n3 <= 6.45, 0,
                                      ifelse(total_n3 > 6.45 & total_n3 <= 8.03, 1,
                                      ifelse(total_n3 > 8.03 & total_n3 <= 13.1, 2, NA))),
                         n3_level_cut = ifelse(total_n3 > 6.8, 1, 0), # cut based on incident analysis non converters
                         # interaction between standardized omegas and SE
                         stnd_totn3_se = stnd_total_n3 * SE_num,
                         stnd_ala_se = stnd_ala * SE_num,
                         stnd_epa_se = stnd_epa * SE_num,
                         stnd_dpa_se = stnd_dpa * SE_num,
                         stnd_dha_se = stnd_dha * SE_num,
                         stnd_totn6_se = stnd_total_n6 * SE_num,
                         stnd_la_se = stnd_la * SE_num,
                         stnd_gla_se = stnd_gla * SE_num,
                         stnd_ara_se = stnd_ara * SE_num,
                         n3_n6_ratio_se = Omega3_6_ratio * SE_num,
                        # dummy coding the n3 categories
                         n3_tert1 = ifelse(tot_n3_cat == 0, 1, 0),
                         n3_tert2 = ifelse(tot_n3_cat == 1, 1, 0),
                         n3_tert3 = ifelse(tot_n3_cat == 2, 1, 0),
                         se_n3tert1 = n3_tert1 * SE_num,
                         se_n3tert2 = n3_tert2 * SE_num,
                         se_n3tert3 = n3_tert3 * SE_num,
                         ia_se_cat = ifelse(ia_ever == 1 & SE_num == 1, 3,
                                     ifelse(ia_ever == 1 & SE_num == 0, 2,
                                     ifelse(ia_ever == 0 & SE_num == 1, 1, 0))))

# check dummy coding of n3 tertiles                    
xtabs(~tot_n3_cat + ia_ever + SE_num, ia_analysis_df)
xtabs(~ia_se_cat + ia_ever + SE_num, ia_analysis_df)

ia_analysis_df$ia_se_cat <- as.factor(ia_analysis_df$ia_se_cat)
glimpse(ia_analysis_df)

# Table 2: association with increasing n-3 biomarker and IA at baseline --------
glimpse(ia_analysis_df)

# total n3 standardized
base_mod <- glm(IA_base ~ stnd_total_n3 + py10 + Sex + n3_supp_impute, 
                family='binomial', data=ia_analysis_df)
summary(base_mod)
exp(cbind(OR = coef(base_mod), confint.default(base_mod)))
# ala
base_mod <- glm(IA_base ~ stnd_ala + py10 + Sex + n3_supp_impute, 
                family='binomial', data=ia_analysis_df)
summary(base_mod)
exp(cbind(OR = coef(base_mod), confint.default(base_mod)))
# epa
base_mod <- glm(IA_base ~ stnd_epa + py10 + Sex + n3_supp_impute, 
                family='binomial', data=ia_analysis_df)

summary(base_mod)
exp(cbind(OR = coef(base_mod), confint.default(base_mod)))
exp(confint(base_mod))
anova(base_mod, test = 'LRT')

# dpa
base_mod <- glm(IA_base ~ stnd_dpa + py10 + Sex + n3_supp_impute, 
                family='binomial', data=ia_analysis_df)
summary(base_mod)
exp(cbind(OR = coef(base_mod), confint.default(base_mod)))
anova(base_mod, test = 'LRT')
# dha
base_mod <- glm(IA_base ~ stnd_dha + py10 + Sex + n3_supp_impute, 
                family='binomial', data=ia_analysis_df)
summary(base_mod)
exp(cbind(OR = coef(base_mod), confint.default(base_mod)))
# epa + dha
base_mod <- glm(IA_base ~ stnd_epa_dha + py10 + Sex + n3_supp_impute, 
                family='binomial'(link='logit'), data=ia_analysis_df)
summary(base_mod)
anova(base_mod, test='Chisq')
anova(base_mod, test = 'LRT')

# this confint uses profiled log likelihood
exp(cbind(OR = coef(base_mod), confint(base_mod)))

# CIs using standard errors
exp(confint.default(base_mod))

?anova.glm

# Table 3: descriptive characteristics for incident IA vs IA-free --------------

incident_ia <- ia_analysis_df %>% filter(ia_cat == 2 | ia_cat == 0)

xtabs(~ IA_ever, incident_ia)

# age 
mod <- lm(Age ~ as.factor(IA_ever), data = incident_ia)
summary(mod)
anova(mod)
# check distributions of age
mu <- group_by(incident_ia, ia_cat) %>% 
  summarise(n = n(), mean_age = mean(Age), std_age = sd(Age), med_age = median(Age))
mu

ggplot(incident_ia, aes(x=Age, group=ia_cat, color=ia_cat)) + 
  geom_density() + 
  geom_vline(data=mu, aes(xintercept=mean_age, color=ia_cat), linetype='dashed')

# Age > 50
tab <- xtabs(~IA_ever + age50, data=incident_ia)
tab
#prop.table(tab) # cell %
prop.table(tab, 1) # row %
#prop.table(tab, 2) # column %
summary(tab)
fisher.test(tab)

# parametric tests with age probably okay

# sex
tab <- xtabs(~ IA_ever + Sex, data=incident_ia)
tab
prop.table(tab) # cell %
prop.table(tab, 1) # row %
prop.table(tab, 2) # column %
summary(tab)
fisher.test(tab)

# race binary (nhw vs other)
tab <- xtabs(~IA_ever + Race_bi, data=incident_ia)
tab
#prop.table(tab) # cell %
prop.table(tab, 1) # row %
#prop.table(tab, 2) # column %
summary(tab)
fisher.test(tab)

# Education >HS (this is probably okay based on looking at education below)
tab <- xtabs(~IA_ever + Educ, data=incident_ia)
tab
# i want row percents the way the table is set up
prop.table(tab, 1) # row %
summary(tab)
fisher.test(tab)

# look at education in more detail
tab <- xtabs(~IA_ever + Education, data=incident_ia)
tab
prop.table(tab) # cell %
prop.table(tab, 1) # row %
prop.table(tab, 2) # column %
fisher.test(tab)

# Income (looks okay to treat as binary)
tab <- xtabs(~IA_ever + Inc, data=incident_ia)
tab
# i want row percents the way the table is set up
prop.table(tab, 1) # row %
summary(tab)
fisher.test(tab)

# BMI (with some imputed)
mod <- lm(BMI_Impute ~ as.factor(IA_ever), data=incident_ia)
summary(mod)
anova(mod)

# checking distribution
mu <- group_by(incident_ia, IA_ever) %>% 
  summarise(n = n(), mean_bmi = mean(BMI_Impute), std_bmi = sd(BMI_Impute), 
            med_bmi = median(BMI_Impute))
mu

ggplot(incident_ia, aes(x=BMI_Impute, color=IA_ever)) + 
  geom_density() + 
  geom_vline(data=mu, aes(xintercept=mean_bmi, color=IA_ever), linetype='dashed')

# check BMI category
tab <- xtabs(~IA_ever + Eversmoke, data=incident_ia)
tab
prop.table(tab, 1) # row %
summary(tab)
fisher.test(tab)

# Ever Smoke (impute)
tab <- xtabs(~IA_ever + EverSmoke_Impute, data=incident_ia)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)

# Current smoke (impute)
tab <- xtabs(~IA_ever + Cursmoke_Impute, data=incident_ia)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)

# packeyears continous
mod <- lm(PackYears~ as.factor(IA_ever), data=incident_ia)
summary(mod)
anova(mod)

xtabs(~ PackYears + Packyears_new, incident_ia)

# checking distribution
mu <- group_by(incident_ia, IA_ever) %>% 
  summarise(n = n(), mean_py = mean(PackYears), med_py = median(PackYears))
mu

ggplot(incident_ia, aes(x=PackYears, color=IA_ever)) + 
  geom_density()

# py10 (looks like I'm missing py data on one control)
tab <- xtabs(~IA_ever + py10, data=incident_ia)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)

# SE pos
tab <- xtabs(~IA_ever + SE_num, data=incident_ia)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)
# SE count
tab <- xtabs(~IA_ever + SEcount, data=incident_ia)
tab
prop.table(tab)
prop.table(tab,1)
summary(tab)
fisher.test(tab)

# CCP2 positivity
tab <- xtabs(~IA_ever + CCP2, data=incident_ia)
tab
prop.table(tab,1)
fisher.test(tab)

# RF IgM Positivity
tab <- xtabs(~IA_ever + RFIgM, data=incident_ia)
tab
prop.table(tab,1)
fisher.test(tab)

# RF IgA Pos
tab <- xtabs(~IA_ever + RFIgA, data=incident_ia)
tab
prop.table(tab,1)
fisher.test(tab)

# RF IgG Pos
tab <- xtabs(~IA_ever + RFIgG, data=incident_ia)
tab
prop.table(tab,1)
fisher.test(tab)
# could look at antibody distributions too

# swollen joint count
tab <- xtabs(~ IA_ever + SwollenWrstMCP, data=incident_ia)
tab
prop.table(tab,1)
fisher.test(tab)

# omega3 supplement use
tab <- xtabs(~ IA_ever + n3_supp_impute, data = incident_ia)
tab
prop.table(tab,1)
fisher.test(tab)

# omega6 supplement use
tab <- xtabs(~ IA_ever + Omega6_bi, data=incident_ia)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)

# vit d supplement use
tab <- xtabs(~ IA_ever + VitaminD_bi, data=incident_ia)
tab
prop.table(tab,1)
fisher.test(tab)
summary(tab)

# anti-oxidant(Antioxidant_bi)
tab <- xtabs(~ IA_ever + Antioxidant_bi, data=incident_ia)
tab
prop.table(tab,1)
fisher.test(tab)

# Multivits_bi
tab <- xtabs(~ IA_ever + Multivits_bi, data=incident_ia)
tab
prop.table(tab,1)
summary(tab)
fisher.test(tab)


# swollen joints
tab <- xtabs(~ IA_ever + RASpecficTenderJtCnt_num, data=incident_ia)
tab
prop.table(tab,1)
fisher.test(tab)
summary(tab)

# Tabel 4: Incident time-varying survival model --------------------------------

# note: survival dataset needs some QCing
# import survival dataset
tv_path <- paste0('C:/Users/RGan/Google Drive/9_Health_PUFA/SAS Data Sets/',
                  'pufasurvival.sas7bdat')

tv_ia <- read_sas(tv_path)

# add some variables 
tv_ia <- tv_ia %>% mutate(age50 = ifelse(Age >= 50, 1, 0))

# se stratum
se_pos_ia <- tv_ia %>% filter(SE_num == 1)
se_neg_ia <- tv_ia %>% filter(SE_num == 0)

glimpse(tv_ia)
xtabs(~IA + RA, tv_ia)

# total_n3
surv_mod <- coxph(Surv(T1, T2, IA== 1) ~ sd_total_n3 + age50 + SE_num + Omega3_bi, 
              data=tv_ia)
summary(surv_mod)


# there are 5 missing values of n3 fatty acids, should not be in dataset


surv_mod <- coxph(Surv(T1, T2, IA== 1) ~ nextvisit_sd_total_n3 + age50 + SE_num + 
                    Omega3_bi, data=tv_ia)
summary(surv_mod)

# se pos stratum
surv_mod <- coxph(Surv(T1, T2, IA== 1) ~ sd_total_n3 + Omega3_bi, 
              data=se_pos_ia )
summary(surv_mod)

# se neg stratum
surv_mod <- coxph(Surv(T1, T2, IA== 1) ~ sd_total_n3 + Omega3_bi, 
              data=se_neg_ia)
summary(surv_mod)

# ala
surv_mod <- coxph(Surv(T1, T2, IA== 1) ~ sd_ala + age50 + SE_num + Omega3_bi, 
              data=tv_ia)
summary(surv_mod)

