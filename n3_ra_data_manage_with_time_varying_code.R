#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Title: Data management of 9Health PUFA data
# Author: Ryan Gan
# Date: 12/4/15
# Date Modified: 12/22/15
# R Version: 3.2.3
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# installing packages ----------------------------------------------------------
# use dplyr and tidyr package
library(dplyr) # for some data management
library(tidyr)
library(data.table) # data.table package has a nice transpose funciton

# setting directories and importing data ---------------------------------------
# setting absolute path when working off google drive
# commenting out either mac or pc depending on where I'm working
 dir <- paste0('C:/Users/RGan/Google Drive/9_Health_PUFA/r_analyses/',
              'n3_ra_9health') # this directory is for my work pc
# mac
# dir <- paste0('/Users/ryangan/Google Drive/9_Health_PUFA/r_analyses/',
#              'n3_ra_9health') # this directory is for my mac


setwd(dir) # set home directory
getwd()
list.files()

# read in horizontal csv files 
# pc path
 long <- paste0('C:/Users/RGan/Google Drive/9_Health_PUFA/r_analyses/',
               'n3_ra_9health/n3_ra_9health_jan2016.csv')
# mac path
# long <- paste0('/Users/ryangan/Google Drive/9_Health_PUFA/r_analyses/',
#               'n3_ra_9health/n3_ra_9health_jan2016.csv')
ninehealth_long <- read.csv(long)

# check the dataset
glimpse(ninehealth_long)
summary(ninehealth_long)

# data management --------------------------------------------------------------

summary(ninehealth_long$DOVisit)
#dovisit is as.factor, I want to convert it to date
ninehealth_long$DOVisit <- as.Date(ninehealth_long$DOVisit, "%d%b%Y")
# dateRADX as.Date
ninehealth_long$DateRADX <- as.Date(ninehealth_long$DateRADX, "%d%b%Y")

# thinking of subsetting to ccp3.1 or ccp3, xtabs to check coding
xtabs(~ninehealth_long$CCP31 + ninehealth_long$CCP31Outcome) # looks good
xtabs(~ninehealth_long$CCP3 + ninehealth_long$CCP3Outcome)

first_vis <- subset(ninehealth_long, Visit_ == 1)
# sampcode doesn't seem to be as useful as it was before; avoid it for now
xtabs(~first_vis$CCP31 ) # only 48 at risk folks
xtabs(~ SampCode, ninehealth_long)

# rename visit number and subset to the population defined by Kevin
pufa_9hlth_long <- rename(ninehealth_long, visit_num = Visit_ ) %>% # rename vis
  filter(IN_Study == 1) # subset to CCP3 or 3.1 pos and IA/RA free

# looks like something is there; can I make a dataset that shifts the outcome 1 row?
# long dataset that contains the outcome observation at the subsequent visit
glimpse(pufa_9hlth_long)

# Jan 18th Edits: I think it would be useful to reduce the dataset, but for now I'm
# going to move on to analytics ----------------------------------------------------
# I should try the alt rectangle copy for efficiently picking out the variables from
# glimpse

xtabs(~ rfneph + CCP31Outcome, pufa_9hlth_long) # checking the antibody pos neg vars
# make a test dataframe with select number of variables; easier to view
pufa_9hlth_clean <- select(pufa_9hlth_long, ID, LABID, DOVisit, visit_num, IN_Study,
                           Age, Gender, Race, BMI, BMI_Impute, Education, IncomeLevel, 
                           RADX, DateRADX, DaysToRADX, ConvertedDuringSERA, CCP2Final, 
                           CCP2, CCP2RefRange, CCP3Outcome, CCP3, CCP31Outcome, CCP31, 
                           RFfinalNum, years_contributed
                           ) %>%
                  mutate(py10 = ifelse(packyears_new >= 10, 1, 0))

# lead command takes the next observation while lag takes the previous observation
test_lead <- group_by(test_df, ID) %>%
             mutate(next_swell_jt = lead(SwollenWristMCP, order_by = visit_num),
                    next_il6 = lead(IL_6, order_by = visit_num),
                    next_tnfa = lead(TNF_a, order_by = visit_num)) 

head(test_lead) # looks like the lead command worked 

# making IA/RA outcome ---------------------------------------------------------
# reduce dataset to CCP3 or CCP3.1 positive
xtabs( ~ CCP3 + CCP31, pufa_9hlth_long) # check who is ccp3 pos
xtabs( ~ IA + IA_base, pufa_9hlth_long)
xtabs( ~ RA + RA_base, pufa_9hlth_long)
# make an IA/RA outcome (IA and RA has some missing, check with elizabeth)
pufa_9hlth_long <- mutate(pufa_9hlth_long, ia_new = ifelse(is.na(IA), 0, IA),
                          ra_new = ifelse(is.na(RA), 0, RA),
                          ia_ra = ifelse(ia_new == 1 | ra_new == 1, 1, 0),
                          py10 = ifelse(Packyears_new >= 10, 1, 0),
                          n3_supp = ifelse(Omega3 == 'Yes', 1, 
                                           ifelse(Omega3 == 'No', 0, NA)))
# checks for data
xtabs( ~ ia_ra + visit_num, pufa_9hlth_long) # 11 people with baseline IA/RA
xtabs( ~ ia_new + ra_new + ia_ra, pufa_9hlth_long) # ia/ra var looks good
xtabs( ~ Omega3, pufa_9hlth_long)
xtabs( ~ Omega3 + n3_supp, pufa_9hlth_long) # check numeric omega3 supp variable
summary(as.factor(pufa_9hlth_long$n3_supp))
summary(as.factor(pufa_9hlth_long$ia_ra))
summary(as.factor(pufa_9hlth_long$IA))

# incident ia_ra dataset -------------------------------------------------------
# I want to see how the data for IA/RA outcome is structured
# subsetting to select variables
glimpse(pufa_9hlth_long)

pufa_9hlth_long <- arrange(pufa_9hlth_long, ID, DOVisit)

# variables I'd like to vary by time
time_vary <- select(pufa_9hlth_long, ID, visit_num, DOVisit, years_contributed, 
                    DateRADX, ia_ra, IA, RA, n3_supp, Omega3_bi, Supp, total_n3, ala, 
                    epa, dpa, dha, epa_dha, total_n3, la, gla, ara, Omega3_6_ratio, 
                    CRPfinalNum, BMI_Impute, CCP2, CCP2_cutoff, RFIgM, SwollenJoint, 
                    TotalSwollenJointCount, TotalTenderJointCount,IN_Study) %>%
                    filter(!is.na(total_n3)) %>% # removing missing n-3 visits
                    group_by(ID) %>% # below are times for time vary dataset
                    mutate(time_bw_vis = ifelse(visit_num == 1, 0, 
                                         (DOVisit - lag(DOVisit, order_by = visit_num))),
                           time_contrib = cumsum(time_bw_vis),
                           start_time = time_contrib,
                           stop_time = lead(time_contrib, order_by = visit_num),
                           lead_ia_ra = lead(ia_ra, order_by = visit_num),
                           lead_tot_n3 = lead(total_n3, order_by = visit_num),
                           tot_n3_chg = lead(total_n3, order_by = visit_num) - total_n3)

# variables I don't want to vary that I will merge back in to time_vary
non_vary_base <- filter(pufa_9hlth_long, visit_num == 1) %>%
                 select(ID, Age, Sex, Race_bi, Educ, Inc, Cursmoke_Impute, EverSmoke_Impute, 
                        Packyears_new, py10, SE_num, SEcount, AnySupplements, 
                        Multivits, Antioxidant, AnyOmega, VitaminD)

# merge non-varying base variables with time varying by id
time_vary2 <- full_join(time_vary, non_vary_base, by = 'ID') %>%
              # time varying age
              group_by(ID) %>%
              mutate(age2 = Age + (time_contrib/365),
                     start_age = age2,
                     stop_age = lead(age2, order_by = visit_num))
                     
# find who had baseline ia_ra
base_ia_ra <- filter(time_vary2, ia_ra == 1 & visit_num == 1) %>%
                  mutate(base_outcome = 1) %>% # indicates outcome at base
                  select(ID, base_outcome) # create an indicator of base out
                  
incident_ia_ra <- anti_join(time_vary2, base_ia_ra, by = 'ID') # select IDs without base ia_ra

xtabs(~ ia_ra + visit_num, incident_ia_ra) # check no base ia_ra

# time varying dataset (need to think about fixed variables like age)
time_vary_ia <- ungroup(incident_ia_ra) %>% # remove the group functions
                filter(!is.na(stop_time)) %>% # remove last observation for everyone with no stop_time
                group_by(ID) %>%
                mutate(visit_convert = ifelse(ia_ra == 1, visit_num, 0)) %>%
                select(ID, visit_convert) %>%
                filter(visit_convert > 0) %>% # merge back in to big dataset
                right_join(incident_ia_ra, by = 'ID') %>%
                mutate(keep = ifelse(is.na(visit_convert) | # creating a keep var
                                       visit_num < visit_convert, 1, 0)) %>%
                filter(!is.na(stop_time) & keep == 1)



# data set: first visit for IA at baseline, IA incident, and no IA -------------

ia_groups_9hlth <- filter(time_vary, ia_ra == 1 & visit_num != 1) %>%
                   group_by(ID) %>%
                   summarise(vist_pos = n()) %>%
                   mutate(incident_outcome = 1) %>%
                   select(ID, incident_outcome) %>%
                   full_join(base_ia_ra, by = 'ID') %>%
                   mutate(ia_group = ifelse(is.na(base_outcome), 2,
                                     ifelse(base_outcome == 1, 1, 0))) %>%
                   select(ID, ia_group) %>%
                   right_join(pufa_9hlth_long, by = 'ID') %>%
                   mutate(ia_cat = ifelse(is.na(ia_group), 0,
                                      ifelse(ia_group == 1, 1,
                                      ifelse(ia_group == 2, 2, NA))),
                          base_ia = ifelse(ia_cat == 1, 1, 0)) %>%
                   select(-ia_group) # remove ia_group

summary(as.factor(ia_groups_9hlth$ia_cat))

# Creating a single observation dataset 
ia_groups_9hlth_first_vis <- filter(ia_groups_9hlth, visit_num == 1)

# what predicts baseline conversion? -------------------------------------------
# I will eventually move this to another section
xtabs(~ base_ia, ia_groups_9hlth_first_vis)

summary(glm(base_ia ~ Sex + Omega3_bi + as.factor(SEcount) + py10 + total_n3, family = 'binomial', data = ia_groups_9hlth_first_vis))

# creation of a wide dataset ---------------------------------------------------
print(tbl_df(pufa_9hlth_long))
# dplyr and tidyr are not as good at big transposes; tidyr dcast seems much better
glimpse(pufa_9hlth_long)

# create a vector of variables I want to transpose by ID and visit_num
transpose_vars <- c("race_bi", "age", "sex", "SE_num", "SEcount", "bmi_new", "ccp2", 
                    "ccp2_cutoff", "ccp31", "ccp31_cutoff","rfigm", "rfigm_cutoff",
                    "rfiga", "rfiga_cutoff", "rfigg", "rfigg_cutoff", "CRPfinalnum",
                    "educ", "Inc", "packyears_new", "supp", "multivits_bi", "omega3_bi", 
                    "vitamind_bi", "antioxidant_bi", "total_n3", "ala", "epa", 
                    "dpa","dha", "f1perc", "RADX", "dateRADX", "SwollenWristMCP", 
                    "TenderWristMCP")

# data.table package might be better for a big transpose
wide <- dcast(setDT(pufa_9hlth_long), ID ~ visit_num, value.var = transpose_vars) 

glimpse(wide)

xtabs(~omega3_bi_1 + SwollenWristMCP_2 +  SwollenWristMCP_1 , wide)

xtabs(~omega3_bi_2 + SwollenWristMCP_3 +  SwollenWristMCP_2 , wide)





