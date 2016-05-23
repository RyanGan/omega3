#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Title: Data management of 9Health PUFA data
# Author: Ryan Gan
# Date: 12/4/15
# Date Modified: 12/22/15
# R Version: 3.2.3
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Note: 5/6/2016

# use dplyr and tidyr package
library(dplyr) # for some data management
library(tidyr)
# library(data.table) # not sure I want data.table package just now

# setting absolute path when working off google drive
# commenting out either mac or pc depending on where I'm working
 dir <- paste0('C:/Users/RGan/Google Drive/9_Health_PUFA/r_analyses/',
              'n3_ra_9health') # this directory is for my work pc

# dir <- paste0(/Users/ryangan/Google Drive/9_Health_PUFA/r_analyses/',
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
#                'n3_ra_9health/n3_ra_9health_jan2016.csv')

ninehealth_long <- read.csv(long)

# check the dataset
glimpse(ninehealth_long)
summary(ninehealth_long)

# data management --------------------------------------------------------------

# dovisit is as.factor, I want to convert it to date
ninehealth_long$DOVISIT <- as.Date(ninehealth_long$DOVisit, "%d%b%Y")
# dateRADX as.Date
ninehealth_long$DateRADX <- as.Date(ninehealth_long$DateRADX, "%d%b%Y")

# thinking of subsetting to ccp3.1 or ccp3, xtabs to check coding
xtabs(~ninehealth_long$CCP31 + ninehealth_long$CCP31Outcome) # looks good
xtabs(~ninehealth_long$CCP3 + ninehealth_long$CCP3Outcome)

first_vis <- subset(ninehealth_long, Visit_ == 1)
xtabs(~first_vis$SampCode, ninehealth_long) # only 48 at risk folks
summary(as.factor(first_vis$CCP31))

# rename visit number
pufa_9hlth_long <- rename(ninehealth_long, visit_num = Visit_) %>% # rename vis
                   filter(IN_Study == 1) # filter to just those who should be in study

# looks like something is there; can I make a dataset that shifts the outcome 1 row?
# long dataset that contains the outcome observation at the subsequent visit
glimpse(pufa_9hlth_long)

# check with IN_Study variable
xtabs(~IN_Study + CCP3, pufa_9hlth_long)

# make a dataset of the baseline visit 

ia_groups_9hlth_first_vis <- pufa_9hlth_long %>% filter(visit_num == 1) %>%
                             mutate(ia_cat = ifelse(IA_ever == 0, 0,
                                             ifelse(IA_base == 1, 1, 
                                      ifelse(IA_base == 0 & IA_ever == 1, 2, NA))))

xtabs(~ IA_base + IA_ever + ia_cat, ia_groups_9hlth_first_vis)

# writing permanent dataset to drive
write_path <- paste0('C:/Users/RGan/Google Drive/9_Health_PUFA/CSV Data Sets/',
                     'ia_baseline_analysis_df.csv')

write.csv(ia_groups_9hlth_first_vis, file = write_path)


# make a test dataframe with select number of variables; easier to view
test_df <- select(pufa_9hlth_long, ID, DOVISIT, visit_num, omega3_bi, total_n3,
                  ala, epa, dpa, dha, SwollenWristMCP, age, sex, race_bi, 
                  bmi_new, educ, Inc, eversmoke, cursmoke, ccp2, rfigm, rfiga, 
                  rfigg, SEcount, SE_num, packyears_new, IL_6, TNF_a) %>%
                  mutate(py10 = ifelse(packyears_new >= 10, 1, 0))

# lead command takes the next observation while lag takes the previous observation
test_lead <- group_by(test_df, ID) %>%
             mutate(next_swell_jt = lead(SwollenWristMCP, order_by = visit_num),
                    next_il6 = lead(IL_6, order_by = visit_num),
                    next_tnfa = lead(TNF_a, order_by = visit_num)) 

head(test_lead) # looks like the lead command worked 


  

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






