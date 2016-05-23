#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Title: Time-varying survival analysis
# Author: Ryan Gan
# Date: 1/20/16
# Date Modified: 1/20/16
# R Version: 3.2.3
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# libraries --------------------------------------------------------------------
library(survival)
library(ggplot2)

# use time_vary_ia dataset
mod1 <- coxph(Surv(start_time, stop_time, lead_ia_ra == 1) ~ total_n3 , 
              data=time_vary_ia)
summary(mod1)

# make interaction term between total_n3 and SE_num
time_vary_ia$se_n3_int <- time_vary_ia$total_n3 * time_vary_ia$SE_num 
time_vary_ia$se_lead_n3 <- time_vary_ia$lead_tot_n3 * time_vary_ia$SE_num
time_vary_ia$se_n3_chg <- time_vary_ia$tot_n3_chg * time_vary_ia$SE_num
time_vary_ia$se_dpa <- time_vary_ia$dpa * time_vary_ia$SE_num
time_vary_ia$Omega3_bi[time_vary_ia$ID == '09-011-00'] <- 0
time_vary_ia$Cursmoke_Impute[time_vary_ia$ID == '09-011-00'] <- 0

# interaction between CCP2 and n3
time_vary_ia$ccp2_n3 <- time_vary_ia$total_n3 * time_vary_ia$CCP2

# interaction mod
mod2 <- coxph(Surv(start_time, stop_time, lead_ia_ra == 1) ~ total_n3 + SE_num + 
                se_n3_int + Omega3_bi , data=time_vary_ia)
summary(mod2)

# dpa
dpa_mod <- coxph(Surv(start_time, stop_time, lead_ia_ra == 1) ~ total_n3 + SE_num + 
                se_n3_int + Omega3_bi , data=time_vary_ia)
summary(mod2)


# what about with lead omega3?
mod3 <- coxph(Surv(start_time, stop_time, lead_ia_ra == 1) ~ total_n3 + Omega3_bi, 
              data=time_vary_ia)
summary(mod3)

# changing omega3 visit to visit?
chg_mod <- coxph(Surv(start_time, stop_time, lead_ia_ra == 1) ~ tot_n3_chg + Omega3_bi, 
              data=time_vary_ia)
summary(chg_mod)

# similar results

# se pos stratum
se_pos_tv <- filter(time_vary_ia, SE_num == 1)
se_neg_tv <- filter(time_vary_ia, SE_num == 0)

# SE negative folks
mod3 <- coxph(Surv(start_time, stop_time, lead_ia_ra == 1) ~ total_n3, data=se_neg_tv)
summary(mod3)

# SE positive folks
mod3 <- coxph(Surv(start_time, stop_time, lead_ia_ra == 1) ~ total_n3 + Omega3_bi, data=se_pos_tv)
summary(mod3)


# how many ccp2 x se
xtabs(~ SEcount + CCP2, time_vary_ia)


mod3 <- coxph(Surv(start_time, stop_time, lead_ia_ra == 1) ~ lead_tot_n3 + CCP2 + ccp2_n3, 
              data=se_neg_tv)
summary(mod3)

# 3 way interaction se_num x ccp2 x n3
time_vary_ia  <- mutate(time_vary_ia, ccp2_totn3_se = lead_tot_n3 * CCP2 * SE_num,
                                      ccp2_totn3 = lead_tot_n3 * CCP2,
                                      ccp2_se = CCP2 * SE_num,
                                      totn3_se = lead_tot_n3 * SE_num)

xtabs(~SE_num + CCP2 + n3_supp + lead_ia_ra, time_vary_ia)
                          

surv_mod <- coxph(Surv(start_time, stop_time, lead_ia_ra == 1) ~ lead_tot_n3 + 
                    CCP2 + SE_num + ccp2_totn3 + ccp2_se + totn3_se + ccp2_totn3_se, 
              data=time_vary_ia)
summary(surv_mod)

# look fortrends in relationship -----------------------------------------------
n3_cut <- quantile(time_vary_ia$epa_dha, c(0, 1/5, 2/5, 3/5, 4/5, 5/5))
n3_cut

time_vary_ia <- mutate(time_vary_ia, n3_level_cat =
                       ifelse(epa_dha >= 1 & epa_dha < 3, 0,
                       ifelse(epa_dha>= 3 & epa_dha < 4, 1,
                       ifelse(epa_dha >= 4 & epa_dha < 5, 2,
                       ifelse(epa_dha >= 5 & epa_dha < 6, 3,
                       ifelse(epa_dha >= 6 & epa_dha < 7, 4,
                       ifelse(epa_dha >= 7 & epa_dha < 8, 5,
                       ifelse(epa_dha >= 8 & epa_dha < 9, 6, 
                       ifelse(epa_dha >= 9 & epa_dha < 10, 7, NA)))))))),
                       custom_epa_dha_cat = ifelse(epa_dha > 7, 1, 0))
                        

summary(as.factor(time_vary_ia$n3_level_cat))
xtabs(~ lead_ia_ra + custom_epa_dha_cat, time_vary_ia)

# group by n3_level_cat
n3_cat_ia <- group_by(time_vary_ia, custom_epa_dha_cat) %>% 
  summarise(denom = n(), count_ia_ra = sum(lead_ia_ra)) %>%
  mutate(prop_ia_ra = count_ia_ra / denom)

ggplot(n3_cat_ia, aes(x=custom_epa_dha_cat, y=prop_ia_ra)) + geom_point()


mod2 <- coxph(Surv(start_time, stop_time, lead_ia_ra == 1) ~ as.factor(n3_level_cat) + SE_num, 
              data=time_vary_ia)
summary(mod2)


# try letting age vary ---------------------------------------------------------
mod1 <- coxph(Surv(start_age, stop_age, lead_ia_ra == 1) ~ lead_tot_n3 , 
              data=time_vary_ia)
summary(mod1)


# subset the df to first observation 
first_obs_df <- filter(time_vary_ia, visit_num == 1) %>%
                mutate(ia_incident = ifelse(is.na(visit_convert), 0, 1),
                       tot_n3_cat = ifelse(total_n3 > 4.5 & total_n3 <= 6.77, 0,
                                    ifelse(total_n3 > 6.77 & total_n3 <= 7.66, 1,
                                    ifelse(total_n3 > 7.66 & total_n3 <= 13.1, 2, NA))))
# should have 15 converters
xtabs(~ia_incident + SE_num, first_obs_df) # yup
# set tertiles of omega3s in those who do not convert
no_convert <- filter(first_obs_df, ia_incident == 0)
tert <- quantile(no_convert$total_n3, c(0, 1/3, 2/3, 3/3))
tert

convert <- filter(first_obs_df, ia_incident == 1)
tert2 <- quantile(convert$total_n3, c(0, 1/3, 2/3, 3/3))
tert2

xtabs(~tot_n3_cat + ia_incident + SE_num, first_obs_df)

summary(lm(total_n3 ~ SE_num, first_obs_df))

# pattern with epa_dha
n3_cut <- quantile(no_convert$epa_dha, c(0, 1/5, 2/5, 3/5, 4/5, 5/5))
n3_cut

# converters
n3_cut2 <- quantile(convert$epa_dha, c(0, 1/5, 2/5, 3/5, 4/5, 5/5))
n3_cut2

first_obs_df <- mutate(first_obs_df, epa_dha_quint =
                       ifelse(epa_dha >= 2 & epa_dha < 3.78, 0,
                       ifelse(epa_dha >= 3.78 & epa_dha < 4.17, 1,
                       ifelse(epa_dha >= 4.17 & epa_dha < 4.68, 2,
                       ifelse(epa_dha >= 4.68 & epa_dha < 5.75, 3,
                       ifelse(epa_dha >= 5.75 & epa_dha < 10, 4, NA))))),
                       custom_epa_dha_cat = ifelse(epa_dha > 4.2, 1, 0))

xtabs(~ epa_dha_quint + ia_incident + SE_num, first_obs_df)

se_pos_df <- filter(first_obs_df, SE_num == 1)
summary(lm(ia_incident ~ as.factor(epa_dha_quint), data = se_pos_df))

# simple change in n3 time to event analysis -----------------------------------
mod <- coxph(Surv(time_contrib, ia_incident == 1) ~ dha_chg + SE_num + Race_bi + Omega3_bi, data=survival_df)
summary(mod)

xtabs(~ CCP2 + SE_num, survival_df)

# ala
summary(glm(ia_incident ~ ala_chg + Omega3_bi + Race_bi + SE_num, family = 'poisson', data=survival_df))
# epa
summary(glm(ia_incident ~ epa_chg + Omega3_bi + Race_bi + SE_num, family = 'poisson', data=survival_df))
# dpa
summary(glm(ia_incident ~ dpa_chg + Omega3_bi + Race_bi + SE_num, family = 'poisson', data=survival_df))
# dha
summary(glm(ia_incident ~ dha_chg + Omega3_bi + Race_bi + SE_num, family = 'poisson', data=survival_df))



# SE subset
se_pos_sub <- filter(survival_df, SE_num == 1)
se_neg_sub <- filter(survival_df, SE_num == 0)

summary(glm(ia_incident ~ dha_chg + CCP2 + Age + Omega3_bi, family = 'poisson', se_pos_sub))

# create interaction term with SE
# not seeing what elizabeth saw with chg
mod <- coxph(Surv(time_contrib, ia_incident == 1) ~ tot_n3_chg + SE_num + se_n3_chg + Omega3_bi + Race_bi, data=survival_df)
summary(mod)

glimpse(survival_df)
# base omega 3 lvls
survival_df$Omega3_bi[survival_df$ID == '09-011-00'] <- 0 # fill missing val
survival_df$se_dpa <- survival_df$dpa * survival_df$SE_num

mod <- coxph(Surv(time_contrib, ia_incident ==1) ~ dpa + SE_num + se_dpa + Race_bi, data = survival_df)
summary(mod)

mod <- coxph(Surv(time_contrib, ia_incident ==1) ~ Omega3_6_ratio + SE_num, data = survival_df)
summary(mod)
