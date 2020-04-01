# MAIC on synthetic dataset

# Read in required libraries
library(haven)
library(dplyr)
library(survival)
library(tidyr)
library(survminer)
library(RCurl)

# Reading in data- get directly from the github URLs below if this doesn't work
# Read in main IPD MAIC dataset (one row per patient)
x <- getURL("https://raw.githubusercontent.com/tim-estima/maic_tutorial/master/maic_data.csv")
df_main <- read.csv(text = x) %>% mutate(biomarker_binary=ifelse(biomarker>600,1,0))

# Read in aggregate dataset being compared to (from ANCHOR study)
y<- getURL('https://raw.githubusercontent.com/tim-estima/maic_tutorial/master/anchor_aggregate_data.csv')
df_anchor_aggregate<- read.csv(text=y)

# Summarise baseline characteristics of MAIN IPD trial at aggregate level
df_main %>%
  select(
    age,
    male,
    metastatic,
    biomarker
  ) %>%
  drop_na() %>%
  summarise(
    age_baseline=round(mean(age),2),
    male=paste0(as.character(round(sum(df_main$male)/length(df_main$male)*100),4),'%'),
    metastatic=paste0(as.character(round(sum(df_main$metastatic)/length(df_main$metastatic)*100),4),'%'),
    biomarker=round(mean(biomarker), 2)
  )

# Display baseline characteristics of ANCHOR trial (binary variables reformatted as percentages)
df_anchor_aggregate %>%
  mutate(
    male=sprintf("%1.0f%%", 100*male),
    metastatic=sprintf("%1.0f%%", 100*metastatic)
  )

df_main$TRT <- relevel(df_main$TRT, 'Estimab')

# Check effect modifiers in MAIN trial
cox_all_covariates_main <-
  coxph(
    Surv(AVAL, CNSR)~
      TRT + age + male + metastatic + biomarker_binary + # prognostic facotors
      TRT*age + TRT*male + TRT*metastatic + TRT*biomarker_binary, # effect modifiers
    data=df_main
  )

summary(cox_all_covariates_main, conf.int=0.9)

km_unadjusted <- survfit(
  Surv(AVAL, CNSR)~ TRT ,
  data=df_main
)

cox_unadjusted <- coxph(
  Surv(AVAL, CNSR)~ TRT ,
  data=df_main
)

# Select vars to include in MAIC

maic_vars <- c(
  'age',
  'male',
  'metastatic',
  'biomarker'
)

df_main_centered <- df_main %>%
  mutate(
    age=age-df_anchor_aggregate$age,
    male=male-df_anchor_aggregate$male,
    metastatic=metastatic-df_anchor_aggregate$metastatic,
    biomarker=biomarker-df_anchor_aggregate$biomarker
  ) %>%
  drop_na() # Standard MAIC can't handle missing values

# Complicated looking function to calculate weights for MAIC by Newton-Raphson

objfn <- function(a1, X){
  sum(exp(X %*% a1))
} # This is a function for Q(b) defined above.
# Gradient function => Derivative of Q(b).
gradfn <- function(a1, X){
  colSums(sweep(X, 1, exp(X %*% a1), "*"))
}
## The following is the in built R function used to optimise Q(b) ##
## using Newton-Raphson techniques ##
print(
  opt1 <- optim(
    par = rep(0,dim(df_main_centered %>% select(maic_vars))[2]) ,
    fn = objfn, gr = gradfn,
    X = as.matrix(df_main_centered %>% select(maic_vars)),
    method = "BFGS"
  )
)

a1 = opt1$par

wt = as.vector(exp(as.matrix(df_main_centered %>% select(maic_vars)) %*% a1))

# Create dataset for analysis by making sure weights are attached to dataset
df_analysis <- df_main %>%
  drop_na() %>%
  mutate(wt = as.vector(exp(as.matrix(df_main_centered %>% select(maic_vars)) %*% a1)))

# Check estimation of weights has worked- miniumum requirement is that the
# baseline characteristics match at the aggregate level
reweighted_anchor_baseline <- as.data.frame(
  df_analysis %>%
    drop_na() %>%
    summarise(
      age=weighted.mean(age, wt),
      male = 100*weighted.mean(male,wt),
      metastatic = 100*weighted.mean(metastatic,wt),
      biomarker=weighted.mean(biomarker, wt)
    )
)

ESS = sum(df_analysis$wt)^2/sum(df_analysis$wt^2)

# Data used for indirect comparison from anchor trial
anchor_pfs_log_hr <- log(0.97)
anchor_pfs_log_hr_lci <- log(0.86)
anchor_pfs_log_hr_uci <- log(1.11)



# Adjusted K-M from MAIN trial

km_adjusted <- survfit(
  Surv(AVAL, CNSR)~ TRT ,
  data=df_analysis,
  weights=wt
)

cox_adjusted <- coxph(
  Surv(AVAL, CNSR)~ TRT ,
  data=df_analysis,
  weights=wt
)

# Rescaled weights and histogram
wt.rs <- (wt / sum(wt))*nrow(df_analysis)

qplot(wt.rs, geom='histogram', xlab='Rescaled weight', binwidth=0.25)

# Anchor log hazard ratio and standard error
anchor_log_hr_pe <- log(0.97)
anchor_log_hr_se <- (log(1.11)- log(0.86))/(2*1.96)

# Main log hazard ratio and standard error- adjusted and unadjusted
main_unadjusted_log_hr_pe <- log(summary(cox_unadjusted)$conf.int[1])
main_unadjusted_log_hr_se <- (log(summary(cox_unadjusted)$conf.int[4])- log(summary(cox_unadjusted)$conf.int[3]))/(2*1.96)

main_adjusted_log_hr_pe <- log(summary(cox_adjusted)$conf.int[1])
main_adjusted_log_hr_se <- (log(summary(cox_adjusted)$conf.int[4])- log(summary(cox_adjusted)$conf.int[3]))/(2*1.96)

#ITC- adjusted and unadjusted
adjusted_itc_pe <- exp(anchor_log_hr_pe - main_adjusted_log_hr_pe)
adjusted_itc_se <- sqrt(anchor_log_hr_se^2)
adjusted_itc_lci <- exp(log(adjusted_itc_pe) - 1.96*adjusted_itc_se)
adjusted_itc_uci <- exp(log(adjusted_itc_pe) + 1.96*adjusted_itc_se)

unadjusted_itc_pe <- exp(anchor_log_hr_pe - main_unadjusted_log_hr_pe)
unadjusted_itc_se <- sqrt(anchor_log_hr_se^2)
unadjusted_itc_lci <- exp(log(unadjusted_itc_pe) - 1.96*unadjusted_itc_se)
unadjusted_itc_uci <- exp(log(unadjusted_itc_pe) + 1.96*unadjusted_itc_se)