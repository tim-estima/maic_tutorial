# MAIC on synthetic dataset

# Read in required libraries
library(haven)
library(dplyr)
library(survival)
library(tidyr)
library(survminer)

# Set baseline directory
base_dir <- 'C:/Users/kvmr654/Desktop/maic/'

# Read in main IPD MAIC dataset (one row per patient)
df_main <- read.csv(paste0(base_dir, 'maic_data.csv')) %>%
  rename(TRT=ï..TRT) %>% # Correct strange naming quirk when data read in
  mutate(biomarker_binary=ifelse(biomarker>600,1,0)) %>%
  select(-X)

# Read in aggregate dataset being compared to (from ANCHOR study)
df_anchor_aggregate<- read.csv(paste0(base_dir, 'anchor_aggregate_data.csv')) %>%
  rename(age=ï..age)

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
} # This is a fuction for Q(b) defined above.
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

# Create dataset for analysis by making sure weights are attached to dataset
df_analysis <- df_main %>%
  drop_na() %>%
  mutate(wt = as.vector(exp(as.matrix(df_main_centered %>% select(maic_vars)) %*% a1)))

# Check estimation of weights has worked- miniumum requirement is that the
# baseline characteristics match at the aggregate level
reweighted_anchor_baseline <- as.data.frame(
  df_analysis %>%
    drop_na() %>%
    bind_cols(as.data.frame(wt)) %>%
    summarise(
      age=weighted.mean(age, wt),
      male = 100*weighted.mean(male,wt),
      metastatic = 100*weighted.mean(metastatic,wt),
      biomarker=weighted.mean(biomarker, wt)
    )
)

ESS = sum(wt)^2/sum(wt^2)

###############Need to add histogram of rescaled weights here

# Data used for indirect comparison from anchor trial
anchor_pfs_log_hr <- log(0.97)
anchor_pfs_log_hr_lci <- log(0.86)
anchor_pfs_log_hr_uci <- log(1.11)

# Adjusted K-M from MAIN trial

km_adjusted <- survfit(
  Surv(AVAL, CNSR)~ TRT , 
  data=df_main,
  weights=wt
)

cox_adjusted <- coxph(
  Surv(AVAL, CNSR)~ TRT , 
  data=df_main,
  weights=wt
)

