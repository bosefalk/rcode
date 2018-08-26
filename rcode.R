## ----------------------------------------------------------------------------

# in case you don't have a dataset ready run this part first:
n = 1000

simulWeib <- function(N, lambda, rho, beta, rateC, grp1)
{
  # at least have 1 covariate be related to survival
  x <- sample(x = c(0, 1), size = N, replace = T, prob = c(grp1, 1 - grp1))
  v <- runif(n = N)
  Tlat <- (-log(v) / (lambda * exp(x * beta))) ^ (1 / rho)
  # censoring times
  C <- rexp(n = N, rate = rateC)
  # follow-up times and event indicators
  srv_allo1 <- pmin(Tlat, C)
  srv_s_allo1 <- as.numeric(Tlat <= C)
  ci_allo1 <- srv_allo1
  ci_temp <- sample(0:1, N, prob = c(0.78, 0.22), replace = T)
  ci_s_allo1 <- ifelse(ci_temp == 1, 1, srv_s_allo1 * 2)
  
  # data set
  data.frame(srv_allo1, srv_s_allo1, ci_allo1, ci_s_allo1, x)
}

sd0 <- simulWeib(1000, 0.025, 1, 0.001, 0.01, 0.6)

sd0$srv_s_allo1 <- ifelse(sd0$srv_allo1 > 48, 0, sd0$srv_s_allo1)
sd0$srv_allo1 <- ifelse(sd0$srv_allo1 > 48, 48, sd0$srv_allo1)

sd0$ci_s_allo1 <- factor(ifelse(sd0$ci_allo1 > 48, 0, sd0$ci_s_allo1), levels = 0:2, labels = c('rfs','relapse','nrm'))
sd0$ci_allo1 <- ifelse(sd0$ci_allo1 > 48, 48, sd0$ci_allo1)

sd0$PATSEX <- factor(sd0$x, levels = 0:1, labels = c('Male','Female')) 
sd0$age_allo1 <- runif(n, min = 18, max = 85)
sd0$intdiagtr_allo1 <- runif(n, min = 0.5, max = 24)
sd0$ric_allo1 <- factor(sample(0:1, 1000, prob = c(0.22, 0.78), replace = T), levels = 0:1, labels = c('Standard','Reduced')) 
miss <- sample(0:1, n, replace = T, prob = c(0.88, 0.12))
sd0$ric_allo1[miss == 1] <- NA

sd0$source_allo1 <- factor(sample(0:1, 1000, prob = c(0.29, 0.71), replace = T), levels = 0:1, labels = c('PB','BM')) 

sd0$DONRL_allo1_1 <- factor(sample(0:1, 1000, prob = c(0.47, 0.53), replace = T), levels = 0:1, labels = c('MRD','MUD')) 
miss <- sample(0:1, n, replace = T, prob = c(0.91, 0.09))
sd0$DONRL_allo1_1[miss == 1] <- NA

## ----------------------------------------------------------------------------

# source(survtables.R)
# source(utils.R)
# source(survplot.R)

# Tables: all are rendered as knitr tables by default

## ----------------------------------------------------------------------------

## 1. Descriptive tables 

## ----------------------------------------------------------------------------

## Frequency table

# Standard frequencies and proportions can be provided by the function freqtab
# freqtab(x, caption = '', digits = 1, order = '', valid.only = F, percent = T, ...)

# with defaults
freqtab(sd$PATSEX, caption = 'Patient sex')

# you have the option to only report valid percent: 
freqtab(sd$PATSEX, caption = 'Patient sex', valid.only = T)

# or not report percentages at all:
freqtab(sd$PATSEX, caption = 'Patient sex', percent = F)

# the number of decimal points can be specified with digits = ...
freqtab(sd$PATSEX, caption = 'Patient sex', digits = 3)

## Continuous variables

# Overview of quantiles for a continuous variable is given by the function examine
# examine(cont, cat = '', colname = '', caption = '', xcol = '', digits = 1)
examine(sd$age_allo1)

# as above, the table name and decimal points can be specified.
# the quantiles can be stratified by any categorical variable in your data frame
examine(sd$age, cat = 'PATSEX')

# can specify how the strata are named:
examine(sd$age, cat = 'PATSEX', xcol = c('M','F'))

## Table one

# This function creates a table one: 
# tableone(data, by = '', margin = 2, cat.names = '', cont.names = '', caption = '', digits = 1, test = F, df = F)

# it assumes the variables are of the correct format and expects a dataframe with only the variables of interest

# variables of interest
sd1 = sd[, c('PATSEX','ric_allo1','source_allo1','DONRL_allo1_1','intdiagtr_allo1','age_allo1')]

# default table
tableone(sd1)

# it might be useful to explicitly name the continuous variables
tableone(sd1, cont.names = c('Interval diag','Age'))

# the same can be done for the categorical variables (mind the order)
tableone(sd1, cat.names = c('Sex','RIC','Source','Donor'),
         cont.names = c('Interval diag tx','Age'))

# if needed, the table can be stratified by any categorical variable
tableone(sd1, by = 'PATSEX', cont.names = c('Interval diag','Age'))

# and you can test for significant differences if you want by specifying test = T
tableone(sd1, by = 'PATSEX', cont.names = c('Interval diag','Age'), test = T)
## note this part does not work with specified cat.names yet

## ----------------------------------------------------------------------------

## 2. Univariate survival outcomes

## ----------------------------------------------------------------------------

# This section requires output from the survival library
# for standard survival outcomes, use survfit(Surv() ~ x, data) as usual
# for competing risks, ensure the status indicator is a factor and the first level is the censoring indicator

#sd$srv_s_allo1 = as.numeric(sd$srv_s_allo1) - 1

# replace these by which ever stratification you want
os <- survfit(Surv(srv_allo1, srv_s_allo1) ~  1, data = sd)
sex <- survfit(Surv(srv_allo1, srv_s_allo1) ~ PATSEX, data = sd)
ric <- survfit(Surv(srv_allo1, srv_s_allo1) ~ ric_allo1, data = sd)

ci <- survfit(Surv(ci_allo1, ci_s_allo1) ~ PATSEX, data = sd0)
ci0 <- survfit(Surv(ci_allo1, ci_s_allo1) ~ 1, data = sd)

## Figures
# survplot(fit, cause = 1, byval = 12, tstart = 0, xlab = '', order = 0, stacked = F, title = '', colormap = 'Set1', bw = F, legendpos = 'topright', lgd = '')

survplot(os)

survplot(ci, cause = 1)
survplot(ci, cause = 2)

# both can be shown side by side
survplot(ci, cause = 1:2)


survplot(ci0, stacked = T, cause = 1:2, bw = T, lgd = c('relapse','nrm'))

## Survival table

# this table provides pretty printing for normal and competing risks outcomes
# survtab(fit, times, cause = 1, digits = 0, caption = '', outcomes = '', median = F, test = T)

# default table: note that times must always be specified
survtab(os, times = c(24, 60))

# stratifying:
survtab(ric, times = c(24, 60))

# say we want to create a type of table with all outcomes
survtab(list(os, ric, sex), times = c(24, 60))

# we can do something similar with competing risks, except now we can specify the causes as well
survtab(ci, times = c(24), cause = 1:2)

# if i do this, it might be better to specify the cause names:
survtab(ci, times = c(24), cause = 1:2, outcomes = c('Relapse','NRM'))

# median can be provided separately but is included in survtab
# function(fit, digits = 1, cause = 1, timescale = 'months', caption = '', internal = F)

# with default settings
med(os)

# works the same with competing risks
med(ci0)

# and with stratification
med(ci, cause = 1)
med(ci, cause = 2)


## ----------------------------------------------------------------------------

## 3. Multivariable analysis table

## ----------------------------------------------------------------------------

## Regression table
# mvatab (fit, caption = '', showevents = TRUE, digits = 2, df = F)
cox = coxph(Surv(srv_allo1, srv_s_allo1) ~ PATSEX + intdiagtr_allo1 + ric_allo1 + age_allo1, data = sd)

mvatab(cox)

mvatab(cox, digits = 3)

## Also works with interactions
cox1 = coxph(Surv(srv_allo1, srv_s_allo1) ~ PATSEX + intdiagtr_allo1*ric_allo1 + age_allo1, data = sd)

mvatab(cox1)

# competing risks work with cause specific hazards
# however it does expect a binary status variable:
sd$ci_s_1 = ifelse(sd$ci_s_allo1 == 'relapse', 1, 0)
sd$ci_s_2 = ifelse(sd$ci_s_allo1 == 'nrm', 1, 0)

csh1 = coxph(Surv(ci_allo1, ci_s_1) ~ PATSEX + intdiagtr_allo1 + ric_allo1 + age_allo1, data = sd)

mvatab(csh1)

# and for NRM
csh2 = coxph(Surv(ci_allo1, ci_s_2) ~ PATSEX + intdiagtr_allo1 + ric_allo1 + age_allo1, data = sd)

mvatab(csh2)
