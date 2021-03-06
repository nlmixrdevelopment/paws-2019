## load the required libraries
library(xpose.nlmixr)
library(nlmixr)
library(RxODE)
library(lattice)
library(data.table)

## specify array of colours for curves
nlmixCOLS <- c("#28466A","#8DB6CD","#B40000")  

#################################################################################
##                                                                             ##
## nlmixr analysis Part 3                                                      ##
##                                                                             ##
## Generating Bayesian feedback estimates                                      ##
##                                                                             ##
#################################################################################

## read in the Warfarin PK-only data set using data.table syntax (fast and efficient!)
PKdata <- fread("warfarin_PK.csv")


## Generating Bayesian feedback estimates
## nlmixr can generate empirical Bayes estimates for Bayesian feedback:
##  individual EBEs for a new data set using existing population parameters
## Useful in a therapeutic drug monitoring setting 
## Or for generating exposure estimates with a particularly nasty model that you 
##  do not want to refit on new data :-)

KA1tr5posthoc <- function() {
  ini({
    # Where initial conditions/variables are specified
    lktr <-  1.18994619   #log ktr (/h)
    lcl  <- -2.01737477   #log Cl (L/h)
    lv   <-  2.06631620   #log V (L)
    prop.err <- 0.07883633#proportional error (SD/mean)
    add.err <- 0.37249666 #additive error (mg/L)
    eta.ktr ~ 0.2532964   #IIV ktr
    eta.cl ~ 0.08073339   #IIV Cl
    eta.v ~ 0.04490733   #IIV V
  })
  model({
    # Where the model is specified
    # The model uses the ini-defined variable names
    cl  <- exp(lcl + eta.cl)
    v   <- exp(lv + eta.v)
    ktr <- exp(lktr + eta.ktr)
    # RxODE-style differential equations are supported
    d/dt(trns1) = -ktr * trns1
    d/dt(trns2) =  ktr * trns1 - ktr * trns2
    d/dt(trns3) =  ktr * trns2 - ktr * trns3
    d/dt(trns4) =  ktr * trns3 - ktr * trns4
    d/dt(trns5) =  ktr * trns4 - ktr * trns5
    d/dt(central) =  ktr * trns5 - (cl/v) * central
    ## Concentration is calculated
    cp = central/v
    # And is assumed to follow proportional and additive error
    cp ~ prop(prop.err) + add(add.err)
  })
}

## specify posthoc as estimation method to only generate EBEs and not do any estimation

fitKA1tr5_Fph <- nlmixr(KA1tr5posthoc, PKdata, est = "posthoc")

##Generate individual smooth profiles using augPred
indivpk <- augPred(fitKA1tr5_Fph)

##and plot using lattice

xyplot(
  values ~ time | id,
  data = indivpk,
  groups = ind,
  type = c("l", "l", "p"),
  col = nlmixCOLS[c(2, 1, 3)],
  cex = c(0.1, 0.1, 1),
  lwd = c(2, 2, 2),
  pch = c(1, 1, 19),
  xlab = "Time (h)\n",
  ylab = "Warfarin (mg/L)",
  layout=c(8,4),
  as.table = TRUE,
  scales = list(alternating = 1),
  main = "Five transit compartment model estimated using Bayesian feedback",
  auto.key = list(
    adj = 1,
    col = nlmixCOLS[c(2, 1, 3)],
    columns = 3,
    space = "bottom",
    rectangles = FALSE,
    points = FALSE
  )
)



#################################################################################
##                                                                             ##
## nlmixr analysis Part 4                                                      ##
##                                                                             ##
## Implement mu-referenced covariates on log-scale                             ##
##                                                                             ##
#################################################################################


Covs <- PKdata[!duplicated(ID)]
table(Covs$SEX)
# 0  1
# 5 27

## One compartment transit model with Sex on V
KAtr1_sexV <- function() {
  ini({
    # Where initial conditions/variables are specified
    lktr <- log(1.15)  #log k transit (/h)
    lcl  <- log(0.135) #log CL (L/h)
    lv   <- log(8)     #log V (L)
    Sex_V <- 0.1       #log Sex on v
    prop.err <- 0.15   #proportional error (SD/mean)
    add.err <- 0.6     #additive error (mg/L)
    eta.ktr ~ 0.5   #IIV ktr
    eta.cl ~ 0.1   #IIV Cl
    eta.v ~ 0.1   #IIV V
  })
  model({
    #Sex on volume
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv + eta.v + Sex_V * SEX)
    ktr <- exp(lktr + eta.ktr)
    # RxODE-style differential equations are supported
    d/dt(depot) = -ktr * depot
    d/dt(central) =  ktr * trans - (cl/v) * central
    d/dt(trans)   =  ktr * depot - ktr * trans
    ## Concentration is calculated
    cp = central/v
    # And is assumed to follow proportional and additive error
    cp ~ prop(prop.err) + add(add.err)
  })
}

fitKAtr1_sexV_F <-
  nlmixr(KAtr1_sexV, PKdata, est = "focei", foceiControl(print = 20))
fitKAtr1_sexV_F
save(fitKAtr1_sexV_F, file = "fitKAtr1_sexV_F.Rdata")
#load(file = "fitKAtr1_sexV_F.Rdata")


## binary categorical effects estimated on log scale can be back-transformed
THETA_SexV <- fitKAtr1_sexV_F$parFixedDf[4, ]
THETA_SexV
#       Estimate        SE     %RSE Back-transformed    CI Lower  CI Upper BSV(CV%) Shrink(SD)%
#Sex_V 0.3937038 0.2440992 62.00072        0.3937038 -0.08472182 0.8721294       NA            

## this results in fold-change estimates
exp(c(
  THETA_SexV$Estimate,
  THETA_SexV$"CI Lower",
  THETA_SexV$"CI Upper"))
# 1.4824614 0.9187678 2.3919989 

## that can be translated in percentage change with proper confidence intervals
paste0(round(100 * (exp(c(
  THETA_SexV$Estimate,
  THETA_SexV$"CI Lower",
  THETA_SexV$"CI Upper")) - 1), 1), "%")
#[1] "48.2%"  "-8.1%"  "139.2%"



#################################################################################
##                                                                             ##
## Implement allometric covariates                                             ##
##                                                                             ##
#################################################################################

## Code is most stable if transformations are carried out in the data file
## instead of in the model code, especially for SAEM
## Using standard R syntax:
## PKdata$logWT70 <- log(PKdata$WT/70)

## Or using data.table syntax
PKdata[,logWT70:=log(WT/70)]

## One compartment transit model with allometric scaling on WT
One.comp.transit.allo <- function() {
  ini({
    # Where initial conditions/variables are specified
    lktr <- log(1.15)  #log k transit (/h)
    lcl  <- log(0.135) #log Cl (L/h)
    lv   <- log(8)     #log V (L)
    ALLC <- fix(0.75)  #allometric exponent cl
    ALLV <- fix(1.00)  #allometric exponent v
    prop.err <- 0.15   #proportional error (SD/mean)
    add.err <- 0.6     #additive error (mg/L)
    eta.ktr ~ 0.5   #IIV ktr
    eta.cl ~ 0.1   #IIV Cl
    eta.v ~ 0.1   #IIV V
  })
  model({
    #Allometric scaling on weight
    cl <- exp(lcl + eta.cl + ALLC * logWT70)
    v  <- exp(lv + eta.v + ALLV * logWT70)
    ktr <- exp(lktr + eta.ktr)
    # RxODE-style differential equations are supported
    d/dt(depot)   = -ktr * depot
    d/dt(central) =  ktr * trans - (cl/v) * central
    d/dt(trans)   =  ktr * depot - ktr * trans
    ## Concentration is calculated
    cp = central/v
    # And is assumed to follow proportional and additive error
    cp ~ prop(prop.err) + add(add.err)
  })
}

fitOne.comp.transit.allo_F <-
  nlmixr(One.comp.transit.allo, 
         PKdata, 
         est = "focei", 
         foceiControl(print = 5))
fitOne.comp.transit.allo_F
save(fitOne.comp.transit.allo_F, file = "fitOne.comp.transit.allo_F.Rdata")
#load(file = "fitOne.comp.transit.allo_F.Rdata")

vpc_ui(
  fitOne.comp.transit.allo_F,    #the nlmixr object
  n = 500,  #number of replicates simulated using estimated parameters and study sampling structure
  show = list(obs_dv = TRUE),    #additional items to show, like the observations
  xlab = "Time (h)",             #x-axis label
  ylab = "Concentration (mg/L)", #y-axis label
  title = "One transit compartment and allometric scaling FOCEI"
)

## do you get a significant drop in OFV by including allometric weight?
load(file="fitOne.comp.transit_F.Rdata")
fitOne.comp.transit.allo_F$OBJF-fitOne.comp.transit_F$OBJF
#[1] -29.49277

#################################################################################
##                                                                             ##
##  Hands-on assignments: nlmixr model development                             ##
##  Run the allometric model without fixing the exponents                      ##
##  Do you get a better fit?                                                   ##
##                                                                             ##
##  fitPK001$OBJF-fitPK002$OBJF                                                ##
##                                                                             ##
#################################################################################

#...

#################################################################################
##                                                                             ##
#################################################################################




