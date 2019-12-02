#################################################################################
##                                                                             ##
##  Solutions to hands-on assignments: nlmixr model development                ##
##  Examine the GOF plots and implement models with                            ##
##  one or transit compartments                                                ##
##                                                                             ##
##  Compare vpcs of alternatives and compare OFVs:                             ##
##                                                                             ##
##  fitPK001$OBJF-fitPK002$OBJF                                                ##
##                                                                             ##
#################################################################################


## load the required libraries
library(xpose.nlmixr)
library(nlmixr)
library(RxODE)
library(lattice)
library(data.table)

## read in the Warfarin PK-only data set using data.table syntax (fast and efficient!)
PKdata <- fread("warfarin_PK.csv")


#################################################################################
##                                                                             ##
##   Update the model with a single effect compartment                         ##
##                                                                             ##
#################################################################################

## One compartment transit model
One.comp.transit <- function() {
  ini({
    # Where initial conditions/variables are specified
    lktr <- log(1.15)  #log k transit (/h)
    lcl  <- log(0.135) #log Cl (L/h)
    lv   <- log(8)     #log V (L)
    prop.err <- 0.15   #proportional error (SD/mean)
    add.err <- 0.6     #additive error (mg/L)
    eta.ktr ~ 0.5   #IIV ktr
    eta.cl ~ 0.1   #IIV Cl
    eta.v ~ 0.1   #IIV V
  })
  model({
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv + eta.v)
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

#################################################################################
##                                                                             ##
##   Run using SAEM                                                            ##
##                                                                             ##
#################################################################################

fitOne.comp.transit_S <-
  nlmixr(One.comp.transit,
         PKdata,
         est = "saem",
         saemControl(print = 100),
         tableControl(cwres = TRUE))
save(fitOne.comp.transit_S,file="fitOne.comp.transit_S.Rdata")
fitOne.comp.transit_S


vpc_ui( 
  fitOne.comp.transit_S,                  #the nlmixr object
  n = 500,                       #number of replicates simulated using estimated parameters and study sampling structure
  show = list(obs_dv = TRUE),    #additional items to show, like the observations
  xlab = "Time (h)",             #x-axis label
  ylab = "Concentration (mg/L)", #y-axis label
  title= "One transit compartment SAEM"
)

#################################################################################
##                                                                             ##
##   Run using FOCEI                                                           ##
##                                                                             ##
#################################################################################

fitOne.comp.transit_F <-
  nlmixr(One.comp.transit, 
         PKdata, 
         est = "focei", 
         foceiControl(print = 5))
fitOne.comp.transit_F
save(fitOne.comp.transit_F, file = "fitOne.comp.transit_F.Rdata")


vpc_ui(
  fitOne.comp.transit_F,    #the nlmixr object
  n = 500,  #number of replicates simulated using estimated parameters and study sampling structure
  show = list(obs_dv = TRUE),    #additional items to show, like the observations
  xlab = "Time (h)",             #x-axis label
  ylab = "Concentration (mg/L)", #y-axis label
  title = "One transit compartment FOCEI"
)

xpdb.1f <- xpose_data_nlmixr(fitOne.comp.transit_F)

#Absolute values of individual weighted residual vs time
IWRES1<-absval_res_vs_idv(xpdb.1f,        #the xpose object
                  res = "IWRES",  #examine absolute values (absval) of individual weighted residuals
                  idv = "TIME",   #as a function of time
                  caption = NULL) #if not NULL provides the directory where this was run

IWRES1



#################################################################################
##                                                                             ##
##   Update the model with five effect compartments                            ##
##                                                                             ##
#################################################################################

## 5 transit compartments
KA1tr5ode <- function() {
  ini({
    # Where initial conditions/variables are specified
    lktr  <- log(1.15) #log transit rate constant (/h)
    lcl  <- log(0.135) #log Cl (L/h)
    lv   <- log(8)     #log V (L)
    prop.err <- 0.15   #proportional error (SD/mean)
    add.err  <- 0.6    #additive error (mg)
    eta.ktr ~ 0.5   #IIV ktr                
    eta.cl ~ 0.1    #IIV cl
    eta.v  ~ 0.1    #IIV v 
  })
  model({
    # Where the model is specified
    ktr <- exp(lktr + eta.ktr)
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv + eta.v)
    ## ODE example
    cc=central/v
    d/dt(depot)= - ktr*depot
    d/dt(central) = ktr*transit5 - cl*cc
    d/dt(transit1)= ktr*(depot - transit1) 
    d/dt(transit2)= ktr*(transit1 - transit2) 
    d/dt(transit3)= ktr*(transit2 - transit3) 
    d/dt(transit4)= ktr*(transit3 - transit4) 
    d/dt(transit5)= ktr*(transit4 - transit5) 
    
    ## where residual error is assumed to follow proportional and additive error
    cc ~ prop(prop.err) + add(add.err)
  })
}
nlmixr(KA1tr5ode)

#################################################################################
##                                                                             ##
##   Run using SAEM                                                            ##
##                                                                             ##
#################################################################################

fitKA1tr5ode_S <-
  nlmixr(KA1tr5ode,
         PKdata,
         est = "saem",
         saemControl(print = 100),
         tableControl(cwres = TRUE))
fitKA1tr5ode_S
save(fitKA1tr5ode_S,file="fitKA1tr5ode_S.Rdata")

vpc_ui( 
  fitKA1tr5ode_S,                #the nlmixr object
  n = 500,                       #number of replicates simulated using estimated parameters and study sampling structure
  show = list(obs_dv = TRUE),    #additional items to show, like the observations
  xlab = "Time (h)",             #x-axis label
  ylab = "Concentration (mg/L)", #y-axis label
  title= "VPC with 5 transit compartments SAEM"
)


fitKA1tr5ode_S$OBJF
#[1] 270.5943
fitOne.comp.transit_F$OBJF
#[1] 321.1425
fitKA1tr5ode_S$OBJF - fitOne.comp.transit_F$OBJF
#[1] -50.5482


#################################################################################
##                                                                             ##
##   Run using FOCEI                                                           ##
##                                                                             ##
#################################################################################
  
fitKA1tr5ode_F <-
  nlmixr(KA1tr5ode,
         PKdata,
         est = "focei",
         foceiControl(print = 20),
         tableControl(cwres = TRUE))
fitKA1tr5ode_F
save(fitKA1tr5ode_F,file="fitKA1tr5ode_F.Rdata")
#load(file="fitKA1tr5ode_F.Rdata")

vpc_ui( 
  fitKA1tr5ode_F,                    #the nlmixr object
  n = 500,                       #number of replicates simulated using estimated parameters and study sampling structure
  show = list(obs_dv = TRUE),    #additional items to show, like the observations
  xlab = "Time (h)",             #x-axis label
  ylab = "Concentration (mg/L)", #y-axis label
  title= "VPC with 5 transit compartments FOCEI"
)



xpdb.3f <- xpose_data_nlmixr(fitKA1tr5ode_F)

#Absolute values of individual weighted residual vs time
IWRES3<-absval_res_vs_idv(xpdb.3f,        #the xpose object
                  res = "IWRES",  #examine absolute values (absval) of individual weighted residuals
                  idv = "TIME",   #as a function of time
                  caption = NULL) #if not NULL provides the directory where this was run
IWRES3



#One transit compartment vs 5 transit compartments:
fitKA1tr5ode_F$OBJF - fitOne.comp.transit_F$OBJF
#[1] -90.81914
