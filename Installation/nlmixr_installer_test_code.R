#Make sure there is only a single library path, and that it points to the nlmixr installer location:
.libPaths()
# should be:
#[1] "C:/R/nlmixr_1.1.1-3/R/library"
#
#if not, set it to the correct directory:
#.libPaths("C:/R/nlmixr_1.1.1-3/R/library")
#or run the following code to remove libPaths that do not point to the installer directory
#.libPaths(.libPaths()[regexpr("nlmixr",.libPaths())!=-1])

## load the required libraries
library(xpose.nlmixr)
library(nlmixr)
library(RxODE)
library(data.table)
library(lattice)
## (or use your favourite graphics package)

#################################################################################
##                                                                             ##
## RxODE simulation Part 1                                                     ##
##                                                                             ##
## Simulation of a single (warfarin concentration) curve with a single dose    ##
##                                                                             ##
#################################################################################

## set up the system of differential equations (ODEs)
odeKA1 <- "
 d/dt(depot)   = -ka*depot;
 d/dt(central) =  ka*depot-(cl/v)*central;
 C1=central/v;
"

## compile the model
modKA1 <- RxODE(model = odeKA1)

## provide the parameter values to be simulated:
Params <- c(
 ka = log(2) / 0.5, # 1/h (aborption half-life of 30 minutes)
 cl = 0.135,        # L/h
 v = 8              # L
)

## create an empty event table that stores both dosing and sampling information :
ev <- eventTable()

## add a dose to the event table:
ev$add.dosing(dose = 500) #mg
              
## add time points to the event table where concentrations will be simulated; these actions are cumulative
ev$add.sampling(seq(0, 120, 0.1))

## Then solve the system
##
## The output from rxSolve is a solved RxODE object,
##  but by making it a data.frame only the simulated values are kept:
Res<-data.frame(rxSolve(modKA1,Params,ev))

## then plot the simulated outcomes in the compartments:
xyplot(depot~time,data=Res,type='b')
xyplot(C1~time,data=Res,type='b')


#################################################################################
##                                                                             ##
## nlmixr analysis Part 1                                                      ##
##                                                                             ##
## Warfarin population PK using SAEM estimation                                ##
##                                                                             ##
#################################################################################

## read in the Warfarin PK-only data set
data(warfarin)
PKdata <- as.data.table(warfarin)
PKdata<-PKdata[dvid=="cp"]

## Define a first order-absorption linear elimination model using a solved system solution

One.comp.KA.solved <- function() {
  ini({
    # Where initial conditions/variables are specified
    lka  <- log(1.15)  #log ka (1/h)
    lcl  <- log(0.135) #log Cl (L/h)
    lv   <- log(8)     #log V (L)
    prop.err <- 0.15
    add.err  <- 0.6
    eta.ka ~ 0.5   #IIV ka                
    eta.cl ~ 0.1   #IIV cl
    eta.v  ~ 0.1   #IIV v 
  })
  model({
    # Where the model is specified
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv + eta.v)
    ka <- exp(lka + eta.ka)
    ## solved system example
    ## where residual error is assumed to follow proportional and additive error
    linCmt() ~ prop(prop.err) + add(add.err)
  })
}


## estimate parameters using nlmixr:
fitOne.comp.KA.solved_S <-
  nlmixr(
    One.comp.KA.solved,           #the model definition
    PKdata,          #the data set
    est = "saem",    #the estimation algorithm (SAEM)
    #the SAEM minimisation options:
    saemControl(nBurn = 200, #200 SAEM burn-in iterations (the default)
                nEm   = 300, #300 EM iterations (the default)
                print = 50), #only print every 50th estimation step (default=1 which gives endless output)
    tableControl(cwres=TRUE) #calculates NONMEM-style conditional weighted residuals for diagnostics
  )
fitOne.comp.KA.solved_S

#################################################################################
##                                                                             ##
## nlmixr analysis Part 2                                                      ##
##                                                                             ##
## Goodness of fit plots using Ben Guiastrennec's xpose                        ##
##                                                                             ##
#################################################################################

## the nlmixr object can be tranformed into an xpose object to allow diagnostics with the new xpose package
xpdb.1s <- xpose_data_nlmixr(fitOne.comp.KA.solved_S)

## dv vs ipred plot:
dv_vs_ipred(xpdb.1s,        #the xpose object
            caption = NULL) #stops printing of the directory where this was run

## CWRES vs time:
res_vs_idv(xpdb.1s,          #the xpose object
           res = "CWRES",    #examine conditional weighted residuals
           idv = "TIME")     #as a function of time


#################################################################################
##                                                                             ##
## nlmixr analysis Part 3                                                      ##
##                                                                             ##
## VPCs using the Ron Keizer's vpc and built-in nlmixr functionality           ##
##                                                                             ##
#################################################################################

## nlmixr comes with its built-in vpc functionality that uses Ron Keizer's vpc package
vpc_ui(
  fitOne.comp.KA.solved_S,                    #the nlmixr object
  n = 500,                       #number of trials simulated using estimated parameters and study sampling structure
  show = list(obs_dv = TRUE),    #additional items to show, like the observations
  xlab = "Time (h)",             #x-axis label
  ylab = "Concentration (mg/L)", #y-axis label
  title= "VPC for first order absorption PopPK model\nSAEM"
)
