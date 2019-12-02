#Make sure there is only a single library path, and that it points to the nlmixr installer location:
.libPaths()
# should be:
#[1] "C:/R/nlmixr_1.1.1-3/R/library"
#
#if not, and nlmixr is not working, set it to the correct directory:
#.libPaths("C:/R/nlmixr_1.1.1-3/R/library")
#or run this command to remove libPaths that do not point to the installer directory
#.libPaths(.libPaths()[regexpr("nlmixr",.libPaths())!=-1])

## load the required libraries
library(RxODE)
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
## the amounts in the depot compartment
xyplot(depot~time,data=Res,type='b')
## the concentrations in the central compartment
xyplot(C1~time,data=Res,type='b')


## Extend the eventTable by adding three infusions to the central compartment
## Remember: updates to the eventTable are cumulative

ev$add.dosing(
  dose = 250,           #mg
  nbr.doses = 3,        #add three doses
  dosing.to = 2,        #add them to the second ODE in the model (=central)
  dosing.interval = 12, #h; set the doses 12 hours apart
  rate = 125,           #mg/h; infuse at a rate of 125 mg/h, resulting in 2-hour infusions
  start.time = 36       #h; have the three doses start at 36h
)

Res<-data.frame(rxSolve(modKA1,Params,ev))

## only the first dose goes into the depot compartment, so nothing changes here:
xyplot(depot~time,data=Res,type='b')

## the concentrations in the central compartment now reflect the three additional infusions: in the compartments:
xyplot(C1~time,data=Res,type='b')




#################################################################################
##                                                                             ##
##  Hands-on assignments: RxODE simulation                                     ##
##  -Simulate the effect of changing parameters                                ##
##  -Simulate the effect of changing doses                                     ##
##                                                                             ##
#################################################################################

#...

#################################################################################
##                                                                             ##
#################################################################################




#################################################################################
##                                                                             ##
## Extending the model with a transit compartment                              ##
##                                                                             ##
#################################################################################

## This requires an update to the ODE equation for the central compartment
##  because amounts now come from trans instead of depot:

odeKA1trans <- "
 d/dt(depot)   = -ka*depot;
 d/dt(central) =  ktr*trans-(cl/v)*central;
 d/dt(trans)   =  ka*depot-ktr*trans; 
 C1=central/v;
"

## compile the model
modKA1trans <- RxODE(model = odeKA1trans)

## provide the extra ktr parameter:
Params2 <- c(
  ka = log(2) / 0.5, # 1/h (aborption half-life of 30 minutes)
  cl = 0.135,        # L/h
  v = 8,             # L
  ktr = log(2) / 5   # 1/h (transit half-life of 5 hours)
)

## the eventTable does not have to change

Res<-data.frame(rxSolve(modKA1trans,Params2,ev))

## Examine the results:
## Only the first dose goes into the depot compartment, so nothing changes here:
xyplot(depot~time,data=Res,type='b')

## Only the first dose gets transferred into the transit compartment:
xyplot(trans~time,data=Res,type='b')

## The concentration in the central compartment then becomes:
xyplot(C1~time,data=Res,type='b')




#################################################################################
##                                                                             ##
##  Hands-on bonus assignment: ODE update                                      ##
##  -Change the system of ODEs and examine the results                         ##
##  For example: extend the original model with five transit compartments      ##
##  and use 4 bolus doses in the 1st compartment                               ##
##  -Add an effect compartment and simulate an effect                          ##
##                                                                             ##
#################################################################################

#...

#################################################################################
##                                                                             ##
#################################################################################
