## load the required libraries
library(xpose.nlmixr)
library(nlmixr)
library(RxODE)
library(lattice)
library(data.table)

#################################################################################
##                                                                             ##
## nlmixr analysis Part 1                                                      ##
##                                                                             ##
## Warfarin population PK using SAEM estimation                                ##
##                                                                             ##
#################################################################################

## read in the Warfarin PK-only data set using data.table syntax (fast and efficient!)
PKdata <- fread("warfarin_PK.csv")

## look at data structure and identify what each column is for
View(PKdata)

## lattice spaghetti plot
xyplot(
  DV ~ TIME,         #plot DV by TIME
  groups = ID,       #make separate curves by ID
  data = PKdata,     #use PKdata
  type = 'l',        #only plot a line"
  lwd = 2,           #make the lines '2' wide
  xlab = "Time (h)", #x-axis label
  ylab = "Warfarin (mg/L)" #y-axis label
)

## Define a first order-absorption linear elimination model using a solved system solution

One.comp.KA.solved <- function() {
  ini({
    # Where initial conditions/variables are specified
    lka  <- log(1.15)  #log ka (1/h)
    lcl  <- log(0.135) #log Cl (L/h)
    lv   <- log(8)     #log V (L)
    prop.err <- 0.15   #proportional error (SD/mean)
    add.err  <- 0.6    #additive error (mg/L)
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

## Check the model and some of the assumptions made by nlmixr 
##  note assumption that AMT goes into CMT=1 is not shown
nlmixr(One.comp.KA.solved) 

## estimate parameters using nlmixr:
fitOne.comp.KA.solved_S <-
  nlmixr(
    One.comp.KA.solved,    #the model definition
    PKdata,                #the data set
    est = "saem",          #the estimation algorithm (SAEM)
                             #the SAEM minimisation options:
    saemControl(nBurn = 200, #200 SAEM burn-in iterations (the default)
                nEm   = 300, #300 EM iterations (the default)
                print = 50),
    #only print every 50th estimation step (default=1 which gives endless output)
    tableControl(cwres = TRUE) #calculates NONMEM-style conditional weighted residuals for diagnostics
  )

## results are stored in the nlmixr object and can be viewed:
fitOne.comp.KA.solved_S

## and saved for future use or reference:
save(fitOne.comp.KA.solved_S, file = "fitOne.comp.KA.solved_S.Rdata")
#load(file = "fitOne.comp.KA.solved_S.Rdata")

## and for SAEM, convergence can be checked using a parameter trace plot:
traceplot(fitOne.comp.KA.solved_S)

## the nlmixr object can be transformed into an xpose object to allow diagnostics with the new xpose package
## the link between nlmixr and xpose is provided by the xpose.nlmixr package
## only xpose_data_nlmixr is from xpose.nlmixr, all further commands are from the xpose package
xpdb.1s <- xpose_data_nlmixr(fitOne.comp.KA.solved_S)

## this can also be used to generate trace plots (parameters vs iterations:)
prm_vs_iteration(xpdb.1s)
## to remove the path to the script from the plot use:
prm_vs_iteration(xpdb.1s,caption=NULL)

## and many types of diagnostic plots (see cheatsheet)

## dv vs cpred plot:
dv_vs_pred(xpdb.1s, 
           caption = NULL)

# by default model typical predictions (PRED) are assigned to CPRED (conditional population predictions):
list_vars(xpdb.1s)
# if you want this to be PRED instead, these can be updated, either using 'standard' syntax:
xpdb.1s<-set_var_types(xpdb.1s,pred = 'PRED')
# or using magrittr piping type code:
# xpdb.1s<-xpdb.1s %>% set_var_types(pred = 'PRED')

## dv vs pred plot:
dv_vs_pred(xpdb.1s, 
           caption = NULL)

## dv vs ipred plot:
dv_vs_ipred(xpdb.1s,        #the xpose object
            caption = NULL) #if not NULL provides the directory where this was run

## CWRES vs time:
res_vs_idv(xpdb.1s,           #the xpose object
           res = "CWRES",     #examine conditional weighted residuals
           idv = "TIME",      #as a function of time
           caption = NULL)    #if not NULL provides the directory where this was run

#Absolute values of individual weighted residual vs time
absval_res_vs_idv(xpdb.1s,        #the xpose object
                  res = "IWRES",  #examine absolute values (absval) of individual weighted residuals
                  idv = "TIME",   #as a function of time
                  caption = NULL) #if not NULL provides the directory where this was run

## nlmixr comes with its built-in vpc functionality that uses Ron Keizer's vpc package
## see the cheatsheet for further options
vpc_ui(
  fitOne.comp.KA.solved_S,        #the nlmixr object
  n = 500,                        #number of trials simulated using estimated 
                                  # parameters and study sampling structure
  show = list(obs_dv = TRUE),     #additional items to show, like the observations
  xlab = "Time (h)",              #x-axis label
  ylab = "Concentration (mg/L)",  #y-axis label
  title = "VPC for first order absorption PopPK model"
)


## or with a log y-axis starting at 0.5
vpc_ui(
  fitOne.comp.KA.solved_S,
  n = 500,
  show = list(obs_dv = TRUE),
  xlab = "Time (h)",
  ylab = "Concentration (mg/L)",
  title = "VPC for first order absorption PopPK model\nwith log y-axis",
  log_y = TRUE,            #to request a log y-axis
  log_y_min = 0.5          #starting at 0.5
)


## Individual fits can be generated using using xpose:
ind_plots(xpdb.1s,caption = NULL,ncol = 4,nrow = 4)
## ...use the arrows in the plot window to examine the earlier curves

## Individual fits can also be generated using augPred (augmented predictions)
## that provides smooth profiles by interpolating the predictions between observations:
plot(augPred(fitOne.comp.KA.solved_S))
## ...use the arrows in the plot window to examine the earlier curves

#or the augPred output can be formatted to your liking for instance using lattice:
indivpk<-augPred(fitOne.comp.KA.solved_S)
## look at individual PK structure
View(indivpk) 
table(indivpk$ind)
#Population Individual   Observed 
#1824       1824        251 

## specify array of colours for curves
nlmixCOLS <- c("#28466A","#8DB6CD","#B40000")  

xyplot(
  values~time|id,          ## plot the variable 'values' by time and make a separate panel for each id
  data=indivpk,            ## data source
  groups=ind,              ## make separate curves for ind that indicates Observed data, 
                           ## Indivdual predictions and Population predictions
  layout=c(8,4),           ## arrange as 8 columns and 4 rows
  type=c("l","l","p"),     ## represent these three by a line, a line and only markers (l=line, p=points)
  col=nlmixCOLS[c(2,1,3)], ## colours for each curve
  cex=c(0.1,0.1,1),        ## character size for the markers
  lwd=c(2,2,0.1),          ## line width of the lines
  pch=19,                  ## use closed circles as marker
  xlab="Time (hr)\n",      ## x-axis label
  ylab="Warfarin (mg/L)",  ## y-axis label
  as.table=TRUE,           ## have the first plot at the top left (otherwise plots start at the lower left corner)
  scales=list(alternating=1),  ## have axis labels at left and bottom (and not alternating)
  main="First order-absorption linear elimination", ## title for plot
  auto.key=list(adj=1,col=nlmixCOLS[c(2,1,3)],columns=3,space="bottom",rectangles=FALSE,points=FALSE) ## key for curves
)


#################################################################################
##                                                                             ##
## nlmixr analysis Part 2                                                      ##
##                                                                             ##
## Warfarin population PK using FOCEi and ODEs                                 ##
##                                                                             ##
#################################################################################


One.comp.KA.ODE <- function() {
  ini({
    # Where initial conditions/variables are specified
    lka  <- log(1.15)  #log ka (1/h)
    lcl  <- log(0.135) #log Cl (L/h)
    lv   <- log(8)     #log V (L)
    prop.err <- 0.15   #proportional error (SD/mean)
    add.err  <- 0.6    #additive error (mg/L)
    eta.ka ~ 0.5   #IIV ka
    eta.cl ~ 0.1   #IIV cl
    eta.v  ~ 0.1   #IIV v
  })
  model({
    # Where the model is specified
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv + eta.v)
    ka <- exp(lka + eta.ka)
    ## ODE example
    d/dt(depot)   = -ka * depot
    d/dt(central) =  ka * depot - (cl/v) * central
    ## Concentration is calculated
    cp = central/v
    ## And is assumed to follow proportional and additive error
    cp ~ prop(prop.err) + add(add.err)
  })
}


## estimate parameters using nlmixr/FOCEI:
fitOne.comp.KA.ODE_F <-
  nlmixr(One.comp.KA.ODE,          #the model definition
         PKdata,                   #the data set
         est = "focei",            #the estimation algorithm (FOCEi)
                                   #FOCEi options:
         foceiControl(print = 5))  #only print every 5th estimation step

## results are stored in the nlmixr object and can be viewed:
fitOne.comp.KA.ODE_F

## and saved for future use or reference:
save(fitOne.comp.KA.ODE_F, file = "fitOne.comp.KA.ODE_F.Rdata")

#################################################################################
##                                                                             ##
##  Hands-on assignment: nlmixr model development                              ##
##  Examine the GOF plots and implement models with                            ##
##  one or transit compartments                                                ##
##                                                                             ##
##  Compare vpcs of alternatives and compare OFVs:                             ##
##                                                                             ##
##  fitPK001$OBJF-fitPK002$OBJF                                                ##
##                                                                             ##
#################################################################################

#...

#################################################################################
##                                                                             ##
#################################################################################
