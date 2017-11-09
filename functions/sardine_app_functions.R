#######################################################################################
## Project: EDF Philippines Sardine Simulation Functions
## Purpose: Functions for Sardine web simulation app
## Date: April 16th, 2017
## Author: Tyler Clavelle
#######################################################################################

# Probability of selection function
gill_sel <- function(l,lopt){
  selectivity <- exp(-(l-lopt)^2/(2*1.747^2))
  return(selectivity)
  }

# Beverton-Holt
bev_holt	<- 	function(alpha,beta,b){
  rec <- b * (alpha + beta * b) ^ (-1)
  return(rec)
}

# Negative log-likelihood fucntion
NLL	<- function(PAR.start,Rec, B){
  ALPHA = PAR.start[1]
  BETA  = PAR.start[2]
  Recruits = bev_holt(ALPHA,BETA,B)
  
  return(sum((Recruits - Rec)^2))
}

# Logit model for estimating probability of capture (selectivity)
Logit	<- function(a,b,x){	1 / (1+exp(-(a+b*x)))}

# Function for turning CV on the log-normal scale into standard deviations you can plug into R functions 
log.norm.cv		<- function(cv){
  sigma2	<-	 log(cv^2 + 1)
  return(sqrt(sigma2))
}

