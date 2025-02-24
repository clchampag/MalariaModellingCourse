library(deSolve)
library(ggplot2)


theta_before=c("r"=1/100,"V"=0.1,"P"=0.9, "Q"=0.3,"eta"=1/(5*365))
theta_control=c("r"=1/100,"V"=0.005,"P"=0.9, "Q"=0.3,"eta"=1/(5*365))
duration_control=6000

##Model from Smith et al. "A sticky situation"
StickyModel<-function(t,x,parms){
  # version from Smith and McKenzie, 2004
  
  ## state variables
  X = x[1]
  I = x[2]
  
  ## parameter values
  r = parms["r"]
  V = parms["V"]
  P = parms["P"]
  Q = parms["Q"]
  eta = parms["eta"]
  
  ##variations
  dX=V*X*(1-P+(P-Q)*I)*(1-X) - r*X
  dI= eta*(X-I)
  res  = c(dX, dI)
  list(res)
}

simulate_sticky=function(parameters, init, max_time){
  
  #parameters
  parms = c(parameters["r"],parameters["V"],
            parameters["P"],parameters["Q"],
            parameters["eta"])
  
  #simulation
  temps <- seq(0,max_time) 
  solveRM <- lsoda(y =init, times=temps, func = StickyModel, 
                   parms = parms)
  solutionRM=as.data.frame(solveRM)
  names(solutionRM)=c("time","X","I")
  
  return(solutionRM)
}


calculate_X_eq=function(parameters){
  r = parameters["r"]
  V = parameters["V"]
  P = parameters["P"]
  Q = parameters["Q"]
  radical=(2*P-1-Q)^2 - 4*(Q-P)*(1-P-r/V)
  numerator=-(2*P-1-Q)-sqrt(radical)
  return(as.numeric(numerator/2/(Q-P)))
}

simul_before=simulate_sticky(parameters = theta_before, init=c(calculate_X_eq(theta_before), calculate_X_eq(theta_before)), max_time = 2000)
simul_control=simulate_sticky(parameters = theta_control, init=as.numeric(simul_before[2000,c(2,3)]), max_time = duration_control)
simul_after=simulate_sticky(parameters = theta_before, init=as.numeric(simul_control[duration_control,c(2,3)]), max_time = 5000)
simul_control$time=simul_control$time+2000
simul_after$time=simul_after$time+2000+duration_control


ggplot(rbind(simul_before,simul_control, simul_after))+
  geom_line(aes(x=time, y=X))+ ylim(0,1)+
  geom_line(aes(x=time, y=I), col="red")+ 
  labs(y="Proportion of\ninfected individuals", x="Time (days)")

