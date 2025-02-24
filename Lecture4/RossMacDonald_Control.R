library(deSolve)
library(ggplot2)

theta_before=c("m"=5,"a"=0.5,"b"=0.09, "c"=1,"p"=0.95,"v"=12, "r"=1/80)
theta_control=c("m"=0.05,"a"=0.5,"b"=0.09, "c"=1,"p"=0.95,"v"=12, "r"=1/80)


##SIR model
RossMacDonald<-function(t,x,parms){
  # version from Smith and McKenzie, 2004
  
  ## state variables
  Hi = x[1]
  Vi = x[2]
  
  ## parameter values
  m = parms["m"]
  a = parms["a"]
  b = parms["b"]
  c = parms["c"]
  p = parms["p"]
  v = parms["v"]
  r = parms["r"]
  g= -log(p)
  
  ##variations
  dHi=m*a*b*(1-Hi)*Vi - r*Hi
  dVi= a*c*Hi*(exp(-g*v)-Vi) - g*Vi
  res  = c(dHi, dVi)
  list(res)
}

simulate_RM=function(parameters, init,max_time){
  
  #parameters
  parms = c(parameters["m"],parameters["a"],
            parameters["b"],parameters["c"],
            parameters["p"],parameters["v"],
            parameters["r"])

  #simulation
  temps <- seq(0,max_time) 
  solveRM <- lsoda(y =init, times=temps, func = RossMacDonald, 
                   parms = parms)
  solutionRM=as.data.frame(solveRM)
  names(solutionRM)=c("time","Hi","Vi")
  
  return(solutionRM)
}

calculate_R0=function(parameters){
  return(c("R0"=as.numeric(parameters["m"]*(parameters["a"]^2)*parameters["b"]*parameters["c"]*(parameters["p"]^parameters["v"])/parameters["r"]/(-log(parameters["p"])))))
}

simul_before=simulate_RM(parameters = theta_before, init=c(0.5, 0.5),max_time = 2000)
simul_control=simulate_RM(parameters = theta_control, init=as.numeric(simul_before[2000,c(2,3)]), max_time = 5000)
simul_after=simulate_RM(parameters = theta_before, init=as.numeric(simul_control[2000,c(2,3)]), max_time = 5000)
simul_control$time=simul_control$time+2000
simul_after$time=simul_after$time+7000


calculate_R0(theta_before)
calculate_R0(theta_control)

ggplot(rbind(simul_before,simul_control, simul_after))+
  geom_line(aes(x=time, y=Hi))+ ylim(0,1)+xlim(100,NA)+
  labs(y="Proportion of\ninfected individuals", x="Time (days)")

