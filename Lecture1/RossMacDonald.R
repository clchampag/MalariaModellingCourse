library(deSolve)
library(ggplot2)

theta_init =c("m"=10,"a"=0.5,"b"=0.09, "c"=1,"p"=0.95,"v"=12, "r"=1/80)
#theta_init =c("m"=10,"a"=0.25,"b"=0.09, "c"=1,"p"=0.95,"v"=12, "r"=1/80)
#theta_init =c("m"=100,"a"=0.5,"b"=0.09, "c"=1,"p"=0.9,"v"=12, "r"=1/80)


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

simulate_RM=function(parameters){
  
  #parameters
  parms = c(parameters["m"],parameters["a"],
            parameters["b"],parameters["c"],
            parameters["p"],parameters["v"],
            parameters["r"])

  #initial conditions
  init <- c(0.05,0)  
  
  #simulation
  temps <- seq(0,500) 
  solveRM <- lsoda(y =init, times=temps, func = RossMacDonald, 
                    parms = parms)
  solutionRM=as.data.frame(solveRM)
  names(solutionRM)=c("time","Hi","Vi")

  return(solutionRM)
}

calculate_R0=function(parameters){
  return(c("R0"=as.numeric(parameters["m"]*(parameters["a"]^2)*parameters["b"]*parameters["c"]*(parameters["p"]^parameters["v"])/parameters["r"]/(-log(parameters["p"])))))
}

simul=simulate_RM(parameters = theta_init)

ggplot(simul)+
  geom_line(aes(x=time, y=Hi))+ ylim(0,1)+
  labs(y="Proportion of\ninfected individuals", x="Time (days)")+
  ggtitle(paste0("R0=",round(calculate_R0(theta_init),2)))
