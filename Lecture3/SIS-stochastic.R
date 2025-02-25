library(deSolve)
library(ggplot2)

lambda=0.05
r=1/200 # recovery rate in human
N=100


theta_init =c("lambda"=lambda,"r"=r, "N"=N)


##SIR model
SIS_deterministic<-function(t,x,parms){
  # version from Smith and McKenzie, 2004
  
  ## state variables
  Hi = x[1]
  
  ## parameter values
  lambda = parms["lambda"]
  r = parms["r"]
  N= parms["N"]
  
  ##variations
  dHi=lambda*(N-Hi)*Hi/N - r*Hi
  res  = c(dHi)
  list(res)
}

simulate_SIS=function(parameters){
  
  #parameters
  parms = c(parameters["lambda"],parameters["r"], parameters["N"])
  
  #initial conditions
  init <- c(0.005*parameters["N"])  
  
  #simulation
  temps <- seq(0,500) 
  solveSIS <- lsoda(y =init, times=temps, func = SIS_deterministic, 
                   parms = parms)
  solutionSIS=as.data.frame(solveSIS)
  names(solutionSIS)=c("time","Hi")
  
  return(solutionSIS)
}


simul=simulate_SIS(parameters = theta_init)

ggplot(simul)+
  geom_line(aes(x=time, y=Hi))+ 
  labs(y="Proportion of\ninfected individuals", x="Time (days)")




simulate_SIS_stochastic=function(lambda, r, N, init, steps=1000){
  
  #initial conditions
  solution=data.frame(I=rep(NA,steps ), time=rep(NA,steps))
  solution[1,]=c(0.005*N,1)

  U1=runif(steps)
  U2=runif(steps)
  
  for(t in 2:steps){
    current_i=solution$I[t-1]
    rate_infect=lambda*current_i*(N-current_i)/N
    rate_recovery=r*current_i
    alpha0=rate_infect+rate_recovery
    if(alpha0>0){
      solution$time[t]=solution$time[t-1]+(1/alpha0)*log(1/U1[t])
      if(rate_infect/alpha0>U2[t]){
        solution$I[t]=solution$I[t-1]+1
      } else {
        solution$I[t]=solution$I[t-1]-1
      }
    } else {
      solution$I[t]=0
      solution$time[t]=solution$time[t-1]+1
    }

  }

  return(solution)
}


simul_sto=simulate_SIS_stochastic(lambda, r, N=N, init=50, steps=1000)
simul_sto2=simulate_SIS_stochastic(lambda, r, N=N, init=50, steps=1000)
simul_sto3=simulate_SIS_stochastic(lambda, r, N=N, init=50, steps=1000)


ggplot(simul)+
  geom_line(aes(x=time, y=Hi))+ 
  geom_line(data=simul_sto,aes(x=time, y=I, col="sto"))+ 
  geom_line(data=simul_sto2,aes(x=time, y=I, col="sto2"))+ 
  geom_line(data=simul_sto3,aes(x=time, y=I, col="sto3"))+ 
  labs(y="Proportion of\ninfected individuals", x="Time (days)")
