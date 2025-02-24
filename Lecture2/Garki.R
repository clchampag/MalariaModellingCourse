library(deSolve)
library(ggplot2)

m=1.45 # relative density mosquito/human
a=0.5 # biting rate
b=0.097 # probability that a bite successfully infects the human
c=1 # probability that a bite successfully infects the mosquito
g =0.05 #  mortality rate for the mosquito
alpha1=0.002 # rate from I to R (loss of infectivity)
alpha2=0.00019 # rate from R to R2 (develop long term immunity)
r=1/80 # recovery rate in human
delta=1/10000 # mortality rate human
N= 15 #duration of the intrinsic incubation period (in human)
v=10 # duration extrinsic incubation period (in mosquito)
q1=1 # infectivity of I
q2=1 # infectivity of R
q3=0.7 # infectivity of R2

theta_init =c(
  "m"=5,"a"=0.5,"b"=0.097,"g"=0.05,"v"=10, 
                "r1"=0.0023, 'r2'=0.023, "alpha1"=0.002, "alpha2"=0.00019, 
              "N"=15, "delta"=1/10000)

####### 
# WARNING
# This is a simplified version of the Garki model for teaching purposes
# (simplified version in continuous time without delays and without seasonality)
# This is not the full and exact Garki model
Garki<-function(t,x,parms){

  ## state variables
  S = x[1]  #x1
  E = x[2]  #x2
  I = x[3]  #y1
  R = x[4]  #y2
  R2 = x[5] #y3
  S2 = x[6] #x3
  E2 = x[7] #x4
  
  
  ## parameter values
  m = parms["m"]
  a = parms["a"]
  #VC = parms["VC"]
  b = parms["b"]
  g = parms["g"]
  v = parms["v"]
  r1 = parms["r1"]
  r2 = parms["r2"]
  alpha1 = parms["alpha1"]
  alpha2 = parms["alpha2"]
  delta = parms["delta"]
  N = parms["N"]
  beta=1/N
  VC=m*a*a*exp(-g*v)/v

  ##variations
  dS=delta+ R*b*(1-exp(-VC*I))/(exp(b*(1-exp(-VC*I))/r1)-1)-b*(1-exp(-VC*I))*S- delta*S
  dE=b*(1-exp(-VC*I))*S - beta*E - delta*E
  dI=beta*E - alpha1*I - delta*I
  dR=alpha1*I - alpha2*R - R*b*(1-exp(-VC*I))/(exp(b*(1-exp(-VC*I))/r1)-1)-  delta*R
  dR2=alpha2*R + beta*E2 -R2*b*(1-exp(-VC*I))/(exp(b*(1-exp(-VC*I))/r2)-1)- delta*R2
  dS2=R2*b*(1-exp(-VC*I))/(exp(b*(1-exp(-VC*I))/r2)-1)- b*(1-exp(-VC*I))*S2 - delta*S2
  dE2=b*(1-exp(-VC*I))*S2- beta*E2- delta*E2
  res  = c(dS, dE, dI, dR, dR2, dS2, dE2)
  list(res)
}

simulate_Garki=function(parameters){
  
  #parameters
  parms = c(parameters["m"],parameters["a"],
            parameters["b"],#parameters["c"],
            parameters["g"],parameters["v"],
            parameters["r1"],parameters["r2"],
            parameters["alpha1"],parameters["alpha2"],
            parameters["N"],parameters["delta"])
  
  #initial conditions
  init <- c(0.99,0, 0.01, 0, 0, 0, 0)  
  
  #simulation
  temps <- seq(0,10000) 
  solveRM <- lsoda(y =init, times=temps, func = Garki, 
                   parms = parms)
  solutionRM=as.data.frame(solveRM)
  names(solutionRM)=c("time","S", "E", "I", "R", "R2", "S2", "E2")
  
  solutionRM$obs_I=solutionRM$I+solutionRM$R+solutionRM$R2* parameters["q3"]
  solutionRM$positive=solutionRM$I+solutionRM$R+solutionRM$R2
  return(solutionRM)
}


simul=simulate_Garki(parameters = theta_init)

VC=theta_init["m"]*theta_init["a"]*theta_init["a"]*exp(-theta_init["g"]*theta_init["v"])/theta_init["v"]
R0=VC*theta_init["b"]/(theta_init["alpha1"]*theta_init["delta"])
ggplot(simul)+
  geom_line(aes(x=time, y=positive))+ ylim(0,1)+
  labs(y="Proportion of\ninfected individuals", x="Time (days)")+
  ggtitle(paste0("VC=",round(VC, 3), "\n R0=", round(R0, 3)))

