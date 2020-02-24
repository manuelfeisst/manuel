#####Input:

  #c_old:   number of responders in the historical control arm
  #t_old:   number of responders in the historical treatment arm
  #nc_old:  number of patients in the historical control arm
  #nt_old:  number of patients in the historical treatment arm
  #nc:      number of patients in the new control arm
  #nt:      number of patients in the new tretment arm
  #delta:   factor that determines the amount of borrowed historical information
  #pi:      true control proportion
  #ES:      Effect size pi_T-pi_C
  #parts:   The number of parts in which the area of the true control proportion is divided (rec.:100)
  #alpha:   The significance niveau
  #power:   The power to detect the effect
  #gamma:   Parameter of the Berger and Boos procedure


#############################################################################
############# General functions
#############################################################################

require(compiler)


#####Indicator function
#############################################################################
Indicator<-function(x,min,max){
  if(min<=x && x<=max){
    y<-1
    
  }
  else{
    y<-0
  }
  return(y)
}

##### Function for calculating the true type I error
#############################################################################
truealpha<-function(c_old,t_old,nc_old,nt_old,nc,nt,delta,pi){
  p<-matrix(c(rep(0,nc*nt)),nrow=nc)
  for (i in 1:(nc-1)){
    for(j in 1:(nt-1)){
      tab<-matrix(c(c_old*delta+i,nc_old*delta+nc-c_old*delta-i,t_old*delta+j,nt_old*delta+nt-t_old*delta-j),nrow=2)
      p[i,j]<-dbinom(i,nc,pi)*dbinom(j,nt,pi)*(1-Indicator(chisq.test(tab,correct=FALSE)$statistic,-1,qchisq(0.95,1)))
    }
  }
  
  g3<-sum(p)
  return(g3)
}
truealpha<-Vectorize(truealpha)
truealpha<-cmpfun(truealpha)



##### Function for calculating the true power
#############################################################################
truepower<-function(c_old,t_old,nc_old,nt_old,nc,nt,delta,pi,ES){
  p<-matrix(c(rep(0,nc*nt)),nrow=nc)
  for (i in 1:(nc-1)){
    for(j in 1:(nt-1)){
      tab<-matrix(c(c_old*delta+i,nc_old*delta+nc-c_old*delta-i,t_old*delta+j,nt_old*delta+nt-t_old*delta-j),nrow=2)
      p[i,j]<-dbinom(i,nc,pi)*dbinom(j,nt,pi+ES)*(1-Indicator(chisq.test(tab,correct=FALSE)$statistic,-1,qchisq(0.95,1)))
    }
  }
  
  g3<-sum(p)
  return(g3)
}
truepower<-Vectorize(truepower)
truepower<-cmpfun(truepower)






#############################################################################
#### Local Procedure
#############################################################################


#### Function that calculates delta* for a fixed true control proportion
#############################################################################
sc<-function(c_old,t_old,nc_old,nt_old,nc,nt,pi,ES,gamma,alpha){
  
  thealpha<-alpha-gamma
  
  
  #Nested intervals procedure
  delta_opt<-0
  delta<-0:2/2
  alpha1<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta[1],pi)
  alpha2<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta[2],pi)
  alpha3<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta[3],pi)
  
  #1. Step
  if(thealpha>alpha3)
  {delta_opt<-1}
  
  else{
    
    if(alpha2>thealpha)
    {delta1<-0
    delta2<-0.25
    delta3<-0.5
    alpha_l<-alpha1
    alpha_r<-alpha2}
    
    else
    {delta1<-0.5
    delta2<-0.75
    delta3<-1
    alpha_l<-alpha2
    alpha_r<-alpha3}
    
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta2,pi)
  
    
    
    #2.Step
    if(alpha_m<thealpha){
      delta_l<-delta2
      delta_m<-round((delta2+delta3)/2,2)
      delta_r<-delta3
    }
    
    else{
      delta_l<-delta1
      delta_m<-round((delta1+delta2)/2,2)
      delta_r<-delta2
    }
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta_m,pi)
    
    
    #3.Step
    if(alpha_m<thealpha){
      delta1<-delta_m
      delta2<-round((delta_m+delta_r)/2,2)
      delta3<-delta_r
    }
    
    else{
      delta1<-delta_l
      delta2<-round((delta_l+delta_m)/2,2)
      delta3<-delta_m
    }
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta2,pi)
    
    #4.Step
    if(alpha_m<thealpha){
      delta_l<-delta2
      delta_m<-round((delta2+delta3)/2,2)
      delta_r<-delta3
    }
    
    else{
      delta_l<-delta1
      delta_m<-round((delta1+delta2)/2,2)
      delta_r<-delta2
    }
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta_m,pi)
    
    #5.Step
    if(alpha_m<thealpha){
      delta1<-delta_m
      delta2<-round((delta_m+delta_r)/2,2)
      delta3<-delta_r
    }
    
    else{
      delta1<-delta_l
      delta2<-round((delta_l+delta_m)/2,2)
      delta3<-delta_m
    }
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta2,pi)
    
    #6.Step
    if(alpha_m<thealpha){
      delta_l<-delta2
      delta_m<-round((delta2+delta3)/2,2)
      delta_r<-delta3
    }
    
    else{
      delta_l<-delta1
      delta_m<-round((delta1+delta2)/2,2)
      delta_r<-delta2
    }
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta_m,pi)
   
    #7.Step
    if(alpha_m<thealpha){
      delta1<-delta_m
      delta2<-round((delta_m+delta_r)/2,2)
      delta3<-delta_r
    }
    
    else{
      delta1<-delta_l
      delta2<-round((delta_l+delta_m)/2,2)
      delta3<-delta_m
    }
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta2,pi)
    
    #8.Step
    if(alpha_m<thealpha){
      delta_l<-delta2
      delta_m<-round((delta2+delta3)/2,2) 
      delta_r<-delta3
    }
    
    else{
      delta_l<-delta1
      delta_m<-(floor(100*(delta1+delta2)/2))/100
      delta_r<-delta2
    }
    
    delta_opt<-max(0,delta_m)
  }
  
    return(list(delta_opt=delta_opt))
 
  
}



### Function that calculates delta* 
#############################################################################
localpro<-function(c_old,t_old,nc_old,nt_old,nc,nt,ES,parts,gamma,alpha){
  cpu<-0
  cpo<-0
  delta<-0
  power0<-0
  power1<-0
  power2<-0
  alpha0<-0
  alpha1<-0
  alpha2<-0
  power<-0
  delta_min<-0
  
  
  #Calculate the Pearson Clopper Confidence Interval
  cpu<-qbeta(gamma/2,c_old,nc_old-c_old+1) # 	pu = BETAINV(??/2;k;n-k+1)
  cpo<-qbeta(1-gamma/2,c_old+1,nc_old-c_old) #po = BETAINV(1-??/2;k+1;n-k) 	
  
  #Seperate the interval in parts
  steps<-0
  step<-(cpo-cpu)/parts
  for(i in 1:parts){
    steps[i]<-cpu+(i)*step
  }
  steps<-c(cpu,steps)
  
  #Calculate delta for every pi in the Confidence interval
  for(i in 1:length(steps)){
    fit<-sc(c_old,t_old,nc_old,nt_old,nc,nt,steps[i],ES,gamma,alpha)
    delta[i]<-fit$delta_opt
  }
  
  delta_max<-min(delta)
  
  alpha<-0
  power<-0
  for(i in 1:length(steps)){
    alpha1[i]<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta_max,steps[i])
    power1[i]<-truepower(c_old,t_old,nc_old,nt_old,nc,nt,delta_max,steps[i],ES)
    alpha0[i]<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,0,steps[i])
    power0[i]<-truepower(c_old,t_old,nc_old,nt_old,nc,nt,0,steps[i],ES)  
  }
  
  return(list(alpha1=alpha1,power1=power1,delta_max=delta_max,delta=delta,power0=power0,
              alpha0=alpha0,steps=steps))
  
  
}


##### Function for the sample size calculation for the local procedure
#############################################################################
samplesizerl<-function(c_old,t_old,nc_old,nt_old,pi,ES,parts,alpha,power,gamma){
  
  #1.SampleSizeCalculation initial
  library(pwr)
  w<-ES.h(pi+ES,pi)
  n<-ceiling(pwr.2p.test(h=w,power=power,sig.level=alpha)$n)
  nc<-n
  nt<-n
  
  #2.Algorithm
  #Calculate Confidence interval
  cpu<-qbeta(gamma/2,c_old,nc_old-c_old+1) # 	pu = BETAINV(??/2;k;n-k+1)
  cpo<-qbeta(1-gamma/2,c_old+1,nc_old-c_old) #po = BETAINV(1-??/2;k+1;n-k) 	
  
  # Steps in the algortihm
  if(pi<cpu || pi>cpo){print("Error: pi is not in the confidence interval")}
  else{
    
    steps<-1
    fit<-localpro(c_old,t_old,nc_old,nt_old,nc,nt,ES,parts,gamma,alpha)
    delta<-fit$delta_max
    pi<-max(fit$steps[which(round(fit$steps,digits=2)==pi)])
    power1<-fit$power1[which(fit$steps==pi)]
    
    while(power1>power){
      power1<-truepower(c_old,t_old,nc_old,nt_old,nc,nt,delta,pi,ES)
      nc<-nc-1
      nt<-nt-1
    }
    nc<-nc+1
    nt<-nt+1  
    
    fit<-localpro(c_old,t_old,nc_old,nt_old,nc,nt,ES,parts,gamma,alpha)
    delta1<-fit$delta_max
    power1<-fit$power1[which(fit$steps==pi)]
    
    while(delta1>delta){
      steps<-steps+1
      while(power1>power){
        power1<-truepower(c_old,t_old,nc_old,nt_old,nc,nt,delta1,pi,ES)
        nc<-nc-1
        nt<-nt-1
      }
      delta<-delta1
      nc<-nc+1
      nt<-nt+1 
      fit<-localpro(c_old,t_old,nc_old,nt_old,nc,nt,ES,parts,gamma,alpha)
      delta1<-fit$delta_max
      
    }
    nc+1
    nt+1
    n_won1<-n-nc
    n_won2<-n-nt
    
    return(list(fit=fit,n_won1=n_won1,n_won2=n_won2,n=n,nc=nc,nt=nt,delta=delta,delta1=delta1,steps=steps))}
}

###Output
#fit:           alpha and power values at each pi
#n_won1;n_won2: saved sample size in the control arm and the treatment arm, respectively
#n:             initial sample size  
#nc;nt:         new reduced sample size
#delta:         delta*
#delta1:        last delta in the algorithm
#steps:         values of pi where delta* was calculated

#############################################################################
##### Global procedure
#############################################################################

### Function that calculates delta* for a fixed true control proportion pi
#############################################################################
sc1<-function(c_old,t_old,nc_old,nt_old,nc,nt,pi,ES,alpha){
  
  thealpha<-alpha
  delta_opt<-0
  
  #Nested intervals procedure 
  delta<-0:2/2
  alpha1<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta[1],pi)
  alpha2<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta[2],pi)
  alpha3<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta[3],pi)
  

  #1. Step
  if(thealpha>alpha3)
  {delta_opt<-1}
  
  else{
    
    if(alpha2>thealpha)
    {delta1<-0
    delta2<-0.25
    delta3<-0.5
    alpha_l<-alpha1
    alpha_r<-alpha2}
    
    else
    {delta1<-0.5
    delta2<-0.75
    delta3<-1
    alpha_l<-alpha2
    alpha_r<-alpha3}
    
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta2,pi)
    
    
    #2.Step
    if(alpha_m<thealpha){
      delta_l<-delta2
      delta_m<-round((delta2+delta3)/2,2)
      delta_r<-delta3
    }
    
    else{
      delta_l<-delta1
      delta_m<-round((delta1+delta2)/2,2)
      delta_r<-delta2
    }
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta_m,pi)
    
    
    #3.Step
    if(alpha_m<thealpha){
      delta1<-delta_m
      delta2<-round((delta_m+delta_r)/2,2)
      delta3<-delta_r
    }
    
    else{
      delta1<-delta_l
      delta2<-round((delta_l+delta_m)/2,2)
      delta3<-delta_m
    }
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta2,pi)

    
    #4.Step
    if(alpha_m<thealpha){
      delta_l<-delta2
      delta_m<-round((delta2+delta3)/2,2)
      delta_r<-delta3
    }
    
    else{
      delta_l<-delta1
      delta_m<-round((delta1+delta2)/2,2)
      delta_r<-delta2
    }
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta_m,pi)
    
    #5.Step
    if(alpha_m<thealpha){
      delta1<-delta_m
      delta2<-round((delta_m+delta_r)/2,2)
      delta3<-delta_r
    }
    
    else{
      delta1<-delta_l
      delta2<-round((delta_l+delta_m)/2,2)
      delta3<-delta_m
    }
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta2,pi)
    
    
    #6.Step
    if(alpha_m<thealpha){
      delta_l<-delta2
      delta_m<-round((delta2+delta3)/2,2)
      delta_r<-delta3
    }
    
    else{
      delta_l<-delta1
      delta_m<-round((delta1+delta2)/2,2)
      delta_r<-delta2
    }
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta_m,pi)
  
    
    #7.Step
    if(alpha_m<thealpha){
      delta1<-delta_m
      delta2<-round((delta_m+delta_r)/2,2)
      delta3<-delta_r
    }
    
    else{
      delta1<-delta_l
      delta2<-round((delta_l+delta_m)/2,2)
      delta3<-delta_m
    }
    
    alpha_m<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta2,pi)
    
    #8.Step
    if(alpha_m<thealpha){
      delta_l<-delta2
      delta_m<-round((delta2+delta3)/2,2) 
      delta_r<-delta3
    }
    
    else{
      delta_l<-delta1
      delta_m<-(floor(100*(delta1+delta2)/2))/100
      delta_r<-delta2
    }
    
    delta_opt<-max(0,delta_m)
  }
  
    return(list(delta_opt=delta_opt))
}



### Function that calculates delta* 
#############################################################################
globalpro<-function(c_old,t_old,nc_old,nt_old,nc,nt,ES,parts,alpha){
  cpu<-0
  cpo<-0
  delta<-0
  power0<-0
  power1<-0
  power2<-0
  alpha0<-0
  alpha1<-0
  alpha2<-0
  power<-0
  delta_min<-0
  
  
  #seperate the values of pi in parts
  steps<-(0:parts/parts)[2:parts]
  
  #Calculate delta for each pi
  for(i in 1:length(steps)){
    fit<-sc1(c_old,t_old,nc_old,nt_old,nc,nt,steps[i],ES,alpha)
    delta[i]<-fit$delta_opt
  }
  
  delta_max<-min(delta)
  
  alpha<-0
  power<-0
  for(i in 1:length(steps)){
    alpha1[i]<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,delta_max,steps[i])
    power1[i]<-truepower(c_old,t_old,nc_old,nt_old,nc,nt,delta_max,steps[i],ES)
    alpha0[i]<-truealpha(c_old,t_old,nc_old,nt_old,nc,nt,0,steps[i])
    power0[i]<-truepower(c_old,t_old,nc_old,nt_old,nc,nt,0,steps[i],ES)  
  }
  
  return(list(alpha1=alpha1,power1=power1,delta_max=delta_max,delta=delta,power0=power0,
              alpha0=alpha0,steps=steps))
  
  
}


##### Function for the sample size calculation for the global procedure
#############################################################################
samplesizer<-function(c_old,t_old,nc_old,nt_old,pi,ES,parts,alpha,power){
  
  #1.SampleSizeCalculation inital
  library(pwr)
  w<-ES.h(pi+ES,pi)
  n<-ceiling(pwr.2p.test(h=w,power=power,sig.level=alpha)$n)
  nc<-n
  nt<-n
  
  #2.Algorithm
  steps<-1
  fit<-globalpro(c_old,t_old,nc_old,nt_old,nc,nt,ES,parts,alpha)
  delta<-fit$delta_max
  power1<-fit$power1[pi*100]
  while(power1>power){
    power1<-truepower(c_old,t_old,nc_old,nt_old,nc,nt,delta,pi,ES)
    nc<-nc-1
    nt<-nt-1
  }
  nc<-nc+1
  nt<-nt+1  
  
  fit<-globalpro(c_old,t_old,nc_old,nt_old,nc,nt,ES,parts,alpha)
  delta1<-fit$delta_max
  power1<-fit$power1[pi*100]
  
  while(delta1>delta){
    steps<-steps+1
    while(power1>power){
      power1<-truepower(c_old,t_old,nc_old,nt_old,nc,nt,delta1,pi,ES)
      nc<-nc-1
      nt<-nt-1
    }
    delta<-delta1
    nc<-nc+1
    nt<-nt+1 
    fit<-globalpro(c_old,t_old,nc_old,nt_old,nc,nt,ES,parts,alpha)
    delta1<-fit$delta_max
    
  }
  nc+1
  nt+1
  n_won1<-n-nc
  n_won2<-n-nt
  
  return(list(fit=fit,n_won1=n_won1,n_won2=n_won2,n=n,nc=nc,nt=nt,delta=delta,delta1=delta1,steps=steps))
}




###Input:
#c_old:   number of responders in the historical control arm
#t_old:   number of responders in the historical treatment arm
#nc_old:  number of patients in the historical control arm
#nt_old:  number of patients in the historical treatment arm
#nc:      number of patients in the new control arm
#nt:      number of patients in the new tretment arm
#delta:   factor that determines the amount of borrowed historical information
#pi:      true control proportion
#ES:      Effect size pi_T-pi_C
#parts:   The number of parts in which the area of the true control proportion is divided (rec.:100)
#alpha:   The significance niveau
#power:   The power to detect the effect
#gamma:   Parameter of the Berger and Boos procedure

###Output
#fit:                   # alpha1 and power1: alpha and power values after borrowing for each pi
                        # delta_max, delta: delta* and the max delta for each pi
                        # power0, alph0: alpha and power values without borrowing for each pi
#n_won1;n_won2: saved sample size in the control arm and the treatment arm, respectively
#n:             initial sample size  
#nc;nt:         new reduced sample size
#delta:         delta*
#delta1:        last delta of the algorithm
#steps:         values of pi where delta* was calculated

#############################################################################
###### Example from the Manuscript (Clincial trial example from section 4)
#############################################################################

##########Lokal procedure:
#############################################################################
fit1<-samplesizerl(c_old = 10,t_old = 16,nc_old = 44,nt_old = 43,pi = 0.23,ES = 0.14,parts = 100,alpha = 0.05,power = 0.8,gamma = 0.0001)

##########Global procedure:
#############################################################################
fit2<-samplesizer(c_old = 10,t_old = 16,nc_old = 44,nt_old = 43,pi = 0.23,ES = 0.14,parts = 100,alpha = 0.05,power = 0.8)



#### Time for each calculation ~ 3h