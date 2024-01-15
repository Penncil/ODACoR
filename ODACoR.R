

############################
#####Helper functions#######
############################
#the weibull distribution
local_log_likelihood<-function(beta,data_site_i,weights_i){
  
  px<-ncol(data_site_i[,-c(1:4)])
  partial_local_likelihood<-0
  
  
  T_local<-data_site_i$t_surv[data_site_i$type==1]
  T_local <- sort(unique(T_local))
  nt_local<- length(T_local)
  
  
  
  
  for(d in 1:nt_local){
    index_risk_local<-which(((data_site_i$type!=1)&(data_site_i$censor==1)&(data_site_i$t_surv<=T_local[d]))|(data_site_i$t_surv>=T_local[d]))
    weights_d<-rep(0,length(index_risk_local))
    weights_d[data_site_i$t_surv[index_risk_local]<T_local[d]]<-weights_i[data_site_i$t_surv==T_local[d]]/weights_i[index_risk_local][data_site_i$t_surv[index_risk_local]<T_local[d]]
    weights_d[data_site_i$t_surv[index_risk_local]>=T_local[d]]=1
    
    ##the covariates times beta in the risk set at time d 
    cov_risk<-as.matrix(data_site_i[index_risk_local,-c(1:4)])
    #  # take the cov_risk matrix times intial value as a vector
    A<-as.vector(exp(cov_risk%*%beta))
    B<-weights_d*A
    partial_local_likelihood<-partial_local_likelihood+(t(unlist(data_site_i[data_site_i$t_surv==T_local[d],-c(1:4)]))%*%beta)-log(sum(B))
    
  }
  
  
  
  partial_local_likelihood=partial_local_likelihood/length(data_site_i[,1])
  
  -partial_local_likelihood
}

#local_log_likelihood(beta,data_site_i,weights_i)
#get the global weights for each site
get_global_weights<-function(data_site_i,fit){
  n_local=length(data_site_i$id.site)
  weights_i=rep(0,length(data_site_i$id.site))
  for(i in 1:n_local){
    weights_i[i]=summary(fit,times=data_site_i$t_surv[i])$surv
  }
  return(weights_i)
  
}


get_local_weights<-function(data_site_i){
  Censor_local<-data_site_i$censor
  Censor_local<-1-Censor_local
  Censor_time_local<-data_site_i$t_surv
  fit = do.call('survfit',list(formula= Surv(Censor_time_local,Censor_local) ~ 1))
  n_local=length(data_site_i$id.site)
  weights_i=rep(0,length(data_site_i$id.site))
  for(i in 1:n_local){
    weights_i[i]=summary(fit,times=data_site_i$t_surv[i])$surv
  }
  return(weights_i)
}
get_local_est<-function(data_site_i){
  tryCatch({
    s1<-summary(crr(ftime =data_site_i$t_surv,fstatus=data_site_i$type,cov1=data_site_i[,-c(1:4)],failcode = 1))
    beta_h=s1$coef[,1]
    var_h=(s1$coef[,3])^2
    return(list(
      beta_h=beta_h,
      var_h=var_h
    ))
    
    
  }, error = function(ex) {
    px <- ncol(data_site_i[,-c(1:4)])
    beta_h=rep(NA,px)
    var_h=rep(NA,px)
    
    return(list(
      beta_h=beta_h,
      var_h=var_h
    ))
    #print(ex$message)
    
  }
  )
  
  
}


#this function is used to calcualte the first and second gradients for local log-likelihood (log-likelihood for local-site data)
first_second_local<-function(data_site_i,beta_bar,weights_i){
  
  px <- ncol(data_site_i[,-c(1:4)])
  #time of type=1 event for local site 
  T_local <- c()
  #time of local site
  T_local<-data_site_i$t_surv[data_site_i$type==1]
  T_local <- sort(unique(T_local))
  nt_local<- length(T_local)
  
  
  
  
  # local likelihood gradient (first derivative) at beta bar
  #covariate_sum for local site i
  cov_sum_i<-rep(0,px)
  
  #component of first at site i
  U_i<-rep(0,nt_local)
  #component of second at site i
  W_i<-matrix(data = 0, nrow =px, ncol = nt_local, byrow = FALSE, dimnames = NULL)
  #component of second at site i
  Z_i<-array(data=rep(0,px),dim=c(px,px,nt_local))
  
  
  cov_sum_i<-colSums(data_site_i[data_site_i$type==1,-c(1:4)])
  for(d in 1:nt_local){
    
    #get the index of the risk set for site i at time d
    index_risk<-which(((!data_site_i$type==1)&(data_site_i$censor==1)&(data_site_i$t_surv<=T_local[d]))|(data_site_i$t_surv>=T_local[d]))
    if(length(index_risk)==0){
      U_i[d]<-0
      W_i[,d]<-rep(0,px)
      Z_i[,,d]<-matrix(data = 0, nrow =px, ncol = px, byrow = FALSE, dimnames = NULL)
      
    }else{
      weights_d<-rep(0,length(index_risk))
      weights_d[data_site_i$t_surv[index_risk]<T_local[d]]<-weights_i[data_site_i$t_surv==T_local[d]]/weights_i[index_risk][data_site_i$t_surv[index_risk]<T_local[d]]
      weights_d[data_site_i$t_surv[index_risk]>=T_local[d]]=1
      
      #the covariates times beta in the risk set at time d 
      cov_risk<-as.matrix(data_site_i[index_risk,-c(1:4)])
      
      
      # take the cov_risk matrix times intial value as a vector
      A<-as.vector(exp(cov_risk%*%beta_bar))
      # take the value times the weights and sum up which is U_i
      B<-A*weights_d
      U_i[d]<-sum(B)
      
      #find the value of W for each d
      W_i[,d]<-W_i[,d]+colSums(weights_d*A*cov_risk)
      #find the value of Z for each d
      sum_1<-0
      for(j in 1:length(cov_risk[,1])){
        sum_1<-sum_1+weights_d[j]*A[j]*(outer(cov_risk[j,],cov_risk[j,]))
      }
      
      Z_i[,,d]<-sum_1
      
      
      
      
      
    }
    
  }
  
  
  
  
  #local likelihood first derivative
  expected_cov_local<-rep(0,px)
  for(i in 1:nt_local){
    expected_cov_local<-expected_cov_local+W_i[,i]/U_i[i]
  }
  ##local likelihood 's first derivative at the beta_bar##
  local_first_bbar<-(cov_sum_i-expected_cov_local)/length(data_site_i[,1])
  
  
  
  
  ##local likelihood 's first derivative at the beta_bar##
  
  sum_1<-0
  for(j in 1:length(W_i[1,])){
    first_component<-(W_i[,j]%*%t(W_i[,j]))/(U_i[j])^2
    second_component<-Z_i[,,j]/U_i[j]
    sum_1<-sum_1+first_component-second_component
  }
  
  local_second_bbar<-sum_1/length(data_site_i[,1])
  
  
  
  
  
  
  return(list(local_first_bbar=local_first_bbar,local_second_bbar=local_second_bbar))
}



get_meta_est<-function(all_site,n_sites){
  tryCatch({
    px<-ncol(all_site[,-c(1:4)])
    bhat_mat=matrix(data=NA,nrow=n_sites,ncol=px)
    vhat_mat=matrix(data=NA,nrow=n_sites,ncol=px)
    for(i in 1:n_sites){
      data_site_i=all_site[all_site$id.site==i,]
      a= get_local_est(data_site_i)
      bhat_mat[i,]=a$beta_h
      vhat_mat[i,]=a$var_h
      
      
    }
    
    betameta = apply(bhat_mat/vhat_mat,2,function(x){sum(x, na.rm = T)})/apply(1/vhat_mat,2,function(x){sum(x, na.rm = T)})
    
    
    return(betameta)
  },error=function(ex){
    #print(ex$message)
    betameta=rep(NA,px)
    return(betameta)
  })
  
  
}


competing_risk_sim<-function(site,n,p,beta,rho,gamma,min_a,max_b){
  #details see the fine and gray's paper
  #1-p is the mass at the infinity
  #also note that type 1 failure has the weibull
  #type 2 weibull distribution 
  rho_1<-rho[1]
  rho_2<-rho[2]
  gamma_1<-gamma[1]
  gamma_2<-gamma[2]
  
  id.site<-rep(site,n)
  t_surv<-rep(0,n)
  type<-rep(0,n)
  censor<-rep(0,n)
  z1<-rbinom(n,size=1,prob=0.5)
  z2<-runif(n,min=-2,max=2)
  beta_11<-beta[1]
  beta_12<-beta[2]
  beta_21<-beta[3]
  beta_22<-beta[4]
  
  for(i in 1:n){
    a<-exp(z1[i]*beta_11+z2[i]*beta_12)
    b<-exp(z1[i]*beta_21+z2[i]*beta_22)
    type[i]<-rbinom(1,size=1,prob=(1-p)^a)+1
    if(type[i]==2){
      u<-runif(1,min=0,max=1)
      c<--log((1-u)^(1/b))
      t_surv[i]<-rho_2*c^(1/gamma_2)
    }else{
      u<-runif(1,min=0,max=1)
      c<--log(1-(1-(1-u*(1-(1-p)^a))^(1/a))/p)
      t_surv[i]<-rho_1*c^(1/gamma_1)
    }
    
    
    t_censor<-runif(1,min=min_a,max=max_b)
    
    if(t_surv[i]<t_censor){
      censor[i]<-1
    }else{
      t_surv[i]<-t_censor
      type[i]<-0
      censor[i]<-0
    }
    
  }
  data_site<-data.frame(
    id.site=id.site,
    t_surv=t_surv,
    type=type,
    censor=censor,
    z1=z1,
    z2=z2
  )
  return(data_site)
}


get_local_UWZ<-function(data_site_i,beta_bar,T_all,fit){
  
  #dimension of covariates
  px <- ncol(data_site_i[,-c(1:4)])
  #weights for each site from the fit, fit is got distributively
  weights_i<-get_global_weights(data_site_i,fit)
  
  nt <- length(T_all)
  
  #weights_all
  weights_all<-rep(0,nt)
  for(i in 1:nt){
    weights_all[i]=summary(fit,times=T_all[i])$surv
  }
  
  
  
  if(sum(data_site_i$type==1)!=0){
    cov_sum_i=colSums(data_site_i[data_site_i$type==1,-c(1:4)])
  }else{
    cov_sum_i=rep(0,px)
  }
  
  #component of first at site i
  U_i<-rep(0,nt)
  #component of second at site i
  W_i<-matrix(data = 0, nrow =px, ncol = nt, byrow = FALSE, dimnames = NULL)
  #component of second at site i
  Z_i<-array(data=rep(0,px),dim=c(px,px,nt))
  
  for(d in 1:nt){
    #get the index of the risk set for site i at time d
    index_risk<-which(((!data_site_i$type==1)&(data_site_i$censor==1)&(data_site_i$t_surv<=T_all[d]))|(data_site_i$t_surv>=T_all[d]))
    if(length(index_risk)==0){
      U_i[d]<-0
      W_i[,d]<-rep(0,px)
      Z_i[,,d]<-matrix(data = 0, nrow =px, ncol = px, byrow = FALSE, dimnames = NULL)
      
    }else{
      
      #the covariates times beta in the risk set at time d 
      cov_risk<-as.matrix(data_site_i[index_risk,-c(1:4)])
      
      weights_d<-rep(0,length(index_risk))
      weights_d[data_site_i$t_surv[index_risk]<T_all[d]]=weights_all[d]/weights_i[index_risk][data_site_i$t_surv[index_risk]<T_all[d]]
      weights_d[data_site_i$t_surv[index_risk]>=T_all[d]]=1
      
      # take the cov_risk matrix times intial value as a vector
      A<-as.vector(exp(cov_risk%*%beta_bar))
      # take the value times the weights and sum up which is U_i
      B<-A*weights_d
      U_i[d]<-sum(B)
      
      #find the value of W for each d
      W_i[,d]<-W_i[,d]+colSums(weights_d*A*cov_risk)
      #find the value of Z for each d
      sum<-0
      for(j in 1:length(cov_risk[,1])){
        sum<-sum+weights_d[j]*A[j]*(outer(cov_risk[j,],cov_risk[j,]))
      }
      
      Z_i[,,d]<-sum
      
      
      
    }
    
  }
  
  
  
  
  
  
  
  
  return(list( 
    cov_sum_i=cov_sum_i,
    U_i=U_i,
    W_i=W_i,
    Z_i=Z_i))
}


d_distribute<-function(all_site,beta_bar,T_all,fit){
  # number of sites
  n_sites<-length(unique(all_site$id.site))
  
  px<-ncol(all_site[,-c(1:4)])
  n_list<-rep(0,n_sites)
  U_list<-list()
  W_list<-list()
  Z_list<-list()
  cov_sum_list<-list()
  for(i in 1:n_sites){
    data_site_i=all_site[all_site$id.site==i,]
    n_list[i]=length(data_site_i$id.site)
    a=get_local_UWZ(data_site_i,beta_bar,T_all,fit)
    cov_sum_list[[i]]=a$cov_sum_i
    U_list[[i]]=a$U_i
    W_list[[i]]=a$W_i
    Z_list[[i]]=a$Z_i
  }
  
  return(
    list( n_list=n_list,
          cov_sum_list=cov_sum_list,
          U_list=U_list,
          W_list=W_list,
          Z_list=Z_list)
  )
  
  
}



d_assemble<-function(n_list,cov_sum_list,U_list,W_list,Z_list){
  #dimension of covariates
  px <- length(cov_sum_list[[1]])
  
  nt <- length(U_list[[1]])
  
  n_sites=length(n_list)
  #component of first at site i
  U<-rep(0,nt)
  cov_sum=rep(0,px)
  #component of second at site i
  W<-matrix(data = 0, nrow =px, ncol = nt, byrow = FALSE, dimnames = NULL)
  #component of second at site i
  Z<-array(data=rep(0,px),dim=c(px,px,nt))
  for(i in 1:n_sites){
    U=U+U_list[[i]]
    W=W+W_list[[i]]
    Z=Z+Z_list[[i]]
    cov_sum=cov_sum+cov_sum_list[[i]]
  }
  
  
  expected_cov<-rep(0,px)
  for(i in 1:nt){
    expected_cov<-expected_cov+W[,i]/U[i]
  }
  ##global likelihood 's first derivative at the beta_bar##
  global_first_bbar<-(cov_sum-expected_cov)/sum(n_list)
  
  
  
  ##global likelihood 's second derivative at the beta_bar##
  
  sum_1<-0
  for(i in 1:length(W[1,])){
    first_component<-(W[,i]%*%t(W[,i]))/(U[i])^2
    second_component<-Z[,,i]/U[i]
    sum_1<-sum_1+first_component-second_component
  }
  global_second_bbar<-sum_1/sum(n_list)
  
  
  return(
    list(
      global_first_bbar=global_first_bbar,
      global_second_bbar=global_second_bbar
    )
  )
  
}



#this function is used to construct the surrogate likelihood 
surrogate_log_likelihood<-function(beta,beta_bar,global_first_bbar,global_second_bbar,local_first_bbar,local_second_bbar,data_site_i,weights_i){
  
  
  
  local_log_likelihood(beta=beta,data_site_i=data_site_i,weights_i=weights_i)-(global_first_bbar-local_first_bbar)%*%beta-(0.5)*(t(beta-beta_bar))%*%(global_second_bbar-local_second_bbar)%*%(beta-beta_bar)
  
  
  
  
}

#this function is used to get the surrogate estimator
get_surrogate<-function(beta_bar,global_first_bbar,global_second_bbar,local_first_bbar,local_second_bbar,data_site_i,weights_i){
  
  tryCatch({ 
    
    sol_sur=optim(par =beta_bar, 
                  fn = surrogate_log_likelihood,
                  beta_bar=beta_bar,
                  global_first_bbar=global_first_bbar,
                  global_second_bbar=global_second_bbar,
                  local_first_bbar=local_first_bbar,
                  local_second_bbar=local_second_bbar,
                  data_site_i=data_site_i,
                  weights_i=weights_i,
                  method ="BFGS",
                  control = list(maxit = 20))
    est_sur=sol_sur$par
    
    return(est_sur)
  },error=function(ex){
    print(ex$message)
    est_sur=rep(NA,length(beta_bar))
    return(est_sur)
  })
  
  
  
}












d_initialize<-function(all_site){
  
  # time for event type of 1
  T_all<-c()
  # time 
  Censor_time <- c()
  #indicator if censor 0, not censor 1
  Censor<-c()
  
  
  # number of sites
  control_sites<-1:n_sites
  
  for(site_i in control_sites){
    
    T_all <- c(T_all, all_site$t_surv[all_site$id.site==site_i&all_site$type==1])
    Censor<-c(Censor,all_site$censor[all_site$id.site==site_i])
    Censor_time<-c(Censor_time,all_site$t_surv[all_site$id.site==site_i])
    
    
  }
  
  
  #get the indicator 1 if censoring 0 if not
  Censor<-1-Censor
  T_all <- sort(unique(T_all))
  #get the estimate of the G/ the survival function
  fit<-survfit(formula = Surv(Censor_time,Censor) ~ 1)
  
  return(
    list(
      T_all=T_all,
      fit=fit
    )
  )
  
  
}


##############################
###### Main function ########
############################
#start from here!
# first install the packages ###
#this one is for the competing risk model
library(cmprsk)

#this one is to handle survival dataset
library(survival)




###########################################
### Part 1: Generate simulated dataset ####
############################################
#Each site corresponding to one specific sub-hazard distribution
#n n_sites beta gamma rho p min_a max_b are used to simulate the competing risk dataset
#does not use these parameters when analysis true dataset, but can be used to show the example
#beta is the true coefficient
#how many sites
n_sites=3
#how many participantes for each site assumed that all sites have the same participants
n=100
# true coefficients for main event (event 1) and competing risk event (event 2)
#notice that always set the event of interest as the event 1.
beta<-c(0.5,-0.5,0.5,-0.5)

#p rho gamma are about to generate the underlying sub-distribution hazard
p<-0.6
rho<-c(100,90)
gamma<-c(10,12)

#min_a max_b are about the parameters of censoring
min_a<-0
max_b<-110



set.seed(1)

all_site<-NULL

for(i in 1:n_sites){
  site_i<-competing_risk_sim(site=i,n=n,p=p,beta=beta,rho=rho,gamma=gamma,min_a=min_a,max_b=max_b)
  all_site<-rbind(all_site,site_i)
}
## try to make the real_data looks like the simulated data##
## the first column should be site id, needs to be number starting from 1,2...
## the second column should be the time of events
## the third column should be event type, put the event of interest as 1 and the competing event as 2
## the fourth column is the censoring. If uncensored, mark it as 1. If censored, mark it as 0. Same as before
## the fifth column to the rest of the columns are about the covariates, the dimension can be more than 2.d=2 is the example

## first try the binary (0,1) and the continuous covariates,
## second if the covariates are categorical for example k levels, you should first make them to k dummy variables and put
## k-1 dummy variable into the model for part 2-4. The part 5 can automatically deal with the categorical covariates. 
## You can first start there.

#####################################################
######################Code to break ties################
###########################################
##codes to break ties... 
##add this chunk of codes before you run the gold standard/surrogate/meta/one-shot estimator
unique_times=unique(all_site$t_surv)
if(length(unique_times)<length(all_site$t_surv)){
  for(i in 1:length(unique_times)){
    index_ties_i=which(all_site$t_surv==unique_times[i])
    small_noise=rnorm(length(index_ties_i),sd=0.01)
    all_site$t_surv[index_ties_i]=all_site$t_surv[index_ties_i]+small_noise
  }
}


###########################################
### Part 2: Calculate meta-analysis estimator ####
############################################
n_sites=length(unique(all_site$id.site))
#meta analysis estimator
beta_bar=get_meta_est(all_site,n_sites)




###########################################
### Part 3: Calculate surrogate estimator ####
############################################


### Initial Step ###
#communicate the surv fit/ time of event of interest/ time of event of censoring
initialize_results<-d_initialize(all_site)
T_all=initialize_results$T_all
fit=initialize_results$fit


##set the log-likelihood for the  local site to be the one with type 1 event#

for(i in 1:n_sites){
  #check if there is type 1 event in that set.
  if(sum(all_site[all_site$id.site==i,]$type==1)>=1){
    data_site_i=all_site[all_site$id.site==i,]
    break
  }
  
}


#then we calculate the summary-level statistics for each site 

b=d_distribute(all_site,beta_bar,T_all,fit)
n_list=b$n_list
cov_sum_list=b$cov_sum_list
U_list=b$U_list
W_list=b$W_list
Z_list=b$Z_list
#then we assemble to get the global first and second gradients
c=d_assemble(n_list,cov_sum_list,U_list,W_list,Z_list)
global_first_bbar=c$global_first_bbar
global_second_bbar=c$global_second_bbar

#also the gradients for the local likelihood
weights_i=get_local_weights(data_site_i)
a<-first_second_local(data_site_i, beta_bar,weights_i)


local_first_bbar<-a$local_first_bbar
local_second_bbar<-a$local_second_bbar
#surrogate
#default: 200 number of iterations

sol_ODACoR_S<- get_surrogate(beta_bar,global_first_bbar,global_second_bbar,local_first_bbar,local_second_bbar,data_site_i,weights_i)


###########################################
### Part 4: Calculate one-step estimator ####
############################################

#Can directly use the summary-level information proposed before##
sol_ODACoR_O<-beta_bar-solve(global_second_bbar)%*%global_first_bbar




