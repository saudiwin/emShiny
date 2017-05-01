# example code from wikipedia

#load library for multivariate normal

require(haven)
require(dplyr)
require(forcats)

anes <- read_dta('anes_timeseries_2016.dta')

# We will use two variables: feeling thermometer for 'Christian Fundamentalists' V162095 and 
# feeling thermometer for 'Muslims' V162106
# Points will be coded by V162171 'Ideology scale' (1-7), collapsed to 1-5
# Points colored by V162079 'Feeling thermometer towards Trump'
# Survey weights V160102

anes_to_use <- select(anes,V162095,V162106,V162171,V162079,V160102) %>% haven::as_factor() %>% 
  mutate_if(is.factor,as.character) %>% mutate_if(is.character,function(x) {
    x[grepl(x = x,pattern="-[0-9]|thought much about this|[Dd]on't know")] <- NA
    return(x)
  }) %>% rename(trump_feeling=V162079,
                chris_feeling=V162095,
                muslim_feeling=V162106,
                ideology=V162171,
                weights=V160102) %>% 
  mutate(ideology=factor(ideology,levels=c('7. Extremely conservative',
                                           '6. Conservative',
                                           '5. Slightly conservative',
                                           '4',
                                           '3. Slightly liberal',
                                           '2. Liberal',
                                           '1. Extremely liberal')),
         ideology=fct_collapse(ideology,C=c('5. Slightly conservative','6. Conservative'),
                               `VC`=c('7. Extremely conservative'),
                               L=c('3. Slightly liberal','2. Liberal'),
                               `VL`='1. Extremely liberal',
                               M='4'))
anes_to_use <- filter(anes_to_use,complete.cases(anes_to_use))

# We are left with 2,963 obs from the original 4,271

anes_matrix <- cbind(as.numeric(anes_to_use$chris_feeling),as.numeric(anes_to_use$muslim_feeling))

#setup grid for plotting
xpts <- seq(from=min(anes_matrix[,1]),to=max(anes_matrix[,1]),length.out=100)
ypts <- seq(from=min(anes_matrix[,2]),to=max(anes_matrix[,2]),length.out=100)

#initial parameter estimates
# correlation matrix puts negative correlations between the two clusters-- e.g., push each other apart
initial_sigma <- rep(list(cov(anes_matrix)),4)
create_theta <- function(tau,mu,sigma,run_data) {
  out_data <- data_frame(tau=tau,mu1=mu[,1],mu2=mu[,2],sigma=sigma)
  return(out_data)
}

calc_dens <- function(params,run_data) {

  all_dens <- sapply(1:nrow(params),function(tau) {
    params$tau[i] * dmvnorm(run_data,mean=c(params$mu1[i],params$mu2[i]),
                            sigma=params$sigma[[i]])
  })
  return(all_dens)
} 

theta <- create_theta(
  tau=c(0.25,0.25,0.25,0.25),
  mu=cbind(c(80,80,20,20),c(80,20,80,20)),
  sigma=initial_sigma,
  run_data=anes_matrix)

#E step: calculates conditional probabilities for latent variables
e_step <- function(theta,run_data) {
  densities <- calc_dens(theta,run_data=run_data)
  stand_dens <- t(apply(densities,1,function(x) x/sum(x)))
  return(stand_dens)
}
  # t(apply(cbind(
  #   theta$tau[1] * dmvnorm(anes_matrix,mean=theta$mu1,sigma=theta$sigma1),
  #   theta$tau[2] * dmvnorm(anes_matrix,mean=theta$mu2,sigma=theta$sigma2),
  #   theta$tau[3] * dmvnorm(anes_matrix,mean=theta$mu3,sigma=theta$sigma3),
  #   theta$tau[4] * dmvnorm(anes_matrix,mean=theta$mu4,sigma=theta$sigma4)
  # ),1,function(x) x/sum(x)))
#M step: calculates the parameter estimates which maximise Q
m_step <- function(stan_dens,run_data=NULL) {
  run_data <- force(run_data)
  all_mus <- apply(stan_dens,2,function(col) {

    this_mu <- apply(run_data,2,weighted.mean,col)
    return(this_mu)
  })
  
  all_sigmas <- lapply(1:ncol(stan_dens),function(col) {
    out_sigma <- cov.wt(run_data,col)$cov
  })
  
  params <- create_theta(tau=apply(stan_dens,2,mean),
                         mu1=all_mus[,1],
                         mu2=all_mus[,2],
                         sigma=all_sigmas)
}

run_em <- function(run_data=NULL,theta,threshold) {
  stan_dens <- e_step(theta=theta,run_data=run_data)
  params <- try(m_step(stan_dens=stan_dens,run_data=run_data))
  if(class(params)=='try-error') {
    print('Cluster collapsed, trying to sumplement with random values.')
    i <- 1
    while(class(params)=='try-error') {
      print(paste0('Resample ',i))
      
      stan_dens <- apply(stan_dens,2, function(x) {
        if(sum(x,na.rm=TRUE)==0) {
          x[is.na(x)] <- 0
          x <- x + abs( rnorm(n=length(x),mean=0,sd=0.01))
        }
        return(x)
      })
      params <- try(m_step(stan_dens))
    }
  return(list(theta=params,stan_dens=stan_dens))
  }
}

generate_plot <- function(stan_dens,params,current_data=NULL,iter) {
  all_data <- list(stan_dens=as_data_frame(stan_dens),params=as_data_frame(params))
  all_data <- lapply(all_data,mutate,iter=paste0('Iteration_',iter))
  if(!is.null(current_data)) {
    current_data$stan_dens <- bind_rows(current_data$stan_dens,all_data$stan_dens)
    current_data$params <- bind_rows(current_data$params,all_data$params)
    return(current_data)
  }
  return(all_data)
}
  
iterate_over <- function(reps,run_data=NULL,theta,threshold,current_data=NULL) {

  if(is.null(current_data)) {
    current_data <-  data_frame()
  }
  for(i in 1:reps) {
    old_tau <- theta$tau
    print(paste0('Now on Iteration: ',i))
    out_em <- run_em(run_data=run_data,theta=theta,threshold=threshold)
    theta <- out_em$theta
    
    if(is.null(current_data)) {
      current_data <- generate_plot(stan_dens=out_em$stan_dens,params=theta,iter=i)
    } else {
      current_data <- generate_plot(stan_dens=out_em$stan_dens,params=theta,current_data=current_data,iter=i)
    }
    
    big_clust <- which(theta$tau==max(theta$tau))
    if((abs(theta$tau[big_clust]-old_tau[big_clust])<threshold)) {
      print('Converged!')
      break
    }
  }
  return(current_data)
}  

gen_data <- iterate_over(reps = 10,run_data=anes_matrix,theta = theta,threshold=.00001)
  # }list(
  # tau= apply(T,2,mean),
  # mu1= apply(anes_matrix,2,weighted.mean,T[,1]),
  # mu2= apply(anes_matrix,2,weighted.mean,T[,2]),
  # mu3= apply(anes_matrix,2,weighted.mean,T[,3]),
  # mu4= apply(anes_matrix,2,weighted.mean,T[,4]),
  # sigma1= cov.wt(anes_matrix,T[,1])$cov,
  # sigma2= cov.wt(anes_matrix,T[,2])$cov,
  # sigma3= cov.wt(anes_matrix,T[,3])$cov,
  # sigma4= cov.wt(anes_matrix,T[,4])$cov)

# #function to plot current data
# plot.em <- function(theta){
#   mixture.contour <- outer(xpts,ypts,function(x,y) {
#     theta$tau[1]*dmvnorm(cbind(x,y),mean=theta$mu1,sigma=theta$sigma1) + theta$tau[2]*dmvnorm(cbind(x,y),mean=theta$mu2,sigma=theta$sigma2)
#   })
#   contour(xpts,ypts,mixture.contour,nlevels=5,drawlabel=FALSE,col="red",xlab="Eruption time (mins)",ylab="Waiting time (mins)",main="Waiting time vs Eruption time of the Old anes_matrix geyser")
#   points(anes_matrix)
# }



# file.remove(list.files(pattern='.png'))
# #plot initial contours
# iter <- 1
# png(filename=paste("em",formatC(iter,width=4,flag="0"),".png",sep=""))
# plot.em(theta)
# dev.off()
# 
# #run EM and plot until stopping
# threshold <- .00001
# for (iter in 2:1000){
#   print(paste0('On iteration ',iter,' with tau1=',theta$tau[1],' and tau2=',theta$tau[2],'and tau3=',theta$tau[3],'and tau4=',theta$tau[4]))
#   old_theta <- theta
#   old_tau <- theta$tau
#   T <- E.step(theta)
#   theta <- try(M.step(T))
#   
#   if(class(theta)=='try-error') {
#     print('Cluster collapsed, trying to sumplement with random values.')
#     i <- 1
#     while(class(theta)=='try-error') {
#       print(paste0('Resample ',i))
#       
#       T <- apply(T,2, function(x) {
#           if(sum(x,na.rm=TRUE)==0) {
#             x[is.na(x)] <- 0
#             x <- x + abs( rnorm(n=length(x),mean=0,sd=0.01))
#           }
#         return(x)
#         })
#       theta <- try(M.step(T))
#     }
#   }
#   #check for convergence in biggest cluster
# 
#   if((iter %% 10)==0) {
#   png(filename=paste("em",formatC(iter,width=4,flag="0"),".png",sep=""))
#   plot.em(theta)
#   dev.off()
#   }
# 
# }