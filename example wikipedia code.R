# example code from wikipedia

#load library for multivariate normal

require(haven)
require(dplyr)
require(forcats)
require(tidyr)
require(ggplot2)

anes <- read_dta('anes_timeseries_2016.dta')

num_clust <- 2
threshold <- .0001

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
initial_sigma <- rep(list(cov(anes_matrix)),num_clust)
create_theta <- function(tau=NULL,mu=NULL,sigma=NULL) {
  out_data <- data_frame(tau=tau,mu1=mu[,1],mu2=mu[,2],sigma=sigma)
  return(out_data)
}

calc_dens <- function(params,run_data=NULL) {
  all_dens <- sapply(1:nrow(params),function(i) {
    params$tau[i] * mvtnorm::dmvnorm(run_data,mean=c(params$mu1[i],params$mu2[i]),
                            sigma=params$sigma[[i]])
  })
  return(all_dens)
} 

theta <- create_theta(
  tau=rep(1/num_clust,num_clust),
  mu=cbind(runif(n = num_clust,min = 10,max=90),runif(n = num_clust,min = 10,max=90)),
  sigma=initial_sigma)

#E step: calculates conditional probabilities for latent variables
e_step <- function(theta,run_data) {
  densities <- calc_dens(theta,run_data=run_data$run_data)
  stan_dens <- t(apply(densities,1,function(x) x/sum(x)))
  return(list(stan_dens=stan_dens,unstan_dens=densities))
}
  # t(apply(cbind(
  #   theta$tau[1] * dmvnorm(anes_matrix,mean=theta$mu1,sigma=theta$sigma1),
  #   theta$tau[2] * dmvnorm(anes_matrix,mean=theta$mu2,sigma=theta$sigma2),
  #   theta$tau[3] * dmvnorm(anes_matrix,mean=theta$mu3,sigma=theta$sigma3),
  #   theta$tau[4] * dmvnorm(anes_matrix,mean=theta$mu4,sigma=theta$sigma4)
  # ),1,function(x) x/sum(x)))
#M step: calculates the parameter estimates which maximise Q
m_step <- function(stan_dens,m_run_data=NULL) {

  all_mus <- apply(stan_dens,2,function(col) {

    this_mu <- apply(m_run_data$run_data,2,weighted.mean,col)
    return(this_mu)
  })
  
  all_sigmas <- lapply(1:ncol(stan_dens),function(col,this_data=NULL) {
    out_sigma <- cov.wt(this_data$run_data,stan_dens[,col])$cov
  },this_data=m_run_data)

  params <- create_theta(tau=apply(stan_dens,2,mean),
                         mu=t(all_mus),
                         sigma=all_sigmas)
}

run_em <- function(this_data=NULL,theta,threshold) {

  all_dens <- e_step(theta=theta,run_data=this_data)
  stan_dens <- all_dens$stan_dens
  params <- try(m_step(stan_dens=stan_dens,m_run_data=this_data))
  if(class(params)[1]=='try-error') {
    print('Cluster collapsed, trying to sumplement with random values.')
    i <- 1
    while(class(params)[1]=='try-error') {
      print(paste0('Resample ',i))
      
      stan_dens <- apply(stan_dens,2, function(x) {
        if(sum(x,na.rm=TRUE)==0) {
          x[is.na(x)] <- 0
          x <- x + abs( rnorm(n=length(x),mean=0,sd=0.01))
        }
        return(x)
      })
      params <- try(m_step(stan_dens,m_run_data=this_data))
    }
  }
  return(list(theta=params,stan_dens=stan_dens,unstan_dens=all_dens$unstan_dens))
}

generate_plot <- function(stan_dens,params,current_data=NULL,iter,new_dens=NULL) {
  params <- mutate(params,new_dens=sum(new_dens))
  all_data <- list(stan_dens=as_data_frame(stan_dens),params=as_data_frame(params))
  all_data <- lapply(all_data,mutate,iter=iter)
  if(!is.null(current_data)) {
    current_data$stan_dens <- bind_rows(current_data$stan_dens,all_data$stan_dens)
    current_data$params <- bind_rows(current_data$params,all_data$params)
    return(current_data)
  }
  return(all_data)
}
  
iterate_over <- function(reps,run_data=NULL,theta,threshold,current_data=NULL,i=1,old_dens=NULL) {
  this_data <- new.env(parent = emptyenv())
  this_data$run_data <- run_data  

  store_data <- new.env(parent=emptyenv())
  store_data$compute <- run_data
  print(paste0('Now on Iteration: ',i))
  if(is.null(old_dens)) old_dens <- 0
  
  out_em <- run_em(theta=theta,threshold=threshold,this_data=this_data)
  theta <- out_em$theta
  new_dens <- out_em$unstan_dens
  current_data <- generate_plot(stan_dens=out_em$stan_dens,params=theta,iter=i,current_data=current_data,new_dens=new_dens)
  big_clust <- which(theta$tau==max(theta$tau))

  i <- i + 1
  
  if(((sum(new_dens)-sum(old_dens))<threshold) && ((sum(new_dens)-sum(old_dens))>0)) {
    print('Converged!')
    return(list(current_data=current_data,new_dens=new_dens))
  } else if(i==reps) {
    print('Finished iterations')
    return(list(current_data=current_data,new_dens=new_dens)) 
    } else {
    iterate_over(reps=reps,run_data=this_data$run_data,theta=theta,threshold=threshold,current_data=current_data,
                 i=i,old_dens=new_dens)
  }
}  

gen_data <- iterate_over(reps = 500,run_data=anes_matrix,theta = theta,threshold=threshold)
plot_data <- mutate(gen_data$current_data$params,cluster=rep(paste0('Cluster_',1:num_clust),n()/num_clust)) %>% 
  gather(param,param_est,tau,mu1,mu2,-cluster,-iter)
plot_data %>% ggplot(aes(y=param_est,x=iter,colour=cluster)) +
  geom_line() + theme_minimal() + facet_wrap(~param,scales = 'free')
plot_data %>% distinct(new_dens,iter) %>% 
  mutate(dens_diff=new_dens - lag(new_dens)) %>% 
           gather(dens_type,dens_val,new_dens,dens_diff) %>% 
           ggplot(aes(y=dens_val,x=iter)) +
  geom_line() + theme_minimal() + facet_wrap(~dens_type,scales='free')

# filter(plot_data,iter>450) %>% distinct(new_dens,iter) %>% 
#   mutate(dens_diff=new_dens - lag(new_dens)) %>% 
#   gather(dens_type,dens_val,new_dens,dens_diff) %>% ggplot(aes(y=dens_val,x=iter)) +
#   geom_line() + theme_minimal() + facet_wrap(~dens_type,scales='free') + 
#   geom_hline(yintercept=threshold,linetype=2,colour='red')


contour_data <- spread(plot_data,key=param,value=param_est)

new_dens <- t(gen_data$new_dens %>% apply(1,function(x) x/sum(x)))
# names(new_dens) <- paste0('Clust_',1:num_clust)
# new_dens <- as_data_frame(new_dens)
# run_data_plot <- as_data_frame(anes_matrix) 
# names(run_data_plot) <- c('Var1','Var2')
# interp_var <- lapply(new_dens,function(y) {
#   model_data <- data_frame(y=y,Var1=anes_matrix[,1],Var2=anes_matrix[,2])
#   all_vals <- expand.grid(0:100,0:100)
#   loess_out <- splines::interpSpline(y~Var1 + Var2,model_data)
#   loess_out_d <- as_data_frame(loess_out) 
#   loess_out_d %>% mutate(Var1=row.names(loess_out)) %>% gather(Var2,density_val,-Var1)
# })
# interp_var <- lapply(interp_var,function(x) x$density_val) %>% as_data_frame
# names(interp_var) <- paste0('Cluster_',1:num_clust)
# out_data <- expand.grid(0:100,0:100)
# out_data <- bind_cols(out_data,interp_var)
# interp_data <- gather(out_data,clust_num,density_val,-Var1,-Var2)
# 
# filter(interp_data,clust_num=='Cluster_2') %>% ggplot(aes(y=Var1,x=Var2,z=density_val)) + geom_contour()
#   geom_point(aes(colour=density_val)) + theme_minimal()
# new_dens_matrix <- as.matrix(new_dens)
check_data <- data_frame(Var1=anes_matrix[,1],Var2=anes_matrix[,2],z=apply(new_dens,1,max),
                         clust_pred=apply(new_dens,1,function(x) which(x==max(x))),
                         trump_feel=as.numeric(anes_to_use$trump_feeling)) %>% 
                        mutate(trump_feel=(trump_feel - min(trump_feel))/(max(trump_feel)-min(trump_feel)))
# 
# ggplot(check_data,aes(y=Var1,x=Var2,z=mean_dens)) + geom_contour() + geom_point(aes(colour=mean_dens))


# Add annotations
  
clust_labels <- gen_data$current_data$params %>% filter(iter==max(iter)) %>% 
  mutate(label=paste0('Clust_',1:num_clust))

all_data <- expand.grid(0:100,0:100)

try_lin <- lm(z~Var1 + Var2 + Var1*Var2 + I(Var1^2) + I(Var2^2) + 
                I(Var1^2)*I(Var2^2) + I(Var1^3) + I(Var2^3) + 
                I(Var1^3)*I(Var2^3) + I(Var2^4) + I(Var1^4) + I(Var1^4)*I(Var2^4) +
                I(Var1^5) + I(Var2^5) + I(Var1^5)*I(Var2^5) +
                I(Var1^6) + I(Var2^6) + I(Var1^6)*I(Var2^6) ,data=check_data)

out_predict <- predict(try_lin,newdata=all_data)

out_res <- predict(try_lin)

#RMSE

mean(sqrt((out_res-check_data$z)^2))
sd(sqrt((out_res-check_data$z)^2))

z_vals <- matrix(out_predict,nrow=length(0:100),ncol=length(0:100))

plot_ly(check_data) %>% 
  add_contour(z=~z_vals,contours=list(coloring='lines')) %>% 
  add_markers(x=~Var1,y=~Var2,symbol=~clust_pred,color=~trump_feel,colors=colorRampPalette(brewer.pal(11,"Spectral"))(100)) %>% 
  add_text(data=clust_labels,y=~mu1,x=~mu2,text=~label,textfont=list(color='white',size=18)) 
