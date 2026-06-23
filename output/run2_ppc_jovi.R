#!/usr/bin/env Rscript

rm(list=ls()) 
cat("\f")

library(ggplot2)
library(cmdstanr)
library(bayesplot)
library(HDInterval)
library(RWiener)
library(posterior)
library(dplyr)
library(tidyr)
library(data.table)
library(loo)
library(parallel)
library(doParallel)
library(foreach)

args <- commandArgs(trailingOnly = TRUE) # get arguments from batch script
run<-as.numeric(args[1]) #break up into different concurrent jobs on hpc
group<-group_id<-as.numeric(args[2]) #1=hc, 2=sz, all=both
modelname <-args[3] # "hddm_m[#]_psychophys"
task<-args[4] # task (gaze or gender or jovi)
seed<-42+run # different runs get different rand seeds
set.seed(seed)

cores <- 10       
options(mc.cores = cores)
setDTthreads(threads=cores)

#############################################################################
# SET DIRECTORIES
#############################################################################

dirname<-'[PATH TO FIT DIRECTORY]' 
datadir<-paste0(dirname,'/data')
scriptsdir<-paste0(dirname,'/scripts')
behfile<-paste0(datadir,"/schizgaze2_jovi_beh.csv")
modeldir<-paste0(dirname,'/output/',modelname,'/',task)
initdir<-paste0(modeldir,'/init_fit')
loglikdir<-paste0(modeldir,'/log_lik')
ppcdir<-paste0(modeldir,'/ppc')
parrecoverdir<-paste0(modeldir,'/par_recover')
finaldir<-paste0(modeldir,'/final_fit')

#############################################################################
# LOAD AND PREPROCESS BEHAVIORAL DATA FOR STAN
#############################################################################

# load preprocessed behavioral data and subset to schizgaze2 only
alldata<-read.csv(behfile)
alldata<-subset(alldata,alldata$Study=="schizgaze2")

# sort subj id's in order (hc's [1000s] will come first, then sz [2000s]) so  when we assign sequential index id's, they will be the same for gender and gaze tasks)
alldata <- alldata[order(alldata$Subj), ]

# remove 2 subjects that didn't have valid data from gaze task
alldata<-subset(alldata,alldata$Subj!="2032")
alldata<-subset(alldata,alldata$Subj!="2042")

#subset to get the appropriate task data
alldata$SubjID<-match(alldata$Subj, unique(alldata$Subj)) #subj IDs to sequential indexes 
alldata$Resp<-ifelse(alldata$Resp=="2",1,2) #recode L (2) as 1 and R (3) as 2
alldata$Group<-match(alldata$Group, unique(alldata$Group))
alldata <- alldata[order(alldata$SubjID, alldata$Resp), ]

minRT<-alldata %>%
  dplyr::group_by(SubjID) %>%
  dplyr::summarise(minRT = min(RT)) %>%
  dplyr::ungroup()

if (max(alldata$RT)>100){alldata$RT<-as.numeric(alldata$RT)/1000} #if rt col is in ms, convert to seconds

#############################################################################
# SET DDM SIMULATION SPECS
#############################################################################

N_choice<-length(unique(alldata$RespBound))
groups<-unique(alldata$Group)
groups_names<-c("HC","SZ")
subjs<-unique(alldata$Subj)
n_obs<-nrow(alldata)
N_jitters<-length(unique(alldata$Jitter))
jitters<-unique(alldata$Jitter)

estimated_quantile_file<-paste(ppcdir,'/quantile_choice_estimates.RData',sep="")
quantiles<-c(.1,.3,.5,.7,.9)

if(file.exists(estimated_quantile_file)){ #if file exists,load it
  load(estimated_quantile_file)
}else{ #if file doesn't exist, run code to estimate quantiles
  
  posterior<-as.data.frame(fread(paste(initdir,'/allgroups_init_fit_samples.csv',sep=""),nThread=cores))
  iterations=nrow(posterior) #number of posterwarmup draws
  
  #setup cluster on computer for parallel processing
  cluster <- makeCluster(cores) 
  registerDoParallel(cluster)
  
  # get subj level variables (group and minRT)
  subj_vars<-alldata %>%
    dplyr::group_by(Group, SubjID) %>%
    dplyr::summarize(minRT = min(RT))
  
  subsample<-1000 #how many draws should we randomly subsample from the full posterior
  
  #create empty array for rt quantiles: #obs x #subsampled iterations x #quantiles
  estim_quantiles_prop <- array(numeric(),c(n_obs,subsample,(length(quantiles)+1)*2))
  
  for (k in 1:n_obs){
    
    #rand subsample
    rand_draws<-sample(seq(from=1,to=iterations,by=1),subsample)
    tmp_posterior<-posterior[rand_draws,]
    
    tmp_subj<-alldata$SubjID[k]
    tmp_group<-alldata$Group[k]
    tmp_jitter<-alldata$Jitter[k]
    tmp_choice<-alldata$Resp[k]
    tmp_rt<-alldata$RT[k]
    tmp_minRT<-subj_vars$minRT[tmp_subj] 
    
    if(task=='gaze'|task=='jovi'){
      tmp_gendername<-''
    }else{
      if(tmp_gender==1){
        tmp_gendername<-'female_'
      }else{
        tmp_gendername<-'male_'
      }
    }
    
    # some versions of stan output use different delimiters. check and update accordingly.
    num_brackets <- sum(grepl("\\[", colnames(tmp_posterior))) + sum(grepl("\\]", colnames(tmp_posterior)))
    if(num_brackets==0){
      start_del<-'.'
      end_del<-'.'
    }else{
      start_del<-'['
      end_del<-']'
    }
    
    #randomly draw group mean, variance values from appropriate posterior
    tmp_mu_grp_alpha_pr<-tmp_posterior[[paste("mu_grp_alpha_pr",start_del,tmp_group,end_del,sep="")]]
    tmp_mu_grp_beta_pr<-tmp_posterior[[paste("mu_grp_beta_pr",start_del,tmp_group,end_del,sep="")]]
    tmp_mu_grp_ndt_pr<-tmp_posterior[[paste("mu_grp_ndt_pr",start_del,tmp_group,end_del,sep="")]]
    tmp_mu_grp_delta_pr<-tmp_posterior[[paste("mu_grp_delta_",tmp_gendername,"pr",start_del,tmp_group,end_del,sep="")]]
    
    tmp_sig_grp_alpha_pr<-tmp_posterior[[paste("sig_grp_alpha_pr",start_del,tmp_group,end_del,sep="")]]
    tmp_sig_grp_beta_pr<-tmp_posterior[[paste("sig_grp_beta_pr",start_del,tmp_group,end_del,sep="")]]
    tmp_sig_grp_delta_pr<-tmp_posterior[[paste("sig_grp_delta_pr",start_del,tmp_group,end_del,sep="")]]
    tmp_sig_grp_ndt_pr<-tmp_posterior[[paste("sig_grp_ndt_pr",start_del,tmp_group,end_del,sep="")]]
    
    # get subj level values
    tmp_sub_alpha_pr<-tmp_posterior[[paste("sub_alpha_pr",start_del,tmp_subj,end_del,sep="")]]
    tmp_sub_beta_pr<-tmp_posterior[[paste("sub_beta_pr",start_del,tmp_subj,end_del,sep="")]]
    tmp_sub_ndt_pr<-tmp_posterior[[paste("sub_ndt_pr",start_del,tmp_subj,end_del,sep="")]]
    tmp_sub_delta_pr<-tmp_posterior[[paste("sub_delta_",tmp_gendername,"pr",start_del,tmp_subj,end_del,sep="")]]
    
    # if b1 parameters exist, get them. otherwise set to 0
    if(length(grep("b1",colnames(tmp_posterior)))>0){
      tmp_mu_grp_b1_pr<-tmp_posterior[[paste("mu_grp_b1_pr",start_del,tmp_group,end_del,sep="")]]
      tmp_sig_grp_b1_pr<-tmp_posterior[[paste("sig_grp_b1_pr",start_del,tmp_group,end_del,sep="")]]
      tmp_sub_b1_pr<-tmp_posterior[[paste("sub_b1_pr",start_del,tmp_subj,end_del,sep="")]]
      #do transform for non-centered parameterization
      tmp_sub_b1 <- tmp_mu_grp_b1_pr + tmp_sig_grp_b1_pr*tmp_sub_b1_pr
    }else{
      tmp_sub_b1 <- rep(0,subsample)
    }
    
    # if sub delta variance parameters exist, get them. otherwise set to 0
    if(length(grep("sig_sub_delta",colnames(posterior)))>0){
      tmp_sig_sub_delta_pr<-tmp_posterior[[paste("sig_sub_delta_pr",start_del,tmp_subj,end_del,sep="")]]
    }else{
      tmp_sig_sub_delta_pr <- rep(0,subsample)
    }
    
    # do transformations for non centered parameterization (scale raw subject par by group mean and variance)
    ### Note: deliberately NOT using (1-beta) or delta (-1*delta) bc qwiener does this internally based on resp specified
    tmp_alpha<-pnorm(tmp_mu_grp_alpha_pr + tmp_sub_alpha_pr*tmp_sig_grp_alpha_pr)*3.9+0.1
    tmp_beta<-pnorm(tmp_mu_grp_beta_pr + tmp_sub_beta_pr*tmp_sig_grp_beta_pr)
    tmp_ndt<-pnorm(tmp_mu_grp_ndt_pr + tmp_sub_ndt_pr*tmp_sig_grp_ndt_pr)*tmp_minRT*0.98
    tmp_delta<-tmp_mu_grp_delta_pr + tmp_sub_delta_pr*tmp_sig_grp_delta_pr
    tmp_delta_cond<-pnorm(tmp_delta+tmp_jitter*tmp_sub_b1)*10-5
    
    #define some variables that we will use to calculate the choice probs below
    sigma<-1 #diffusion coefficient. fixed at 1
    tinc<-.01 #time increment (in seconds)
    
    # estimated quantiles
    qc<-foreach(l=1:length(tmp_alpha),.combine='rbind',.packages = c("RWiener","dplyr")) %dopar% {
      
      ## calc absolute start point
      z<-tmp_alpha[l]*tmp_beta[l]
      
      # calculate prob of each choice 
      tmp_choice_prob_upper <- (exp(2*tmp_delta_cond[l]*tmp_alpha[l]/sigma^2) - exp(2*tmp_delta_cond[l]*(tmp_alpha[l]-z)/sigma^2))/(exp(2*tmp_delta_cond[l]*tmp_alpha[l]/sigma^2)-1)
      tmp_choice_prob_lower <- (exp(-2*tmp_delta_cond[l]*tmp_alpha[l]/sigma^2) - exp(-2*tmp_delta_cond[l]*z/sigma^2))/(exp(-2*tmp_delta_cond[l]*tmp_alpha[l]/sigma^2)-1)
      
      #################################################
      # 1) DO UPPER BOUNDARY RESPONSE
      #################################################
      
      m<-1 #iteration tracker
      
      t<-c(tinc+tmp_ndt[l]) #ndt plus some tolerance
      looking<-1 #status of routine within while statement
      pt<-c(NA)
      
      #loop thru all samples, calculating + fixing defective cdf, and terminating when prob from fixed CDF exceeds .9 + .01
      while (looking==1){
        
        #calculate DEFECTIVE cdf for the response given on this trial (note: im deliberately making rt here tinc+tmp_ndt bc pwiener subtracts out the ndt from the rt for us)
        pt[m]<-pwiener(t[m],tmp_alpha[l],tmp_ndt[l],tmp_beta[l],tmp_delta_cond[l],"upper")
        
        #fix the defective CDF (by dividing current prob value by the prob of upper bound choice)
        pt[m]<-(pt[m])/tmp_choice_prob_upper
        test<-as.numeric(pt>max(quantiles)+.01)
        looking<-as.numeric(sum(test)<2)
        m<-m+1
        t[m]<-t[m-1]+tinc
      }
      qc1<-c(NA)
      for (n in 1:length(quantiles)){
        o <- last(which(pt<=quantiles[n]))
        if (length(o)>0){#if it's not empty
          qc1[n] <- (t[o+1] - t[o])*(quantiles[n]-pt[o])/(pt[o+1]-pt[o]) + t[o]
        }else{
          qc1[n] <- tmp_ndt[l]
        }
      }
      
      #append prob of choice as final index in vector
      qc1[length(quantiles)+1]<-tmp_choice_prob_upper
      
      #################################################
      # 2) DO LOWER BOUNDARY RESPONSE
      #################################################
      
      m<-1 #iteration tracker
      
      t<-c(tinc+tmp_ndt[l]) #ndt plus some tolerance
      looking<-1 #status of routine within while statement
      pt<-c(NA)
      
      #loop thru all samples, calculating + fixing defective cdf, and terminating when prob from fixed CDF exceeds .9 + .01
      while (looking==1){
        
        #calculate DEFECTIVE cdf for the response given on this trial (note: im deliberately making rt here tinc+tmp_ndt bc pwiener subtracts out the ndt from the rt for us)
        pt[m]<-pwiener(t[m],tmp_alpha[l],tmp_ndt[l],tmp_beta[l],tmp_delta_cond[l],"lower")
        
        #fix the defective CDF (by dividing current prob value by the prob of upper bound choice)
        pt[m]<-(pt[m])/tmp_choice_prob_lower
        test<-as.numeric(pt>max(quantiles)+.01)
        looking<-as.numeric(sum(test)<2)
        m<-m+1
        t[m]<-t[m-1]+tinc
      }
      qc2<-c(NA)
      for (n in 1:length(quantiles)){
        o <- last(which(pt<=quantiles[n]))
        if (length(o)>0){#if it's not empty
          qc2[n] <- (t[o+1] - t[o])*(quantiles[n]-pt[o])/(pt[o+1]-pt[o]) + t[o]
        }else{
          qc2[n] <- tmp_ndt[l]
        }
      }
      
      #append prob of choice as final index in vector
      qc2[length(quantiles)+1]<-tmp_choice_prob_lower
      
      qc<-c(qc1,qc2)
      
      return(qc)
    }
    
    estim_quantiles_prop[k,,]<-qc
    print(k)
    
    #save quantile estimates every 1000 observations (in case things crash)
    if(as.numeric(endsWith(as.character(unlist(k)), "000"))==1){
      save(estim_quantiles_prop, file = estimated_quantile_file)
      invisible(gc())
    }
  }
  
  # save full estimated quantile data to .RData file
  save(estim_quantiles_prop, file = estimated_quantile_file)
  
}


