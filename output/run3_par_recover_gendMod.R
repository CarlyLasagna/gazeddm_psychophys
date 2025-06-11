#!/usr/bin/env Rscript

# This script runs parameter recovery on the winning model for the gender condition (Model 3)

library(ggplot2)
library(cmdstanr)
library(bayesplot)
library(HDInterval)
library(RWiener)
library(posterior)
library(dplyr)
library(tidyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE) # get arguments from batch script
run<-as.numeric(args[1]) #break up into different concurrent jobs on hpc
group<-group_id<-as.numeric(args[2]) #1=hc, 2=sz, all=both
modelname <-args[3] # "hddm_m[#]_psychophys"
task<-args[4] # task (gaze or gender)
seed<-42+run # different runs get different rand seeds
set.seed(seed)

#############################################################################
# SET DIRECTORIES
#############################################################################

dirname<-'/N/slate/clasagn/gazeddm_psychophys'
datadir<-paste0(dirname,'/data')
scriptsdir<-paste0(dirname,'/scripts')
behfile<-paste0(datadir,"/schizgaze12_gaze_beh.csv")
modeldir<-paste0(dirname,'/output/',modelname,'/',task)
initdir<-paste0(modeldir,'/init_fit')
loglikdir<-paste0(modeldir,'/log_lik')
ppcdir<-paste0(modeldir,'/ppc')
parrecoverdir<-paste0(modeldir,'/par_recover')
finaldir<-paste0(modeldir,'/final_fit')

#############################################################################
# SET MCMC SAMPLER SPECS FOR STAN
#######################################################f#####################

warmup <- 1000        # number of warmup samples (1000)
cores <- 8            # number of cores to use (default=1) (8)
iter <- 500           # number of samples per chains (total postwarmup draws = (iter-warmup)*chains) (500)
chains <- 8           # number of chains (8)
adapt_delta <- 0.95   # default=.95 (target mean proposal acceptance probability during adaptation period)
max_treedepth <- 12   # default=10 (when max_treedepth reached, sampler terminates prematurely)
stepsize <- 1         # default=1 (discretization interval)
options(mc.cores = cores)
setDTthreads(threads=cores)

#############################################################################
# LOAD AND PREPROCESS BEHAVIORAL DATA FOR STAN
#############################################################################

# load preprocessed behavioral data and subset to schizgaze2 only
alldata<-read.csv(behfile)
alldata<-subset(alldata,alldata$Study=="schizgaze2")

# note on behavioral data: gaze angles in this data are originally coded as 2,3,4,5,6,7,8,9,10 (where 2 is most deviated gaze and 10 is most direct)
# sort subj id's in order (hc's [1000s] will come first, then sz [2000s]) so  when we assign sequential index id's, they will be the same for gender and gaze tasks)

#subset to get the appropriate task data
if(task=='gaze'){
  alldata<-subset(alldata,alldata$Task=="Eyes")
  alldata$SubjID<-match(alldata$Subj, unique(alldata$Subj)) #subj IDs to sequential indexes 
  alldata$Resp<-ifelse(alldata$Resp=="Y",1,2) #recode "yes" as 1 and "no as 2
  alldata$Gender<-ifelse(alldata$Gender=="F",1,2) #recode female gender of stim as 1 and male as 2
  alldata$Group<-match(alldata$Group, unique(alldata$Group))
  alldata <- alldata[order(alldata$SubjID, alldata$Resp), ]
}else if(task=='gender'){
  alldata<-subset(alldata,alldata$Task=="GenderID")
  alldata$SubjID<-match(alldata$Subj, unique(alldata$Subj)) #subj IDs to sequential indexes 
  alldata$Resp<-ifelse(alldata$Resp=="F",1,2) #recode male as 1 and female as 2
  alldata$Gender<-ifelse(alldata$Gender=="F",1,2) #recode female gender of stim as 1 and male as 2
  alldata$Group<-match(alldata$Group, unique(alldata$Group))
  alldata <- alldata[order(alldata$SubjID, alldata$Resp), ]
}

minRT<-alldata %>%
  dplyr::group_by(SubjID) %>%
  dplyr::summarise(minRT = min(RT)) %>%
  dplyr::ungroup()

# zscore gaze angle
alldata$z_GazeAngle<-(alldata$GazeAngle-mean(alldata$GazeAngle))/sd(alldata$GazeAngle)

if (max(alldata$RT)>100){alldata$RT<-as.numeric(alldata$RT)/1000} #if rt col is in ms, convert to seconds

#############################################################################
# LOAD IN TRUE POSTERIOR SAMPLES (USED TO DEFINE SAMPLES WE WILL SIMULATE FROM)
#############################################################################

posterior<-as.data.frame(fread(paste0(initdir,'/allgroups_init_fit_samples.csv'),nThread=cores))

# extract the names, means, and sds, of all parameters
all_pars<-colnames(posterior)
all_means<-rowMeans(posterior)
all_sds<-c()

for(i in 1:ncol(posterior)){
  all_sds[i]<-sd(posterior[,i])
}

#############################################################################
# SET DDM SIMULATION SPECS 
#############################################################################

N_trials<-108 #total trials in task
N_subj<-35*2 #35 per group
N_choice<-2
N_gend<-2

N_trials_per_gend<-N_trials/N_gend

subjs<-seq(from=1,to=N_subj,by=1)
choices<-seq(from=1,to=N_choice,by=1)
gends <-unique(alldata$Gender)
gendnames<-c('female','male')

if (file.exists(parrecoverdir) == FALSE){dir.create(parrecoverdir)} #if output directory doesn't exist, create it

# get simulation start and end points (safety check in case job fails on cluster)
simidx_filename<-paste(parrecoverdir,"/sim_start_point.csv",sep="")

if (file.exists(simidx_filename) == FALSE){ #if sim start/end file doesn't exist, create it
  startend_data<-data.frame(run=c(1,2,3,4,5),
                            simstart=c(1,11,21,31,41),
                            simend=c(10,20,30,40,50))
  write.csv(startend_data,simidx_filename,row.names=F)
} 

simidxdata<-read.csv(simidx_filename,row.names = NULL)
simstart<-simidxdata$simstart[run]
simend<-simidxdata$simend[run]

subj_vars<-data.frame(Group=c(rep(1,N_subj/2),rep(2,N_subj/2)),
                      SubjID=seq(from=1,to=N_subj,by=1),
                      Trials=rep(N_trials,N_subj))

# simulate using vals from group parameters (m=mean of grp level mu, sd = standard dev of group level mean dist)
for (i in simstart:simend){
  
  msg<-paste("par_recov_mod:",modelname,"group",group, "run",run,"sim",i,sep=" ")
  print(msg)
  
  for (j in 1:N_subj){
    
    idx<-which(subj_vars$SubjID==j)
    tmp_group<-subj_vars$Group[idx]
    
    # define all group-level posteriors / NOTE: all values (except NDT) are on natural scale. 
    alpha_pr<- mean(posterior[[paste("mu_grp_alpha_pr[",tmp_group,"]",sep="")]]) 
    beta_pr<- mean(posterior[[paste("mu_grp_beta_pr[",tmp_group,"]",sep="")]]) 
    ndt_pr<- 0.2 # (in seconds) we fix it here for parameter recovery
    
    #draw random generating value for this subject from group level posterior (except delta which we do below per gend cond)
    rand_alpha<-pnorm(rnorm(1,mean=alpha_pr,sd=0.1))*3.9+0.1 
    rand_beta<-pnorm(rnorm(1,mean=beta_pr,sd=0.1))
    rand_ndt<-ndt_pr #fixed at 0.2
    
    for(k in 1:N_gend){
      
      tmp_gend<-gends[k]
      tmp_subj<-j
      tmp_gendname<-gendnames[k]
      
      #define delta posterior for this gender condition and draw rand subject drift rate
      delta_pr<- mean(posterior[[paste("mu_grp_delta_",tmp_gendname,"_pr[",tmp_group,"]",sep="")]]) # delta 
      rand_delta<-pnorm(rnorm(1,mean=delta_pr,sd=0.1))*10-5
      
      choice<-rt<-c()
      
      for (m in 1:N_trials_per_gend){
        trial_result<-rwiener(1,rand_alpha,rand_ndt,rand_beta,rand_delta)
        choice[m] <- trial_result$resp[1]
        rt[m] <- trial_result$q[1]
      }
      
      #trial level simulated choice/RT data
      tmp_simdata<-data.frame(sim=i,
                              subj=j,
                              group=tmp_group,
                              gend=tmp_gend,
                              choice=choice,
                              rt=rt)
      
      # group and subject level generating parameters
      tmp_gendata<-data.frame(sim=i,
                              subj=j,
                              group=tmp_group,
                              gend=tmp_gend,
                              grp_genAlpha=pnorm(alpha_pr)*3.9+.1,
                              grp_genBeta=pnorm(beta_pr),
                              grp_genDelta=pnorm(delta_pr)*10-5, 
                              grp_genNDT=ndt_pr, #already in secs so no transform needed
                              sub_genAlpha=rand_alpha,
                              sub_genBeta=rand_beta,
                              sub_genDelta=rand_delta,
                              sub_genNDT=rand_ndt)
      
      if(j==1&k==1){
        all_simdata<-tmp_simdata
        all_gendata<-tmp_gendata
      }else{
        all_simdata<-rbind(all_simdata,tmp_simdata)
        all_gendata<-rbind(all_gendata,tmp_gendata)
      }
    }
  }
  
  isimpath <- paste(parrecoverdir,"/",i,"/",sep="") 
  if (file.exists(isimpath) == FALSE){dir.create(isimpath)} #if output directory doesn't exist, create it
  
  setwd(isimpath)
  
  #save simulated generating values 
  genoutfile <-paste(isimpath,"simgendata_",i,".csv", sep="")
  write.csv(all_gendata,file=genoutfile, row.names = FALSE)
  
  #save simulated behavioral dataset 
  simoutfile <- paste(isimpath,"simtrialdata_",i,".csv", sep="")
  write.csv(all_simdata,file=simoutfile, row.names = FALSE)
  
  # calculate minRT for each subject
  minRT<-all_simdata %>%
    dplyr::group_by(subj) %>%
    dplyr::summarise(minRT = min(rt)) %>%
    dplyr::ungroup()
  
  all_simdata$level<-1
  
  # prep data for stan 
  data_stan<-list(
    N_obs=nrow(all_simdata), # number of observations [integer]
    N_subj=length(unique(all_simdata$subj)),# Number of subjects [integer]
    N_groups=length(unique(all_simdata$group)),# Number of groups [integer]
    N_levels=length(unique(all_simdata$level)), #number of stimulus strength levels
    N_choice=length(unique(all_simdata$choice)), #number of choices [integer]
    RT=all_simdata$rt, # RT in seconds for each observation [vector of doubles of  length 'N_obs']
    choice=all_simdata$choice, # choice for each observation [integer vector of length 'N_obs']
    gender=all_simdata$gend, # gender for each observation [integer vector of length 'N_obs'] 1=female, 2=male
    subj=match(all_simdata$subj, unique(all_simdata$subj)), # subject id for each observation [integer vector length 'N_obs']
    group=all_simdata$group, # group id for each observation 
    level=all_simdata$level, # zscored signal strength (already zscored) for each observation 
    minRT=minRT$minRT, #min rt in seconds for each subject
    rtBound=0.0001) #lower bound on RT in seconds
  
  # define name of model you want to run
  stanmodelname=paste(scriptsdir,"/",modelname,".stan",sep="")
  
  # run the model
  print("Running model in Stan") 
  model <- cmdstan_model(stanmodelname)
  
  fit <- model$sample(data=data_stan, 
                      iter_warmup=warmup, 
                      iter_sampling=iter,
                      init=0, 
                      chains=chains, 
                      parallel_chains=chains, #num cores
                      save_warmup = FALSE,
                      adapt_delta=adapt_delta,
                      max_treedepth=max_treedepth,
                      step_size=stepsize,
                      refresh=100,
                      seed=seed+1,
                      output_dir=isimpath,
                      diagnostics = c("divergences", "treedepth", "ebfmi"))
  
  draws<-fit$draws()
  all_simmcmc <- data.frame(as_draws_df(fit$draws()))
  
  sampler_pars<-fit$sampler_diagnostics(format = "df")
  fit_summary<-fit$summary()
  
  this_genData<-subset(all_gendata,all_gendata$sim==i)
  
  ############################################################################
  # group-level summary of fit diagnostics, group-level generating/recovered values
  ############################################################################
  
  tmp_grp_sim_summary<-data.frame(sim=i,
                                  run=run,
                                  divergences=sum(sampler_pars$divergent__),
                                  alpha_rhat1=fit_summary$rhat[which(fit_summary$variable=="mu_grp_alpha_pr[1]")],
                                  beta_rhat1=fit_summary$rhat[which(fit_summary$variable=="mu_grp_beta_pr[1]")],
                                  ndt_rhat1=NA,
                                  delta_female_rhat1=fit_summary$rhat[which(fit_summary$variable=="mu_grp_delta_female_pr[1]")],
                                  delta_male_rhat1=fit_summary$rhat[which(fit_summary$variable=="mu_grp_delta_male_pr[1]")],
                                  alpha_rhat2=fit_summary$rhat[which(fit_summary$variable=="mu_grp_alpha_pr[2]")],
                                  beta_rhat2=fit_summary$rhat[which(fit_summary$variable=="mu_grp_beta_pr[2]")],
                                  ndt_rhat2=NA,
                                  delta_female_rhat2=fit_summary$rhat[which(fit_summary$variable=="mu_grp_delta_female_pr[2]")],
                                  delta_male_rhat2=fit_summary$rhat[which(fit_summary$variable=="mu_grp_delta_male_pr[2]")],
                                  alpha_ess1=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_alpha_pr[1]")],
                                  beta_ess1=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_beta_pr[1]")],
                                  ndt_ess1=NA,
                                  delta_female_ess1=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_delta_female_pr[1]")],
                                  delta_male_ess1=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_delta_male_pr[1]")],
                                  alpha_ess2=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_alpha_pr[2]")],
                                  beta_ess2=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_beta_pr[2]")],
                                  ndt_ess2=NA,
                                  delta_female_ess2=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_delta_female_pr[2]")],
                                  delta_male_ess2=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_delta_male_pr[2]")],
                                  alpha_grp_gen1=unique(this_genData$grp_genAlpha[which(this_genData$group==1)]), # get group-level generating values
                                  beta_grp_gen1=unique(this_genData$grp_genBeta[which(this_genData$group==1)]),
                                  ndt_grp_gen1=0.2,
                                  delta_female_grp_gen1=unique(this_genData$grp_genDelta[which(this_genData$group==1&this_genData$gend==1)]),
                                  delta_male_grp_gen1=unique(this_genData$grp_genDelta[which(this_genData$group==1&this_genData$gend==2)]),
                                  alpha_grp_gen2=unique(this_genData$grp_genAlpha[which(this_genData$group==2)]), # get group-level generating values
                                  beta_grp_gen2=unique(this_genData$grp_genBeta[which(this_genData$group==2)]),
                                  ndt_grp_gen2=0.2,
                                  delta_female_grp_gen2=unique(this_genData$grp_genDelta[which(this_genData$group==2&this_genData$gend==1)]),
                                  delta_male_grp_gen2=unique(this_genData$grp_genDelta[which(this_genData$group==2&this_genData$gend==2)]),
                                  alpha_grp_sim_mean1=mean(all_simmcmc$mu_alpha.1), # calculate means of posteriors on sim data
                                  beta_grp_sim_mean1=mean(all_simmcmc$mu_beta.1),
                                  ndt_grp_sim_mean1=0.2,
                                  delta_female_grp_sim_mean1=mean(all_simmcmc$mu_delta_female.1),
                                  delta_male_grp_sim_mean1=mean(all_simmcmc$mu_delta_male.1),
                                  alpha_grp_sim_mean2=mean(all_simmcmc$mu_alpha.2),
                                  beta_grp_sim_mean2=mean(all_simmcmc$mu_beta.2),
                                  ndt_grp_sim_mean2=0.2,
                                  delta_female_grp_sim_mean2=mean(all_simmcmc$mu_delta_female.2),
                                  delta_male_grp_sim_mean2=mean(all_simmcmc$mu_delta_male.2),
                                  alpha_grp_sim_hdi_lo1=hdi(all_simmcmc$mu_alpha.1)[1], #calculate hdis of posteriors on sim data
                                  alpha_grp_sim_hdi_hi1=hdi(all_simmcmc$mu_alpha.1)[2],
                                  beta_grp_sim_hdi_lo1=hdi(all_simmcmc$mu_beta.1)[1],
                                  beta_grp_sim_hdi_hi1=hdi(all_simmcmc$mu_beta.1)[2],
                                  ndt_grp_sim_hdi_lo1=0.2,
                                  ndt_grp_sim_hdi_hi1=0.2,
                                  delta_female_grp_sim_hdi_lo1=hdi(all_simmcmc$mu_delta_female.1)[1],
                                  delta_female_grp_sim_hdi_hi1=hdi(all_simmcmc$mu_delta_female.1)[2],
                                  delta_male_grp_sim_hdi_lo1=hdi(all_simmcmc$mu_delta_male.1)[1],
                                  delta_male_grp_sim_hdi_hi1=hdi(all_simmcmc$mu_delta_male.1)[2],
                                  alpha_grp_sim_hdi_lo2=hdi(all_simmcmc$mu_alpha.2)[1], #calculate hdis of posteriors on sim data
                                  alpha_grp_sim_hdi_hi2=hdi(all_simmcmc$mu_alpha.2)[2],
                                  beta_grp_sim_hdi_lo2=hdi(all_simmcmc$mu_beta.2)[1],
                                  beta_grp_sim_hdi_hi2=hdi(all_simmcmc$mu_beta.2)[2],
                                  ndt_grp_sim_hdi_lo2=0.2,
                                  ndt_grp_sim_hdi_hi2=0.2,
                                  delta_female_grp_sim_hdi_lo2=hdi(all_simmcmc$mu_delta_female.2)[1],
                                  delta_female_grp_sim_hdi_hi2=hdi(all_simmcmc$mu_delta_female.2)[2],
                                  delta_male_grp_sim_hdi_lo2=hdi(all_simmcmc$mu_delta_male.2)[1],
                                  delta_male_grp_sim_hdi_hi2=hdi(all_simmcmc$mu_delta_male.2)[2])
  
  # note whether generating group value was recovered
  tmp_grp_sim_summary$alpha_recover1<-ifelse(between(tmp_grp_sim_summary$alpha_grp_gen1,tmp_grp_sim_summary$alpha_grp_sim_hdi_lo1,tmp_grp_sim_summary$alpha_grp_sim_hdi_hi1)==TRUE,1,0)
  tmp_grp_sim_summary$beta_recover1<-ifelse(between(tmp_grp_sim_summary$beta_grp_gen1,tmp_grp_sim_summary$beta_grp_sim_hdi_lo1,tmp_grp_sim_summary$beta_grp_sim_hdi_hi1)==TRUE,1,0)
  tmp_grp_sim_summary$ndt_recover1<-NA
  tmp_grp_sim_summary$delta_female_recover1<-ifelse(between(tmp_grp_sim_summary$delta_female_grp_gen1,tmp_grp_sim_summary$delta_female_grp_sim_hdi_lo1,tmp_grp_sim_summary$delta_female_grp_sim_hdi_hi1)==TRUE,1,0) 
  tmp_grp_sim_summary$delta_male_recover1<-ifelse(between(tmp_grp_sim_summary$delta_male_grp_gen1,tmp_grp_sim_summary$delta_male_grp_sim_hdi_lo1,tmp_grp_sim_summary$delta_male_grp_sim_hdi_hi1)==TRUE,1,0) 
  
  tmp_grp_sim_summary$alpha_recover2<-ifelse(between(tmp_grp_sim_summary$alpha_grp_gen2,tmp_grp_sim_summary$alpha_grp_sim_hdi_lo2,tmp_grp_sim_summary$alpha_grp_sim_hdi_hi2)==TRUE,1,0)
  tmp_grp_sim_summary$beta_recover2<-ifelse(between(tmp_grp_sim_summary$beta_grp_gen2,tmp_grp_sim_summary$beta_grp_sim_hdi_lo2,tmp_grp_sim_summary$beta_grp_sim_hdi_hi2)==TRUE,1,0)
  tmp_grp_sim_summary$ndt_recover2<-NA
  tmp_grp_sim_summary$delta_female_recover2<-ifelse(between(tmp_grp_sim_summary$delta_female_grp_gen2,tmp_grp_sim_summary$delta_female_grp_sim_hdi_lo2,tmp_grp_sim_summary$delta_female_grp_sim_hdi_hi2)==TRUE,1,0) 
  tmp_grp_sim_summary$delta_male_recover2<-ifelse(between(tmp_grp_sim_summary$delta_male_grp_gen2,tmp_grp_sim_summary$delta_male_grp_sim_hdi_lo2,tmp_grp_sim_summary$delta_male_grp_sim_hdi_hi2)==TRUE,1,0) 
  
  #save a copy of the group level par recovery summary in the parent directory after each iteration
  grpsummaryfile<-paste(parrecoverdir,"/parRecovery_groupLvl_summary",".csv",sep="")
  
  if (file.exists(grpsummaryfile) == FALSE){#if file doesn't exist, then save
    write.csv(tmp_grp_sim_summary,grpsummaryfile,row.names = F)
  }else{#if it does exist, then open it, append to it, and resave
    all_grp_sim_summary<-read.csv(grpsummaryfile,row.names=NULL)
    all_grp_sim_summary<-rbind(all_grp_sim_summary,tmp_grp_sim_summary)
    write.csv(all_grp_sim_summary,grpsummaryfile,row.names=FALSE)
  }
  
  ############################################################################
  # subj-level summary of subj-level generating/recovered values
  ############################################################################
  for (j in 1:N_subj){
    
    tmp_genvals<-subset(this_genData,this_genData$subj==j)
    tmp_sub_sim_summary<-data.frame(sim=i,
                                    group=group,
                                    run=run,
                                    subj=j,
                                    alpha_sub_gen=unique(tmp_genvals$sub_genAlpha), # group-level generating values
                                    beta_sub_gen=unique(tmp_genvals$sub_genBeta),
                                    ndt_sub_gen=0.2,
                                    delta_female_sub_gen=unique(tmp_genvals$sub_genDelta[which(tmp_genvals$gend==1)]),
                                    delta_male_sub_gen=unique(tmp_genvals$sub_genDelta[which(tmp_genvals$gend==2)]),
                                    alpha_sub_sim_mean=mean(all_simmcmc[[paste("sub_alpha.",j,".",sep="")]]), # calculate means of posteriors on sim data
                                    beta_sub_sim_mean=mean(all_simmcmc[[paste("sub_beta.",j,".",sep="")]]),
                                    ndt_sub_sim_mean=0.2,
                                    delta_female_sub_sim_mean=mean(pnorm(all_simmcmc[[paste("sub_delta_female.",j,".",sep="")]])*10-5),
                                    delta_male_sub_sim_mean=mean(pnorm(all_simmcmc[[paste("sub_delta_male.",j,".",sep="")]])*10-5),
                                    alpha_sub_sim_hdi_lo=hdi(all_simmcmc[[paste("sub_alpha.",j,".",sep="")]])[1], #calculate hdis of posteriors on sim data
                                    alpha_sub_sim_hdi_hi=hdi(all_simmcmc[[paste("sub_alpha.",j,".",sep="")]])[2],
                                    beta_sub_sim_hdi_lo=hdi(all_simmcmc[[paste("sub_beta.",j,".",sep="")]])[1],
                                    beta_sub_sim_hdi_hi=hdi(all_simmcmc[[paste("sub_beta.",j,".",sep="")]])[2],
                                    ndt_sub_sim_hdi_lo=0.2,
                                    ndt_sub_sim_hdi_hi=0.2,
                                    delta_female_sub_sim_hdi_lo=hdi(pnorm(all_simmcmc[[paste("sub_delta_female.",j,".",sep="")]])*10-5)[1],
                                    delta_female_sub_sim_hdi_hi=hdi(pnorm(all_simmcmc[[paste("sub_delta_female.",j,".",sep="")]])*10-5)[2],
                                    delta_male_sub_sim_hdi_lo=hdi(pnorm(all_simmcmc[[paste("sub_delta_male.",j,".",sep="")]])*10-5)[1],
                                    delta_male_sub_sim_hdi_hi=hdi(pnorm(all_simmcmc[[paste("sub_delta_male.",j,".",sep="")]])*10-5)[2])
    
    # note whether generating subj value was recovered
    tmp_sub_sim_summary$alpha_recover<-ifelse(between(tmp_sub_sim_summary$alpha_sub_gen,tmp_sub_sim_summary$alpha_sub_sim_hdi_lo,tmp_sub_sim_summary$alpha_sub_sim_hdi_hi)==TRUE,1,0)
    tmp_sub_sim_summary$beta_recover<-ifelse(between(tmp_sub_sim_summary$beta_sub_gen,tmp_sub_sim_summary$beta_sub_sim_hdi_lo,tmp_sub_sim_summary$beta_sub_sim_hdi_hi)==TRUE,1,0)
    tmp_sub_sim_summary$ndt_recover<-NA
    tmp_sub_sim_summary$delta_female_recover<-ifelse(between(tmp_sub_sim_summary$delta_female_sub_gen,tmp_sub_sim_summary$delta_female_sub_sim_hdi_lo,tmp_sub_sim_summary$delta_female_sub_sim_hdi_hi)==TRUE,1,0)
    tmp_sub_sim_summary$delta_male_recover<-ifelse(between(tmp_sub_sim_summary$delta_male_sub_gen,tmp_sub_sim_summary$delta_male_sub_sim_hdi_lo,tmp_sub_sim_summary$delta_male_sub_sim_hdi_hi)==TRUE,1,0)
    
    if (j==1){
      sub_sim_summary<-tmp_sub_sim_summary
    }else{
      sub_sim_summary<-rbind(sub_sim_summary,tmp_sub_sim_summary)
    }
  }
  
  #save a copy of the subj level par recovery summary in the parent directory after each iteration
  subsummaryfile<-paste(parrecoverdir,"/parRecovery_subjLvl_summary",".csv",sep="")
  
  if (file.exists(subsummaryfile) == FALSE){#if file doesn't exist, then save
    write.csv(sub_sim_summary,subsummaryfile,row.names = F)
  }else{#if it does exist, then open it, append to it, and resave
    all_sub_sim_summary<-read.csv(subsummaryfile,row.names=NULL)
    all_sub_sim_summary<-rbind(all_sub_sim_summary,sub_sim_summary)
    write.csv(all_sub_sim_summary,subsummaryfile,row.names=FALSE)
  }
  
  #update the simstart value (helps pickup where job left off in case of node fail on cluster). re-read file (in case other runs have completed)
  simidxdata<-read.csv(simidx_filename,row.names = NULL)
  simidxdata$simstart[run]<-i+1
  write.csv(simidxdata,simidx_filename,row.names=FALSE)
  
  rm(list = c("fit", "draws","all_simmcmc"))
  
}

