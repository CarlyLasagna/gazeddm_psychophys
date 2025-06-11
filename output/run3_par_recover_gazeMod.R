#!/usr/bin/env Rscript

# This script runs parameter recovery on the winning model for the gaze condition (Model 2)

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

warmup <- 1000        # number of warmup samples
cores <- 8            # number of cores to use (default=1)
iter <- 500           # number of samples per chains (total postwarmup draws = (iter-warmup)*chains)
chains <- 8           # number of chains
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
  alldata$Resp<-ifelse(alldata$Resp=="F",1,2) #recode female as 1 and male as 2
  alldata$Gender<-ifelse(alldata$Gender=="F",1,2) #recode female gender of stim as 1 and male as 2
  alldata$Group<-match(alldata$Group, unique(alldata$Group))
  alldata <- alldata[order(alldata$SubjID, alldata$Resp), ]
}

minRT<-alldata %>%
  dplyr::group_by(SubjID) %>%
  dplyr::summarise(minRT = min(RT)) %>%
  dplyr::ungroup()

alldata$z_GazeAngle<-(alldata$GazeAngle-mean(alldata$GazeAngle))/sd(alldata$GazeAngle) # zscore gaze angle
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
N_level<-9

N_trials_per_level<-N_trials/N_level

subjs<-seq(from=1,to=N_subj,by=1)
choices<-seq(from=1,to=N_choice,by=1)
levels <-unique(alldata$z_GazeAngle)
levels<-levels[order(levels)]

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

# get average # yes responses at each signal strength
alldata$RespProp<-ifelse(alldata$RespBound==2,0,alldata$RespBound)
group_vars_yes <- alldata %>%
  dplyr::group_by(Group,z_GazeAngle) %>%
  dplyr::summarize(mean_prop = mean(RespProp)) 
group_vars_yes$counts<-round(group_vars_yes$mean_prop*N_trials_per_level,0)
group_vars_yes$Resp<-1
group_vars_no<-group_vars_yes
group_vars_no$Resp<-2
group_vars_no$counts<- N_trials_per_level-group_vars_no$counts
group_vars<-rbind(group_vars_yes,group_vars_no)

subj_vars<-data.frame(Group=c(rep(1,N_subj/2),rep(2,N_subj/2)),
                      SubjID=seq(from=1,to=N_subj,by=1),
                      Trials=rep(N_trials,N_subj))

# get the counts of all trials for each subj, level, resp
cond_vars <- alldata %>%
  dplyr::group_by(Group,SubjID, z_GazeAngle, Resp) %>%
  dplyr::summarize(Trials = n(), .groups = 'drop') %>%   # Ungroup the data automatically
  tidyr::complete(Group,SubjID,z_GazeAngle, Resp, fill = list(n = 0))

cond_vars$Trials[is.na(cond_vars$Trials)] <- 0

# simulate using vals from group parameters (m=mean of grp level mu, sd = standard dev of group level mean dist)
for (i in simstart:simend){
  
  msg<-paste("par_recov_mod:",modelname,"group",group, "run",run,"sim",i,sep=" ")
  print(msg)
  
  for (j in 1:N_subj){
    
    # get this subject's group
    idx<-which(subj_vars$SubjID==j)
    tmp_group<-subj_vars$Group[idx]
    
    # define all group-level posteriors / NOTE: all values (except NDT) are on natural scale. 
    alpha_pr<- mean(posterior[[paste("mu_grp_alpha_pr[",tmp_group,"]",sep="")]]) 
    beta_pr<- mean(posterior[[paste("mu_grp_beta_pr[",tmp_group,"]",sep="")]]) 
    ndt_pr<- 0.2 # (in seconds) we fix it here for parameter recovery
    b1_pr<- mean(posterior[[paste("mu_grp_b1_pr[",tmp_group,"]",sep="")]])
    delta_pr<- mean(posterior[[paste("mu_grp_delta_pr[",tmp_group,"]",sep="")]]) # delta bias
    
    #draw random generating value for this subject from group level posterior
    rand_alpha<-pnorm(rnorm(1,mean=alpha_pr,sd=0.1))*3.9+0.1 
    rand_beta<-pnorm(rnorm(1,mean=beta_pr,sd=0.1))
    rand_ndt<-ndt_pr #fixed at 0.2
    rand_b1<-rnorm(1,mean=b1_pr,sd=0.05) 
    rand_delta_bias <- rnorm(1,mean=delta_pr,sd=0.1) 
    
      for (l in 1:length(levels)){
        
        tmp_trials<-N_trials_per_level
        
        if(tmp_trials>0){ #if subject has at least 1 trial for this level/choice 
          
          tmp_subj<-j
          tmp_level<-levels[l]
          rand_delta <- pnorm(rand_delta_bias+rand_b1*tmp_level)*10-5 
          
          choice<-rt<-c()
          
          for (m in 1:tmp_trials){
            trial_result<-rwiener(1,rand_alpha,rand_ndt,rand_beta,rand_delta)
            choice[m] <- trial_result$resp[1]
            rt[m] <- trial_result$q[1]
          }
          
          #trial level simulated choice/RT data
          tmp_simdata<-data.frame(sim=i,
                                  subj=j,
                                  group=tmp_group,
                                  level=tmp_level,
                                  choice=choice,
                                  rt=rt)
          
          # group and subject level generating parameters
          tmp_gendata<-data.frame(sim=i,
                                  subj=j,
                                  group=tmp_group,
                                  level=tmp_level,
                                  grp_genAlpha=pnorm(alpha_pr)*3.9+.1,
                                  grp_genBeta=pnorm(beta_pr),
                                  grp_genDelta_bias=pnorm(delta_pr)*10-5, 
                                  grp_genNDT=ndt_pr, #already in secs so no transform needed
                                  grp_genB1=b1_pr,
                                  sub_genAlpha=rand_alpha,
                                  sub_genBeta=rand_beta,
                                  sub_genDelta_bias=pnorm(rand_delta_bias)*10-5,
                                  sub_genDelta_transform=rand_delta,
                                  sub_genNDT=rand_ndt,
                                  sub_genB1=rand_b1)
  
          if(j==1&l==1){
            all_simdata<-tmp_simdata
            all_gendata<-tmp_gendata
          }else{
            all_simdata<-rbind(all_simdata,tmp_simdata)
            all_gendata<-rbind(all_gendata,tmp_gendata)
          }
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
  
  # prep data for stan 
  data_stan<-list(
    N_obs=nrow(all_simdata), # number of observations [integer]
    N_subj=length(unique(all_simdata$subj)),# Number of subjects [integer]
    N_groups=length(unique(all_simdata$group)),# Number of groups [integer]
    N_levels=length(unique(all_simdata$level)), #number of stimulus strength levels
    N_choice=length(unique(all_simdata$choice)), #number of choices [integer]
    RT=all_simdata$rt, # RT in seconds for each observation [vector of doubles of  length 'N_obs']
    choice=all_simdata$choice, # choice for each observation [integer vector of length 'N_obs']
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
                                  delta_rhat1=fit_summary$rhat[which(fit_summary$variable=="mu_grp_delta_pr[1]")],
                                  b1_rhat1=fit_summary$rhat[which(fit_summary$variable=="mu_grp_b1_pr[1]")],
                                  alpha_rhat2=fit_summary$rhat[which(fit_summary$variable=="mu_grp_alpha_pr[2]")],
                                  beta_rhat2=fit_summary$rhat[which(fit_summary$variable=="mu_grp_beta_pr[2]")],
                                  ndt_rhat2=NA,
                                  delta_rhat2=fit_summary$rhat[which(fit_summary$variable=="mu_grp_delta_pr[2]")],
                                  b1_rhat2=fit_summary$rhat[which(fit_summary$variable=="mu_grp_b1_pr[2]")],
                                  alpha_ess1=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_alpha_pr[1]")],
                                  beta_ess1=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_beta_pr[1]")],
                                  ndt_ess1=NA,
                                  delta_ess1=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_delta_pr[1]")],
                                  b1_ess1=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_b1_pr[1]")],
                                  alpha_ess2=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_alpha_pr[2]")],
                                  beta_ess2=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_beta_pr[2]")],
                                  ndt_ess2=NA,
                                  delta_ess2=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_delta_pr[2]")],
                                  b1_ess2=fit_summary$ess_bulk[which(fit_summary$variable=="mu_grp_b1_pr[2]")],
                                  alpha_grp_gen1=unique(this_genData$grp_genAlpha[which(this_genData$group==1)]), # get group-level generating values
                                  beta_grp_gen1=unique(this_genData$grp_genBeta[which(this_genData$group==1)]),
                                  ndt_grp_gen1=0.2,
                                  delta_grp_gen1=unique(this_genData$grp_genDelta_bias[which(this_genData$group==1)]),
                                  b1_grp_gen1=unique(this_genData$grp_genB1[which(this_genData$group==1)]),
                                  alpha_grp_gen2=unique(this_genData$grp_genAlpha[which(this_genData$group==2)]), # get group-level generating values
                                  beta_grp_gen2=unique(this_genData$grp_genBeta[which(this_genData$group==2)]),
                                  ndt_grp_gen2=0.2,
                                  delta_grp_gen2=unique(this_genData$grp_genDelta_bias[which(this_genData$group==2)]),
                                  b1_grp_gen2=unique(this_genData$grp_genB1[which(this_genData$group==2)]),
                                  alpha_grp_sim_mean1=mean(all_simmcmc$mu_alpha.1), # calculate means of posteriors on sim data
                                  beta_grp_sim_mean1=mean(all_simmcmc$mu_beta.1),
                                  ndt_grp_sim_mean1=0.2,
                                  delta_grp_sim_mean1=mean(all_simmcmc$mu_delta.1),
                                  b1_grp_sim_mean1=mean(all_simmcmc$mu_grp_b1_pr.1),
                                  alpha_grp_sim_mean2=mean(all_simmcmc$mu_alpha.2),
                                  beta_grp_sim_mean2=mean(all_simmcmc$mu_beta.2),
                                  ndt_grp_sim_mean2=0.2,
                                  delta_grp_sim_mean2=mean(all_simmcmc$mu_delta.2),
                                  b1_grp_sim_mean2=mean(all_simmcmc$mu_grp_b1_pr.2),
                                  alpha_grp_sim_hdi_lo1=hdi(all_simmcmc$mu_alpha.1)[1], #calculate hdis of posteriors on sim data
                                  alpha_grp_sim_hdi_hi1=hdi(all_simmcmc$mu_alpha.1)[2],
                                  beta_grp_sim_hdi_lo1=hdi(all_simmcmc$mu_beta.1)[1],
                                  beta_grp_sim_hdi_hi1=hdi(all_simmcmc$mu_beta.1)[2],
                                  ndt_grp_sim_hdi_lo1=0.2,
                                  ndt_grp_sim_hdi_hi1=0.2,
                                  delta_grp_sim_hdi_lo1=hdi(all_simmcmc$mu_delta.1)[1],
                                  delta_grp_sim_hdi_hi1=hdi(all_simmcmc$mu_delta.1)[2],
                                  b1_grp_sim_hdi_lo1=hdi(all_simmcmc$mu_grp_b1_pr.1)[1],
                                  b1_grp_sim_hdi_hi1=hdi(all_simmcmc$mu_grp_b1_pr.1)[2],
                                  alpha_grp_sim_hdi_lo2=hdi(all_simmcmc$mu_alpha.2)[1], #calculate hdis of posteriors on sim data
                                  alpha_grp_sim_hdi_hi2=hdi(all_simmcmc$mu_alpha.2)[2],
                                  beta_grp_sim_hdi_lo2=hdi(all_simmcmc$mu_beta.2)[1],
                                  beta_grp_sim_hdi_hi2=hdi(all_simmcmc$mu_beta.2)[2],
                                  ndt_grp_sim_hdi_lo2=0.2,
                                  ndt_grp_sim_hdi_hi2=0.2,
                                  delta_grp_sim_hdi_lo2=hdi(all_simmcmc$mu_delta.2)[1],
                                  delta_grp_sim_hdi_hi2=hdi(all_simmcmc$mu_delta.2)[2],
                                  b1_grp_sim_hdi_lo2=hdi(all_simmcmc$mu_grp_b1_pr.2)[1],
                                  b1_grp_sim_hdi_hi2=hdi(all_simmcmc$mu_grp_b1_pr.2)[2])
  
  # note whether generating group value was recovered
  tmp_grp_sim_summary$alpha_recover1<-ifelse(between(tmp_grp_sim_summary$alpha_grp_gen1,tmp_grp_sim_summary$alpha_grp_sim_hdi_lo1,tmp_grp_sim_summary$alpha_grp_sim_hdi_hi1)==TRUE,1,0)
  tmp_grp_sim_summary$beta_recover1<-ifelse(between(tmp_grp_sim_summary$beta_grp_gen1,tmp_grp_sim_summary$beta_grp_sim_hdi_lo1,tmp_grp_sim_summary$beta_grp_sim_hdi_hi1)==TRUE,1,0)
  tmp_grp_sim_summary$ndt_recover1<-NA
  tmp_grp_sim_summary$delta_recover1<-ifelse(between(tmp_grp_sim_summary$delta_grp_gen1,tmp_grp_sim_summary$delta_grp_sim_hdi_lo1,tmp_grp_sim_summary$delta_grp_sim_hdi_hi1)==TRUE,1,0) 
  tmp_grp_sim_summary$b1_recover1<-ifelse(between(tmp_grp_sim_summary$b1_grp_gen1,tmp_grp_sim_summary$b1_grp_sim_hdi_lo1,tmp_grp_sim_summary$b1_grp_sim_hdi_hi1)==TRUE,1,0) 
  
  tmp_grp_sim_summary$alpha_recover2<-ifelse(between(tmp_grp_sim_summary$alpha_grp_gen2,tmp_grp_sim_summary$alpha_grp_sim_hdi_lo2,tmp_grp_sim_summary$alpha_grp_sim_hdi_hi2)==TRUE,1,0)
  tmp_grp_sim_summary$beta_recover2<-ifelse(between(tmp_grp_sim_summary$beta_grp_gen2,tmp_grp_sim_summary$beta_grp_sim_hdi_lo2,tmp_grp_sim_summary$beta_grp_sim_hdi_hi2)==TRUE,1,0)
  tmp_grp_sim_summary$ndt_recover2<-NA
  tmp_grp_sim_summary$delta_recover2<-ifelse(between(tmp_grp_sim_summary$delta_grp_gen2,tmp_grp_sim_summary$delta_grp_sim_hdi_lo2,tmp_grp_sim_summary$delta_grp_sim_hdi_hi2)==TRUE,1,0) 
  tmp_grp_sim_summary$b1_recover2<-ifelse(between(tmp_grp_sim_summary$b1_grp_gen2,tmp_grp_sim_summary$b1_grp_sim_hdi_lo2,tmp_grp_sim_summary$b1_grp_sim_hdi_hi2)==TRUE,1,0) 
  
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
                                    delta_sub_gen=unique(tmp_genvals$sub_genDelta_bias),
                                    b1_sub_gen=unique(tmp_genvals$sub_genB1),
                                    alpha_sub_sim_mean=mean(all_simmcmc[[paste("sub_alpha.",j,".",sep="")]]), # calculate means of posteriors on sim data
                                    beta_sub_sim_mean=mean(all_simmcmc[[paste("sub_beta.",j,".",sep="")]]),
                                    ndt_sub_sim_mean=0.2,
                                    delta_sub_sim_mean=mean(pnorm(all_simmcmc[[paste("sub_delta.",j,".",sep="")]])*10-5),
                                    b1_sub_sim_mean=mean(all_simmcmc[[paste("sub_b1.",j,".",sep="")]]),
                                    alpha_sub_sim_hdi_lo=hdi(all_simmcmc[[paste("sub_alpha.",j,".",sep="")]])[1], #calculate hdis of posteriors on sim data
                                    alpha_sub_sim_hdi_hi=hdi(all_simmcmc[[paste("sub_alpha.",j,".",sep="")]])[2],
                                    beta_sub_sim_hdi_lo=hdi(all_simmcmc[[paste("sub_beta.",j,".",sep="")]])[1],
                                    beta_sub_sim_hdi_hi=hdi(all_simmcmc[[paste("sub_beta.",j,".",sep="")]])[2],
                                    ndt_sub_sim_hdi_lo=0.2,
                                    ndt_sub_sim_hdi_hi=0.2,
                                    delta_sub_sim_hdi_lo=hdi(pnorm(all_simmcmc[[paste("sub_delta.",j,".",sep="")]])*10-5)[1],
                                    delta_sub_sim_hdi_hi=hdi(pnorm(all_simmcmc[[paste("sub_delta.",j,".",sep="")]])*10-5)[2],
                                    b1_sub_sim_hdi_lo=hdi(all_simmcmc[[paste("sub_b1.",j,".",sep="")]])[1],
                                    b1_sub_sim_hdi_hi=hdi(all_simmcmc[[paste("sub_b1.",j,".",sep="")]])[2])
    
    # note whether generating subj value was recovered
    tmp_sub_sim_summary$alpha_recover<-ifelse(between(tmp_sub_sim_summary$alpha_sub_gen,tmp_sub_sim_summary$alpha_sub_sim_hdi_lo,tmp_sub_sim_summary$alpha_sub_sim_hdi_hi)==TRUE,1,0)
    tmp_sub_sim_summary$beta_recover<-ifelse(between(tmp_sub_sim_summary$beta_sub_gen,tmp_sub_sim_summary$beta_sub_sim_hdi_lo,tmp_sub_sim_summary$beta_sub_sim_hdi_hi)==TRUE,1,0)
    tmp_sub_sim_summary$ndt_recover<-NA
    tmp_sub_sim_summary$delta_recover<-ifelse(between(tmp_sub_sim_summary$delta_sub_gen,tmp_sub_sim_summary$delta_sub_sim_hdi_lo,tmp_sub_sim_summary$delta_sub_sim_hdi_hi)==TRUE,1,0)
    tmp_sub_sim_summary$b1_recover<-ifelse(between(tmp_sub_sim_summary$b1_sub_gen,tmp_sub_sim_summary$b1_sub_sim_hdi_lo,tmp_sub_sim_summary$b1_sub_sim_hdi_hi)==TRUE,1,0)
    
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

