#!/usr/bin/env Rscript

library(readxl)
library(dplyr)
library(ggplot2)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(data.table)
library(loo)

# run group modelname task
args <- commandArgs(trailingOnly = TRUE) # get arguments from batch script
run<-as.numeric(args[1]) #break up into different concurrent jobs on hpc
group<-as.numeric(args[2]) #1=hc, 2=sz, all=both
modelname <-args[3] # "hddm_m[#]_psychophys"
task<-args[4] # task (gaze or gender or jovi)
seed<-42+run # different runs get different rand seeds
set.seed(seed)

#############################################################################
# SET DIRECTORIES
#############################################################################

dirname<-'[PATH TO FIT DIRECTORY]' 
datadir<-paste0(dirname,'/data')
scriptsdir<-paste0(dirname,'/scripts')
behfile<-paste0(datadir,"/schizgaze2_jovi_beh_all.csv")
outdir<-paste0(dirname,'/output')
modeldir<-paste0(dirname,'/output/',modelname,'/',task)
initdir<-paste0(modeldir,'/init_fit')
loglikdir<-paste0(modeldir,'/log_lik')
ppcdir<-paste0(modeldir,'/ppc')
parrecoverdir<-paste0(modeldir,'/par_recover')
finaldir<-paste0(modeldir,'/final_fit')

#############################################################################
# SET MCMC SAMPLER SPECS FOR STAN
#############################################################################

warmup <- 1500        # number of warmup samples 
cores <- 10           # number of cores to use (default=1) 
iter <- 1000          # number of samples per chain (total postwarmup draws = (iter-warmup)*chains) 
chains <- 10          # number of chains
adapt_delta <- 0.95   # default=.95 (target mean proposal acceptance probability during adaptation period)
max_treedepth <- 10   # default=10 (when max_treedepth reached, sampler terminates prematurely)
stepsize <- 1         # default=1 (discretization interval)

options(mc.cores = cores)
setDTthreads(threads=cores)

#############################################################################
# LOAD AND PREPROCESS BEHAVIORAL DATA FOR STAN
#############################################################################

# load preprocessed behavioral data and subset to schizgaze2 only
alldata<-read.csv(behfile)
alldata<-subset(alldata,alldata$Study=="schizgaze2")

# sort subj id's in order (hc's [1000s] will come first, then sz [2000s])
# (do this so that when we assign sequential index id's, they will be the same for gender and gaze tasks)
alldata <- alldata[order(alldata$Subj), ]

# remove 2 subjects that didn't have valid data from gaze task
alldata<-subset(alldata,alldata$Subj!="2032")
alldata<-subset(alldata,alldata$Subj!="2042")

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
}else if(task=='jovi'){
  alldata$SubjID<-match(alldata$Subj, unique(alldata$Subj)) #subj IDs to sequential indexes 
  alldata$Resp<-ifelse(alldata$Resp=="2",1,2) #recode L as 1 and R as 2
  alldata$Group<-match(alldata$Group, unique(alldata$Group))
  alldata <- alldata[order(alldata$SubjID, alldata$Resp), ]
}

minRT<-alldata %>%
  dplyr::group_by(SubjID) %>%
  dplyr::summarise(minRT = min(RT)) %>%
  dplyr::ungroup()

# prep data for stan
data_stan<-list(
  N_obs=nrow(alldata), # number of observations 
  N_subj=length(unique(alldata$Subj)),# Number of subjects 
  N_groups=length(unique(alldata$Group)),# Number of groups 
  N_jitters=length(unique(alldata$Jitter)), #number of jitter/stimulus strength levels
  N_choice=length(unique(alldata$Resp)), #number of choices 
  RT=alldata$RT, # RT in seconds for each observation 
  choice=alldata$Resp, # choice for each observation; 1=L,2=R
  subj=match(alldata$SubjID, unique(alldata$SubjID)), # subject id for each observation
  group=alldata$Group,   # group id for each observation
  jitter=alldata$Jitter, # jitter for each trial (already zscored)
  minRT=minRT$minRT, # min rt in seconds for each subject
  censored=alldata$censored, #1=trial not censored, 0=trial is censored  (col already exists in data)
  rtBound=0.0001) #lower bound on RT in seconds


data_stan<-list(
  N_obs=nrow(alldata), # number of observations [integer]
  N_subj=length(unique(alldata$Subj)),# Number of subjects [integer]
  N_groups=length(unique(alldata$Group)),# Number of groups [integer]
  N_jitters=length(unique(alldata$Jitter)), #number of stimulus strength levels
  N_choice=length(unique(alldata$Resp)), #number of choices [integer]
  RT=alldata$RT, # RT in seconds for each observation [vector of doubles of  length 'N_obs']
  choice=alldata$Resp, # choice for each observation [integer vector of length 'N_obs']; 1=L,2=R
  subj=match(alldata$SubjID, unique(alldata$SubjID)), # subject id for each observation [integer vector length 'N_obs']
  group=alldata$Group, # group id for each observation
  jitter=alldata$Jitter, # jitter is already zscored
  minRT=minRT$minRT, #min rt in seconds for each subject
  
  rtBound=0.0001) #lower bound on RT in seconds

stanmodelname=paste(scriptsdir,"/",modelname,".stan",sep="")

# run the model
print("Running model in Stan") 
model <- cmdstan_model(stanmodelname)

fit <- model$sample(data=data_stan, 
                    iter_warmup=warmup, 
                    iter_sampling=iter,
                    chains=chains, 
                    parallel_chains=cores, #num cores
                    save_warmup = FALSE,
                    adapt_delta=adapt_delta,
                    max_treedepth=max_treedepth,
                    step_size=stepsize,
                    init=0,
                    refresh=100,
                    seed=seed,
                    output_dir=initdir,
                    diagnostics = c("divergences", "treedepth", "ebfmi"))

#############################################################################
# EXTRACT/SAVE SAMPLE DATA 
#############################################################################

posterior_df <- as_draws_df(fit$draws())

# get parameter samples and save
#samples_cols <- posterior_df[, -grep("log_lik", names(posterior_df))]
samples_cols <- posterior_df
sample_file<-paste0(initdir,'/allgroups_init_fit_samples.csv')
write.csv(samples_cols,sample_file,row.names = FALSE)
