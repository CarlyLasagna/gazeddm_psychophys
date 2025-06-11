#!/usr/bin/env Rscript

# This supplemental model re-runs the winning gender model (model 3)
# separately for male and female subjects to ensure that results were not the by-product of participant gender.

library(readxl)
library(dplyr)
library(ggplot2)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(data.table)
library(loo)

args <- commandArgs(trailingOnly = TRUE) # get arguments from batch script
run<-as.numeric(args[1]) #break up into different concurrent jobs on hpc
group<-as.numeric(args[2]) #1=hc, 2=sz, all=both
modelname <-args[3] # "hddm_m[#]_psychophys"
task<-args[4] # task (gaze or gender)
gend_group<-args[5] # which gender group to run (female_only or male_only)
seed<-42+run # different runs get different rand seeds
set.seed(seed)

#############################################################################
# SET DIRECTORIES
#############################################################################

dirname<-'/N/slate/clasagn/gazeddm_psychophys' 
datadir<-paste0(dirname,'/data')
scriptsdir<-paste0(dirname,'/scripts')
behfile<-paste0(datadir,"/schizgaze12_gaze_beh.csv")
corrfile<-paste0(datadir,"/schizgaze2_preproc_correlates.csv")
outdir<-paste0(dirname,'/output')
modeldir<-paste0(dirname,'/output/',modelname,'/',task)
initdir<-paste0(modeldir,'/init_fit')
supp_bygend_dir<-paste0(modeldir,'/supplement_bygend/',gend_group)
loglikdir<-paste0(modeldir,'/log_lik')
ppcdir<-paste0(modeldir,'/ppc')
parrecoverdir<-paste0(modeldir,'/par_recover')
finaldir<-paste0(modeldir,'/final_fit')

#############################################################################
# SET MCMC SAMPLER SPECS FOR STAN
#############################################################################

warmup <- 1500        # number of warmup samples (1500)
cores <- 10           # number of cores to use (default=1) (10)
iter <- 1000          # number of samples per chain (total postwarmup draws = (iter-warmup)*chains) (1000)
chains <- 10          # number of chains  (10)
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

alldata<-subset(alldata,alldata$Task=="GenderID")
alldata$SubjID<-match(alldata$Subj, unique(alldata$Subj)) #subj IDs to sequential indexes 
alldata$Resp<-ifelse(alldata$Resp=="F",1,2) #recode female as 1 and male as 2
alldata$Gender<-ifelse(alldata$Gender=="F",1,2) #recode female gender of stim as 1 and male as 2
alldata$Group<-match(alldata$Group, unique(alldata$Group))
alldata <- alldata[order(alldata$SubjID, alldata$Resp), ]

#########################################################################
#########################################################################
# NEW! For supplement byGend analysis. 
## Add each subject's gender to task data, subset to current participant 
## gender group ("gend_group") and run model separate by sub gender,
#########################################################################

## Load schizgaze2 gender data
demodata<-read.csv(corrfile)
demodata$gender<-ifelse(demodata$gender=="Female","1",demodata$gender)
demodata$gender<-ifelse(demodata$gender=="Male","2",demodata$gender)

alldata$Subj_Gender<-NA

## Add each subject's gender to task data file
for(i in 1:nrow(alldata)){
  temp_sub<-alldata$Subj[i]
  temp_sub_row<-which(demodata$subj==temp_sub)
  if(length(temp_sub_row)>0){
    alldata$Subj_Gender[i]<-demodata$gender[temp_sub_row]
  }
}

# subset to include only current gender group
if (gend_group=='female_only'){
  alldata<-subset(alldata,alldata$Subj_Gender==1)
}else if(gend_group=='male_only'){
  alldata<-subset(alldata,alldata$Subj_Gender==2)
}else{
  alldata<-NA #otherwise make data NA so it throws an error
}

minRT<-alldata %>%
  dplyr::group_by(SubjID) %>%
  dplyr::summarise(minRT = min(RT)) %>%
  dplyr::ungroup()

# prep data for stan 
data_stan<-list(
  N_obs=nrow(alldata), # number of observations [integer]
  N_subj=length(unique(alldata$Subj)),# Number of subjects [integer]
  N_groups=length(unique(alldata$Group)),# Number of groups [integer]
  N_levels=length(unique(alldata$GazeAngle)), #number of stimulus strength levels
  N_choice=length(unique(alldata$Resp)), #number of choices [integer]
  RT=alldata$RT/1000, # RT in seconds for each observation [vector of doubles of  length 'N_obs']
  choice=alldata$Resp, # choice for each observation [integer vector of length 'N_obs']
  gender=alldata$Gender, # gender for each observation [integer vector of length 'N_obs'] 1=female, 2=male
  subj=match(alldata$SubjID, unique(alldata$SubjID)), # subject id for each observation [integer vector length 'N_obs']
  group=alldata$Group, # group id for each observation
  level=(alldata$GazeAngle-mean(alldata$GazeAngle))/sd(alldata$GazeAngle), # zscored signal strength for each observation (will transform back to original units post-fitting)
  minRT=minRT$minRT/1000, #min rt in seconds for each subject
  rtBound=0.0001) #lower bound on RT in seconds

stanmodelname=paste(scriptsdir,"/",modelname,".stan",sep="")

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
                    output_dir=supp_bygend_dir,
                    diagnostics = c("divergences", "treedepth", "ebfmi"))

#############################################################################
# EXTRACT/SAVE LOGLIK AND SAMPLE DATA FROM INITIAL FIT
#############################################################################

posterior_df <- as_draws_df(fit$draws())

# get parameter samples and save
samples_cols <- posterior_df[, -grep("log_lik", names(posterior_df))]
sample_file<-paste0(supp_bygend_dir,'/allgroups_init_fit_samples.csv')
write.csv(samples_cols,sample_file,row.names = FALSE)
