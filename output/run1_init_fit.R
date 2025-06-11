#!/usr/bin/env Rscript

# This script runs initial DDMs in Stan for the gaze/gender conditions, saves the outputs, and does loo calculations by group

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
seed<-42+run # different runs get different rand seeds
set.seed(seed)

#############################################################################
# SET DIRECTORIES
#############################################################################

dirname<-'/N/slate/clasagn/gazeddm_psychophys' ## IU cluster
datadir<-paste0(dirname,'/data')
scriptsdir<-paste0(dirname,'/scripts')
behfile<-paste0(datadir,"/schizgaze12_gaze_beh.csv")
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

# sort subj id's in order (hc's [1000s] will come first, then sz [2000s]). (do this so that when we assign sequential index id's, they will be the same for gender and gaze tasks)
alldata <- alldata[order(alldata$Subj), ]

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

# prep data for stan 
data_stan<-list(
  N_obs=nrow(alldata), # Number of observations [integer]
  N_subj=length(unique(alldata$Subj)),# Number of subjects [integer]
  N_groups=length(unique(alldata$Group)),# Number of groups [integer]
  N_levels=length(unique(alldata$GazeAngle)), # Number of stimulus strength levels [integer]
  N_choice=length(unique(alldata$Resp)), # Number of choices [integer]
  RT=alldata$RT/1000, # RT in seconds for each observation [vector of doubles of  length 'N_obs']
  choice=alldata$Resp, # Choice for each observation [integer vector of length 'N_obs']
  gender=alldata$Gender, # Gender for each observation [integer vector of length 'N_obs'] 1=female, 2=male
  subj=match(alldata$SubjID, unique(alldata$SubjID)), # Subject id for each observation [integer vector length 'N_obs']
  group=alldata$Group, # Group id for each observation [integer vector length 'N_obs']
  level=(alldata$GazeAngle-mean(alldata$GazeAngle))/sd(alldata$GazeAngle), # zscored signal strength for each observation [vector of reals; length = N_obs]
  minRT=minRT$minRT/1000, # Min rt in seconds for each subject [vector of reals; length = N_subj]
  rtBound=0.0001) # Lower bound on RT in seconds [real]

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
# EXTRACT/SAVE LOGLIK AND SAMPLE DATA FROM INITIAL FIT
#############################################################################

posterior_df <- as_draws_df(fit$draws())

# get parameter samples and save
samples_cols <- posterior_df[, -grep("log_lik", names(posterior_df))]
sample_file<-paste0(initdir,'/allgroups_init_fit_samples.csv')
write.csv(samples_cols,sample_file,row.names = FALSE)

# get loglik for each trial and save
log_lik_cols <- posterior_df[, grep(".chain|log_lik", names(posterior_df))]
loglik_file<-paste0(loglikdir,'/allgroups_init_fit_rawloglik.csv')
write.csv(log_lik_cols,loglik_file,row.names = FALSE)

#############################################################################
# CALCULATE LOO
#############################################################################

# get loglik data
chain <- posterior_df[, grep(".chain", names(posterior_df))]
log_lik_cols <- posterior_df[, grep("log_lik", names(posterior_df))]

# CALC LOO FROM LOGLIK FOR ALL SUBJECTS COMBINED
temp<-as.matrix(log_lik_cols)
reff<-relative_eff(exp(temp),chain_id=chain$.chain,cores = getOption("mc.cores", cores))

loo_1 <- loo(temp, r_eff = reff, cores = cores)
loo_2<-list(loo=loo_1,chain=chain$.chain)
filename<-paste0(loglikdir,'/loglik_all.RData')
save(loo_2,file = filename)

# CALC LOO FROM LOGLIK FOR ALL SEPARATE GROUPS
for (i in 1:length(unique(alldata$Group))){
  temp_grp_ids<-which(alldata$Group == i)
  temp<-as.matrix(log_lik_cols[,temp_grp_ids])
  reff<-relative_eff(exp(temp),chain_id=chain$.chain,cores = getOption("mc.cores", cores))
  loo_1 <- loo(temp, r_eff = reff, cores = cores)
  loo_2<-list(loo=loo_1,chain=chain$.chain)
  filename<-paste0(loglikdir,'/loglik_group',i,'.RData')
  save(loo_2,file = filename)
  rm(temp,temp_grp_ids,loo_1,loo_2,reff)
  gc()
}
