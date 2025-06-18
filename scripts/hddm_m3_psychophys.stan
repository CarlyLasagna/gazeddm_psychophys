// PSYCHOPHYS STAN MODEL 3 (Gender Models Only)
// This program defines a hierarchical drift diffusion model (Ratcliff 2008), 
// with a non-centered parameterization. Female responses are 
// modeled as upper boundary responses and Male responses are modeled as
// lower boundary responses. In this version, drift rate (delta) is NOT
// modulated by the signal strength of eye gaze of the stimuli, 
// In this version, the magnitude of drift rates for stimuli of diff genders
// are allowed to vary. 

data {
  int N_obs;               // number of observations [single integer]                           
  int N_subj;              // number of subjects [single integer]
  int N_levels;            // number of stimuli signal strengths (i.e., gaze angles) [single integer]
  array[N_obs] real level;  // self referential signal strength for each trial (9 levels in .1 increments from .2 to 1; zscored)
  int N_choice;            // number of choice alternatives [single integer]
  int N_groups;            // number of diagnostic groups [single integer]
  array[N_obs] real RT;     // RT in seconds for each trial [numeric vector; length N_obs]
  array[N_obs] int subj;    // subj id for each trial [integer vector; length N_obs]
  array[N_obs] int choice;  // response for each trial [integer vector; length N_obs]
  array[N_obs] int gender;  // gender id of stim for each trial [integer vector; length N_obs];1=male;2=female
  array[N_obs] int group;   // group id for each trial [integer vector; length N_obs]
  array[N_subj] real minRT; // minimum RT for each subject [vector of reals; length N_subj]
  real rtBound;             // lowest rt allowed [single number] 
}
transformed data {
  array[N_subj] int subj_group;   // gives a vector of group id's at the subject-level (as opposed to the observation level as in 'group' above)
  for (i in 1:N_obs){
    subj_group[subj[i]]=group[i];
  }
}

parameters { 
  
  // GROUP-level parameters
  vector[N_groups] mu_grp_alpha_pr;           // threshold sep. group mean
  vector[N_groups] mu_grp_beta_pr;            // start point group mean
  vector[N_groups] mu_grp_delta_female_pr;    // drift rate group mean, gender FEMALE conditions
  vector[N_groups] mu_grp_delta_male_pr;      // drift rate group mean, gender MALE conditions
  vector[N_groups] mu_grp_ndt_pr;             // non-decision time group mean 
  
  vector<lower=0>[N_groups] sig_grp_alpha_pr; // threshold sep. group SD 
  vector<lower=0>[N_groups] sig_grp_beta_pr;  // start point group SD
  vector<lower=0>[N_groups] sig_grp_delta_pr; // drift rate group SD
  vector<lower=0>[N_groups] sig_grp_ndt_pr;   // non-decision time group SD
  
  //SUBJECT-level parameters
  vector[N_subj] sub_alpha_pr;         // threshold sep. subject mean
  vector[N_subj] sub_beta_pr;          // start point subject mean
  vector[N_subj] sub_delta_female_pr;  // drift rate subject mean, female YES conditions
  vector[N_subj] sub_delta_male_pr;    // drift rate subject mean, male NO conditions
  vector[N_subj] sub_ndt_pr;           // non-decision time subject mean in sec
}

transformed parameters { 
  
  // SUBJECT-level transformed pars for non-centered parameterization
  vector<lower=0.1,upper=4>[N_subj] sub_alpha;        // threshold sep. TRANSFORMED subject mean
  vector<lower=0,upper=1>[N_subj] sub_beta;           // start point TRANSFORMED subject mean
  vector[N_subj] sub_delta_female;                    // drift rate TRANSFORMED subject mean (not yet scaled!)
  vector[N_subj] sub_delta_male;                      // drift rate TRANSFORMED subject mean (not yet scaled!)
  vector<lower=rtBound,upper=max(minRT)*.98>[N_subj] sub_ndt;   // non-decision time in sec TRANSFORMED subject mean
  
  for (i in 1:N_subj) { 
    sub_alpha[i] = 0.1 + 3.9 * Phi(mu_grp_alpha_pr[subj_group[i]] + sig_grp_alpha_pr[subj_group[i]] * sub_alpha_pr[i]);
    sub_beta[i] = Phi(mu_grp_beta_pr[subj_group[i]] + sig_grp_beta_pr[subj_group[i]] * sub_beta_pr[i]);
    sub_ndt[i] = ((minRT[i]*.98 - rtBound) * Phi(mu_grp_ndt_pr[subj_group[i]] + sig_grp_ndt_pr[subj_group[i]] * sub_ndt_pr[i]))+rtBound;
    sub_delta_female[i] = mu_grp_delta_female_pr[subj_group[i]] + sig_grp_delta_pr[subj_group[i]] * sub_delta_female_pr[i];
    sub_delta_male[i] = mu_grp_delta_male_pr[subj_group[i]] + sig_grp_delta_pr[subj_group[i]] * sub_delta_male_pr[i];
  }
}

model {
  
  // GROUP-level hyperpriors
  mu_grp_alpha_pr ~ normal(0, 1);   // prior on threshold sep group mean
  mu_grp_beta_pr ~ normal(0, 1);    // prior on start point group mean
  mu_grp_delta_female_pr ~ normal(0, 1);   // prior on drift rate group mean, female YES conditions
  mu_grp_delta_male_pr ~ normal(0, 1);   // prior on drift rate group mean, male NO conditions
  mu_grp_ndt_pr ~ normal(0, 1);     // prior on NDT group mean
  
  sig_grp_alpha_pr ~ normal(0, .2); // prior on threshold sep group SD 
  sig_grp_beta_pr ~ normal(0, .2);  // prior on start point group SD
  sig_grp_delta_pr ~ normal(0, .2); // prior on drift rate group SD
  sig_grp_ndt_pr ~ normal(0, .2);   // prior on NDT group SD
  
  //SUBJECT-level priors
  sub_alpha_pr ~ normal(0, 1);           // prior on untransformed threshold sep subj mean
  sub_beta_pr  ~ normal(0, 1);           // prior on untransformed start point subj mean
  sub_delta_female_pr  ~ normal(0, 1);   // prior on untransformed drift rate subj mean, female conditions
  sub_delta_male_pr  ~ normal(0, 1);     // prior on untransformed drift rate subj mean, male conditions
  sub_ndt_pr  ~ normal(0, 1);            // prior on untransformed NDT subj mean
  
  // loop through observations
  for (i in 1:N_obs){ 
    real drift; // trial level-drift rate nuisance parameter
    if(gender[i]==1){ //if STIM is  MALE 
      if(choice[i]==1){ // if choice =  MALE (correct)
        drift = Phi(sub_delta_male[subj[i]])*10-5;
        RT[i] ~ wiener(sub_alpha[subj[i]], sub_ndt[subj[i]], sub_beta[subj[i]], drift);
      }else{// if choice = FEMALE (incorrect)
        drift = Phi(sub_delta_male[subj[i]])*10-5;
        RT[i] ~ wiener(sub_alpha[subj[i]], sub_ndt[subj[i]], 1-sub_beta[subj[i]], -drift);
      }
    } else { //if STIM is FEMALE 
      if(choice[i]==1){ //if choice =  MALE (incorrect)
        drift = Phi(sub_delta_female[subj[i]])*10-5;
        RT[i] ~ wiener(sub_alpha[subj[i]], sub_ndt[subj[i]], sub_beta[subj[i]], drift);
      }else{//if choice = FEMALE (correct)
        drift = Phi(sub_delta_female[subj[i]])*10-5;
        RT[i] ~ wiener(sub_alpha[subj[i]], sub_ndt[subj[i]], 1-sub_beta[subj[i]], -drift);
      }
    }
  }
}
generated quantities {
  vector[N_obs] log_lik = rep_vector(0, N_obs); // log liklihood for each observation
  
  // GROUP-level transformed parameters
  vector<lower=0,upper=4>[N_groups] mu_alpha = 0.1 + 3.9*Phi(mu_grp_alpha_pr);                // threshold sep group mean
  vector<lower=0,upper=1>[N_groups] mu_beta = Phi(mu_grp_beta_pr);                            // start point group mean
  vector<lower=-5,upper=5>[N_groups] mu_delta_female = -5 + 10*Phi(mu_grp_delta_female_pr);   // drift rate group mean, female YES conditions
  vector<lower=-5,upper=5>[N_groups] mu_delta_male = -5 + 10*Phi(mu_grp_delta_male_pr);       // drift rate group mean, male NO conditions
  vector<lower=0, upper=0.98>[N_groups] mu_ndt = Phi(mu_grp_ndt_pr);                          // NDT group, proportion
  vector<lower=-10,upper=10>[N_groups] mu_delta_bias;
  
  for (i in 1:N_groups){
    mu_delta_bias[i] = mu_delta_male[i]-(-1*mu_delta_female[i]);  
  }
  
  {
    for (i in 1:N_obs){
    real drift; // trial level-drift rate nuisance parameter
    if(gender[i]==1){ //if STIM is MALE
      if(choice[i]==1){ // if choice = MALE (correct)
        drift = Phi(sub_delta_male[subj[i]])*10-5;
        log_lik[i] += wiener_lpdf(RT[i]|sub_alpha[subj[i]], sub_ndt[subj[i]], sub_beta[subj[i]], drift);
      }else{// if choice = FEMALE (incorrect)
        drift = Phi(sub_delta_male[subj[i]])*10-5;
        log_lik[i] += wiener_lpdf(RT[i]|sub_alpha[subj[i]], sub_ndt[subj[i]], 1-sub_beta[subj[i]], -drift);
      }
    } else { //if STIM is FEMALE
      if(choice[i]==1){ //if choice =  MALE (incorrect)
        drift = Phi(sub_delta_female[subj[i]])*10-5;
        log_lik[i] += wiener_lpdf(RT[i]|sub_alpha[subj[i]], sub_ndt[subj[i]], sub_beta[subj[i]], drift);
      }else{//if choice = FEMALE (correct)
        drift = Phi(sub_delta_female[subj[i]])*10-5;
        log_lik[i] += wiener_lpdf(RT[i]|sub_alpha[subj[i]], sub_ndt[subj[i]], 1-sub_beta[subj[i]], -drift);
      }
    }
  }
  }
}

