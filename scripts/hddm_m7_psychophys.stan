// PSYCHOPHYS STAN MODEL 7 (JOVI/Visual Integration Models Only)
// This program defines a hierarchical drift diffusion model 
// with a non-centered parameterization. Left responses are 
// modeled as upper boundary responses. Right responses are modeled as
// lower boundary responses. In this version, drift rate (delta) is NOT
// modulated by the jitter level of the stimuli, which has larger magnitude as 
// as signal gets stronger (i.e., jitter decreases) and is negative for "Right"
// stimuli and positive for "Left" stimuli. This model accounts for 
// right censoring of responses/RTs in the JOVI task.
// Note that jitter is zscored 

functions{
  
  real partial_sum_fullddm(real rt, real a, real t0, real w, real v, int resp, int censored){
      
      real sv=0;
      real sw=0;
      real st=0;
      real right_bound=2.5;
      
      real ans=0;

        if(censored==1){//not censored data
          if(resp==1){ //upper bound
            ans+=wiener_lpdf(rt|a,t0,w,v,sv,sw,st);
          }else{//lower bound
            ans+=wiener_lpdf(rt|a,t0,1-w,-v,sv,sw,st);
          }
    }else{ //censored data
      ans+=log_sum_exp(wiener_lccdf_unnorm(right_bound,a,t0,w,v,sv,sw,st),wiener_lccdf_unnorm(right_bound,a,t0,1-w,-v,sv,sw,st));
      }
    
    return ans;
  }
}

data {
  int N_obs;                 // number of observations                         
  int N_subj;                // number of subjects
  int N_jitters;             // number of jitter signal strengths 
  array[N_obs] real jitter;  // jitter level for each trial (zscored)
  int N_choice;              // number of choice alternatives
  int N_groups;              // number of diagnostic groups 
  array[N_obs] real RT;      // RT in seconds for each trial
  array[N_obs] int subj;     // subj id for each trial 
  array[N_obs] int choice;   // response for each trial 
  array[N_obs] int group;    // group id for each trial
  array[N_subj] real minRT;  // minimum RT for each subject
  array[N_obs] int censored; // is observation censored (0) or not (1)
  real rtBound;              // lowest rt allowed 
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
  vector[N_groups] mu_grp_delta_pr;           // drift rate group mean
  vector[N_groups] mu_grp_ndt_pr;             // non-decision time group mean in sec
  
  vector<lower=0>[N_groups] sig_grp_alpha_pr; // threshold sep. group SD 
  vector<lower=0>[N_groups] sig_grp_beta_pr;  // start point group SD
  vector<lower=0>[N_groups] sig_grp_delta_pr; // drift rate group SD
  vector<lower=0>[N_groups] sig_grp_ndt_pr;   // non-decision time group SD
  
  //SUBJECT-level parameters
  vector[N_subj] sub_alpha_pr;  // threshold sep. subject mean
  vector[N_subj] sub_beta_pr;   // start point subject mean
  vector[N_subj] sub_delta_pr;  // drift rate subject mean
  vector[N_subj] sub_ndt_pr;    // non-decision time subject mean in sec
  
  //parameter indexing effect of signal strength/jitter on drift rate
  vector[N_groups] mu_grp_b1_pr;           // group level hyperparamter
  vector<lower=0>[N_groups] sig_grp_b1_pr; // group level between subject variance
  vector[N_subj] sub_b1_pr;                // subject level parameter
}

transformed parameters { 
  
  // SUBJECT-level transformed pars for non-centered parameterization
  vector<lower=0.1,upper=4>[N_subj] sub_alpha;                   // threshold sep. TRANSFORMED subject mean
  vector<lower=0,upper=1>[N_subj] sub_beta;                      // start point TRANSFORMED subject mean
  vector[N_subj] sub_delta;                                      // drift rate TRANSFORMED subject mean (not yet scaled!)
  vector<lower=rtBound,upper=max(minRT)*0.98>[N_subj] sub_ndt;   // non-decision time in sec TRANSFORMED subject mean
  vector[N_subj] sub_b1;                                         // effect of jitter on drift rate
  
  for (i in 1:N_subj) { 
    sub_alpha[i] = 0.1 + 3.9 * Phi(mu_grp_alpha_pr[subj_group[i]] + sig_grp_alpha_pr[subj_group[i]] * sub_alpha_pr[i]);
    sub_beta[i] = Phi(mu_grp_beta_pr[subj_group[i]] + sig_grp_beta_pr[subj_group[i]] * sub_beta_pr[i]);
    sub_ndt[i] = ((minRT[i]*0.98 - rtBound) * Phi(mu_grp_ndt_pr[subj_group[i]] + sig_grp_ndt_pr[subj_group[i]] * sub_ndt_pr[i]))+rtBound;
    sub_b1[i] = mu_grp_b1_pr[subj_group[i]]+sig_grp_b1_pr[subj_group[i]]*sub_b1_pr[i];
    sub_delta[i] = mu_grp_delta_pr[subj_group[i]] + sig_grp_delta_pr[subj_group[i]] * sub_delta_pr[i];
  }
}

model {
  
  // GROUP-level hyperpriors
  mu_grp_alpha_pr ~ normal(0, 1);   // prior on threshold sep group mean
  mu_grp_beta_pr ~ normal(0, 1);    // prior on start point group mean
  mu_grp_delta_pr ~ normal(0, 1);   // prior on drift rate group mean
  mu_grp_ndt_pr ~ normal(0, 1);     // prior on NDT group mean
  mu_grp_b1_pr ~ normal(0,1);       // prior on modulation of drift rate by signal strength (group mean)
  
  sig_grp_alpha_pr ~ normal(0, .2); // prior on threshold sep group SD
  sig_grp_beta_pr ~ normal(0, .2);  // prior on start point group SD
  sig_grp_delta_pr ~ normal(0, .2); // prior on drift rate group SD
  sig_grp_ndt_pr ~ normal(0, .2);   // prior on NDT group SD
  
  //SUBJECT-level priors
  sub_alpha_pr ~ normal(0, 1);      // prior on untransformed threshold sep subj mean
  sub_beta_pr  ~ normal(0, 1);      // prior on untransformed start point subj mean
  sub_delta_pr  ~ normal(0, 1);     // prior on untransformed drift rate subj mean
  sub_ndt_pr  ~ normal(0, 1);       // prior on untransformed NDT subj mean
  sub_b1_pr ~ normal(0,1);          // prior on modulation of drift rate by signal strength (subj mean)
  sig_grp_b1_pr ~ normal(0,.2);     // prior on modulation of drift rate by signal strength (group variance)
  
  // loop through observations
  for (i in 1:N_obs){ 
    real drift; // trial level-drift rate parameter
    drift = Phi(sub_delta[subj[i]]+sub_b1[subj[i]]*jitter[i])*10-5;
    target+=partial_sum_fullddm(RT[i],sub_alpha[subj[i]], sub_ndt[subj[i]],sub_beta[subj[i]],drift,choice[i],censored[i]);
  }
}

generated quantities {
 
  // GROUP-level transformed parameters
  vector<lower=0.1,upper=4>[N_groups] mu_alpha = 0.1 + 3.9*Phi(mu_grp_alpha_pr);// threshold sep group mean
  vector<lower=0,upper=1>[N_groups] mu_beta = Phi(mu_grp_beta_pr);              // start point group mean
  vector<lower=-5,upper=5>[N_groups] mu_delta = -5 + 10*Phi(mu_grp_delta_pr);   // drift rate group mean
  vector<lower=0, upper=0.98>[N_groups] mu_ndt = Phi(mu_grp_ndt_pr);            // NDT group, proportion
  
}
