##############################################################################
# create summary digest data based on choices/quantiles per subject per
# signal strength/task condition per response. 
# 
# Note: 
# if a subject has no  trials for a given resp/stimulus combination, then
# n_trials will be set to 0 and all quantiles will be NA.
# For the visual integration task only, subject will have a row designating 
# censored data for a given stimulus strength only if censoring was present
# there are no rows depicting the number of trials for which censoring
# did NOT occur
##############################################################################

library(dplyr)

datadir<-"[PATH TO DATA]"

########################################################
# 1) SUMMARIZE GAZE TASK DATA 
########################################################

outfile<-paste(datadir,"/schizgaze12_gaze_beh.csv",sep="")
data<-read.csv(outfile)

data<-subset(data,data$Study=="schizgaze2")

# sort subj id's in order (hc's [1000s] will come first, then sz [2000s]). (do this so that when we assign sequential index id's, they will be the same for gender and gaze tasks)
data <- data[order(data$Subj), ]

#subset to get the appropriate task data
data<-subset(data,data$Task=="Eyes")
data$SubjID<-match(data$Subj, unique(data$Subj)) #subj IDs to sequential indexes 
data$Resp<-ifelse(data$Resp=="Y",1,2) #recode "yes" as 1 and "no as 2
data$Gender<-ifelse(data$Gender=="F",1,2) #recode female gender of stim as 1 and male as 2
data$Group<-match(data$Group, unique(data$Group))
data <- data[order(data$SubjID, data$Resp), ]
data$SignalStrength<-(data$GazeAngle-mean(data$GazeAngle))/sd(data$GazeAngle) # zscored signal strength 
if (max(data$RT)>100){data$RT<-as.numeric(data$RT)/1000} #if rt col is in ms, convert to seconds

summary<-data %>%
  group_by(SubjID,Group,Resp,SignalStrength) %>%
  summarise(
    n_trials=n(),
    q1 = quantile(RT, 0.1, na.rm = TRUE),
    q3 = quantile(RT, 0.3, na.rm = TRUE),
    q5 = quantile(RT, 0.5, na.rm = TRUE), # median
    q7 = quantile(RT, 0.7, na.rm = TRUE),
    q9 = quantile(RT, 0.9, na.rm = TRUE),
    .groups = "drop")

signals<-unique(summary$SignalStrength)
signals<-signals[order(signals)]

subjects<-unique(summary$SubjID)
subjects<-subjects[order(subjects)]

resps<-unique(summary$Resp)

summary<-as.data.frame(summary)

for(i in 1:length(subjects)){
  tmp_all<-subset(summary,summary$SubjID==subjects[i])
  for(j in 1:length(signals)){
    for(k in 1:length(resps)){
      tmp<-subset(summary,summary$SubjID==subjects[i]&summary$SignalStrength==signals[j]&summary$Resp==resps[k])
      if(nrow(tmp)==0){
        row_id<-nrow(summary)+1
        summary[row_id,]<-NA
        summary$SubjID[row_id]<-subjects[i]
        summary$Group[row_id]<-unique(tmp_all$Group)
        summary$Resp[row_id]<-resps[k]
        summary$SignalStrength[row_id]<-signals[j]
        summary$n_trials[row_id]<-0
      }
    }
  }
}

summary$GroupName<-ifelse(summary$Group==1,"Control","Schizophrenia")
summary$RespName<-ifelse(summary$Resp=="1","LookingAtMe","NotLookingAtMe")
summary <- summary[order(summary$SubjID, summary$SignalStrength), ]
write.csv(summary,"[PATH TO DATA]/gazeTask_summary.csv",row.names = F)

########################################################
# 2) SUMMARIZE GENDER TASK DATA 
########################################################

outfile<-paste(datadir,"/schizgaze12_gaze_beh.csv",sep="")
data<-read.csv(outfile)

data<-subset(data,data$Study=="schizgaze2")

# sort subj id's in order (hc's [1000s] will come first, then sz [2000s]). (do this so that when we assign sequential index id's, they will be the same for gender and gaze tasks)
data <- data[order(data$Subj), ]

#subset to get the appropriate task data
data<-subset(data,data$Task=="GenderID")
data$SubjID<-match(data$Subj, unique(data$Subj)) #subj IDs to sequential indexes
data$Resp<-ifelse(data$Resp=="F",1,2) #recode female as 1 and male as 2
data$Gender<-ifelse(data$Gender=="F",1,2) #recode female gender of stim as 1 and male as 2
data$Group<-match(data$Group, unique(data$Group))
data <- data[order(data$SubjID, data$Resp), ]
if (max(data$RT)>100){data$RT<-as.numeric(data$RT)/1000} #if rt col is in ms, convert to seconds

summary<-data %>%
  group_by(SubjID,Group,Gender,Resp) %>%
  summarise(
    n_trials=n(),
    q1 = quantile(RT, 0.1, na.rm = TRUE),
    q3 = quantile(RT, 0.3, na.rm = TRUE),
    q5 = quantile(RT, 0.5, na.rm = TRUE), # median
    q7 = quantile(RT, 0.7, na.rm = TRUE),
    q9 = quantile(RT, 0.9, na.rm = TRUE),
    .groups = "drop")

genders<-unique(summary$Gender)
genders<-genders[order(genders)]

subjects<-unique(summary$SubjID)
subjects<-subjects[order(subjects)]

resps<-unique(summary$Resp)

summary<-as.data.frame(summary)

for(i in 1:length(subjects)){
  tmp_all<-subset(summary,summary$SubjID==subjects[i])
  for(j in 1:length(genders)){
    for(k in 1:length(resps)){
      tmp<-subset(summary,summary$SubjID==subjects[i]&summary$Gender==genders[j]&summary$Resp==resps[k])
      if(nrow(tmp)==0){
        row_id<-nrow(summary)+1
        summary[row_id,]<-NA
        summary$SubjID[row_id]<-subjects[i]
        summary$Group[row_id]<-unique(tmp_all$Group)
        summary$Resp[row_id]<-resps[k]
        summary$Gender[row_id]<-genders[j]
        summary$n_trials[row_id]<-0
      }
    }
  }
}

summary$GroupName<-ifelse(summary$Group==1,"Control","Schizophrenia")
summary$RespName<-ifelse(summary$Resp=="1","Female","Male")
summary$GenderName<-ifelse(summary$Gender=="1","Female","Male")
summary$Accuracy<-ifelse(summary$Resp==summary$Gender,1,0)

summary <- summary[order(summary$SubjID, summary$Gender), ]
write.csv(summary,"[PATH TO DATA]/genderTask_summary.csv",row.names = F)

########################################################
# 3) SUMMARIZE VISUAL INTEGRATION (JOVI) TASK DATA 
########################################################

outfile<-paste(datadir,"/schizgaze2_jovi_beh_all.csv",sep="")
data<-read.csv(outfile)

data<-subset(data,data$Study=="schizgaze2")

# sort subj id's in order (hc's [1000s] will come first, then sz [2000s]). (do this so that when we assign sequential index id's, they will be the same for gender and gaze tasks)
data <- data[order(data$Subj), ]

#flip responses for subject that mixed up response buttons
for(i in 1:nrow(data)){
  if(data$Subj[i]=="2057"){
    data$Resp[i]<-ifelse(data$Resp[i]==2,3,2)
    data$RespBound[i]<-ifelse(data$RespBound[i]==1,2,1)
    data$Acc[i]<-ifelse(data$Direction[i]==data$RespBound[i],1,0)
  }
}

data$SubjID<-match(data$Subj, unique(data$Subj)) #subj IDs to sequential indexes for fitting
data$Resp<-ifelse(data$Resp=="2",1,2) #recode L response (2) as 1 and R response (3) as 2
data$Group<-match(data$Group, unique(data$Group))
data <- data[order(data$SubjID, data$Resp), ]
if (max(data$RT)>100){data$RT<-as.numeric(data$RT)/1000} #if rt col is in ms, convert to seconds
data$Resp<-ifelse(data$censored==0,NA,data$Resp)
data$RespBound<-ifelse(data$censored==0,NA,data$RespBound)
data$Acc<-ifelse(data$censored==0,NA,data$Acc)
data$RT<-ifelse(data$censored==0,NA,data$RT)

summary<-data %>%
  group_by(SubjID,Group,Jitter,censored,Resp) %>%
  summarise(
    n_trials=n(),
    q1 = quantile(RT, 0.1, na.rm = TRUE),
    q3 = quantile(RT, 0.3, na.rm = TRUE),
    q5 = quantile(RT, 0.5, na.rm = TRUE), # median
    q7 = quantile(RT, 0.7, na.rm = TRUE),
    q9 = quantile(RT, 0.9, na.rm = TRUE),
    .groups = "drop")

jitters<-unique(summary$Jitter)
jitters<-jitters[order(jitters)]

subjects<-unique(summary$SubjID)
subjects<-subjects[order(subjects)]

resps<-unique(summary$Resp)[-1]
censors<-unique(summary$censored) #0=censored, 1=not censored

summary<-as.data.frame(summary)

for(i in 1:length(subjects)){
  tmp_all<-subset(summary,summary$SubjID==subjects[i])
  for(j in 1:length(jitters)){
    for(k in 1:length(resps)){
        tmp<-subset(summary,summary$SubjID==subjects[i]&summary$Jitter==jitters[j]&summary$Resp==resps[k])
        if(nrow(tmp)==0){
          row_id<-nrow(summary)+1
          summary[row_id,]<-NA
          summary$SubjID[row_id]<-subjects[i]
          summary$Group[row_id]<-unique(tmp_all$Group)
          summary$Resp[row_id]<-resps[k]
          summary$Jitter[row_id]<-jitters[j]
          summary$n_trials[row_id]<-0
        }
    }
  }
}

summary$GroupName<-ifelse(summary$Group==1,"Control","Schizophrenia")
summary$RespName<-ifelse(summary$Resp=="1","Right","Left")
summary$DirectionName<-ifelse(summary$Jitter<0,"Left","Right")
summary$CensorName<-ifelse(summary$censored==0,"Censored","NotCensored")
summary$Accuracy<-ifelse(summary$RespName==summary$DirectionName,1,0)

summary <- summary[order(summary$SubjID, summary$censored,summary$Jitter),]

write.csv(summary,"[PATH TO DATA]/joviTask_summary.csv",row.names = F)
