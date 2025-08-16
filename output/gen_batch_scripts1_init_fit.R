# generate and run batch scripts for initial model fit

input_models <- c('3')
runs <- c(1)
type<-c('init_fit')#par_recovery or #init_fit, ppc, par_recover, or final_fit. 
task<-c('gaze','gender') #gaze or gender
groups<-c('all')
run<-c(1) #1=run batch scripts after generating; 0=1 dont run batch scripts
fitpath<-'[PATH TO FIT DIRECTORY]'

for (i in input_models) {
  for (j in runs) {
    for (k in groups){
      for (l in type){
        for(m in task){
          
          setwd(current_dir<-paste0(fitpath,'/output/hddm_m',i,'_psychophys/',m,'/',l,'/'))
          script_name <- paste0("batch",i,"_",l, "_g",k,"_r", j,".sh")
          script_content <- paste0("#!/bin/bash",
                                   "\n#SBATCH -J INm", i, "r", j, "_",m,l,
                                   "\n#SBATCH -A rXXXXX", # allocation ID on cluster
                                   "\n#SBATCH -o jobname_%j.txt",
                                   "\n#SBATCH -e jobname_%j.err",
                                   "\n#SBATCH --nodes=1",
                                   "\n#SBATCH --ntasks=10",
                                   "\n#SBATCH --cpus-per-task=1",
                                   "\n#SBATCH --time=0-06:00:00",
                                   "\n#SBATCH --mail-user=XXXXXXX", #email address
                                   "\n#SBATCH --mail-type=BEGIN,FAIL,END",
                                   "\n#SBATCH --mem=75G",
                                   "\n#SBATCH --partition=general",
                                   "\n\n#Set up environment",
                                   "\npwd; hostname; date",
                                   "\necho \"Running init fit m", i," group ",k," ",l," for m", i, '_run',j, " on IU Quartz HPC\"",
                                   "\n\nmodule load r/4.3",
                                   "\nexport LD_LIBRARY_PATH=\"/geode2/soft/hps/rhel8/gcc/12.2.0/lib64/:$LD_LIBRARY_PATH\"",
                                   "\n\n#Define variables",
                                   "\nrun=", j, 
                                   "\ngroup=",k,
                                   "\nmodelname='hddm_m", i, "_psychophys'",
                                   "\ntask=",m,
                                   "\n\n#Commands to run",
                                   "\nRscript /[PATH TO RUN FILE]/run1_init_fit.R $run $group $modelname $task")
          
          # Write to file
          writeLines(script_content, script_name)
          
          if(run==1){ 
            # sbatch run job script
            system(paste0("sbatch ",current_dir,script_name),wait=FALSE)
          }
        }
      }
    }
  }
}

