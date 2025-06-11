# generate and run batch scripts for supplemental analysis running gender models separately by participant gender

input_models <- c('3')
runs <- c(1)
type<-c('supplement_bygend')#par_recovery or #init_fit, ppc, par_recover, or final_fit
task<-c('gender')
groups<-c('all')
gend_groups<-c('female_only','male_only')
run<-c(1) #1=run batch scripts after generating; 0=1 dont run batch scripts
fitpath<-'/N/slate/clasagn/gazeddm_psychophys'

for (i in input_models) {
  for (j in runs) {
    for (k in groups){
      for (l in type){
        for(m in task){
          for(n in gend_groups){
          
            setwd(current_dir<-paste0(fitpath,'/output/hddm_m',i,'_psychophys/',m,'/',l,'/',n,'/'))
            script_name <- paste0("batch",i,"_",l, "_g",k,"_r", j,n,".sh")
            script_content <- paste0("#!/bin/bash",
                                     "\n#SBATCH -J SUm", i, "r", j, "_",m,l,n,
                                     "\n#SBATCH -A rXXXXX", #cluster allocation ID
                                     "\n#SBATCH -o jobname_%j.txt",
                                     "\n#SBATCH -e jobname_%j.err",
                                     "\n#SBATCH --nodes=1",
                                     "\n#SBATCH --ntasks=10",
                                     "\n#SBATCH --cpus-per-task=1",
                                     "\n#SBATCH --time=0-4:00:00",
                                     "\n#SBATCH --mail-user=XXXXXXX", #email
                                     "\n#SBATCH --mail-type=BEGIN,FAIL,END",
                                     "\n#SBATCH --mem=100G",
                                     "\n#SBATCH --partition=general",
                                     "\n\n#Set up environment",
                                     "\npwd; hostname; date",
                                     "\necho \"Running init fit m", i," group ",k," ",l," for m", i, '_run',j,n, " on IU Quartz HPC\"",
                                     "\n\nmodule load r/4.3",
                                     "\nexport LD_LIBRARY_PATH=\"/geode2/soft/hps/rhel8/gcc/12.2.0/lib64/:$LD_LIBRARY_PATH\"",
                                     "\n\n#Define variables",
                                     "\nrun=", j, 
                                     "\ngroup=",k,
                                     "\nmodelname='hddm_m", i, "_psychophys'",
                                     "\ntask=",m,
                                     "\ngend_group=",n,
                                     "\n\n#Commands to run",
                                     "\nRscript /N/slate/clasagn/gazeddm_psychophys/output/run4_supplement_bygend.R $run $group $modelname $task $gend_group")
            
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
}

