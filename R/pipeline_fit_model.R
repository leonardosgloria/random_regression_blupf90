pipeline_fit_model <- function(method="blupf90+",Gibbs_option=NULL,slurm_option=NULL,
                               parameter_file="renf90.par",dense,n_threads=2,blupf90_folder=blupf90_folder){
if(dense==T){dense_blupf90=" --dense"}else{dense_blupf90=""}
  if(gsub("([0-9]|\\.)","",version$os)=="linux-gnu"){
    S_OP <- "Linux"
  } else if(gsub("([0-9]|\\.)","",version$os)=="darwin"){
    S_OP <- "Mac_OSX"}else{
      S_OP <- "Windows"
    }

    if(method=="blupf90+"){


      if(S_OP=="Windows"){
        aireml <- paste0(blupf90_folder,"/blupf90+.exe ", parameter_file, dense_blupf90)
      }else{
        aireml <- paste0("ulimit -s unlimited && export OMP_NUM_THREADS=",n_threads,"&&",blupf90_folder,"/blupf90+ ", parameter_file, dense_blupf90)
      }
      system(aireml)
    }
    if(method=="gibbsf90+"){
      if(S_OP=="Windows"){
      gibbs <- paste0(blupf90_folder,
                      paste("/gibbsf90+.exe ", parameter_file,
                            "--samples",Gibbs_option$samples,
                            "--burnin",Gibbs_option$burnin,
                            "--thin",Gibbs_option$thin)
                      )
      }else{
        gibbs <- paste0(blupf90_folder,
                        paste("/gibbsf90+ ", parameter_file,
                              "--samples",Gibbs_option$samples,
                              "--burnin",Gibbs_option$burnin,
                              "--thin",Gibbs_option$thin)
        )
      }
      system(gibbs)
    }

    if(method=="SLURM"){
    ###SLURM FILE
    sink(file = "slurm.txt")
    cat("#!/bin/sh -l","\n")
    cat("# FILENAME: slurm.txt","\n")
    cat("#SBATCH ",slurm_options$node_type," ",slurm_options$account,"\n",sep = "")
    cat("#SBATCH ","--nodes=",slurm_options$nodes,"\n",sep = "")
    cat("#SBATCH ","--ntasks=",slurm_options$ntasks,"\n",sep = "")
    cat("#SBATCH ","--mem=",slurm_options$memory,"\n",sep = "")
    cat("#SBATCH ","--time=",slurm_options$time_analysis,"\n",sep = "")
    cat("#SBATCH ","--job-name",basename(getwd()),"\n",sep = " ")
    cat("#SBATCH -o slurm_output.out","\n",sep = "")
    cat("#SBATCH -e slurm_error.out","\n",sep = "")
    #cat(# cd $SLURM_SUBMIT_DIR)
    cat("export OMP_NUM_THREADS=",slurm_options$n_threads_BLUPF90,"\n",sep = "")
    cat("./slurm.sh")
    sink()

    #####.sh file

    if(is.null(Gibbs_option)==T){
      gibbs_opt <- NULL
    }else{
      gibbs_opt <-paste("--rounds",Gibbs_option$samples,
                        "--burnin",Gibbs_option$burnin,
                        "--thin",Gibbs_option$thin)}

    sink(file = "slurm.sh")
    cat("#!/bin/sh -l","\n")
    cat("ulimit -s unlimited","\n")
    cat("echo ",parameter_file," |","\n")
    cat(paste0(blupf90_folder,"/",slurm_options$programf90," "))
    cat(gibbs_opt)
    sink()

    system("chmod 777 slurm.sh")
    system("sbatch slurm.txt")
  }
}
