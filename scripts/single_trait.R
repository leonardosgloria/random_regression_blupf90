if (!require("pacman")) install.packages("pacman")
pacman::p_load("orthopolynom", "splines","dplyr","tidyr","stringr",
               "data.table","readr","purrr","gtools")

download_BLUPF90(update=T) # Download all blupf90 software to your (folder), the script use the binary from this folder

R_script_folder <- "./R"
files <- list.files(R_script_folder , pattern = "[.][Rr]$", full.names = TRUE, recursive = TRUE)
invisible(lapply(files, function(f) try(source(f, chdir = TRUE), silent = TRUE)))

########
#DATASET
########

pheno_RRM <- data.table::fread("pheno_RRM.txt",h=T,data.table = F)
datarenum1 <- pheno_RRM

datarenum1$Geno <- factor(datarenum1$Geno)
datarenum1$IntBlk <- factor(datarenum1$IntBlk)
datarenum1$Block <- factor(datarenum1$Block)
datarenum1$int <- factor(datarenum1$int)
############
#COV YLD
############
SP <-
  data.frame(Block=datarenum1$Block,Row=datarenum1$Row,Col=datarenum1$Col)
MAP <-
  NAM::NNsrc(SP,2.2,1)

datarenum1$YLD <-
  replace(datarenum1$YLD, which(datarenum1$YLD < 0), NA)

datarenum1$covYLD <-
  NAM::NNcov(MAP,datarenum1$YLD)

########################################################
model <- AdjCC ~ IntBlk + ped|RRM|Geno

residual_start1 <-  matrix(c(0.047))

VC_start1 <-  matrix(c(0.2404,0.5354E-01, -0.2300,
                       0.5354E-01,0.1589E-01,-0.5394E-01,
                       -0.2300,-0.5394E-01,0.2220))

fiels_output_pos1 <- NULL
weights_object1 <- NULL
missing_values1 <- -99
RRM_option1 <- list(poly = 2,
                   Timevar = "Day",
                   Pmin=17,
                   Pmax=73)

genotype_file1 <- "genotype_BLUPF90.txt"
ped_name1 <- NULL
PED_DEPTH1 <- 3
fiels_output_pos1 <- NULL
weights_object1<- NULL


model_single <- blup(datarenum=datarenum1,formula = model,fields_output=fiels_output_pos1,weights_object=weights_object1,
               residual_start=residual_start1,VCA_RRM=VC_start1,ped_name=ped_name1,
               PED_DEPTH=PED_DEPTH1,genotype_file=genotype_file1,
               missing_values=missing_values1,
               RRM_option=RRM_option1,het_res_variance=NULL,
               fit_option=list(yams=T,solution_mean=T,
                               VCE=T,sol_se=T,Inbreeding=T,
                               alpha_size=30,
                               EM_REML=1,
                               maxrounds=3000,
                               alpha_beta=c(0.95,0.05),
                               tunedG=0,
                               conv_crit=1e-10),
              run_model=list(n_threads=20)
              ,keep_files=T)
            #,extra.option.blup="OPTION GBLUP")

fread("h2.txt") %>% filter(Time_var>=RRM_option1$Pmin&Time_var<=RRM_option1$Pmax) # narrow-sense Heritability
fread("VC_Time_var.txt") %>% filter(Time_var>=RRM_option1$Pmin&Time_var<=RRM_option1$Pmax) %>%
  select(Time_var,Additive_variance=Geno)# Additive Genetic variance







