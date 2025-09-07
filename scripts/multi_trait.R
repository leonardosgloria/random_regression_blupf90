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
model <- list(AdjCC ~   int/RRM + IntBlk + covYLD + ped|RRM|Geno,
              YLD ~   IntBlk + covYLD + ped|Geno
              )

residual_start1 <-matrix(c(0.0400,       0.0000
                   ,0.0000,       90.000  ))

VC_start1 <- matrix(c(0.2528, 4.025,0.4715E-01,0.000,-0.2369,0.000,
                            4.025,64.42,0.7586,0.000,-3.772,0.000,
                            0.4715E-01,0.7586,0.1352E-01,0.000,-0.4739E-01,0.000,
                            0.000,0.000,0.000,0.000,0.000,0.000,
                            -0.2369,-3.772,-0.4739E-01,0.000,0.2246,0.000,
                            0.000,0.000,0.000,0.000,0.000,0.000 ))

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


model_multi <- blup(datarenum=datarenum1,formula = model,fields_output=fiels_output_pos1,weights_object=weights_object1,
               residual_start=residual_start1,VCA_RRM=VC_start1,ped_name=ped_name1,
               PED_DEPTH=PED_DEPTH1,genotype_file=genotype_file1,
               missing_values=missing_values1,
               RRM_option=RRM_option1,het_res_variance=NULL,
               fit_option=list(yams=T,solution_mean=T,
                               VCE=T,sol_se=T,Inbreeding=T,
                               alpha_size=30,
                               EM_REML=1000,
                               maxrounds=3000,
                               alpha_beta=c(0.95,0.05),
                               tunedG=0,
                               conv_crit=1e-10),
              run_model=list(n_threads=20)
              ,keep_files=T)
##############
# OUTPUT
##############            
model_multi$h2_multi$AdjCC # Variance components for Adjusted canopy coverage
model_multi$h2_multi$YLD %>% distinct(RESIDUAL,.keep_all = T) # Variance components for yield
model_multi$gen_cor # Genetic correlation between the traits





