pipeline_cross_validation <- function(datarenum,formula,group_name_cv,LOO_cv,
                                      prop,nrep,eff_name,RRM_option,keep_folder=F,
                                      parameter_file="renf90.par",logname="blupf90.log",
                                      genotype_file,n_threads,method="blupf90+",blupf90_folder=blupf90_folder){
  orig <- getwd()
  fold_LR <- c()

#################
### VC names####
#################

  ter = terms(formula)
  labs = labels(ter)

  variab_form <-
    labs %>% gsub("\\+ ", " ",.) %>% strsplit(., split="  ", fixed=F) %>%
    unlist()

  fixed_terms <-
    variab_form %>% str_detect(.,"[/ | ]", negate = T) %>% which(T) %>%
    variab_form[.]

  ped_terms1 <-
    variab_form %>% str_detect(.,"ped ") %>% which(TRUE) %>%
    variab_form[.] %>%
    gsub( "ped ", "",.) %>% gsub( " | ", "",.)

  ped_terms2 <-
    ped_terms1  %>% str_detect(.,"RRM", negate = T) %>% which(TRUE) %>%
    ped_terms1[.]

  ped_terms <-
    ped_terms2  %>% str_detect(.,"/", negate = T) %>% which(TRUE) %>%
    ped_terms2[.]%>% gsub( "\\|", "",.)

  diag_terms1 <-
    variab_form %>% str_detect(.,"1 ") %>% which(TRUE) %>%
    variab_form[.] %>%
    gsub( "1 \\| ", "",.)

  diag_terms2 <- diag_terms1%>%
    str_detect(.,"RRM | ", negate = T) %>% which(TRUE) %>% diag_terms1[.]

  diag_terms <- diag_terms2%>%
    str_detect(.,"/", negate = T) %>% which(TRUE) %>% diag_terms2[.]

  nested_term1 <-
    variab_form %>% str_detect(.,"/", negate = F) %>% which(T)%>%
    variab_form[.]

  nested_term2 <-
    nested_term1 %>% str_detect(.," | ", negate = T) %>% which(T)%>%
    nested_term1[.]

  nestedF_term <- nested_term2 %>%  str_detect(.,"/RRM", negate = T) %>% which(T)%>%
    nested_term2[.] %>% strsplit(., split="/", fixed=F)

  nested_diag_term <-
    nested_term1 %>% str_detect(.,"1 \\| ", negate = F) %>% which(T)%>%
    nested_term1[.] %>% gsub( "1 \\| ", "",.) %>%
    strsplit(., split="/", fixed=F)

  nested_ped_term <-
    nested_term1 %>% str_detect(.,"ped \\| ", negate = F) %>% which(T)%>%
    nested_term1[.] %>% gsub( "ped \\| ", "",.) %>%
    strsplit(., split="/", fixed=F)

  RRM_diag_terms1 <-
    variab_form %>% str_detect(.,"RRM") %>% which(TRUE) %>%
    variab_form[.] %>% gsub( "1 \\| ", "",.)

  RRM_diag_terms2 <-
    RRM_diag_terms1 %>%str_detect(.,"ped", negate = T) %>% which(TRUE) %>%
    RRM_diag_terms1[.]

  RRM_diag_terms <- RRM_diag_terms2 %>% str_detect(.,"RRM ", negate = F) %>% which(TRUE) %>%
    RRM_diag_terms2[.] %>% gsub( "RRM \\| ", "",.)

  RRM_nestedF_terms <-
    RRM_diag_terms2 %>% str_detect(.,"RRM ", negate = T) %>% which(TRUE) %>%
    RRM_diag_terms2[.] %>% gsub( "/RRM", "",.)


  RRM_ped_terms <-
    RRM_diag_terms1 %>%str_detect(.,"ped", negate = F) %>% which(TRUE) %>%
    RRM_diag_terms1[.] %>% gsub( "RRM \\| ", "",.) %>%
    gsub( "ped ", "",.) %>% gsub( "\\| ", "",.)

  trait_terms <-
    ter %>% .[[2]] %>% as.character() %>%
    gsub( "\\|", NA,.) %>% na.exclude() %>%
    c()


  VC_names <- c(diag_terms,ped_terms,RRM_ped_terms,RRM_diag_terms)
  ### Find position of Genetic variance components inside postmean

  pos_VC <- read_lines(logname) %>%
    str_detect(.,"effect") %>% which(TRUE) %>% +1


####
 if(is.null(LOO_cv)==T){

   for (nreply in 1:nrep) {
    setwd(orig)
     dir.create(paste0("fold",nreply))
####FULL EBV#########################
  pos_random <- readLines(parameter_file) %>%grep(" RANDOM_GROUP", .)+1

  VC_ <- list()
  for (r in 1:length(pos_random)) {
    VC_[[r]] <-unlist(strsplit(readLines(parameter_file)[pos_random[r]], split=" ", fixed=TRUE)) %>%
      as.numeric() %>% na.exclude() %>% c()
  }

  vc_matrix <- list()
  for (vc in 1:length(VC_)) {
    #name_vc <- paste0(VC_names[vc],".temp")
    vc_matrix [[vc]]<- read_lines(logname)[seq(pos_VC[vc],pos_VC[vc]+length(VC_[[vc]])-1)] %>%
      strsplit(.,split = " ", fixed=TRUE) %>% unlist() %>%
      as.numeric() %>% na.exclude() %>%
      matrix(ncol = length(VC_[[vc]]),nrow = length(VC_[[vc]]))
    #write.table(vc_matrix,file = name_vc,row.names = F,col.names = F)
  }
  names(vc_matrix) <- VC_names

  genetic_eff_pos <- VC_[[which(names(vc_matrix)==eff_name)]]

  if(method=="blupf90+"){solut="solutions.orig"}else{
    solut="final_solutions.orig"
  }

ncol_sol <- ncol(data.table::fread(solut,skip = 1,nrows = 2))
  if(ncol_sol<6){
  col_sol <- c("trait","effect_","level","id","solution")
  }else{
  col_sol <- c("trait","effect_","level","id","solution","s.e.")
  }

full_ebv <- data.table::fread(solut,skip = 1,col.names=col_sol) %>%
    filter(effect_%in%genetic_eff_pos)

########################################################
#copy all required files to the k_fold folder
########################################################
########################################################
#copy all required files to the k_fold folder
########################################################
file.copy(list.files(orig, "*.ped"), paste0(orig,"/fold",nreply))

file.copy(list.files(orig, parameter_file), paste0(orig,"/fold",nreply))

file.copy(list.files(orig, "renf90.inb"), paste0(orig,"/fold",nreply))

file.copy(list.files(orig, "renf90.fields"), paste0(orig,"/fold",nreply))
file.copy(list.files(orig, "renf90.tables"), paste0(orig,"/fold",nreply))
#Copy Genotype data if exist

  if(is.null(genotype_file)==F){
     file.copy(list.files(orig, genotype_file), paste0(orig,"/fold",nreply))
     }
##########################################
"%ni%" <- Negate("%in%") #opposite of %in%

##########################################################
   data_fold <-
     data.table::fread("renf90.dat",
           data.table = F,
           col.names = c(read.table("cols_data_renumf90.txt",h=T)$names))

#filter_fold <- list(list())

   for (grpby in group_name_cv) {

     fold_group <- levels(factor(data_fold[,grpby]))[sample(c(1:length(levels(factor(data_fold[,grpby])))),
                                                                    ceiling(length(levels(factor(data_fold[,grpby])))*prop))]

#     filter_fold[[nreply]][grpby] <- list(levels(factor(data_fold[,grpby]))[sample(c(1:length(levels(factor(data_fold[,grpby])))),
#                                                                    ceiling(length(levels(factor(data_fold[,grpby])))*prop))])
     data_fold <- data_fold[data_fold[,grpby]%ni%fold_group,]
   }
filter_fold <- data_fold[,eff_name]

cat(paste("############################# ",nreply, "Fold folder created ########################################"," \n"))

setwd(paste0(orig,"/fold",nreply))

data.table::fwrite(data_fold,"renf90.dat",sep = " ", col.names = F)
#END create all kfold folders
##################################
##Calculating EBV by folder#######
##################################
cat(paste("############################# Calculating EBV for the",nreply,"Fold ########################################"," \n"))

    pipeline_fit_model(method="blupf90+",
                       Gibbs_option= NULL,
                       slurm_option=NULL,
                       parameter_file=parameter_file,dense=F,n_threads,blupf90_folder)


ncol_sol <- ncol(data.table::fread(solut,skip = 1,nrows = 2))
    if(ncol_sol<6){
    col_sol <- c("trait","effect_","level","id","solution")
    }else{
    col_sol <- c("trait","effect_","level","id","solution","s.e.")
    }
          #IF model is RR
          if(is.null(RRM_option)==F){
##################################
######PIPELINE CROSSVALIDATION RRM
##################################

    Timemin=min(datarenum[,RRM_option$Timevar])
    Timemax=max(datarenum[,RRM_option$Timevar])
    Time_var <-seq(Timemin,Timemax)

          fi1 <- as.matrix(read.table("../fi.txt",h=F))

          fold_sol <-
            data.table::fread(solut,skip = 1,col.names=col_sol) %>%
              select(effect_,id,solution) %>%
              filter(effect_%in%genetic_eff_pos) %>%
              filter(id%ni%filter_fold) %>%
              pivot_wider(names_from = effect_, values_from = solution,values_fn = list) %>%
              unnest(cols = paste(genetic_eff_pos)) %>% data.frame()

          fold_sol_id <- fold_sol$id

          fold_sol<- fold_sol %>% .[,-1]%>%
              as.matrix(.)%*%t(fi1) %>% data.frame(fold_sol_id,.)

          full_ebv <- full_ebv %>%
            select(effect_,id,solution) %>%
              filter(effect_%in%genetic_eff_pos) %>%
              filter(id%ni%filter_fold) %>%
              pivot_wider(names_from = effect_, values_from = solution,values_fn = list) %>%
              unnest(cols = paste(genetic_eff_pos)) %>% data.frame()

          full_ebv_id <- full_ebv$id

          full_ebv<- full_ebv %>% .[,-1]%>%
              as.matrix(.)%*%t(fi1) %>% data.frame(full_ebv_id,.)

          accu_k <-
            mapply(cor,full_ebv[,-1],fold_sol[,-1])

            x <-
              full_ebv[,-1] %>% colMeans() %>%
              cbind(colMeans(fold_sol[,-1]))
          bias_f <- function(x){
                x[1]-x[2]
              }
          bias_k <- apply(x, 1,bias_f)

          dispersion_k <-
            mapply(cov,full_ebv[,-1],fold_sol[,-1])/apply(fold_sol[,-1], 2, var)

          fold_LR <- data.frame(accu_k,bias_k,dispersion_k) %>% mutate(Time_var)%>%
                rbind(fold_LR)

          #ELSE non-RRM
          }else{
          fold_sol <-
            data.table::fread(solut,skip = 1,col.names=col_sol) %>%
              select(effect_,id,solution) %>%
              filter(effect_%in%genetic_eff_pos) %>%
              filter(id%ni%filter_fold) %>% select(solution)

          full_ebv <- full_ebv %>%
            select(effect_,id,solution) %>%
              filter(effect_%in%genetic_eff_pos) %>%
              filter(id%ni%filter_fold) %>% select(solution)

          fold_LR <- data.frame(full_ebv,fold_sol) %>%
            summarise(accu_k = cor(solution,solution.1),
                      bias_k = mean(solution) - mean(solution.1),
                      dispersion_k = cov(solution,solution.1)/var(solution.1)) %>%
             rbind(fold_LR)
          write.table(fold_LR,"fold_LR.txt",row.names = F)
                }
  }#end loop rep
  setwd(orig)
  if(is.null(RRM_option)==F){
  fold_LR %>% group_by(Time_var) %>% summarise(accur_m=mean(accu_k,na.rm = T),accur_sd=sd(accu_k,na.rm = T),
                                                bias_m=mean(bias_k,na.rm = T),bias_sd=sd(bias_k,na.rm = T),
                                                dispersion_m=mean(dispersion_k,na.rm = T),dispersion_sd=sd(dispersion_k,na.rm = T)) %>%
  write.table("cross_validation_result.txt",row.names = F,quote = F)
  }else{
  fold_LR %>% summarise(accur_m=mean(accu_k,na.rm = T),accur_sd=sd(accu_k,na.rm = T),
                        bias_m=mean(bias_k,na.rm = T),bias_sd=sd(bias_k,na.rm = T),
                        dispersion_m=mean(dispersion_k,na.rm = T),dispersion_sd=sd(dispersion_k,na.rm = T)) %>%
  write.table("cross_validation_result.txt",row.names = F,quote = F)
  }
#ELSE providing the effect
}else{
  data_fold <-
     data.table::fread("renf90.dat",
           data.table = F,
           col.names = c(read.table("cols_data_renumf90.txt",h=T)$names))

  LOO_cv_lvs <- c(levels(factor(data_fold[,LOO_cv])))

  for (nreply in LOO_cv_lvs) {
   setwd(orig)

    data_fold <-
      data.table::fread("renf90.dat",
            data.table = F,
            col.names = c(read.table("cols_data_renumf90.txt",h=T)$names))

    #levels(factor(data_fold[,LOO_cv]))[nreply]

    data_fold <- data_fold[data_fold[,LOO_cv]%ni%nreply,]

    dir.create(paste0("fold",nreply))
####FULL EBV#########################
  pos_random <- readLines(parameter_file) %>%grep(" RANDOM_GROUP", .)+1

  VC_ <- list()
  for (r in 1:length(pos_random)) {
    VC_[[r]] <-unlist(strsplit(readLines(parameter_file)[pos_random[r]], split=" ", fixed=TRUE)) %>%
      as.numeric() %>% na.exclude() %>% c()
  }

  vc_matrix <- list()
  for (vc in 1:length(VC_)) {
    #name_vc <- paste0(VC_names[vc],".temp")
    vc_matrix [[vc]]<- read_lines(logname)[seq(pos_VC[vc],pos_VC[vc]+length(VC_[[vc]])-1)] %>%
      strsplit(.,split = " ", fixed=TRUE) %>% unlist() %>%
      as.numeric() %>% na.exclude() %>%
      matrix(ncol = length(VC_[[vc]]),nrow = length(VC_[[vc]]))
    #write.table(vc_matrix,file = name_vc,row.names = F,col.names = F)
  }
  names(vc_matrix) <- VC_names

  genetic_eff_pos <- VC_[[which(names(vc_matrix)==eff_name)]]

ncol_sol <- ncol(data.table::fread(solut,skip = 1,nrows = 2))
  if(ncol_sol<5){
  col_sol <- c("trait","effect_","id","solution")
  }else{
  col_sol <- c("trait","effect_","id","solution","s.e.")
  }

full_ebv <- data.table::fread(solut,skip = 1,col.names=col_sol) %>%
    filter(effect_%in%genetic_eff_pos)

########################################################
#copy all required files to the k_fold folder
########################################################
  file.copy(list.files(orig, "*.ped"), paste0(orig,"/fold",nreply))

  file.copy(list.files(orig, parameter_file), paste0(orig,"/fold",nreply))

  file.copy(list.files(orig, "renf90.inb"), paste0(orig,"/fold",nreply))

  #Copy Genotype data if exist

  if(is.null(genotype_file)==F){
     file.copy(list.files(orig, genotype_file), paste0(orig,"/fold",nreply))
  }
cat(paste("############################# ",nreply, "Fold folder created ########################################"," \n"))
  setwd(paste0(orig,"/fold",nreply))

  data.table::fwrite(data_fold,"renf90.dat",sep = " ", col.names = F)


##################################
##Calculating EBV by folder#######
##################################
cat(paste("############################# Calculating EBV for the",nreply,"Fold ########################################"," \n"))

    pipeline_fit_model(method=method,
                       Gibbs_option= Gibbs_option,
                       slurm_option=slurm_options,
                       parameter_file=parameter_file,dense=dense,n_threads = n_threads)

ncol_sol <- ncol(data.table::fread(solut,skip = 1,nrows = 2))
  if(ncol_sol<5){
  col_sol <- c("trait","effect_","id","solution")
  }else{
  col_sol <- c("trait","effect_","id","solution","s.e.")
  }
          #IF model is RR
          if(is.null(RRM_option)==F){
          fi1 <- as.matrix(read.table("../fi.txt",h=F))

          fold_sol <-
            data.table::fread(solut,skip = 1,col.names=col_sol) %>%
              select(effect_,id,solution) %>%
              filter(effect_%in%genetic_eff_pos) %>%
              filter(id%ni%filter_fold) %>%
              pivot_wider(names_from = effect_, values_from = solution,values_fn = list) %>%
              unnest(cols = paste(genetic_eff_pos)) %>% data.frame()

          fold_sol_id <- fold_sol$id

          fold_sol<- fold_sol %>% .[,-1]%>%
              as.matrix(.)%*%t(fi1) %>% data.frame(fold_sol_id,.)

          full_ebv <- full_ebv %>%
            select(effect_,id,solution) %>%
              filter(effect_%in%genetic_eff_pos) %>%
              filter(id%ni%filter_fold) %>%
              pivot_wider(names_from = effect_, values_from = solution,values_fn = list) %>%
              unnest(cols = paste(genetic_eff_pos)) %>% data.frame()

          full_ebv_id <- full_ebv$id

          full_ebv<- full_ebv %>% .[,-1]%>%
              as.matrix(.)%*%t(fi1) %>% data.frame(full_ebv_id,.)

          accu_k <-
            mapply(cor,full_ebv,fold_sol) %>% .[-1]

            x <-
              full_ebv %>% colMeans() %>%
              cbind(colMeans(fold_sol))
          bias_f <- function(x){
                x[1]-x[2]
              }
          bias_k <- apply(x, 1,bias_f) %>% .[-1]

          dispersion_k <-
            mapply(cov,full_ebv,fold_sol)[-1]/apply(fold_sol, 2, var)[-1]

          fold_LR <- data.frame(accu_k,bias_k,dispersion_k) %>% mutate(Time_var)%>%
                rbind(fold_LR)

          write.table(fold_LR,"fold_LR.txt",row.names = F)
          #ELSE non-RRM
          }else{
          fold_sol <-
            data.table::fread(solut,skip = 1,col.names=col_sol) %>%
              select(effect_,id,solution) %>%
              filter(effect_%in%genetic_eff_pos) %>%
              filter(id%ni%filter_fold) %>% select(solution)

          full_ebv <- full_ebv %>%
            select(effect_,id,solution) %>%
              filter(effect_%in%genetic_eff_pos) %>%
              filter(id%ni%filter_fold) %>% select(solution)

          fold_LR <- data.frame(full_ebv,fold_sol) %>%
            summarise(accu_k = cor(solution,solution.1),
                      bias_k = mean(solution) - mean(solution.1),
                      dispersion_k = cov(solution,solution.1)/var(solution.1)) %>%
             rbind(fold_LR)
          write.table(fold_LR,"fold_LR.txt",row.names = F)

                }
  }#end loop rep
  setwd(orig)
  if(is.null(RRM_option)==F){
  fold_LR %>% group_by(Time_var) %>% summarise(accur_m=mean(accu_k,na.rm = T),accur_sd=sd(accu_k,na.rm = T),
                                                bias_m=mean(bias_k,na.rm = T),bias_sd=sd(bias_k,na.rm = T),
                                                dispersion_m=mean(dispersion_k,na.rm = T),dispersion_sd=sd(dispersion_k,na.rm = T)) %>%
  write.table("cross_validation_result.txt",row.names = F,quote = F)
  }else{
  fold_LR %>% summarise(accur_m=mean(accu_k,na.rm = T),accur_sd=sd(accu_k,na.rm = T),
                        bias_m=mean(bias_k,na.rm = T),bias_sd=sd(bias_k,na.rm = T),
                        dispersion_m=mean(dispersion_k,na.rm = T),dispersion_sd=sd(dispersion_k,na.rm = T)) %>%
  write.table("cross_validation_result.txt",row.names = F,quote = F)
  }
 }
#Remove fold folders
    if(keep_folder==F){
      fold_name_folder <- list.dirs() %>% gsub("[.-/]+","",.) %>% .[!.==""]
      unlink(fold_name_folder, recursive = TRUE)
    }
}
