
###############################FUNCTION#######################################################
pipeline_renum_single <- function(datarenum,formula,fields_output,weights_object,residual_start,VCA_RRM,VCD_RRM,VCA,VCD,
                            ped_name,PED_DEPTH,genotype_file,missing_values,RRM_option,het_res_variance=NULL,
                            fit_option,extra.option.blup,blupf90_folder){

##################################

ter = terms(formula)
labs = labels(ter)

variab_form <-
  labs %>% gsub("\\+ ", " ",.) %>% strsplit(., split="  ", fixed=F) %>%
  unlist()

fixed_terms <-
  variab_form %>% stringr::str_detect(.,"[/ | ]", negate = T) %>% which(T) %>%
  variab_form[.]

ped_terms1 <-
  variab_form %>% stringr::str_detect(.,"ped ") %>% which(TRUE) %>%
    variab_form[.] %>%
     gsub( "ped ", "",.) %>% gsub( " | ", "",.)

ped_terms2 <-
  ped_terms1  %>% stringr::str_detect(.,"RRM", negate = T) %>% which(TRUE) %>%
  ped_terms1[.]

ped_terms <-
  ped_terms2  %>% stringr::str_detect(.,"/", negate = T) %>% which(TRUE) %>%
  ped_terms2[.]%>% gsub( "\\|", "",.)

diag_terms1 <-
  variab_form %>% stringr::str_detect(.,"1 ") %>% which(TRUE) %>%
  variab_form[.] %>%
  gsub( "1 \\| ", "",.)

diag_terms2 <- diag_terms1%>%
  stringr::str_detect(.,"RRM | ", negate = T) %>% which(TRUE) %>% diag_terms1[.]

diag_terms <- diag_terms2%>%
  stringr::str_detect(.,"/", negate = T) %>% which(TRUE) %>% diag_terms2[.]

nested_term1 <-
variab_form %>% stringr::str_detect(.,"/", negate = F) %>% which(T)%>%
  variab_form[.]

nested_term2 <-
  nested_term1 %>% stringr::str_detect(.," | ", negate = T) %>% which(T)%>%
    nested_term1[.]

nestedF_term <- nested_term2 %>%  stringr::str_detect(.,"/RRM", negate = T) %>% which(T)%>%
  nested_term2[.] %>% strsplit(., split="/", fixed=F)

nested_diag_term <-
  nested_term1 %>% stringr::str_detect(.,"1 \\| ", negate = F) %>% which(T)%>%
  nested_term1[.] %>% gsub( "1 \\| ", "",.) %>%
  strsplit(., split="/", fixed=F)

nested_ped_term <-
  nested_term1 %>% stringr::str_detect(.,"ped \\| ", negate = F) %>% which(T)%>%
  nested_term1[.] %>% gsub( "ped \\| ", "",.) %>%
  strsplit(., split="/", fixed=F)

RRM_diag_terms1 <-
variab_form %>% stringr::str_detect(.,"RRM") %>% which(TRUE) %>%
  variab_form[.] %>% gsub( "1 \\| ", "",.)

RRM_diag_terms2 <-
  RRM_diag_terms1 %>%stringr::str_detect(.,"ped", negate = T) %>% which(TRUE) %>%
  RRM_diag_terms1[.]

RRM_diag_terms <- RRM_diag_terms2 %>% stringr::str_detect(.,"RRM ", negate = F) %>% which(TRUE) %>%
  RRM_diag_terms2[.] %>% gsub( "RRM \\| ", "",.)

RRM_nestedF_terms <-
  RRM_diag_terms2 %>% stringr::str_detect(.,"RRM ", negate = T) %>% which(TRUE) %>%
  RRM_diag_terms2[.] %>% gsub( "/RRM", "",.)


RRM_ped_terms <-
  RRM_diag_terms1 %>%stringr::str_detect(.,"ped", negate = F) %>% which(TRUE) %>%
  RRM_diag_terms1[.] %>% gsub( "RRM \\| ", "",.) %>%
    gsub( "ped ", "",.) %>% gsub( "\\| ", "",.)

trait_terms <-
 ter %>% .[[2]] %>% as.character() %>%
  gsub( "\\|", NA,.) %>% na.exclude() %>%
   c()
############################################
#####Get variable position and class########
variables_loc <-
  c(trait_terms,fixed_terms,diag_terms,ped_terms,
    RRM_ped_terms,RRM_diag_terms,
    nestedF_term,nested_diag_term,nested_ped_term,RRM_nestedF_terms) #need check nested!!

variables_loc_names=variables_loc %>% unlist() %>% unique()

variables_loc <-  variables_loc_names

 if(length(nested_diag_term)==0||length(nested_ped_term)==0)
  {nested_VC <- NULL}else{
   nested_VC <-c(nested_diag_term,nested_ped_term) %>% unlist() %>%
    matrix(ncol=2, byrow = T) %>% .[,2]
   }
names_VC <- variables_loc_names %>% data.frame() %>%
  dplyr::filter(.%in%c(diag_terms,ped_terms,
                nested_VC,
                RRM_ped_terms,RRM_diag_terms)) %>% .$.

names_VC_RRM <<- variables_loc_names %>% data.frame() %>%
  dplyr::filter(.%in%c(RRM_ped_terms,RRM_diag_terms)) %>% .$. %>% paste0(.,".temp")

for (i in 1:length(variables_loc)) {

  ifelse(i==1,variables_pos <- which(colnames(datarenum) == variables_loc[i]),
         variables_pos[i] <- which(colnames(datarenum) == variables_loc[i]))
}
variables_pos <- variables_pos %>%data.frame(terms_v=variables_loc_names,pos=.)
variables_class <-
  lapply(datarenum, class) %>% unlist()%>% data.frame() %>%
  data.frame(terms_v=rownames(.),class_v=.[,1]) %>% dplyr::select(terms_v,class_v)
##########PEDIGREE############################
if(is.null(ped_name)==F){
  ped_name=ped_name
  }else{ped_name="pedBLUPF90.txt"}
############SNP##############################
if(is.null(genotype_file)==T){
  snp_gen=c("")
}else{
  snp_gen <- paste("SNP_FILE","\n", genotype_file,"\n")}


#####Creating Legendre if RRM is true
if(length(RRM_ped_terms)>0||length(RRM_diag_terms)>0){
  #creating RRM based on the RRM column

  #  eval(call("=",paste0(as.name(Timevar),"_renum"),seq(min(datarenum[,Timevar]),
  #                                                      max(datarenum[,Timevar]))))
  Timemin=min(datarenum[,RRM_option$Timevar])
  Timemax=max(datarenum[,RRM_option$Timevar])
  #eval(call("=",as.name(Timevar),seq(Timemin,Timemax)))  # create dim Timemax range
  ##########################
  #function to bind columns
  ##########################
  cbind_poly <- function(x){
    cb1=rep(1,length(x))
    for (i in 1:RRM_option$poly) {
      cb1=cbind(cb1,x^i)
    }
    return(cb1)
  }
  ##########################
  #dir.create(paste0(getwd(),"/poly",RRM_option$poly))
  # find the files that you want and copy the files to the new folder
  # if(is.null(genotype_file)==F){
  #   file.copy(list.files(getwd(), genotype_file), paste0(getwd(),"/poly",RRM_option$poly))}
  # #file.copy(list.files(getwd(), ped_name), paste0(getwd(),"/poly",poly))
  # set as working directory the new folder
  # set as working directory the new folder
  # if(ped_name!="pedBLUPF90.txt"){
  #   file.copy(list.files(getwd(), ped_name), paste0(getwd(),"/poly",RRM_option$poly))}
  #
  # setwd(paste0(getwd(),"/poly",RRM_option$poly))

  ######################################################################################################
  lambda1 <-
    orthopolynom::legendre.polynomials(RRM_option$poly, normalized=T) %>% # calculate the polynomials normalized
    orthopolynom::polynomial.coefficients() %>% # extract the polynomial coefficients
    sapply(., '[', seq(max(sapply(., length)))) %>% # matrix with the polynomial coeff as upper triangular
    as.data.frame() %>%
    dplyr::mutate_all(~replace(., is.na(.), 0)) # substituting NA to 0

  beta = function(x,from=-1,to=1) from+(to-from)*((x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)))

  fi1 <-
    c(-1+2*(seq(Timemin,Timemax)-Timemin)/(Timemax-Timemin)) %>% #Timemax (x variables to regress for) and standardizing [-1;1]
    cbind_poly() %*%  #%*% # creating the variables for each polynomial order
    as.matrix(lambda1) %>% #multipling for the lambda ()
    data.frame(Timevar=eval(call("=",as.name(RRM_option$Timevar),seq(Timemin,Timemax))),.)

  fi1[,3:ncol(fi1)] <- apply(fi1[,3:ncol(fi1)],2,beta,from=0,to=1)
  fi1[,2] = 1

  colnames(fi1) <- c(as.name(RRM_option$Timevar),2:ncol(fi1))

  colnames(datarenum[RRM_option$Timevar]) <- RRM_option$Timevar

  data_final <- dplyr::left_join(datarenum,fi1, by=RRM_option$Timevar,keep=F)

  colnames(data_final) <- c(colnames(datarenum),paste0("poly",c(0:RRM_option$poly)))

  fi1<- fi1[,-1]

  write.table(fi1,"fi.txt",col.names = F,row.names = F)

  poly_pos <- data.frame(col=colnames(data_final)) %>% data.frame(col=.,pos=rownames(.)) %>%
    dplyr::filter(col%in%grep("poly", names(data_final), value = TRUE)) %>%
    .$pos %>% as.numeric()

  ################################
  ####COV start values############
  ################################
  if(class(VCA_RRM)[1]=="list"){
    cov_rrm_add <- matrix(VCA_RRM$COV,ncol = length(poly_pos),nrow = length(poly_pos))
    diag(cov_rrm_add) <- c(rep(VCA_RRM$VAR,length(poly_pos)))
  }else{
    cov_rrm_add <- VCA_RRM

  }
  if(class(VCD_RRM)[1]=="list"){
    cov_rrm_diag <- matrix(VCD_RRM$COV,ncol = length(poly_pos),nrow = length(poly_pos))
    diag(cov_rrm_diag) <- c(rep(VCD_RRM$VAR,length(poly_pos)))
  }else{
    cov_rrm_diag <- VCD_RRM

  }
}else{data_final=datarenum}


#####FIXED CLASS##############################
fixed_terms<- fixed_terms %>% data.frame(terms_v=.)
for (i in 1:nrow(fixed_terms)) {

  fcla <-fixed_terms %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,2]
  ifelse(i==1,class_fixed <- ifelse(fcla=='factor','cross alpha','cov'),
       class_fixed[i] <- ifelse(fcla=='factor','cross alpha','cov'))
}
#####Fixed position#########################################
fixed_pos <-
  fixed_terms %>%  dplyr::left_join(variables_pos,by="terms_v")%>% .[,2]
##############################################
fixed_card <-ifelse(is.na(class_fixed)==TRUE,fixed_card <- c(""),
                    paste(c("EFFECT\n"),
                  fixed_pos, class_fixed,"\n"))
####DIAGONAL CLASS########################
diag_terms<- diag_terms %>% data.frame(terms_v=.)

for (i in 1:nrow(diag_terms)) {
  dcla <- diag_terms %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,2]
  ifelse(i==1,class_diag <- ifelse(dcla=='factor','cross alpha','cov'),
         class_diag[i] <- ifelse(dcla=='factor','cross alpha','cov'))
}
#####Diagonal position#########################################
diag_pos <-
  diag_terms %>%  dplyr::left_join(variables_pos,by="terms_v")%>% .[,2]
##############################################
diag_card <- ifelse(is.na(class_diag)==TRUE,diag_card <- c(""),
  paste(c("EFFECT\n"),
          diag_pos, class_diag,"\n",
          "RANDOM\n","diagonal\n",
          "(CO)VARIANCES","\n",
          VCD,"\n"))
####PED CLASS#########################
ped_terms<- ped_terms %>% data.frame(terms_v=.)


for (i in 1:length(ped_terms)) {
  pedcla <- ped_terms %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,2]
  ifelse(i==1,class_ped <- ifelse(pedcla=='factor','cross alpha','cov'),
         class_ped[i] <- ifelse(pedcla=='factor','cross alpha','cov'))
}
#####Diagonal position#########################################
ped_pos <-
  ped_terms %>%  dplyr::left_join(variables_pos,by="terms_v")%>% .[,2]
#######################################
ped_card <-ifelse(is.na(class_ped)==TRUE,ped_card <- c(""),
       paste(c("EFFECT\n"),
        ped_pos, class_ped,"\n",
        "RANDOM\n","animal\n","FILE\n",ped_name,"\n",
        snp_gen,
        "(CO)VARIANCES","\n",
        VCA,"\n"))
####RRM DIAG CLASS#########################
RRM_diag_terms<- RRM_diag_terms %>% data.frame(terms_v=.)

for (i in 1:length(RRM_diag_terms)) {
  rrmdcla <- RRM_diag_terms %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,2]
  ifelse(i==1,class_RRM_diag <- ifelse(rrmdcla=='factor','cross alpha','cov'),
         class_RRM_diag[i] <- ifelse(rrmdcla=='factor','cross alpha','cov'))
}
#####RRM Diagonal position#################
RRM_diag_pos <-
  RRM_diag_terms %>%  dplyr::left_join(variables_pos,by="terms_v")%>% .[,2]

#########################################
#Create RRM diagonal card for all effects
#########################################
if(is.na(class_RRM_diag)==TRUE){RRM_diag_card <- c("")}else{
RRM_diag_card <-
	paste(c("EFFECT\n"),
        RRM_diag_pos, class_RRM_diag,"\n",
        "RANDOM\n",
      "diagonal\n",
        "RANDOM_REGRESSION\n",
        "data\n",
        "RR_POSITION\n",
      toString(poly_pos),"\n",
        "(CO)VARIANCES\n",
        toString(cov_rrm_diag),"\n"
      ) %>% gsub( ", ", " ",.)

}


####RRM PED CLASS#########################
RRM_ped_terms<- RRM_ped_terms %>% data.frame(terms_v=.)

for (i in 1:length(RRM_ped_terms)) {
  rrmpcla <- RRM_ped_terms %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,2]
  ifelse(i==1,class_RRM_ped <- ifelse(rrmpcla=='factor','cross alpha','cov'),
         class_RRM_ped[i] <- ifelse(rrmpcla=='factor','cross alpha','cov'))
}
#####RRM PED position#########################################
RRM_ped_pos <-
  RRM_ped_terms %>%  dplyr::left_join(variables_pos,by="terms_v")%>% .[,2]

#######################################

if(is.na(class_RRM_ped)==TRUE){RRM_ped_card <- c("")}else{
RRM_ped_card <-
	paste("EFFECT\n",
        RRM_ped_pos, class_RRM_ped,"\n",
        "RANDOM\n",
        "animal\n",
        "FILE\n",
        ped_name,"\n",
        snp_gen,
        "PED_DEPTH\n",
        PED_DEPTH,"\n",
        "INBREEDING\n",
        "pedigree\n",
        "RANDOM_REGRESSION\n",
        "data\n",
        "RR_POSITION\n",
        toString(poly_pos),"\n",
        "(CO)VARIANCES\n",
        toString(cov_rrm_add),"\n"
        ) %>% gsub( ", ", " ",.)}

#####NESTED variables FIXED class#################

nestedF_term<-
  nestedF_term %>% unlist() %>% data.frame(terms_v=.)

class_nestedF <- ifelse(is.na(nestedF_term)==T,class_nestedF <- c(),
                        class_nestedF <-
for (i in 1:nrow(nestedF_term)) {
  nestedF_class <-nestedF_term %>%
    dplyr::left_join(variables_class,by="terms_v") %>% .[i,2]

  ifelse(i==1,class_nestedF <- ifelse(nestedF_class=='factor','cross alpha','cov'),
         class_nestedF[i] <- ifelse(nestedF_class=='factor','cross alpha','cov'))
}

)
class_nestedF <- class_nestedF %>% matrix(ncol=2, byrow = T) %>%  data.frame()
#####NESTED variables FIXED position#################
nestedF_pos <-ifelse(is.na(class_nestedF)==T,nestedF_pos <- c(),
  nestedF_term %>% unlist() %>% data.frame(terms_v=.)%>%  dplyr::left_join(variables_pos,by="terms_v") %>% .[,2] %>%
  matrix(ncol=2, byrow = T) %>%  data.frame()
)
##########NESTED variables FIXED card####################################

nestedF_card <- c()
ifelse(is.na(nestedF_pos)==TRUE,nestedF_card <- c(""),
for (nv in 1:nrow(nestedF_pos)) {
  nestedF_card[nv] <-ifelse(is.na(class_nestedF)==TRUE,nestedF_card <- c(""),
                        paste(c("EFFECT\n"),
                              nestedF_pos[nv,2], class_nestedF[nv,2],"\n",
                              c("NESTED\n"),
                              nestedF_pos[nv,1], class_nestedF[nv,1],"\n")) %>% .[1,1]

}
)
#####NESTED variables Random diag class#################

if(length(nested_diag_term)==0){
                              class_nested_diag <- c("")
                              }else{
                                nested_diag_term<-
                                  nested_diag_term %>% unlist() %>% data.frame(terms_v=.)
                          for (i in 1:nrow(nested_diag_term)) {
                           nested_diag_class <-nested_diag_term %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,2]
                           ifelse(i==1,class_nested_diag <- ifelse(nested_diag_class=='factor','cross alpha','cov'),
                           class_nested_diag[i] <- ifelse(nested_diag_class=='factor','cross alpha','cov'))
                                    }
}
class_nested_diag <- class_nested_diag %>% matrix(ncol=2, byrow = T) %>%  data.frame()
#####NESTED variables Random diag position#################
nested_diag_pos <- if(length(nested_diag_term)==0){
                      nested_diag_pos <- c("")
                    }else{
                    nested_diag_term %>% unlist() %>% data.frame(terms_v=.)%>%
                        dplyr::left_join(variables_pos,by="terms_v") %>% .[,2] %>%
                    matrix(ncol=2, byrow = T) %>%  data.frame()
}
#####NESTED variables Random diag card#########################################
nested_diag_card <- c()
if(length(nested_diag_term)==0){
                          nested_diag_card <- c("")
                          }else{
                          for (nv in 1:nrow(nested_diag_pos)) {
                          nested_diag_card[nv] <-ifelse(is.na(class_nested_diag)==TRUE,nested_diag_card <- c(""),
                          paste(c("EFFECT\n"),
                                  nested_diag_pos[nv,2], class_nested_diag[nv,2],"\n",
                                  c("NESTED\n"),
                                  nested_diag_pos[nv,1], class_nested_diag[nv,1],"\n",
                                  "RANDOM\n","diagonal\n")) %>% .[1,1]
                                    }
}
#####NESTED variables Random PED class#################

nested_ped_term<-
  nested_ped_term %>% unlist() %>% data.frame(terms_v=.)

if(length(nested_ped_term)==0){
                                class_nested_ped <- c("")
                                }else{
                                for (i in 1:nrow(nested_ped_term)) {
                                nested_ped_class <-nested_ped_term %>%
                                  dplyr::left_join(variables_class,by="terms_v") %>% .[i,2]
                              ifelse(i==1,class_nested_ped <- ifelse(nested_ped_class=='factor','cross alpha','cov'),
         class_nested_ped[i] <- ifelse(nested_ped_class=='factor','cross alpha','cov'))
                                  }
}
class_nested_ped <- class_nested_ped %>% matrix(ncol=2, byrow = T) %>%  data.frame()
#####NESTED variables Random PED position#################
nested_ped_pos <-ifelse(length(nested_ped_term)==0,class_nested_ped <- c(""),
  nested_ped_term %>% unlist() %>% data.frame(terms_v=.)%>%  dplyr::left_join(variables_pos,by="terms_v") %>% .[,2] %>%
  matrix(ncol=2, byrow = T) %>%  data.frame()
)
#####NESTED variables Random PED card#########################################
nested_ped_card <- c()

if(length(nested_ped_term)==0){
                                class_nested_ped <- c("")
                                }else{
                                for (nv in 1:nrow(nested_ped_pos)) {
                          nested_ped_card[nv] <-ifelse(is.na(class_nested_ped)==TRUE,nested_ped_card <- c(""),
                                paste(c("EFFECT\n"),
                                      nested_ped_pos[nv,2], class_nested_ped[nv,2],"\n",
                                      c("NESTED\n"),
                                      nested_ped_pos[nv,1], class_nested_ped[nv,1],"\n",
                                      "RANDOM\n","animal\n","FILE\n",ped_name,"\n")) %>% .[1,1]

                                }
}
##########NESTED variables RRM fixed class#################

RRM_nestedF_terms<- RRM_nestedF_terms %>% data.frame(terms_v=.)
colnames(RRM_nestedF_terms) <- c("terms_v")

class_RRM_nestedF <- c()
if(length(RRM_nestedF_terms)==0){
    class_RRM_nestedF<- c("")
    }else{
    for (i in 1:nrow(RRM_nestedF_terms)) {
      RRM_nestedF_class <-RRM_nestedF_terms %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,2]
      ifelse(i==1,class_RRM_nestedF <- ifelse(RRM_nestedF_class=='factor','cross alpha','cov'),
     class_RRM_nestedF[i] <- ifelse(RRM_nestedF_class=='factor','cross alpha','cov'))
        }
}
#####NESTED variables RRM fixed position#################
RRM_nestedF_pos <-if(length(RRM_nestedF_terms)==0){
                 RRM_nestedF_pos<- c("")
                }else{
                RRM_nestedF_terms %>% unlist() %>% data.frame(terms_v=.)%>%
                    dplyr::left_join(variables_pos,by="terms_v") %>% .[,2]
                }
#####NESTED variables RRM fixed card#########################################

if(length(RRM_nestedF_pos)==0){
                            RRM_nestedF_card<- c("")
                            }else{
                              RRM_nestedF_card <- c()
                            for (nv in 1:length(RRM_nestedF_pos)) {
                            ifelse(length(RRM_nestedF_class)==0,RRM_nestedF_card1 <- c(""),
                            RRM_nestedF_card1 <-paste(c("EFFECT\n"),
                                  poly_pos, "cov","\n",
                                  c("NESTED\n"),
                                  RRM_nestedF_pos[nv], class_RRM_nestedF[nv],"\n"))

                              ifelse(is.na(class_RRM_nestedF)==TRUE,RRM_nestedF_card <- c(""),
                                     RRM_nestedF_card <-rbind(RRM_nestedF_card,RRM_nestedF_card1))
                                                                  }
                            RRM_nestedF_card <- RRM_nestedF_card %>%
                              matrix(ncol=length(RRM_nestedF_pos), byrow = T) %>% c()
}
###########TRAIT POSITION##########
trait_terms<- trait_terms %>% data.frame(terms_v=.)

trait_pos <-
  trait_terms %>%  dplyr::left_join(variables_pos,by="terms_v")%>% .[,2]
##################################
#CARD_data##########
if(is.null(genotype_file)==F)
  {geno_card <- paste("OPTION SNP_file", genotype_file)
  }else{geno_card <-c("")}

######################
#Residual Start matrix
######################
if(class(residual_start)[1]=="list"){
  residual_VC <- matrix(residual_start$COV,ncol = ncol(trait_terms),
                        nrow = ncol(trait_terms))
  diag(residual_VC) <- c(rep(residual_start$VAR,ncol(trait_terms)))
}else{
  residual_VC <- residual_start

}


###################################################
#fields_output_pos including  Heterogeneous variance
###################################################
#creating heterogeneous variance vector
if(is.null(het_res_variance)==F){
  het_var <- factor(data_final[,het_res_variance$het_levels]) %>% levels() %>% data.frame()
  colnames(het_var) <- c(het_res_variance$het_levels)
  het_var <-
    het_var %>% mutate(lev=rownames(.)) %>%
    mutate(h_var=cut(as.numeric(lev),het_res_variance$n_res_var,labels = F)) %>% dplyr::select(het_res_variance$het_levels,h_var)

  het_var <- my_dummy(het_var, select_columns = 'h_var')

  data_final<-
    data_final %>% merge.data.frame(.,het_var, by=het_res_variance$het_levels,all.x = T) %>%
    .[,c(colnames(datarenum),colnames(.[,seq(ncol(datarenum)+1,ncol(.))]))]

  data_final<- data_final[,c(colnames(datarenum),colnames(data_final[,seq(ncol(datarenum)+1,ncol(data_final))]))]
}


subj_field <- c(RRM_ped_pos,ped_pos)
het_levels_field <- colnames(data_final) %>%
  stringr::str_detect(.,"h_var", negate = F) %>% which(T)

het_eff_field <- which(colnames(data_final) == het_res_variance$het_levels)

fields_output <- which(names(datarenum)%in%fields_output)

fields_output_pos <-
  c(het_eff_field,fields_output,het_levels_field,subj_field) %>%
  unlist() %>% as.vector()
#################
###RENUM CARD
#################
sink(file = "renum.par")
cat("DATAFILE\n")
cat("data_renumf90.txt","\n")
cat("TRAITS\n")
cat(trait_pos,"\n")
cat("FIELDS_PASSED TO OUTPUT", "\n")
cat(fields_output_pos,"\n")
cat("WEIGHT(S)")
cat("\n")
cat(weights_object)
cat("\n")
cat("RESIDUAL VARIANCE\n")
cat(residual_VC,"\n")
cat(fixed_card)
cat(RRM_nestedF_card)
cat(diag_card)
cat(ped_card)
cat(RRM_ped_card)
cat(RRM_diag_card)
cat("OPTION missing",missing_values)
sink()


#######Fake pedigree#####
if(ped_name=="pedBLUPF90.txt"){
  pos_ID <- readr::read_lines("renum.par") %>%
    stringr::str_detect(.,"animal") %>% which(TRUE)-2

  id<-
    readLines("renum.par",warn = F)[pos_ID]%>%
    gsub(c("([a-z]+).*$"),"",x=.)%>% .[1] %>%
    strsplit(.,split = " ", fixed=TRUE) %>% unlist() %>%
    as.numeric() %>% unique() %>% .[-1] %>%
    data_final[,.] %>% levels()

  if(is.null(genotype_file)==F){
    c(id,data.table::fread(genotype_file,h=F,select = 1,data.table = F)$V1) %>%unique()%>%
      data.frame(ID=.,sire=rep(0,length(.)),dam=rep(0,length(.))) %>%
      write.table(.,"pedBLUPF90.txt",row.names = F, col.names = F, quote = F)
  }else{
    id %>%
      data.frame(ID=.,sire=rep(0,length(.)),dam=rep(0,length(.))) %>%
      write.table(.,"pedBLUPF90.txt",row.names = F, col.names = F, quote = F)
  }

}
# else{
#   setwd("../")
#   file.copy(list.files(getwd(), ped_name), paste0(getwd(),"/poly",RRM_option$poly))
#   setwd(paste0(getwd(),"/poly",RRM_option$poly))
# }
#Saving final dataset as data_renumf90.txt
data.table::fwrite(data_final,"data_renumf90.txt",sep=" ",row.names = F,col.names = F,quote = F,na=missing_values)

####Run RENUMF90
if(gsub("([0-9]|\\.)","",version$os)=="linux-gnu"){
  S_OP <- "Linux"
} else if(gsub("([0-9]|\\.)","",version$os)=="darwin"){
  S_OP <- "Mac_OSX"}else{
    S_OP <- "Windows"
  }
if(S_OP=="Windows"){
  renum <- paste0(blupf90_folder,"/renumf90.exe renum.par")
}else{
  renum <- paste0(blupf90_folder,"/renumf90 renum.par")
}

system(renum)

######################################################
####Find and add heterogeneous variance on renf90.par
######################################################
if(is.null(het_res_variance)==F){
  d_renf90 <-
    data.table::fread("renf90.dat",h=F,data.table = F)

  colnames(d_renf90) <- c(seq(1,ncol(d_renf90)-(het_res_variance$n_res_var+2)),
                          c("h_var",paste0("h_var_",c(1:het_res_variance$n_res_var)),"ori_ID"))


  hetres_pos <- colnames(d_renf90) %>% stringr::str_detect(.,"h_var_", negate = F) %>% which(T)
  hetres_pol <- c(het_res_variance$hetres_int,
                  rep(het_res_variance$hetres_int,length(hetres_pos)-1))

  if(het_res_variance$method=="blupf90+"){

    h_pos <- toString(c("OPTION hetres_pos", hetres_pos)) %>%
      gsub( ",", "",.)

    h_pol <- toString(c("OPTION hetres_pol", hetres_pol)) %>%
      gsub( ",", "",.)

    write(h_pos, file = "renf90.par", append = T)
    write(h_pol, file = "renf90.par", append = T)
  }else{
    het_col <- which(colnames(d_renf90) == "h_var")
    het_nlev <- d_renf90[,"h_var"] %>% factor() %>% levels() %>% length()
    h_pol <- toString(c("OPTION hetres_int", het_col, het_nlev)) %>%
      gsub( ",", "",.)

    write(h_pol, file = "renf90.par", append = T)
      }
  }
#Including OPTION on renf90.par

OPT_default <- list(yams=F,solution_mean=T,VCE=T,sol_se=T,Inbreeding=T,alpha_size=30,EM_REML=10,maxrounds=300,alpha_beta=c(0.95,0.05),tunedG=0,
  conv_crit=1e-10)

opt_list <- c(OPT_default[names(OPT_default)[!(names(OPT_default) %in% names(fit_option))]],fit_option)

sm <- ifelse(opt_list$solution_mean==T,paste("OPTION solution mean","\n"),"")
yams <- ifelse(opt_list$yams==T,paste("OPTION use_yams","\n"),"")
VCE <- ifelse(opt_list$VCE==T,paste("OPTION method VCE","\n"),"")
sol_se <- ifelse(opt_list$sol_se==T,paste("OPTION sol se","\n"),"")
imb <- ifelse(opt_list$Inbreeding==T,paste("OPTION Inbreeding","\n"),"")
alps <- ifelse(opt_list$alpha_size>0,paste("OPTION alpha_size",opt_list$alpha_size,"\n"),"")
EM_reml <- ifelse(opt_list$EM_REML>0,paste("OPTION EM-REML",(opt_list$EM_REML),"\n"),"")
max_r <- ifelse(opt_list$maxrounds>0,paste("OPTION maxrounds",opt_list$maxrounds,"\n"),"")
alph_b <- ifelse(opt_list$alpha_beta>0,paste("OPTION AlphaBeta",opt_list$alpha_beta[1],opt_list$alpha_beta[2],"\n"),"")
tuG <- ifelse(opt_list$tunedG!=0,paste("OPTION tunedG",opt_list$tunedG,"\n"),paste("OPTION tunedG",0,"\n"))
conv_crit <- ifelse(opt_list$conv_crit!=1e-10,paste("OPTION conv_crit",opt_list$conv_crit,"\n"),paste("OPTION tunedG",1e-10,"\n"))

oriID <- "OPTION origID"
opt_blupf90 <- unique(paste0(sm,VCE,yams,sol_se,imb,alps,
                             EM_reml,max_r,alph_b,tuG,conv_crit,oriID))
##############################################################################
write(opt_blupf90, file = "renf90.par", append = T)
#Extra OPTION within paremeter file
write(extra.option.blup, file = "renf90.par", append = T)
######

data.table::fread("renf90.fields",select = c("field","origfield"),data.table = F) %>%
  mutate(names=names(data_final)[.$origfield]) %>%
  write.table(.,"cols_data_renumf90.txt",row.names = F,quote = T)
 }

