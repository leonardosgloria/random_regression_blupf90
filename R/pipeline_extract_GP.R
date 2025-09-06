##########################################
#Function to extract genetics parameters
##########################################
pipeline_extract_vc <- function(datarenum,formula,logname,parms_card,method,RRM_option){


  ####Getting the number of traits

  number_traits <- readLines(parms_card) %>%
    .[grep("NUMBER_OF_TRAITS", .)+1] %>%
    as.numeric()

  ####Finding random group in renf90.par

  pos_random <- readLines(parms_card) %>%grep(" RANDOM_GROUP", .)+1

  VC_ <- list()
  for (r in 1:length(pos_random)) {
    VC_[[r]] <-unlist(strsplit(readLines(parms_card)[pos_random[r]], split=" ", fixed=TRUE)) %>%
      as.numeric() %>% na.exclude() %>% c()
  }
  #################
  #### VC names####
  #################

  ter = terms(as.formula(formula))
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

  if(is.null(RRM_option)==F){
    Timemin=min(datarenum[,RRM_option$Timevar])
    Timemax=max(datarenum[,RRM_option$Timevar])
    Time_var <-seq(Timemin,Timemax)
    #### polynomial columns
    pos_poli_ebv <-  read_lines(parms_card) %>%
      .[grep("add",.)-2] %>% gsub(c("([a-z]+).*$"),"",x=.)%>% .[1] %>%
      strsplit(.,split = "    ", fixed=TRUE) %>% unlist() %>%
      as.numeric() %>% na.exclude() %>% as.numeric()

    len_VC_ <- VC_ %>% map(.,length) %>% unlist()

    names_VC_RRM <- c(RRM_ped_terms,RRM_diag_terms)

    ######################
    #legendre coeff matrix
    ######################
    fi1 <- as.matrix(read.table("fi.txt",h=F))
    ######

    RRMeff <-
      vc_matrix[names_VC_RRM] %>%
      map_dfr(~ fi1%*%as.matrix(.)%*%t(fi1)) %>%
      map_dfr(~ diag(.)) %>% data.frame()

    '%ni%' <- Negate('%in%') # non within

    non_RR <- vc_matrix[VC_names%ni%names_VC_RRM] %>% rep(.,nrow(RRMeff)) %>% unlist()

    if(is.null(non_RR)==F){
      NonRRM_VC <-
        vc_matrix[VC_names%ni%names_VC_RRM] %>% rep(.,nrow(RRMeff)) %>% unlist() %>%
        matrix(ncol=length(vc_matrix[VC_names%ni%names_VC_RRM]),byrow = T) %>%
        data.frame()
      colnames(NonRRM_VC) <- names(vc_matrix[VC_names%ni%names_VC_RRM])
    }else{
      NonRRM_VC <-NULL
    }
  }else{
    RRMeff <-as.data.frame(1)
  }

  ### Residual variance components
  if(method=="blupf90+"){res_name="Res"}else{
    res_name="R"
  }

  pos_res <- read_lines(logname) %>%
    str_detect(.,res_name) %>% which(TRUE) %>% +1


  res_matrix <- read_lines(logname)[seq(pos_res,number_traits+pos_res-1)] %>%
    strsplit(.,split = " ", fixed=TRUE) %>% unlist() %>%
    as.numeric() %>% na.exclude() %>%
    matrix(ncol = number_traits,nrow = number_traits) %>%
    rep(.,nrow(RRMeff))

  res_out <- read_lines(logname)[seq(pos_res,number_traits+pos_res-1)] %>%
    strsplit(.,split = " ", fixed=TRUE) %>% unlist() %>%
    as.numeric() %>% na.exclude() %>%
    matrix(ncol = number_traits,nrow = number_traits)

  #Heritability for RRM
  if(is.null(RRM_option)==F){
  h2_to_plot <-
  if(is.null(NonRRM_VC)==F){
    data.frame(RRMeff,NonRRM_VC,res_matrix) %>%
      mutate(den=rowSums(.)) %>%
      mutate(h2=RRMeff[,RRM_ped_terms]/den,Time_var=Time_var) %>% select(Time_var,h2)
  }else{
    data.frame(RRMeff,res_matrix) %>%
      mutate(den=rowSums(.)) %>%
      mutate(h2=RRMeff[,RRM_ped_terms]/den,Time_var=Time_var) %>% select(Time_var,h2)
  }

genetic_param_Time_var <-
  if(is.null(NonRRM_VC)==F){
    data.frame(RRMeff,NonRRM_VC,res_matrix) %>%
      mutate(den=rowSums(.)) %>%
      mutate(h2=RRMeff[,RRM_ped_terms]/den,Time_var=Time_var)
  }else{
    data.frame(RRMeff,res_matrix) %>%
      mutate(den=rowSums(.)) %>%
      mutate(h2=RRMeff[,RRM_ped_terms]/den,Time_var=Time_var)
  }

 write.table(h2_to_plot, "h2.txt", row.names = F, col.names = T, quote = F)
 write.table(genetic_param_Time_var, "VC_Time_var.txt", row.names = F, col.names = T, quote = F)

############################
### GEBV for RRM############
############################
if(method=="blupf90+"){solut="solutions.orig"}else{
  solut="final_solutions"
}


b <- data.table::fread(solut,skip=1,data.table = F) %>%
  filter(V2 %in%pos_poli_ebv) %>%
  select(c(V2,V5))%>% data.frame() %>%
  pivot_wider(names_from = V2, values_from = V5,values_fn = list) %>%
  unnest(cols = paste(pos_poli_ebv)) %>%
  data.frame()

ped_name <-
  read_lines(parms_card) %>%
  str_detect(.,"renadd") %>% which(TRUE) %>%
  read_lines(parms_card)[.] %>%
  gsub(" ","",.)

id_pos_renfdat <-
  data.table::fread(ped_name,skip=1,
                    data.table = F,nrows = 1) %>% ncol()

get_original_id <-
  data.table::fread(ped_name,data.table = F,select = id_pos_renfdat)

EBV_Time_var <-
  as.matrix(b)%*%t(fi1) %>%
  data.frame(ID=get_original_id,.) %>% #check the ID
  pivot_longer(cols = starts_with("X"),values_to = "EBV",names_to = RRM_option$Timevar) %>%
  data.frame()

EBV_Time_var[,RRM_option$Timevar] <- EBV_Time_var[,RRM_option$Timevar] %>% gsub("\\X","",.) %>% as.numeric()


colnames(EBV_Time_var) <- c("ID",RRM_option$Timevar,"EBV")

EBV_t <- EBV_Time_var %>% aggregate(EBV ~ ID, ., sum)

data.table::fwrite(EBV_Time_var, "EBV_Time_var.txt", row.names = F, col.names = T,sep = " ")
data.table::fwrite(EBV_t, "EBV_t.txt", row.names = F, col.names = T,sep = " ")
###################
#Theorical accuracy
###################
ncol_sol <- data.table::fread(solut,skip=1,data.table = F,nrows = 1) %>% ncol()
if(ncol_sol==6){
  d_add <- diag(as.matrix(vc_matrix[RRM_ped_terms][[1]])) %>%
    data.frame(varg=.) %>%
    mutate(effect_=pos_poli_ebv) # get Varg

  inbreeding_ <-
    read.table("renf90.inb",col.names = c("ori_id","inbreeding_","level"))

    col_sol <- c("trait","effect_","level","original_id","solution","s.e.")

  th_acc <-
    data.table::fread(solut,skip=1,data.table = F,col.names = col_sol) %>% #reading solutions output
    filter(effect_%in%pos_poli_ebv) %>% # filtering coef solutions
    select(effect_,level,solution,s.e.) %>%
    right_join(.,d_add,by="effect_") %>% # merging with Varg
    right_join(.,inbreeding_,by="level") %>% # merging with ID inbreeding
    #mutate(th_acc=((1-(s.e.^2))/((1+inbreeding_)*varg))^0.5 )
    mutate(th_acc=(1-(s.e.^2/((1+inbreeding_)*varg)))^0.5) #calculating Theorical accuracy

  write.table(th_acc, "ACC_th.txt", row.names = F, col.names = T,quote = F)
              }
  }
##############
#OUTPUT
##############
VCA_out <- vc_matrix[ped_terms] %>% unlist() %>% list(GENETIC=.)
VCD_out <- vc_matrix[diag_terms] %>% unlist()%>% list(RANDOM=.)
residual_out <- res_out%>% data.frame(RESIDUAL=.)
RRM_random <-  vc_matrix[RRM_diag_terms] %>% unlist()%>% list(RRM_RANDOM=.)
RRM_genetic <- vc_matrix[RRM_ped_terms] %>% unlist()%>% list(RRM_GENETIC=.)
Variance_Components_model <- c(VCD_out,VCA_out,RRM_random,RRM_genetic,residual_out)
return(Variance_Components_model)

}
