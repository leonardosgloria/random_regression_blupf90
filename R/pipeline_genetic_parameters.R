##########################################
#Function to calculate genetics parameters
##########################################
pipeline_genetic_parameters <- function(form,logname,parms_card,DAPmin,DAPmax,method){
  
  #######################################
  ################Calculate heritability
  ###########################################################################
  DAP <-seq(DAPmin,DAPmax)
  ####Getting the number of traits
  
  
  number_traits <- readLines(parms_card) %>%
    .[grep("NUMBER_OF_TRAITS", .)+1] %>% 
    as.numeric()
  
  #### polynomial columns
  pos_poli_ebv <-  read_lines(parms_card) %>%
    .[grep("add",.)-2] %>% gsub(c("([a-z]+).*$"),"",x=.)%>% .[1] %>% 
    strsplit(.,split = "    ", fixed=TRUE) %>% unlist() %>% 
    as.numeric() %>% na.exclude() %>% as.numeric()
  
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
  
  ter = terms(as.formula(form))
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
  
  #temp <- list.files(pattern = "*.temp")
  #myfiles <- lapply(temp,read.table)
  #names(myfiles) <- paste0(VC_names,".temp") %>% as.factor()
  
  #pos_VC_a <- read_lines(logname) %>% 
  #  str_detect(.,paste("effect", pos_poli_ebv[1])) %>% which(TRUE) %>% +1
  
  #VC_a=myfiles[paste0(which(pos_VC%in%pos_VC_a),".temp")]
  
  ## polynomial columns
  pos_poli_ebv <-  read_lines(parms_card) %>%
    .[grep("add",.)-2] %>% gsub(c("([a-z]+).*$"),"",x=.)%>% .[1] %>% 
    strsplit(.,split = "    ", fixed=TRUE) %>% unlist() %>% 
    as.numeric() %>% na.exclude() %>% as.numeric()
  
  len_VC_ <- VC_ %>% map(.,length) %>% unlist()
  
  names_VC_RRM <- c(RRM_ped_terms,RRM_diag_terms)
  ######################
  #legendre coeff matrix
  ######################
  cbind_poly <- function(x){
    cb1=rep(1,length(x))
    for (i in 1:(length(pos_poli_ebv)-1)) {
      cb1=cbind(cb1,x^i)
    }
    return(cb1)
    
  }
  
  lambda1 <- 
    legendre.polynomials(length(pos_poli_ebv)-1, normalized=T) %>% # calculate the polynomials normalized
    polynomial.coefficients() %>% # extract the polynomial coefficients
    sapply(., '[', seq(max(sapply(., length)))) %>% # matrix with the polynomial coeff as upper triangular
    as.data.frame() %>% 
    mutate_all(~replace(., is.na(.), 0)) # substituting NA to 0
  
  #fi1 <-
  #  c(-1+2*(seq(DAPmin,DAPmax)-DAPmin)/(DAPmax-DAPmin)) %>% #DAP (x variables to regress for) and standardizing [-1;1]
  #  cbind_poly() %*%  #%*% # creating the variables for each polynomial order
  #  as.matrix(lambda1) #multipling for the lambda ()
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
  
  ###################################################
  ## creating and Saving Residual variance components
  ###################################################
  
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
  #write.table(res_matrix,file = "res_coef_vc.temp",row.names = F,col.names = F)
  
  
  h2_to_plot <- 
    if(is.null(NonRRM_VC)==F){
      data.frame(RRMeff,NonRRM_VC,res_matrix) %>%
        mutate(den=rowSums(.)) %>% 
        mutate(h2=RRMeff[,RRM_ped_terms]/den,DAP=DAP) %>% select(DAP,h2)
    }else{
      data.frame(RRMeff,res_matrix) %>%
        mutate(den=rowSums(.)) %>% 
        mutate(h2=RRMeff[,RRM_ped_terms]/den,DAP=DAP) %>% select(DAP,h2)
    }
  
  genetic_param_DAP <-
    if(is.null(NonRRM_VC)==F){
      data.frame(RRMeff,NonRRM_VC,res_matrix) %>%
        mutate(den=rowSums(.)) %>% 
        mutate(h2=RRMeff[,RRM_ped_terms]/den,DAP=DAP)
    }else{
      data.frame(RRMeff,res_matrix) %>%
        mutate(den=rowSums(.)) %>% 
        mutate(h2=RRMeff[,RRM_ped_terms]/den,DAP=DAP)
    }
  
  write.table(h2_to_plot, "h2.txt", row.names = F, col.names = T, quote = F)
  write.table(genetic_param_DAP, "VC_DAP.txt", row.names = F, col.names = T, quote = F)
  
  ####################
  ### GEBV############
  ####################
  if(method=="blupf90+"){solut="solutions"}else{
    solut="final_solutions"
  }
  
  comp = data.table::fread(solut,skip=1,data.table = F)
  
  b <- comp %>% filter(V2 %in%pos_poli_ebv) %>%
    select(c(V2,V4))%>% data.frame() %>% 
    pivot_wider(names_from = V2, values_from = V4,values_fn = list) %>%
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
    fread(ped_name,data.table = F,select = id_pos_renfdat)
  
  EBV_DAP <- 
    as.matrix(b)%*%t(fi1) %>% 
    data.frame(ID=get_original_id,.) %>% #check the ID 
    pivot_longer(cols = starts_with("X"),values_to = "EBV",names_to = "DAP") %>% 
    data.frame()
  
  EBV_DAP$DAP <- EBV_DAP$DAP %>% gsub("\\X","",.) %>% as.numeric()
  colnames(EBV_DAP) <- c("ID","DAP","EBV")
  
  EBV_t <- EBV_DAP %>% aggregate(EBV ~ ID, ., sum)
  
  data.table::fwrite(EBV_DAP, "EBV_DAP.txt", row.names = F, col.names = T,sep = " ")
  data.table::fwrite(EBV_t, "EBV_t.txt", row.names = F, col.names = T,sep = " ")
  
  ###################
  #Theorical accuracy
  ###################
  
  d_add <- diag(as.matrix(vc_matrix[RRM_ped_terms][[1]])) %>%
    data.frame(varg=.) %>%
    mutate(effect_=pos_poli_ebv) # get Varg
  
  inbreeding_ <-
    read.table("renf90.inb",col.names = c("ori_id","inbreeding_","level")) 
  
  th_acc <-
    read.table(solut,skip = 1,col.names=c("trait","effect_","level","solution","s.e.")) %>% #reading solutions output
    filter(effect_%in%pos_poli_ebv) %>% # filtering coef solutions
    select(effect_,level,solution,s.e.) %>%  
    right_join(.,d_add,by="effect_") %>% # merging with Varg
    right_join(.,inbreeding_,by="level") %>% # merging with ID inbreeding 
    mutate(th_acc=(1-(s.e.^2/((1+inbreeding_)*varg)))^0.5) #calculating Theorical accuracy
  
  write.table(th_acc, "ACC_th.txt", row.names = F, col.names = T,quote = F)
}
