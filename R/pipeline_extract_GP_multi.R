##########################################
#Function to extract genetics parameters
##########################################
pipeline_extract_vc_multi <- function(datarenum,formula,logname="blupf90.log",parms_card,method,RRM_option=NULL){

  ####Getting the number of traits

  number_traits <- readLines(parms_card) %>%
    .[grep("NUMBER_OF_TRAITS", .)+1] %>%
    as.numeric()

  for (mod in 1:length(formula)){

    form <- formula[[mod]]

    ter = terms(form)
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

    names_VC_RRM <-  variables_loc_names %>% data.frame() %>%
      dplyr::filter(.%in%c(RRM_ped_terms,RRM_diag_terms)) %>% .$. %>% paste0(.,".temp")

    for (i in 1:length(variables_loc)) {

      ifelse(i==1,variables_pos <- which(colnames(datarenum) == variables_loc[i]),
             variables_pos[i] <- which(colnames(datarenum) == variables_loc[i]))
    }
    variables_pos <- variables_pos %>%data.frame(terms_v=variables_loc_names,pos=.)
    variables_class <-
      lapply(datarenum, class) %>% unlist()%>% data.frame() %>%
      data.frame(terms_v=rownames(.),class_v=.[,1]) %>% .[,-1]
    ###########TRAIT POSITION##########
    trait_terms<- trait_terms %>% data.frame(terms_v=.)

    trait_pos <-
      trait_terms %>%  dplyr::left_join(variables_pos,by="terms_v")%>% .[,2]

    if(mod==1){
      trait_pos_f <- trait_pos}else{
        trait_pos_f <- cbind(trait_pos_f,trait_pos)
      }
    ntraits <- ncol(trait_pos_f) #number of traits

    #####FIXED terms##############################
    fixed_terms<- fixed_terms %>% data.frame(terms_v=.) %>%
      dplyr::arrange(terms_v)

    #####Fixed position#########################################
    fixed_pos <-
      fixed_terms %>%  dplyr::left_join(variables_pos,by="terms_v")%>% .[,2] %>%
      data.frame(terms_v=.)
    ####DIAGONAL terms########################
    diag_terms<- diag_terms %>% data.frame(terms_v=.)

    #####Diagonal position#########################################
    diag_pos <-
      diag_terms %>%  dplyr::left_join(variables_pos,by="terms_v")%>% .[,2] %>%
      data.frame(terms_v=.)

    ####PED terms#########################
    ped_terms<- ped_terms %>% data.frame(terms_v=.)
    #####PED position#########################################
    ped_pos <-
      ped_terms %>%  dplyr::left_join(variables_pos,by="terms_v")%>% .[,2] %>%
      data.frame(terms_v=.)
    ##########NESTED variables RRM fixed class#################

    RRM_nestedF_terms<- RRM_nestedF_terms %>% data.frame(terms_v=.)
    colnames(RRM_nestedF_terms) <- c("terms_v")

    #####NESTED variables RRM fixed position#################
    RRM_nestedF_pos <- if(length(RRM_nestedF_terms)==0){
      RRM_nestedF_pos<- c("")
    }else{
      RRM_nestedF_terms %>% unlist() %>% data.frame(terms_v=.)%>%
        dplyr::left_join(variables_pos,by="terms_v") %>% .[,2] %>%
        data.frame(terms_v=.)
    }
    ####RRM DIAGONAL terms########################
    RRM_diag_terms<- RRM_diag_terms %>% data.frame(terms_v=.)

    #####RRM Diagonal position#########################################
    RRM_diag_pos <-
      RRM_diag_terms %>%  dplyr::left_join(variables_pos,by="terms_v")%>% .[,2] %>%
      data.frame(terms_v=.)
    ####RRM pedigree terms########################
    RRM_ped_terms<- RRM_ped_terms %>% data.frame(terms_v=.)

    #####RRM pedigree position#########################################
    RRM_ped_pos <-
      RRM_ped_terms %>%  dplyr::left_join(variables_pos,by="terms_v")%>% .[,2] %>%
      data.frame(terms_v=.)

    #####Function terms#########################################

    unique_terms_eff <- function(eff,eff_n){

      unique_terms <- eff%>%
        #dplyr::na_if(.,0) %>%
        apply(., 1, na.exclude) %>%
        lapply(., unique,sep=' ', na.rm=T) %>% unlist() %>% data.frame(terms_v=.)

      if(nrow(unique_terms)==0){
        unique_terms <-  unique_terms %>%
          matrix(.,nrow = 1,ncol =1 ) %>%
          data.frame(terms_v=.)}else{
            unique_terms <- unique_terms %>%
              as.matrix(.,nrow = nrow(unique_terms),ncol =1 )%>%
              data.frame(terms_v=.)}

      unique_terms <- apply(unique_terms,2,function(y) sapply(y,function(x) ifelse(is.null(x),0, x)))

      unique_terms <- unique_terms %>% as.data.frame() %>%
        #dplyr::na_if(.,0) %>%
        apply(., 1, na.exclude) %>%
        lapply(., unique,sep=' ', na.rm=T) %>% unlist() %>% data.frame(terms_v=.)

      if(nrow(unique_terms)==0){
        unique_terms <-  unique_terms %>%
          matrix(.,nrow = 1,ncol =1 ) %>%
          data.frame(terms_v=.)}else{
            unique_terms <- unique_terms %>%
              as.matrix(.,nrow = nrow(unique_terms),ncol =1 )%>%
              data.frame(terms_v=.)}

      names(unique_terms)[1] <- c('terms_v')

      unique_terms$terms_v <- unique_terms$terms_v %>% as.character()

      unique_terms <- unique_terms %>%
        dplyr::right_join(.,tidyr::as_tibble(eval(parse(text=paste0(eff_n,'_terms')))),by='terms_v',keep=T)

      unique_terms[,2] <- as.character(unique_terms[,2])

      eff <-
        eff %>%
        dplyr::full_join(.,data.frame(terms_v=unique_terms[,2]),by='terms_v',keep=T) %>%
        data.frame(.)


      eff <-
        eff %>% replace(is.na(.),0)
    }
    #####Function position#########################################
    unique_pos_eff <- function(eff,eff_n){

      unique_pos <- eff%>%
        #dplyr::na_if(.,0) %>%
        apply(., 1, na.exclude) %>%
        lapply(., unique,sep=' ', na.rm=T) %>% unlist() %>% data.frame(terms_v=.)

      unique_pos <- unique_pos %>%
        #dplyr::na_if(.,0) %>%
        apply(., 1, na.exclude) %>%
        lapply(., unique,sep=' ', na.rm=T) %>% unlist() %>% data.frame(terms_v=.)


      if(nrow(unique_pos)==0){
        unique_pos <-  unique_pos %>%
          matrix(.,nrow = 1,ncol =1 ) %>%
          data.frame(terms_v=.)}else{
            unique_pos <- unique_pos %>%
              as.matrix(.,nrow = nrow(unique_pos),ncol =1 )%>%
              data.frame(terms_v=.)}

      unique_pos <- apply(unique_pos,2,function(y) lapply(y,function(x) ifelse(is.null(x),0, x)))

      names(unique_pos)[1] <- c('terms_v')

      unique_pos$terms_v <- unique_pos$terms_v %>% unlist() %>% as.integer()

      unique_pos <- unique_pos %>% as.data.frame() %>%
        dplyr::right_join(.,tidyr::as_tibble(eval(parse(text=paste0(eff_n,'_pos')))),by='terms_v',keep=T) %>%
        data.frame(terms_v=.)

      unique_pos[,2] <- as.integer(unique_pos[,2])

      eff <-
        eff %>%
        dplyr::full_join(.,data.frame(terms_v=unique_pos[,2]),by='terms_v',keep=T) %>%
        data.frame(.)

      eff <-
        eff %>% replace(is.na(.),0) %>% unique()
      return(eff)
    }

    ########Function join eff ###################################
    join_row_eff <- function(A){
      for(i in 1:nrow(A)){
        paste("i:",i)
        for(j in 1:nrow(A)){
          paste("j:",j)
          if(i != j){
            if(sum((A[i,] %in% A[j,]) & A[i,] != 0) > 0){
              A[i,] = A[i,] + A[j,]
              A[j,] = A[j,] - A[j,]
            }
          }
          else{
            c("-----")
          }
        }
      }
      A <-  A[which(rowSums(A) > 0),]
    }
    ########Function get_terms ###################################
    get_terms <- function(eff){
      fx <- eff %>%
        #dplyr::na_if(.,0) %>%
        apply(., 1, na.exclude) %>% lapply(., unique,sep=' ', na.rm=T) %>%
        unlist() %>% data.frame(terms_v=.) %>% unique()
      colnames(datarenum)[c(fx[,1])] %>% noquote(.) %>% data.frame(terms_v=as.character(.))
    }



    ###########################################################
    if(mod==1){
      #####Fixed terms#########################################

      fixed_terms_temp <-
        tidyr::as_tibble(fixed_terms) %>%
        dplyr::full_join(.,tidyr::as_tibble(fixed_terms),by='terms_v',keep=F) %>%
        dplyr::mutate_all(~replace(., .== NA, 0)) %>%
        data.frame(.)

      fixed_terms_f <- fixed_terms_temp %>% data.frame(terms_v=.)
      #####Fixed position#########################################
      fixed_pos_temp <-
        tidyr::as_tibble(fixed_pos) %>%
        dplyr::full_join(.,tidyr::as_tibble(fixed_pos),by='terms_v',keep=F) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(),~tidyr::replace_na(.,0))) %>%
        data.frame(.)

      fixed_pos_f <- fixed_pos_temp %>% data.frame(terms_v=.)

      #####Diagonal terms#########################################

      diag_terms_temp <-
        tidyr::as_tibble(diag_terms) %>%
        dplyr::full_join(.,tidyr::as_tibble(diag_terms),by='terms_v',keep=F) %>%
        dplyr::mutate_all(~replace(., .== NA, 0))  %>%
        data.frame(.)

      diag_terms_f <- diag_terms_temp %>% data.frame(terms_v=.)
      #####Diagonal position#########################################
      diag_pos_temp <-
        tidyr::as_tibble(diag_pos) %>%
        dplyr::full_join(.,tidyr::as_tibble(diag_pos),by='terms_v',keep=F) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(),~tidyr::replace_na(.x,0))) %>%
        data.frame(.)

      diag_pos_f <- diag_pos_temp %>% data.frame(terms_v=.)
      #####Ped terms#########################################

      ped_terms_temp <-
        tidyr::as_tibble(ped_terms) %>%
        dplyr::full_join(.,tidyr::as_tibble(ped_terms),by='terms_v',keep=F) %>%
        dplyr::mutate_all(~replace(., .== NA, 0))  %>%
        data.frame(.)

      ped_terms_f <- ped_terms_temp %>% data.frame(terms_v=.)
      #####Ped position#########################################
      ped_pos_temp <-
        tidyr::as_tibble(ped_pos) %>%
        dplyr::full_join(.,tidyr::as_tibble(ped_pos),by='terms_v',keep=F) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(),~tidyr::replace_na(.x,0))) %>%
        data.frame(.)

      ped_pos_f <- ped_pos_temp %>% data.frame(terms_v=.)
      #####NESTED variables RRM fixed terms######################################
      RRM_nestedF_terms_temp <-
        tidyr::as_tibble(RRM_nestedF_terms) %>%
        dplyr::full_join(.,tidyr::as_tibble(RRM_nestedF_terms),by='terms_v',keep=F) %>%
        dplyr::mutate_all(~replace(., .== NA, 0))  %>%
        data.frame(.)
      #####NESTED variables RRM fixed position###################################
      RRM_nestedF_pos_temp <-
        tidyr::as_tibble(RRM_nestedF_pos) %>%
        dplyr::full_join(.,tidyr::as_tibble(RRM_nestedF_pos),by='terms_v',keep=F) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(),~tidyr::replace_na(.x,0))) %>%
        data.frame(.)
      RRM_nestedF_pos_f <- RRM_nestedF_pos_temp %>% data.frame(terms_v=.)
      #####RRM Diagonal terms#########################################

      RRM_diag_terms_temp <-
        tidyr::as_tibble(RRM_diag_terms) %>%
        dplyr::full_join(.,tidyr::as_tibble(RRM_diag_terms),by='terms_v',keep=F) %>%
        dplyr::mutate_all(~replace(., .== NA, 0))  %>%
        data.frame(.)

      RRM_diag_terms_f <- RRM_diag_terms_temp %>% data.frame(terms_v=.)
      #####RRM Diagonal position#########################################
      RRM_diag_pos_temp <-
        tidyr::as_tibble(RRM_diag_pos) %>%
        dplyr::full_join(.,tidyr::as_tibble(RRM_diag_pos),by='terms_v',keep=F) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(),~tidyr::replace_na(.x,0))) %>%
        data.frame(.)

      RRM_diag_pos_f <- RRM_diag_pos_temp %>% data.frame(terms_v=.)
      #####RRM Pedigree terms#########################################

      RRM_ped_terms_temp <-
        tidyr::as_tibble(RRM_ped_terms) %>%
        dplyr::full_join(.,tidyr::as_tibble(RRM_ped_terms),by='terms_v',keep=F) %>%
        dplyr::mutate_all(~replace(., .== NA, 0))  %>%
        data.frame(.)

      RRM_ped_terms_f <- RRM_ped_terms_temp %>% data.frame(terms_v=.)
      #####RRM Pedigree position#########################################
      RRM_ped_pos_temp <-
        tidyr::as_tibble(RRM_ped_pos) %>%
        dplyr::full_join(.,tidyr::as_tibble(RRM_ped_pos),by='terms_v',keep=F) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(),~tidyr::replace_na(.x,0))) %>%
        data.frame(.)

      RRM_ped_pos_f <- RRM_ped_pos_temp %>% data.frame(terms_v=.)
    }else{

      #####Fixed terms#########################################
      fixed_terms_temp <- unique_terms_eff(fixed_terms_temp,eff_n='fixed')

      #####fixed position#########################################

      fixed_pos_temp <- unique_pos_eff(fixed_pos_temp,eff_n='fixed')

      #####Diagonal terms#########################################
      diag_terms_temp  <- unique_terms_eff(diag_terms_temp,eff_n='diag')

      #####Diagonal position#########################################
      diag_pos_temp <- unique_pos_eff(diag_pos_temp,eff_n='diag')

      #####ped terms#########################################
      ped_terms_temp  <- unique_terms_eff(ped_terms_temp,eff_n='ped')

      #############################PED POS#############################
      ped_pos_temp <- unique_pos_eff(ped_pos_temp,eff_n='ped')

      #####NESTED variables RRM fixed terms######################################
      RRM_nestedF_terms_temp <-  unique_terms_eff(RRM_nestedF_terms_temp,eff_n='RRM_nestedF')

      #####NESTED variables RRM fixed position###################################
      RRM_nestedF_pos_temp <-   unique_pos_eff(RRM_nestedF_pos_temp,eff_n='RRM_nestedF')

      #####RRM Diagonal terms#########################################
      RRM_diag_terms_temp  <- unique_terms_eff(RRM_diag_terms_temp,eff_n='RRM_diag')

      #####RRM Diagonal position#########################################
      RRM_diag_pos_temp <- unique_pos_eff(RRM_diag_pos_temp,eff_n='RRM_diag')

      #####RRM pedigree terms#########################################
      RRM_ped_terms_temp  <- unique_terms_eff(RRM_ped_terms_temp,eff_n='RRM_ped')

      #####RRM pedigree position#########################################
      RRM_ped_pos_temp <- unique_pos_eff(RRM_ped_pos_temp,eff_n='RRM_ped')

    }
    ############################################
    names(fixed_terms_temp)[1] <- c('terms_v')
    names(fixed_pos_temp)[1] <- c('terms_v')
    names(diag_terms_temp)[1] <- c('terms_v')
    names(diag_pos_temp)[1] <- c('terms_v')
    names(ped_terms_temp)[1] <- c('terms_v')
    names(ped_pos_temp)[1] <- c('terms_v')
    names(RRM_nestedF_terms_temp)[1] <- c('terms_v')
    names(RRM_nestedF_pos_temp)[1] <- c('terms_v')
    names(RRM_diag_terms_temp)[1] <- c('terms_v')
    names(RRM_diag_pos_temp)[1] <- c('terms_v')
    names(RRM_ped_terms_temp)[1] <- c('terms_v')
    names(RRM_ped_pos_temp)[1] <- c('terms_v')
  }
  #############Getting final temp######################
  fixed_pos_f <- fixed_pos_temp
  fixed_terms_f <- fixed_terms_temp
  diag_pos_f <- diag_pos_temp
  diag_terms_f <- diag_terms_temp
  ped_pos_f <- ped_pos_temp
  ped_terms_f <- ped_terms_temp
  RRM_nestedF_pos_f <- RRM_nestedF_pos_temp
  RRM_nestedF_terms_f <- RRM_nestedF_terms_temp
  RRM_diag_pos_f <- RRM_diag_pos_temp
  RRM_diag_terms_f <- RRM_diag_terms_temp
  RRM_ped_pos_f <- RRM_ped_pos_temp
  RRM_ped_terms_f <- RRM_ped_terms_temp

  ####Finding random group in renf90.par

  pos_random <- readLines(parms_card) %>%grep(" RANDOM_GROUP", .)+1

  VC_ <- list()
  for (r in 1:length(pos_random)) {
    VC_[[r]] <-unlist(strsplit(readLines(parms_card)[pos_random[r]], split=" ", fixed=TRUE)) %>%
      as.numeric() %>% na.exclude() %>% c()
  }

  diag_terms <- unique(unlist(diag_terms_f)[unlist(diag_terms_f)!=0])
  ped_terms <- unique(unlist(ped_terms_f)[unlist(ped_terms_f)!=0])
  RRM_ped_terms <- unique(unlist(RRM_ped_terms_f)[unlist(RRM_ped_terms_f)!=0])
  RRM_diag_terms <- unique(unlist(RRM_diag_terms_f)[unlist(RRM_diag_terms_f)!=0])

  if(length(RRM_ped_terms)>0){ped_terms <- NULL}

  VC_names <- c(diag_terms,ped_terms,RRM_ped_terms,RRM_diag_terms)
  ### Find position of Genetic variance components inside postmean

  pos_VC <- read_lines(logname) %>%
    str_detect(.,"effect") %>% which(TRUE) %>% +1

  vc_matrix <- list()
  for (vc in 1:length(VC_)) {
    #name_vc <- paste0(VC_names[vc],".temp")
    vc_matrix [[vc]]<- read_lines(logname)[seq(pos_VC[vc],pos_VC[vc]+length(VC_[[vc]])*number_traits-1)] %>%
      strsplit(.,split = " ", fixed=TRUE) %>% unlist() %>%
      as.numeric() %>% na.exclude() %>%
      matrix(ncol = length(VC_[[vc]])*number_traits,nrow = length(VC_[[vc]])*number_traits,byrow = T)
    #write.table(vc_matrix,file = name_vc,row.names = F,col.names = F)
  }
  names(vc_matrix) <- VC_names
  names(VC_)<- VC_names
  ### Residual variance components
  if(method=="blupf90+"){res_name="Res"}else{
    res_name="R"
  }

  pos_res <- read_lines(logname) %>%
    str_detect(.,res_name) %>% which(TRUE) %>% +1

  res_out <- read_lines(logname)[seq(pos_res,number_traits+pos_res-1)] %>%
    strsplit(.,split = " ", fixed=TRUE) %>% unlist() %>%
    as.numeric() %>% na.exclude() %>%
    matrix(ncol = number_traits,nrow = number_traits)
  ##############
  #OUTPUT
  ##############
  if(is.null(ped_terms)==F){
    VCA_out <- vc_matrix[ped_terms] %>% unlist() %>%
      matrix(ncol = number_traits,nrow = number_traits) %>% list(GENETIC=.)
  }else{VCA_out <- vc_matrix[ped_terms] %>% unlist()%>% list(GENETIC=.)}

  if(is.null(ped_terms)==F&length(diag_terms)>0){
  VCD_out <- vc_matrix[diag_terms] %>% unlist()%>%
    matrix(ncol = number_traits,nrow = number_traits)%>% list(RANDOM=.)
  }else{VCD_out <- vc_matrix[diag_terms] %>% unlist()%>% list(RANDOM=.)}


  residual_out <- res_out%>% unlist %>% list(RESIDUAL=.)

  if(length(RRM_diag_terms)>0){
  RRM_random <-  vc_matrix[RRM_diag_terms] %>% unlist  %>% matrix(.,ncol = length(VC_[[RRM_diag_terms]])*number_traits) %>%list(RRM_RANDOM=.)
  }else{RRM_random <- list(RRM_RANDOM=NULL)}
  if(length(RRM_ped_terms)>0){
  RRM_genetic <- vc_matrix[RRM_ped_terms] %>% unlist  %>%
    matrix(.,ncol = length(VC_[[RRM_ped_terms]])*number_traits) %>%  list(RRM_GENETIC=.)
  }else{RRM_genetic <- list(RRM_genetic=NULL)}
  Variance_Components_model <- c(VCD_out,VCA_out,RRM_genetic,RRM_random,residual_out)
  return(Variance_Components_model)
}

