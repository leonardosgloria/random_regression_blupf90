pipeline_renum_multi <- function(datarenum,formula,fields_output=NULL,weights_object=NULL,residual_start,VCA_RRM,VCD_RRM,VCA,VCD,
                           ped_name,PED_DEPTH=3,genotype_file,missing_values,RRM_option,het_res_variance=NULL,
                           fit_option,extra.option.blup,blupf90_folder){
  ##################################


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
    ##########PEDIGREE############################
    if(is.null(ped_name)==F){
      ped_name=ped_name

    }else{ped_name="pedBLUPF90.txt"}
    ############SNP##############################
    if(is.null(genotype_file)==T){
      snp_gen=c("")
    }else{
      snp_gen <- paste("SNP_FILE","\n", genotype_file,"\n")}

    if(is.null(genotype_file)==F){geno_card <- paste("OPTION SNP_file", genotype_file)
    }else{geno_card <-c("")}
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
      fixed_terms_temp <- unique_terms_eff(eff=fixed_terms_temp,eff_n='fixed')

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


  ################
  #END of the LOOP
  ################

  #####Creating Legendre if RRM is true
  if(nrow(RRM_diag_pos_f)>0||nrow(RRM_ped_pos_f>0)){
    #creating RRM based on the RRM column

    #  eval(call("=",paste0(as.name(Timevar),"_renum"),seq(min(datarenum[,Timevar]),
    #                                                      max(datarenum[,Timevar]))))
    Timemin=min(datarenum[,RRM_option$Timevar])
    Timemax=max(datarenum[,RRM_option$Timevar])

    if(is.null(RRM_option$Pmin)==F){Timemin=RRM_option$Pmin}
    if(is.null(RRM_option$Pmax)==F){Timemax=RRM_option$Pmax}
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
    # dir.create(paste0(getwd(),"/poly",RRM_option$poly))
    # # find the files that you want and copy the files to the new folder
    # if(is.null(genotype_file)==F){
    #   file.copy(list.files(getwd(), genotype_file), paste0(getwd(),"/poly",RRM_option$poly))}
    # # set as working directory the new folder
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

    data_final_poly <-
      data_final %>% select(starts_with("poly")) %>% dplyr::mutate_all(~replace_na(.,0))
    data_final_poly$poly0 <- 1

    data_final[,c(seq(ncol(data_final)-RRM_option$poly,ncol(data_final)))] <- data_final_poly

    fi1<- fi1[,-1]

    write.table(fi1,"fi.txt",col.names = F,row.names = F)

    poly_pos <- data.frame(col=colnames(data_final)) %>% data.frame(col=.,pos=rownames(.)) %>%
      dplyr::filter(col%in%grep("poly", names(data_final), value = TRUE)) %>%
      .$pos %>% as.numeric()
  }else{data_final=datarenum}


  #######PED all classes#####################################################
  if(nrow(ped_pos_f)==0){ped_pos_f <-ped_pos_f }else{
    ped_pos_f <-join_row_eff(ped_pos_f) %>% .[which(rowSums(.) > 0),]
  }

  if(nrow(ped_pos_f)==0){get_ped_class <-ped_pos_f }else{
    get_ped_class <-get_terms(ped_pos_f)
  }

  if(nrow(get_ped_class)==0){
    ped_card <- c("")
  }else{
    for (i in 1:nrow(get_ped_class)) {

      fcla <-get_ped_class %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,3]
      ifelse(i==1,class_ped <- ifelse(fcla=='factor','cross alpha','cov'),
             class_ped[i] <- ifelse(fcla=='factor','cross alpha','cov'))
    }

    ################################
    ####COV start values############
    ################################

    if(class(VCA)[1]=="list"){
      cov_add <- matrix(VCA$COV,ncol = ncol(trait_pos_f),
                            nrow = ncol(trait_pos_f))
      diag(cov_add) <- c(rep(VCA$VAR,ncol(trait_pos_f)))
    }else{
      cov_add <-VCA
    }

    ######PED CARD RENUM########################################


    pos_class_ped <- function(x){
      paste("EFFECT\n",toString(x),"\n",
            "RANDOM\n","animal\n","FILE\n",ped_name,"\n",
            snp_gen,
            "(CO)VARIANCES","\n",
            toString(cov_add),"\n") %>%
        gsub("\\,","",.)
    }

    ped_card <-apply(unique(data.frame(ped_pos_f,class_ped)), 1, pos_class_ped)
  }
  #######Fixed all classes#####################################################

  if(nrow(fixed_pos_f)==0){fixed_pos_f <-fixed_pos_f }else{
    fixed_pos_f <-join_row_eff(fixed_pos_f) %>% .[which(rowSums(.) > 0),]
  }

  if(nrow(fixed_pos_f)==0){get_fixed_class <-fixed_pos_f }else{
    get_fixed_class <-get_terms(fixed_pos_f)
  }
  if(nrow(get_fixed_class)==0){
    fixed_card <- c("")
  }else{
    for (i in 1:nrow(get_fixed_class)) {

      fcla <-get_fixed_class %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,3]
      ifelse(i==1,class_fixed <- ifelse(fcla=='factor','cross alpha','cov'),
             class_fixed[i] <- ifelse(fcla=='factor','cross alpha','cov'))
    }

    ######FIXED CARD RENUM########################################
    pos_class_fixed <- function(x){
      paste("EFFECT\n",toString(x),"\n") %>%
        gsub("\\,","",.)
    }

    fixed_card <-  apply(unique(data.frame(fixed_pos_f,class_fixed)), 1, pos_class_fixed)
  }
  #######Diagonal all classes#####################################################

  if(nrow(diag_pos_f)==0){diag_pos_f <-diag_pos_f }else{
    diag_pos_f <-join_row_eff(diag_pos_f) %>% .[which(rowSums(.) > 0),]
  }

  if(nrow(diag_pos_f)==0){get_diag_class <-diag_pos_f }else{
    get_diag_class <-get_terms(diag_pos_f)
  }

  if(nrow(get_diag_class)==0){
    diag_card <- c("")
  }else{

    for (i in 1:nrow(get_diag_class)) {

      fcla <-data.frame(get_diag_class) %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,3]
      ifelse(i==1,class_diag <- ifelse(fcla=='factor','cross alpha','cov'),
             class_diag[i] <- ifelse(fcla=='factor','cross alpha','cov'))
    }
    ################################
    ####COV start values############
    ################################
    if(class(VCD)[1]=="list"){
      cov_diag <- matrix(VCD$COV,ncol = ncol(trait_pos_f),
                        nrow = ncol(trait_pos_f))
      diag(cov_diag) <- c(rep(VCD$VAR,ncol(trait_pos_f)))
    }else{
      cov_diag <-VCD
    }
    ######diagonal CARD RENUM########################################
    pos_class_diag <- function(x){
      paste("EFFECT\n",toString(x),"\n",
            "RANDOM\n","diagonal\n",
            "(CO)VARIANCES","\n",
            toString(cov_diag),"\n") %>%
        gsub("\\,","",.)
    }

    diag_card <- apply(unique(data.frame(diag_pos_f,class_diag)), 1, pos_class_diag)
  }


  #######RRM_nestedF all classes#####################################################

  if(nrow(RRM_nestedF_pos_f)==0){RRM_nestedF_pos_f <-RRM_nestedF_pos_f }else{
    RRM_nestedF_pos_f <-join_row_eff(RRM_nestedF_pos_f) %>% .[which(rowSums(.) > 0),]
  }

  if(nrow(RRM_nestedF_pos_f)==0){get_RRM_nestedF_class <-RRM_nestedF_pos_f }else{
    get_RRM_nestedF_class <-get_terms(RRM_nestedF_pos_f)
  }

  if(nrow(get_RRM_nestedF_class)==0){
    RRM_nestedF_card <- c("")
  }else{

    for (i in 1:nrow(get_RRM_nestedF_class)) {

      fcla <-data.frame(get_RRM_nestedF_class) %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,3]
      ifelse(i==1,class_RRM_nestedF <- ifelse(fcla=='factor','cross alpha','cov'),
             class_RRM_nestedF[i] <- ifelse(fcla=='factor','cross alpha','cov'))
    }

    if(nrow(RRM_nestedF_pos_f)==0){
      return(c(""))
    }else{
      RRM_nestedF_card <- c()
      for (eff in 1:nrow(RRM_nestedF_pos_f)) {

        for (nv in 1:length(poly_pos)) {
          ifelse(length(class_RRM_nestedF)==0,RRM_nestedF_card1 <- c(""),
                 RRM_nestedF_card1 <-paste(c("EFFECT\n"),
                                           toString(cbind(replicate(ncol(RRM_nestedF_pos_f),poly_pos))[nv,]), "cov","\n",
                                           c("NESTED\n"),
                                           toString(data.frame(RRM_nestedF_pos_f,class_RRM_nestedF)[eff,]),"\n"))

          ifelse(is.na(class_RRM_nestedF)==TRUE,RRM_nestedF_card <- c(""),
                 RRM_nestedF_card <-rbind(RRM_nestedF_card,RRM_nestedF_card1))
        }
        RRM_nestedF_card <- RRM_nestedF_card %>% gsub("\\,","",.)

      }
    }
  }
    #######RRM Diagonal all classes#####################################################
    RRM_diag_card_func <- function(x){
      paste(c("EFFECT\n"),
            toString(x),"\n",
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

    if(nrow(RRM_diag_pos_f)==0){RRM_diag_pos_f <-RRM_diag_pos_f }else{
      RRM_diag_pos_f <-join_row_eff(RRM_diag_pos_f) %>% .[which(rowSums(.) > 0),]
    }

    if(nrow(RRM_diag_pos_f)==0){get_RRM_diag_class <-RRM_diag_pos_f }else{
      get_RRM_diag_class <-get_terms(RRM_diag_pos_f)
    }

    if(nrow(get_RRM_diag_class)==0){
      RRM_diag_card <- c("")
    }else{

      for (i in 1:nrow(get_RRM_diag_class)) {

        fcla <-data.frame(get_RRM_diag_class) %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,3]
        ifelse(i==1,class_RRM_diag <- ifelse(fcla=='factor','cross alpha','cov'),
               class_RRM_diag[i] <- ifelse(fcla=='factor','cross alpha','cov'))
      }
      ################################
      ####COV start values############
      ################################
      if(class(VCD_RRM)[1]=="list"){
        cov_rrm_diag <- matrix(VCD_RRM$COV,ncol = length(poly_pos)*ntraits,
                              nrow = length(poly_pos)*ntraits)
        diag(cov_rrm_diag) <- c(rep(VCD_RRM$VAR,length(poly_pos)*ntraits))

        ####COV start values making zero who is not RRM############

        rrmvec <- RRM_diag_pos_f[1,] %>% gsub("[1-9]","1",.) %>% as.numeric()
        rrm_v <-
          c(rep(1,ntraits),rep(rrmvec,length(poly_pos)-1))
        cov_rrm_diag <- t(rrm_v*cov_rrm_diag)*rrm_v
      }else{
        cov_rrm_diag <- VCD_RRM

      }
      ######RRM diagonal CARD RENUM########################################
      for (rrmeff in 1:nrow(RRM_diag_pos_f)) {
        RRM_diag_pos_f[rrmeff,] <-RRM_diag_pos_f[rrmeff,][RRM_diag_pos_f[rrmeff,]!=0] %>%
          unique() %>% rep(.,ncol(RRM_diag_pos_f))
      }
      RRM_diag_card <- apply(unique(data.frame(RRM_diag_pos_f,class_RRM_diag)), 1, RRM_diag_card_func)
    }
    #######RRM ped all classes#####################################################
    RRM_ped_card_func <- function(x){
      paste("EFFECT\n",
            toString(x),"\n",
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
      ) %>% gsub( ", ", " ",.)

    }

    if(nrow(RRM_ped_pos_f)==0){RRM_ped_pos_f <-RRM_ped_pos_f }else{
      RRM_ped_pos_f <-join_row_eff(RRM_ped_pos_f) %>% .[which(rowSums(.) > 0),]
    }

    if(nrow(RRM_ped_pos_f)==0){get_RRM_ped_class <-RRM_ped_pos_f }else{
      get_RRM_ped_class <-get_terms(RRM_ped_pos_f)
    }

    if(is.null(RRM_option)==T){
      RRM_ped_card <- c("")
    }else{

      for (i in 1:nrow(get_RRM_ped_class)) {

        fcla <-data.frame(get_RRM_ped_class) %>% dplyr::left_join(variables_class,by="terms_v") %>% .[i,3]
        ifelse(i==1,class_RRM_ped <- ifelse(fcla=='factor','cross alpha','cov'),
               class_RRM_ped[i] <- ifelse(fcla=='factor','cross alpha','cov'))
      }
      ################################
      ####COV start values############
      ################################
      if(class(VCA_RRM)[1]=="list"){
        cov_rrm_add <- matrix(VCA_RRM$COV,ncol = length(poly_pos)*ntraits,
                              nrow = length(poly_pos)*ntraits)
        diag(cov_rrm_add) <- c(rep(VCA_RRM$VAR,length(poly_pos)*ntraits))

        ####COV start values making zero who is not RRM############

        rrmvec <- RRM_ped_pos_f[1,] %>% gsub("[1-9]","1",.) %>% as.numeric()
        rrm_v <-
          c(rep(1,ntraits),rep(rrmvec,length(poly_pos)-1))
        cov_rrm_add <- t(rrm_v*cov_rrm_add)*rrm_v
      }else{
        cov_rrm_add <- VCA_RRM

      }
      ######RRM ped CARD RENUM########################################
      for (rrmeff in 1:nrow(RRM_ped_pos_f)) {
        RRM_ped_pos_f[rrmeff,] <- RRM_ped_pos_f[rrmeff,][RRM_ped_pos_f[rrmeff,]!=0] %>%
          unique() %>% rep(.,ncol(RRM_ped_pos_f))
      }

      RRM_ped_card <- apply(unique(data.frame(RRM_ped_pos_f,class_RRM_ped)), 1, RRM_ped_card_func)
    }
    ###########################RRM multitrait
    if (is.null(RRM_option)==F){
      ped_card= NULL
    }


  ######################
  #Residual Start matrix
  ######################
  if(class(residual_start)[1]=="list"){
    residual_VC <- matrix(residual_start$COV,ncol = ncol(trait_pos_f),
                        nrow = ncol(trait_pos_f))
    diag(residual_VC) <- c(rep(residual_start$VAR,ncol(trait_pos_f)))
  }else{
    residual_VC <-residual_start
  }
  ###################################################
  #fields_output_pos including  Heterogeneous variance
  ###################################################
  #creating heterogeneous variance vector
  if(is.null(het_res_variance)==F){
    het_var <- factor(data_final[,het_res_variance$het_levels]) %>% levels() %>% data.frame()
    colnames(het_var) <- c(het_res_variance$het_levels)
    het_var <-
      het_var %>% dplyr::mutate(lev=rownames(.)) %>%
      dplyr::mutate(h_var=cut(as.numeric(lev),het_res_variance$n_res_var,labels = F)) %>% dplyr::select(het_res_variance$het_levels,h_var)

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

  fields_output_pos <-c(het_eff_field,fields_output,het_levels_field,subj_field) %>%
    unlist() %>% as.vector()
  #################
  ###RENUM CARD
  #################
  sink(file = "renum.par")
  cat("DATAFILE\n")
  cat("data_renumf90.txt","\n")
  cat("TRAITS\n")
  cat(trait_pos_f,"\n")
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
  system(renum) # you should check "renum.par" file in advance

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

OPT_default <- list(yams=F,solution_mean=T,VCE=T,sol_se=T,Inbreeding=T,alpha_size=30,EM_REML=100,maxrounds=300,alpha_beta=c(0.95,0.05),tunedG=0,
                      conv_crit=1e-10)
#OPTION conv_crit 1e-12
opt_list <- c(OPT_default[names(OPT_default)[!(names(OPT_default) %in% names(fit_option))]],fit_option)

sm <- ifelse(opt_list$solution_mean==T,paste("OPTION solution mean","\n"),"")
yams <- ifelse(opt_list$yams==T,paste("OPTION use_yams","\n"),"")
VCE <- ifelse(opt_list$VCE==T,paste("OPTION method VCE","\n"),"")
sol_se <- ifelse(opt_list$sol_se==T,paste("OPTION sol se","\n"),"")
imb <- ifelse(opt_list$Inbreeding==T,paste("OPTION Inbreeding","\n"),"")
alps <- ifelse(opt_list$alpha_size>0,paste("OPTION alpha_size",opt_list$alpha_size,"\n"),"")
EM_reml <- ifelse(opt_list$EM_REML!=0,paste("OPTION EM-REML",(opt_list$EM_REML),"\n"),"")
max_r <- ifelse(opt_list$maxrounds>0,paste("OPTION maxrounds",opt_list$maxrounds,"\n"),"")
alph_b <- ifelse(opt_list$alpha_beta>0,paste("OPTION AlphaBeta",opt_list$alpha_beta[1],opt_list$alpha_beta[2],"\n"),"")
tuG <- ifelse(opt_list$tunedG!=0,paste("OPTION tunedG",opt_list$tunedG,"\n"),paste("OPTION tunedG",0,"\n"))
conv_crit <- ifelse(opt_list$conv_crit!=1e-10,paste("OPTION conv_crit",opt_list$conv_crit,"\n"),paste("OPTION conv_crit",1e-10,"\n"))

oriID <- "OPTION origID"
opt_blupf90 <- unique(paste0(sm,VCE,yams,sol_se,imb,alps,
                               EM_reml,max_r,alph_b,tuG,conv_crit,oriID))
write(opt_blupf90, file = "renf90.par", append = T)

#Extra OPTION within paremeter file
write(extra.option.blup, file = "renf90.par", append = T)
######
  data.table::fread("renf90.fields",select = c("field","origfield"),data.table = F) %>%
    mutate(names=names(data_final)[.$origfield]) %>%
    write.table(.,"cols_data_renumf90.txt",row.names = F,quote = T)
}
#################
