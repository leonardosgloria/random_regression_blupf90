h2_multi <- function(datarenum,formula,logname,parms_card,
                     method,RRM_option){

fi = as.matrix(read.table("fi.txt", header = F))              # desired matrix, Legendre pol for x (dim) variable
DAP <- seq(RRM_option$Pmin,RRM_option$Pmax,by=1)
###############################
number_traits <- readLines(parms_card) %>%
    .[grep("NUMBER_OF_TRAITS", .)+1] %>%
    as.numeric()

  pos_random <- readLines(parms_card) %>%grep(" RANDOM_GROUP", .)+1

  VC_ <- list()
  for (r in 1:length(pos_random)) {
    VC_[[r]] <-unlist(strsplit(readLines(parms_card)[pos_random[r]], split=" ", fixed=TRUE)) %>%
      as.numeric() %>% na.exclude() %>% c()
  }

  ### Find position of Genetic variance components inside airemlf90.log
  pos_VC <- read_lines(logname) %>%
    str_detect(.,"Genetic") %>% which(TRUE) %>% +1


  ####################################################################
  ######Obtaining the variances components for the genetic parameters
  ###############################################################co#####
  vc_names <- c("add")
  comb_trait <- gtools::combinations(number_traits,2,c(1:number_traits))


  for (vc in 1:length(vc_names)) {
    name_vc <- paste0(vc_names[vc])
    vc_matrix <- read_lines(logname)[seq(pos_VC[vc],pos_VC[vc]+length(VC_[[vc]])*number_traits-1)] %>%
      strsplit(.,split = " ", fixed=TRUE) %>% unlist() %>%
      as.numeric() %>% na.exclude() %>%
      matrix(ncol = length(VC_[[vc]])*number_traits ,nrow = length(VC_[[vc]])*number_traits)
    write.table(vc_matrix,file = name_vc,row.names = F,col.names = F)

  }


############
#Additive RR
############
#add = as.matrix(read.table("RRM_GENETIC", header = F))   # varcomp matrix coeff additive genetic
add <-
  pipeline_extract_vc_multi(datarenum=datarenum,formula,
                            logname=logname,parms_card=parms_card,method,RRM_option=RRM_option) %>%
  .$RRM_GENETIC


  for (i in 1:number_traits) {
    trait_vc <- add[,seq(i,ncol(add),by=2)]
    ifelse(i==1,vc_order_col <- trait_vc,vc_order_col <- cbind(vc_order_col,trait_vc))
  }
  vc_order_col <- t(vc_order_col)
  for (i in 1:number_traits) {
    trait_vc <- vc_order_col[,seq(i,ncol(vc_order_col),by=2)]
    ifelse(i==1,vc_order_row <- trait_vc,vc_order_row <- cbind(vc_order_row,trait_vc))
  }

  vc_order_row <- t(vc_order_row)



  for (add_vc_t in 1:number_traits) {
    block_max <- add_vc_t*length(VC_[[vc]])
    block_min <- add_vc_t*length(VC_[[vc]])-length(VC_[[vc]])+1
    add_vc_trait <- vc_order_row[block_min:block_max,block_min:block_max]
    ifelse(add_vc_t==1,vc_trait <- list(add_vc_trait),vc_trait[add_vc_t] <- list(add_vc_trait))
  }

###############
#nonAdditive RR
###############
#nadd = as.matrix(read.table("RRM_RANDOM", header = F))   # varcomp matrix coeff additive genetic
nadd <-
    pipeline_extract_vc_multi(datarenum=datarenum,formula=formula,
                              logname=logname,parms_card=parms_card,method=method,RRM_option=RRM_option) %>%
    .$RRM_RANDOM

  naddvc_trait <- NULL
if(is.null(nadd)==F){
  for (i in 1:number_traits) {
    trait_vc <- nadd[,seq(i,ncol(nadd),by=2)]
    ifelse(i==1,vc_order_col <- trait_vc,vc_order_col <- cbind(vc_order_col,trait_vc))
  }
  vc_order_col <- t(vc_order_col)
  for (i in 1:number_traits) {
    trait_vc <- vc_order_col[,seq(i,ncol(vc_order_col),by=2)]
    ifelse(i==1,vc_order_row <- trait_vc,vc_order_row <- cbind(vc_order_row,trait_vc))
  }

  vc_order_row <- t(vc_order_row)

  for (nadd_vc_t in 1:number_traits) {
    block_max <- nadd_vc_t*length(VC_[[vc]])
    block_min <- nadd_vc_t*length(VC_[[vc]])-length(VC_[[vc]])+1
    nadd_vc_trait <- vc_order_row[block_min:block_max,block_min:block_max]
    ifelse(nadd_vc_t==1,naddvc_trait <- list(nadd_vc_trait),naddvc_trait[nadd_vc_t] <- list(nadd_vc_trait))
  }
}
#############################################
  gp <-
    pipeline_extract_vc_multi(datarenum=datarenum,formula=formula,
                              logname=logname,parms_card=parms_card,method=method,RRM_option=RRM_option)

  gp$RRM_GENETIC <- NULL
  gp$RRM_RANDOM <- NULL

  gp$RESIDUAL %>% diag() %>% data.frame(trait=seq(1,length(.)),residual=.) %>% split(.$trait) %>%
    .[[1]] %>%
    select(residual) %>% c(.,rep(0,ncol(vc_trait[[1]])*ncol(vc_trait[[1]])-1)) %>%
    unlist() %>% matrix(.,ncol=ncol(vc_trait[[1]]))

  gpl <- list()
  gptrait <- list()
  i=1
  ngpl <- NULL

  for (gpti in 1:number_traits) {
    i=1
    for (gpi in 1:length(gp)) {

      if(length(gp[[gpi]]) == 0) {
        gpl[[gpi]]=NULL
        i=i
      }else{
        gpl[[i]]=
          gp[[gpi]] %>% diag() %>% data.frame(trait=seq(1,length(.)),residual=.) %>% split(.$trait) %>%
          .[[gpti]] %>%
          select(residual) %>% c(.,rep(0,ncol(vc_trait[[1]])*ncol(vc_trait[[1]])-1)) %>%
          unlist() %>% matrix(.,ncol=ncol(vc_trait[[1]]))


        ngpl[i] <- names(gp)[gpi]
        i=i+1
      }
    }
    gptrait[[gpti]] <-gpl
    names(gptrait[[gpti]]) <- ngpl
  }

  gpc <- NULL
  vct <- list()
  for (gpti in 1:number_traits){
    gpc <- NULL
    for (gpi in 1:length(ngpl)){
      gpc1 <- diag(fi%*%gptrait[[gpti]][[gpi]]%*%t(fi))
      gpc <- cbind(gpc,gpc1)

    }
    colnames(gpc) <- ngpl
    gRRM <- fi%*%as.matrix(vc_trait[[gpti]])%*%t(fi) %>%
      diag()

    dRRM <-NULL

    if(is.null(naddvc_trait)==F){
    dRRM <- fi%*%as.matrix(naddvc_trait[[gpti]])%*%t(fi) %>%
      diag()
    }
    # do the same with random RR
    vct[[gpti]] <- cbind(gpc,gRRM,dRRM) %>% data.frame()

    vct[[gpti]] <- vct[[gpti]] %>% mutate(h2=gRRM/apply(., 1, sum),DAP)
  }
#get trait names
trait_terms <- NULL
  for (mod in 1:length(formula)){
    form <- formula[[mod]]

    ter = terms(form)
    labs = labels(ter)

    trait_terms[mod] <-
      ter %>% .[[2]] %>% as.character() %>%
      gsub( "\\|", NA,.) %>% na.exclude() %>%
      c()
  }
names(vct) <- trait_terms
return(vct)
}
