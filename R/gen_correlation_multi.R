
RRM_cor <- function(logname,vc_names,parms_card,DAPmin,DAPmax){
  #vc_names <- c("add")
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
  #vc_names <- c("add")
  comb_trait <- gtools::combinations(number_traits,2,c(1:number_traits))

  vc_names <- vc_names %>% gsub("RESIDUAL",NA,.) %>% na.exclude()

  for (vc in 1:length(vc_names)) {
    name_vc <- paste0(vc_names[vc])
    vc_matrix <- read_lines(logname)[seq(pos_VC[vc],pos_VC[vc]+length(VC_[[vc]])*number_traits-1)] %>%
      strsplit(.,split = " ", fixed=TRUE) %>% unlist() %>%
      as.numeric() %>% na.exclude() %>%
      matrix(ncol = length(VC_[[vc]])*number_traits ,nrow = length(VC_[[vc]])*number_traits)
    write.table(vc_matrix,file = name_vc,row.names = F,col.names = F)

  }

  add = as.matrix(read.table(vc_names[1], header = F))   # varcomp matrix coeff additive genetic


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

  fi = as.matrix(read.table("fi.txt", header = F))              # desired matrix, Legendre pol for x (dim) variable

  for (add_vc_t in 1:number_traits) {
    block_max <- add_vc_t*length(VC_[[vc]])
    block_min <- add_vc_t*length(VC_[[vc]])-length(VC_[[vc]])+1
    add_vc_trait <- vc_order_row[block_min:block_max,block_min:block_max]
    ifelse(add_vc_t==1,vc_trait <- list(add_vc_trait),vc_trait[add_vc_t] <- list(add_vc_trait))
  }

  for (add_cov_t in 1:nrow(comb_trait)) {

    block1_max <- add_cov_t*length(VC_[[vc]])
    block1_min <- add_cov_t*length(VC_[[vc]])-length(VC_[[vc]])+1
    block2_max <- (add_cov_t+1)*length(VC_[[vc]])
    block2_min <- (add_cov_t+1)*length(VC_[[vc]])-length(VC_[[vc]])+1

    add_cov_trait <- vc_order_row[block1_min:block1_max,block2_min:block2_max]
    ifelse(add_cov_t==1,cov_trait <- list(add_cov_trait),cov_trait[add_cov_t] <- list(add_cov_trait))
  }

  for (v_fi in 1:number_traits) {
    vfi1 <- fi%*%vc_trait[[v_fi]]%*%t(fi)%>% diag()
    ifelse(v_fi==1,vfi <- vfi1,vfi <- cbind(vfi,vfi1))
  }

  for (cov_fi in 1:nrow(comb_trait)) {
    covfi1 <- fi%*%cov_trait[[cov_fi]]%*%t(fi)%>% diag()
    ifelse(cov_fi==1,covfi <- covfi1,covfi <- cbind(covfi,covfi1))
  }
  covfi <- as.matrix(covfi)



  for (n_comb in 1:nrow(comb_trait)) {

    vf_data <- data.frame(name_trait1=vfi[,comb_trait[n_comb,1]],name_trait2=vfi[,comb_trait[n_comb,2]])

    covfi_data <- covfi[,n_comb]

    corg12_ <- covfi_data/sqrt(vf_data[,1]*vf_data[,2])

    ifelse(n_comb==1,corg12 <- corg12_,corg12 <- cbind(corg12,corg12_))

    ifelse(n_comb==1,cor_name <- paste0(comb_trait[n_comb,],collapse = "_") %>% paste0("cor_",.),
           cor_name[n_comb] <- paste0(comb_trait[n_comb,],collapse = "_") %>% paste0("cor_",.))
  }

  DAP <- seq(DAPmin,DAPmax,by=1)
  corg12_ <- corg12_ %>% data.frame(correlation=.,DAP=DAP)
  write.table(corg12_,"Gen_correlation.txt",row.names = F,quote = F)
  return(corg12_)

}
