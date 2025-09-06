BLUPF90_genotype <- function(file_genotype,ID_Colunm,BLUPF90_format=T){
  if(class(file_genotype)=="character"){
    if(BLUPF90_format==F){
    ID_subject <- data.table::fread(file_genotype,data.table = T,select = ID_Colunm)
    genotype_blupf90 <- data.table::fread(file_genotype,data.table = T,h=F)
    genotype_blupf90[,ID_Colunm]=NULL
    genotype_blupf90<- genotype_blupf90 %>% data.table::as.data.table() %>%
      data.table::setnafill(., fill=5) %>%
      apply(.,1,paste,collapse='') %>% data.frame(ID_subject,SNP=.)
    blup.write(x=genotype_blupf90,file = "genotype_BLUPF90.txt",rownames = F,colnames = F)
    }else{
    genotype_blupf90 <- data.table::fread(file_genotype,data.table = F,h=F)
    blup.write(x=genotype_blupf90,file = "genotype_BLUPF90.txt",rownames = F,colnames = F)
    }
  }else{
    if(BLUPF90_format==F){
    ID_subject <- file_genotype[,ID_Colunm]
    file_genotype[,ID_Colunm]=NULL
    genotype_blupf90<- file_genotype %>%data.table::setnafill(., fill=5) %>%
      apply(.,1,paste,collapse='') %>% data.frame(ID_subject,SNP=.)
    blup.write(x=genotype_blupf90,file = "genotype_BLUPF90.txt",rownames = F,colnames = F)
    }else{
      blup.write(x=data.frame(file_genotype),file = "genotype_BLUPF90.txt",rownames = F,colnames = F)
        }
      }
}
