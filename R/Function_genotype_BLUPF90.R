#Convert txt genotype file to BLUPF90 genotype file
# ID first column
genotype_BLUPF90_ <- function(geno_file){

sink(file = "geno_BLUPF90.py")
cat("
geno = open('snp.txt', 'wt')
id_sub = open('ID.txt', 'wt')") 
cat("
with open(",glue::double_quote(geno_file),") as a_file:
  for line in a_file:
    id_sub.write(line.split()[0]+",paste0(glue::double_quote('\n')),")","\n")
cat("    geno.write(''.join(line.split()[1:])+",paste0(glue::double_quote('\n')),")","\n")
cat("geno.close()
id_sub.close()
")
sink()

system("python3 geno_BLUPF90.py")

system("paste -d ' ' ID.txt snp.txt > output.txt")

geno <- data.table::fread("output.txt",nrows = 1,skip = 1)

IDc <- data.table::fread("output.txt",select = 1)

sSNP <-
  geno %>% data.frame() %>% apply(.,1,str_count) %>% .[2,1]

sID <-
  IDc %>% data.frame() %>% apply(.,1,str_count) %>% max()

aw <- glue::double_quote(paste0("%-",sID,"s %-",sSNP,"s\n"))

awkv <- 
  glue::glue("awk '{printf [aw] ,$1, $2}' output.txt > genotype_BLUPF90.txt",
             sep = "\"",.open = "[",.close = "]")
system(awkv)

file.remove("ID.txt")
file.remove("snp.txt")
file.remove("output.txt")
file.remove("geno_BLUPF90.py")
}
