# .onLoad <- function(libname, pkgname) {
# #  download_BLUPF90(dest_folder = paste0(.libPaths()[1],"/blupf90"))
#   print("download blupf90 software has finished")
# }
.onAttach <- function(libname, pkgname) {
 packageStartupMessage("✔ Package 'blup' loaded. To download BLUPF90 binaries, run: download_BLUPF90()")
}
