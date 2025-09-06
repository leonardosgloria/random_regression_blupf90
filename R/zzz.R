## R/zzz.R

.onAttach <- function(libname, pkgname) {
  # Destination folder where we expect the binaries
  dest <- file.path(.libPaths()[1], "blupf90")
  
  # Determine OS and the list of expected filenames
  sys  <- Sys.info()[["sysname"]]
  os   <- switch(sys,
                 Linux   = "Linux",
                 Darwin  = "Mac_OSX",
                 Windows = "Windows",
                 stop("Unsupported platform: ", sys)
  )
  
  expected <- if (os == "Windows") {
    c(
      "blupf90+.exe", "gibbsf90+.exe", "idsolf90.exe", "inbupgf90.exe",
      "postgibbsf90.exe", "postGSf90.exe", "predf90.exe", "predictf90.exe",
      "preGSf90.exe", "qcf90.exe", "renumf90.exe", "seekparentf90.exe",
      "validationf90.exe"
    )
  } else {
    c(
      "blupf90+",  "gibbsf90+",  "idsolf90",   "inbupgf90",
      "postgibbsf90", "postGSf90", "predf90",    "predictf90",
      "preGSf90",  "qcf90",      "renumf90",    "rrmebvf90",
      "seekparentf90","validationf90"
    )
  }
  
  # Check which are missing
  missing <- expected[!file.exists(file.path(dest, expected))]
  
  if (length(missing)) {
    packageStartupMessage(
      "⚠️  BLUPF90 binaries not all found in '", dest, "'.\n",
      "    Missing: ", paste(missing, collapse = ", "), "\n",
      "    Please run: download_BLUPF90()"
    )
  }
}
