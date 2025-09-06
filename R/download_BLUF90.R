#' download_BLUPF90
#'
#' Description:
#' This function downloads all the BLUPF90 software files from their official repository
#' (http://nce.ads.uga.edu/html/projects/programs/) and saves them to the specified
#' destination folder. The BLUPF90 software is used for statistical analysis
#' and genetic evaluation of animal and plant breeding data.
#' By default, the function will save the software files in the R user folder,
#' but you can provide a different destination folder using the dest_folder parameter.
#' If the update parameter is set to TRUE, the function will replace any existing
#' files in the destination folder with the latest versions from the repository.
#' If update is set to FALSE, the function will skip the download for files that already exist in the destination folder.
#'
#'
#'
#' Usage:
#' \code{download_BLUPF90(dest_folder = NULL)}
#'
#' Arguments:
#' @param dest_folder (optional) A character string specifying the destination folder where the BLUPF90 files will be saved.
#' @param update (optional): Specifies whether to update the existing BLUPF90 software files if they already exist in the destination folder. If set to TRUE, the function will download and replace any existing files. If set to FALSE, the function will skip the download if the files already exist. The default value is FALSE.
#'
#' @return NULL
#' @export
#' @details
#' The \code{download_BLUPF90} function downloads BLUPF90 software files from the appropriate URL based on the operating system. The function identifies the operating system (Linux, Mac_OSX, or Windows) and constructs the URL accordingly. It retrieves the list of available BLUPF90 files from the URL and compares them to the files in the local destination folder. It then downloads the missing files from the URL to the specified destination folder.
#'
#'
#' @examples
#' \code{
#' # Download BLUPF90 files to the default destination folder
#' download_BLUPF90()
#'
#' # Download BLUPF90 files to a specific destination folder
#' download_BLUPF90(dest_folder = "~/blupf90_files")
#' }
#'
#' Note:
#' The function requires RCurl and httr package to be installed.
#'
#' Please ensure that you have proper permissions to access the destination folder and
#' that you have an active internet connection during the execution of this function.
#'
download_BLUPF90 <- function(dest_folder = NULL, update = FALSE) {
  ## 1. Prepare dest_folder
  if (is.null(dest_folder)) {
    dest_folder <- file.path(.libPaths()[1], "blupf90")
  }
  if (!dir.exists(dest_folder)) {
    dir.create(dest_folder, recursive = TRUE, showWarnings = FALSE)
  }
  
  ## 2. Detect OS
  sys  <- Sys.info()[["sysname"]]
  os   <- switch(sys,
                 Linux   = "Linux",
                 Darwin  = "Mac_OSX",
                 Windows = "Windows",
                 stop("Unsupported platform: ", sys)
  )
  
  ## 3. Helpers
  has_cmd <- function(x) nzchar(Sys.which(x))
  base_url <- paste0(
    "https://nce.ads.uga.edu/html/projects/programs/",
    os, "/64bit/"
  )
  
  ## 4. Fetch HTML listing (skip SSL checks)
  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp), add = TRUE)
  
  fetch_args <- function(url, out) {
    if (os == "Windows" && !has_cmd("curl")) {
      # PowerShell fallback
      cmd <- sprintf(
        "Invoke-WebRequest -Uri %s -OutFile %s -UseBasicParsing -SkipCertificateCheck",
        shQuote(url), shQuote(out)
      )
      system2("powershell", c("-Command", cmd))
    } else if (has_cmd("curl")) {
      system2("curl", c("-k", "-fsSL", url, "-o", shQuote(out)))
    } else if (os != "Windows" && has_cmd("wget")) {
      system2("wget", c("--no-check-certificate", "-qO", shQuote(out), url))
    } else {
      stop("Need 'curl' or 'wget' (or on Windows, PowerShell).")
    }
  }
  
  fetch_args(base_url, tmp)
  html <- readLines(tmp, warn = FALSE)
  
  ## 5. Extract all f90-related filenames
  if (os == "Windows") {
    lines <- grep('href="[^"]+\\.exe"', html, value = TRUE)
    files <- unique(gsub('.*href="([^"]+\\.exe)".*', "\\1", lines))
  } else {
    lines <- grep('href="[^"]*f90[^"]*"', html, value = TRUE)
    files <- unique(gsub('.*href="([^"]+)".*', "\\1", lines))
    # drop parent-directory, any trailing slash
    files <- files[!files %in% c("../", "") & grepl("f90", files)]
  }
  
  if (!length(files)) {
    stop("No BLUPF90 files found at ", base_url)
  }
  
  ## 6. Download each one
  out_paths <- character()
  for (fname in files) {
    url  <- paste0(base_url, fname)
    dest <- file.path(dest_folder, fname)
    if (file.exists(dest) && !update) next
    
    message("Downloading ", fname)
    fetch_args(url, dest)
    if (os != "Windows") Sys.chmod(dest, "0755")
    out_paths <- c(out_paths, dest)
  }
  
  invisible(out_paths)
}