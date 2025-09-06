blup.write <- function (x, file = "", quote = FALSE, sep = " ",
          na = "", rownames = FALSE, colnames = TRUE, rowCol = NULL,
          justify = "left", formatInfo = FALSE, quoteInfo = TRUE,
          width = NULL, eol = "\n",
          scientific = TRUE, ...)
{
  dapply <- function(x, FUN, ..., simplify = TRUE) {
    if (is.data.frame(x))
      return(sapply(x, FUN, ..., simplify = simplify))
    else if (is.matrix(x))
      return(apply(x, 2, FUN, ...))
    else stop("x must be a data.frame or a matrix")
  }
  if (!(is.data.frame(x) || is.matrix(x)))
    stop("'x' must be a data.frame or matrix")
  if (length(na) > 1)
    stop("only single value can be defined for 'na'")
  if (!scientific) {
    option.scipen <- getOption("scipen")
    on.exit(function() options(scipen = option.scipen))
    options(scipen = 100)
  }
  if (rownames) {
    x <- as.data.frame(x)
    x <- cbind(rownames(x), x)
    rowColVal <- ifelse(!is.null(rowCol), rowCol, "row")
    colnames(x)[1] <- rowColVal
  }
  colnamesMy <- colnames(x)
  if (length(colnamesMy) == 0)
    colnamesMy <- paste("V", 1:ncol(x), sep = "")
  nRow <- nrow(x)
  nCol <- length(colnamesMy)
  widthNULL <- is.null(width)
  if (!widthNULL && length(width) != nCol) {
    warning("recycling 'width'")
    widthOld <- width
    width <- integer(length = nCol)
    width[] <- widthOld
  }
  retFormat <- data.frame(colname = colnamesMy, nlevels = 0,
                          position = 0, width = 0, digits = 0, exp = 0, stringsAsFactors = FALSE)
  isNum <- dapply(x, is.numeric)
  isNum <- isNum & !(dapply(x, inherits, what = "Date") |
                       dapply(x, inherits, what = "POSIXt"))
  isFac <- dapply(x, is.factor)
  if (any(isFac))
    x[, isFac] <- sapply(x[, isFac, drop = FALSE], as.character)
  tmp <- dapply(x, format.info, ..., simplify = FALSE)
  if (is.matrix(x))
    tmp <- as.data.frame(tmp)
  tmp1 <- sapply(tmp, length)
  tmp <- t(as.data.frame(tmp))
  retFormat$width <- tmp[, 1]
  if (any(isNum)) {
    test <- tmp1 > 1
    if (any(test)) {
      retFormat[test, c("digits", "exp")] <- tmp[test,
                                                 c(2, 3)]
      test2 <- tmp[test, 3] > 0
      if (any(test2))
        retFormat[test, ][test2, "exp"] <- retFormat[test,
        ][test2, "exp"] + 1
    }
  }
  y <- x
  for (i in 1:nCol) {
    if (widthNULL) {
      tmp <- NULL
    }
    else {
      tmp <- width[i]
    }
    test <- is.na(y[, i])
    x2 <- character(length = nRow)
    x2[!test] <- format(y[!test, i], justify = justify,
                        width = tmp, ...)
    x2[test] <- na
    x[, i] <- x2
    tmp2 <- format.info(x2, ...)[1]
    if (tmp2 != retFormat[i, "width"]) {
      retFormat[i, "width"] <- tmp2
      x[, i] <- format(x[, i], justify = ifelse(isNum[i],
                                                "right", justify), width = tmp, ...)
    }
    if (nchar(na) < retFormat[i, "width"]) {
      x[test, i] <- format(na, justify = ifelse(isNum[i],
                                                "right", justify), width = retFormat[i, "width"],
                           ...)
    }
  }
  if (any(!isNum)) {
    retFormat[!isNum, "nlevels"] <- dapply(x[, !isNum, drop = FALSE],
                                           function(z) length(unique(z)))
  }
  if (!widthNULL) {
    test <- retFormat$width > width
    if (any(test)) {
      tmpCol <- paste(colnamesMy[test], collapse = ", ")
      tmpWidth <- paste(width[test], collapse = ", ")
      tmpNeed <- paste(retFormat$width[test], collapse = ", ")
      stop(paste("'width' (", tmpWidth, ") was too small for columns: ",
                 tmpCol, "\n 'width' should be at least (", tmpNeed,
                 ")", sep = ""))
    }
  }
  if (colnames) {
    if (rownames && is.null(rowCol))
      colnamesMy <- colnamesMy[-1]
    data.table::fwrite(t(as.matrix(colnamesMy)), file = file,
                quote = quote, sep = sep, eol = eol, na = na, row.names = FALSE,
                col.names = FALSE)
  }
  data.table::fwrite(x = x, file = file,
              quote = quote, sep = sep, eol = eol, na = na, row.names = FALSE,
              col.names = FALSE)
  if (formatInfo) {
    retFormat$position[1] <- ifelse(quote, ifelse(quoteInfo,
                                                  1, 2), 1)
    if (ifelse(quote, quoteInfo, FALSE))
      retFormat$width <- retFormat$width + 2
    N <- nrow(retFormat)
    if (N > 1) {
      for (i in 2:N) {
        retFormat$position[i] <- retFormat$position[i -
                                                      1] + retFormat$width[i - 1] + nchar(x = sep,
                                                                                          type = "chars") + ifelse(quote, ifelse(quoteInfo,
                                                                                                                                 0, 1), 0)
      }
    }
    if (rownames && is.null(rowCol)) {
      retFormat <- retFormat[-1, ]
      rownames(retFormat) <- 1:(N - 1)
    }
    return(retFormat)
  }
}
