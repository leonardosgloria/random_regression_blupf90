my_dummy <- function (.data, select_columns = NULL, remove_first_dummy = FALSE, 
          remove_most_frequent_dummy = FALSE, ignore_na = FALSE, split = NULL, 
          remove_selected_columns = FALSE) 
{
  stopifnot(is.null(select_columns) || is.character(select_columns), 
            select_columns != "", is.logical(remove_first_dummy), 
            length(remove_first_dummy) == 1, is.logical(remove_selected_columns))
  if (remove_first_dummy == TRUE & remove_most_frequent_dummy == 
      TRUE) {
    stop("Select either 'remove_first_dummy' or 'remove_most_frequent_dummy'\n         to proceed.")
  }
  if (is.vector(.data)) {
    .data <- data.frame(.data = .data, stringsAsFactors = FALSE)
  }
  if (!is.null(select_columns)) {
    char_cols <- select_columns
    cols_not_in_data <- char_cols[!char_cols %in% names(.data)]
    char_cols <- char_cols[!char_cols %in% cols_not_in_data]
    if (length(char_cols) == 0) {
      stop("select_columns is/are not in data. Please check data and spelling.")
    }
  }
  else if (ncol(.data) == 1) {
    char_cols <- names(.data)
  }
  else {
    char_cols <- sapply(.data, class)
    char_cols <- char_cols[char_cols %in% c("factor", "character")]
    char_cols <- names(char_cols)
  }
  if (length(char_cols) == 0 && is.null(select_columns)) {
    stop(paste0("No character or factor columns found. ", 
                "Please use select_columns to choose columns."))
  }
  if (!is.null(select_columns) && length(cols_not_in_data) > 
      0) {
    warning(paste0("NOTE: The following select_columns input(s) ", 
                   "is not a column in data.\n"), paste0(names(cols_not_in_data), 
                                                         "\t"))
  }
  for (col_name in char_cols) {
    if (is.factor(.data[[col_name]])) {
      unique_vals <- levels(.data[[col_name]])
      if (any(is.na(.data[[col_name]]))) {
        unique_vals <- c(unique_vals, NA)
      }
    }
    else {
      unique_vals <- unique(.data[[col_name]])
      unique_vals <- stringr::str_sort(unique_vals, na_last = TRUE, 
                                       locale = "en_US", numeric = TRUE)
    }
    unique_vals <- as.character(unique_vals)
    if (!is.null(split)) {
      unique_vals <- unique(trimws(unlist(strsplit(unique_vals, 
                                                   split = split))))
    }
    if (ignore_na) {
      unique_vals <- unique_vals[!is.na(unique_vals)]
    }
    if (remove_most_frequent_dummy) {
      vals <- as.character(.data[[col_name]])
      vals <- data.frame(sort(table(vals), decreasing = TRUE), 
                         stringsAsFactors = FALSE)
      top_vals <- vals[vals$Freq %in% max(vals$Freq), ]
      other_vals <- vals$vals[!vals$Freq %in% max(vals$Freq)]
      other_vals <- as.character(other_vals)
      top_vals <- top_vals[stringr::str_order(top_vals$vals, 
                                              na_last = TRUE, locale = "en_US", numeric = TRUE), 
      ]
      if (nrow(top_vals) == 1) {
        top_vals <- NULL
      }
      else {
        top_vals <- as.character(top_vals$vals[2:nrow(top_vals)])
      }
      unique_vals <- c(top_vals, other_vals)
      unique_vals <- stringr::str_sort(unique_vals, na_last = TRUE, 
                                       locale = "en_US", numeric = TRUE)
    }
    if (remove_first_dummy) {
      unique_vals <- unique_vals[-1]
    }
    data.table::alloc.col(.data, ncol(.data) + length(unique_vals))
    .data[, paste0(col_name, "_", unique_vals)] <- 0L
    for (unique_value in unique_vals) {
      data.table::set(.data, i = which(data.table::chmatch(as.character(.data[[col_name]]), 
                                                           unique_value, nomatch = 0) == 1L), j = paste0(col_name, 
                                                                                                         "_", unique_value), value = 1L)
      if (!is.na(unique_value)) {
        data.table::set(.data, i = which(is.na(.data[[col_name]])), 
                        j = paste0(col_name, "_", unique_value), value = NA)
      }
      if (!is.null(split)) {
        max_split_length <- max(sapply(strsplit(as.character(.data[[col_name]]), 
                                                split = split), length))
        for (split_length in 1:max_split_length) {
          data.table::set(.data, i = which(data.table::chmatch(as.character(trimws(sapply(strsplit(as.character(.data[[col_name]]), 
                                                                                                   split = split), `[`, split_length))), unique_value, 
                                                               nomatch = 0) == 1L), j = paste0(col_name, 
                                                                                               "_", unique_value), value = 1L)
        }
        if (is.na(unique_value)) {
          .data[[paste0(col_name, "_", unique_value)]][which(!is.na(.data[[col_name]]))] <- 0
        }
      }
    }
  }
  if (remove_selected_columns) {
    .data <- .data[-which(names(.data) %in% char_cols)]
  }
  
  return(.data)
}
