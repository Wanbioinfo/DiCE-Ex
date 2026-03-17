#' Read a File of Any Supported Format into R
#'
#' `read_any()` is a lightweight wrapper that automatically detects a file's
#' extension and loads the data using the appropriate reader. It supports common
#' formats used in DiCE input workflows, including RDS, Excel files, CSV, and
#' tab-delimited text files (.tsv, .txt).
#'
#'
#' @param path Character string giving the file path.
#'
#' @return A data frame, tibble, matrix, or object read from file. 
#' 
#' @details
#' Supported extensions (case-insensitive):
#' \itemize{
#'   \item \code{"rds"}
#'   \item \code{"xls"}, \code{"xlsx"}
#'   \item \code{"csv"}
#'   \item \code{"tsv"}, \code{"txt"}
#' }
#'
#' Unsupported extensions trigger:  
#' \code{stop("Unsupported file type: <ext>")}.
#' @noRd
read_any <- function(path) {
  ext <- tools::file_ext(path)
  
  ext <- tolower(ext)
  
  if (ext == "rds") {
    return(readRDS(path))
  }
  
  if (ext %in% c("xls", "xlsx")) {
    return(read_excel(path))
  }
  
  if (ext %in% c("csv")) {
    return(read.csv(path, stringsAsFactors = FALSE))
  }
  
  if (ext %in% c("tsv", "txt")) {
    return(fread(path, sep = "\t", data.table = FALSE))
  }
  
  stop(paste("Unsupported file type:", ext))
}
