#' Save data
#'
#' @name saveData
#'
#' @details Basic utilities for saving data in \code{phyloseq}.
#'
#' \code{saveAbundanceData} save the otu_table.
#'
#' \code{saveTaxonomyData} save the taxa_table.
#'
#' \code{saveSampleData} save the sample_data.
#'
#' @param x \code{phyloseq} object
#'
#' @param path_file path to save the file along with name.
#'                  Default is NULL. E.g. "data/output/otu_table.csv"
#'
#' @param file_type One of "delim","csv","csv2","tsv","excel_csv","excel_csv2"
#'                  supported by \code{readr} package
#'
#' @param ... Arguments to pass to readr functions
#'
#' @examples
#' \dontrun{
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' saveAbundanceData(FuentesIliGutData, path_file = "otu.tsv", file_type = "tsv")
#' saveTaxonomyData(FuentesIliGutData, path_file = "taxa.tsv", file_type = "tsv")
#' saveSampleData(FuentesIliGutData, path_file = "sample.tsv", file_type = "tsv")
#' }
#' @return Store file locally specified by user.
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/microsud/biomeUtils}
#'
#'
NULL

#' @rdname saveData
#' @aliases saveAbundanceData
#' @export
saveAbundanceData <- function(x,
                              path_file = NULL,
                              file_type = c("delim",
                                            "csv",
                                            "csv2",
                                            "tsv",
                                            "excel_csv",
                                            "excel_csv2"), ...) {

  if(is.null(path_file) || is.na(path_file)){
    stop("Provide a path and/or file name")
  }

  tib.ab <- getAbundanceTibble(x)

  switch(file_type,
         delim = readr::write_delim(tib.ab, file = path_file, ...),
         csv = readr::write_csv(tib.ab, file = path_file, ...) ,
         csv2 = readr::write_csv2(tib.ab, file = path_file, ...),
         tsv = readr::write_tsv(tib.ab, file = path_file, ...),
         excel_csv = readr::write_excel_csv(tib.ab, file = path_file, ...),
         excel_csv2 = readr::write_excel_csv2(tib.ab, file = path_file, ...)
  )
}



# Taxonomy
#' @rdname saveData
#' @aliases saveTaxonomyData
#' @export
saveTaxonomyData <- function(x,
                             path_file = NULL,
                             file_type = c("delim",
                                           "csv",
                                           "csv2",
                                           "tsv",
                                           "excel_csv",
                                           "excel_csv2"), ...) {

  if(is.null(path_file) || is.na(path_file)){
    stop("Provide a path and/or file name")
  }

  tib.tx <- getTaxaTibble(x)

  switch(file_type,
         delim = readr::write_delim(tib.tx, file = path_file, ...),
         csv = readr::write_csv(tib.tx, file = path_file, ...) ,
         csv2 = readr::write_csv2(tib.tx, file = path_file, ...),
         tsv = readr::write_tsv(tib.tx, file = path_file, ...),
         excel_csv = readr::write_excel_csv(tib.tx, file = path_file, ...),
         excel_csv2 = readr::write_excel_csv2(tib.tx, file = path_file, ...)
  )
}

# SampleData
#' @rdname saveData
#' @aliases saveSampleData
#' @export
saveSampleData <- function(x,
                           path_file = NULL,
                           file_type = c("delim",
                                         "csv",
                                         "csv2",
                                         "tsv",
                                         "excel_csv",
                                         "excel_csv2"), ...) {

  if(is.null(path_file) || is.na(path_file)){
    stop("Provide a path and/or file name")
  }

  tib.sm <- getSampleTibble(x)

  switch(file_type,
         delim = readr::write_delim(tib.sm, file = path_file, ...),
         csv = readr::write_csv(tib.sm, file = path_file, ...) ,
         csv2 = readr::write_csv2(tib.sm, file = path_file, ...),
         tsv = readr::write_tsv(tib.sm, file = path_file, ...),
         excel_csv = readr::write_excel_csv(tib.sm, file = path_file, ...),
         excel_csv2 = readr::write_excel_csv2(tib.sm, file = path_file, ...)
  )
}
