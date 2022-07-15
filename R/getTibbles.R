#' Get tibble
#'
#' @name getTibble
#'
#' @details Convert different \code{phyloseq} slots into tibbles.
#'
#' \code{getAbundanceTibble} gets the otu_table in tibble format.
#'
#' \code{getTaxaTibble} gets the taxa_table in tibble format.
#'
#' \code{getSampleTibble} gets the sample_data in tibble format.
#'
#' @param x \code{phyloseq} object
#'
#' @param column_id A character
#'
#' @param select_rows Rows to return in output.
#'
#' @param select_cols Columns to return in output
#'
#' @examples
#' library(biomeUtils)
#' data("FuentesIliGutData")
#' otu_tib <- getAbundanceTibble(FuentesIliGutData)
#' tax_tib <- getTaxaTibble(FuentesIliGutData)
#' meta_tib <- getSampleTibble(FuentesIliGutData)
#' @return A tibble
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Shetty SA (2020). Utilities for microbiome analytics.
#' \url{https://github.com/RIVM-IIV-Microbiome/biomeUtils}
#'
#'
NULL


#' @rdname getTibble
#' @aliases getAbundanceTibble
#' @importFrom tibble rownames_to_column
#' @importFrom microbiome abundances
#' @importFrom dplyr as_tibble %>%
#' @importFrom rlang .data sym
#' @export
getAbundanceTibble <- function(x,
                               column_id = "FeatureID",
                               select_rows = NULL,
                               select_cols = NULL) {
  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  # column_id_1 <- column_id
  sel_rows <- .select_rows_abundance(x, select_rows)
  sel_cols <- .select_cols_abundance(x, select_cols)

  cols <- c(column_id, sel_cols)
  id.col <- sym(column_id)
  tib_dat <- abundances(x) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rownames_to_column(column_id) %>%
    as_tibble() %>%
    filter(.data[[column_id]] %in% sel_rows)
  # filter(.data[[var]] %in% sel_rows)

  tib_dat <- tib_dat[, cols]
  return(tib_dat)
}


#' @rdname getTibble
#' @aliases getTaxaTibble
#' @importFrom tibble rownames_to_column
#' @importFrom phyloseq tax_table
#' @importFrom dplyr as_tibble %>%
#' @importFrom methods as
#' @importFrom rlang .data sym
#'
#' @export
getTaxaTibble <- function(x,
                          column_id = "FeatureID",
                          select_rows = NULL,
                          select_cols = NULL) {
  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  # column_id_1 <- column_id
  sel_rows <- .select_rows_taxonomy(x, select_rows)
  sel_cols <- .select_cols_taxonomy(x, select_cols)

  cols <- c(column_id, sel_cols)
  tib_dat <- tax_table(x) %>%
    as("matrix") %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rownames_to_column(column_id) %>%
    as_tibble() %>%
    filter(.data[[column_id]] %in% sel_rows)

  tib_dat <- tib_dat[, cols]

  return(tib_dat)
}


#' @rdname getTibble
#' @aliases getSampleTibble
#' @importFrom tibble rownames_to_column
#' @importFrom microbiome meta
#' @importFrom dplyr as_tibble %>%
#' @importFrom rlang .data sym
#'
#' @export
getSampleTibble <- function(x,
                            column_id = "SampleID",
                            select_rows = NULL,
                            select_cols = NULL) {
  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  # column_id_1 <- column_id
  # column_id_1 <- column_id
  sel_rows <- .select_rows_sample(x, select_rows)
  sel_cols <- .select_cols_sample(x, select_cols)


  column_id_nw <- .check_col_sam_data(x, column_id)

  if (column_id_nw == column_id) {
    cols <- c(column_id, sel_cols)
  } else {
    cols <- c(column_id_nw, sel_cols)
  }


  tib_dat <- meta(x) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rownames_to_column(column_id_nw) %>%
    as_tibble() %>%
    filter(.data[[column_id_nw]] %in% sel_rows)

  if(is.null(select_cols) && any(!colnames(tib_dat) %in% cols) || is.null(select_cols) && any(!cols %in% colnames(tib_dat))){
    warning(paste0("Dropped following columns with merging issues ", setdiff(colnames(tib_dat), cols)))
    tib_dat <- tib_dat[, colnames(tib_dat) %in% cols]

  } else {
    tib_dat <- tib_dat[, cols]
  }


  if (!column_id_nw == column_id) {
    tib_dat <- tib_dat %>%
      select(-column_id_nw)
  }
  return(tib_dat)
}


###################### For getSampleTibble ######################
#' @importFrom phyloseq sample_variables

.select_cols_sample <- function(x, select_cols) {
  if (is.null(select_cols) || any(is.na(select_cols))) {
    return(sample_variables(x))
  } else if (any(select_cols %in% sample_variables(x))) {
    return(select_cols)
  } else if (!is.null(select_cols) && !any(is.na(select_cols)) && !any(select_cols %in% sample_variables(x))) {
    stop("For getSampleTibble, please provide valid sample_variables(x) in 'select_cols' ")
  }
}

#' @importFrom phyloseq sample_names
.select_rows_sample <- function(x, select_rows) {
  if (is.null(select_rows) || any(is.na(select_rows))) {
    return(sample_names(x))
  } else if (any(select_rows %in% sample_names(x))) {
    return(select_rows)
  } else if (!is.null(select_rows) && !any(is.na(select_rows)) && !any(select_rows %in% sample_names(x))) {
    stop("For getSampleTibble, please provide valid sample_names(x) in 'select_rows' ")
  }
}

###################### For getTaxaTibble ######################
#' @importFrom phyloseq rank_names
.select_cols_taxonomy <- function(x, select_cols) {
  if (is.null(select_cols) || any(is.na(select_cols))) {
    return(rank_names(x))
  } else if (any(select_cols %in% rank_names(x))) {
    return(select_cols)
  } else if (!is.null(select_cols) && !any(is.na(select_cols)) && !any(select_cols %in% rank_names(x))) {
    stop("For getTaxaTibble, provide valid rank_names(x) in 'select_cols' ")
  }
}

#' @importFrom phyloseq taxa_names
.select_rows_taxonomy <- function(x, select_rows) {

  # )
  if (is.null(select_rows) || any(is.na(select_rows))) {
    return(taxa_names(x))
  } else if (any(select_rows %in% taxa_names(x))) {
    return(select_rows)
  } else if (!is.null(select_rows) && !any(is.na(select_rows)) && !any(select_rows %in% taxa_names(x))) {
    stop("For getTaxaTibble, please provide valid taxa_names(x) in 'select_rows' ")
  }
}



###################### For getAbundanceTibble ######################
#' @importFrom phyloseq sample_names
.select_cols_abundance <- function(x, select_cols) {

  #
  if (is.null(select_cols) || any(is.na(select_cols))) {
    return(sample_names(x))
  } else if (any(select_cols %in% sample_names(x))) {
    return(select_cols)
  } else if (!is.null(select_cols) && !any(is.na(select_cols)) && !any(select_cols %in% sample_names(x))) {
    stop("For getAbundanceTibble, please provide valid sample_names(x) in 'select_cols' ")
  }
}

#' @importFrom phyloseq taxa_names
.select_rows_abundance <- function(x, select_rows) {

  # )
  if (is.null(select_rows) || any(is.na(select_rows))) {
    return(taxa_names(x))
  } else if (any(select_rows %in% taxa_names(x))) {
    return(select_rows)
  } else if (!is.null(select_rows) && !any(is.na(select_rows)) && !any(select_rows %in% taxa_names(x))) {
    stop("For getAbundanceTibble, please provide valid taxa_names(x) in 'select_rows' ")
  }
}

##################################################################

.check_col_sam_data <- function(x, column_id) {

  #global_var
  column_id_nw <- NULL

  if (!(column_id %in% sample_variables(x))){

    column_id_nw <- column_id
    return(column_id_nw)

  } else if (column_id %in% sample_variables(x)) {

    column_id_nw <- paste0(column_id,"XX")
    #message(paste0(column_id, "is already present, hence creating a
     #              new column variable", column_id_nw))
    return(column_id_nw)
  }

}

##################################################################

