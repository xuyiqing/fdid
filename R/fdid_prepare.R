#' Prepare Data for Factorial Difference-in-Differences Analysis
#'
#' Prepares a dataset for factorial difference-in-differences (FDID) analysis by reshaping the data into a wide format,
#' averaging time-varying covariates, and renaming columns for consistency in subsequent analysis.
#'
#' @param data A data frame containing the dataset to be processed.
#' @param Y_label A string specifying the column name of the outcome variable.
#' @param X_labels A character vector specifying the column names of the time-varying covariates.
#' @param G_label A string specifying the column name of the group variable (e.g., treatment vs. control).
#' @param unit_label A string specifying the column name of the unit identifier (e.g., individual or entity).
#' @param time_label A string specifying the column name of the time variable.
#' @param cluster_label An optional string specifying the column name of the clustering variable. Default is `NULL`.
#'
#' @return A data frame in wide format with the following:
#' - Outcome variable pivoted to wide format with time columns.
#' - Time-varying covariates averaged across time.
#' - Columns renamed:
#'   - Unit identifier -> `unit`
#'   - Covariates -> `x1`, `x2`, ...
#'   - Group variable -> `G`
#'   - Clustering variable (if provided) -> `c`
#'
#' @examples
#' data <- data.frame(
#'   id = rep(1:3, each = 4),
#'   time = rep(1:4, times = 3),
#'   outcome = rnorm(12),
#'   covar1 = runif(12),
#'   covar2 = runif(12),
#'   group = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1)
#' )
#' fdid_data <- fdid_prepare(
#'   data = data,
#'   Y_label = "outcome",
#'   X_labels = c("covar1", "covar2"),
#'   G_label = "group",
#'   unit_label = "id",
#'   time_label = "time"
#' )
#' head(fdid_data)
#'
#' @import dplyr
#' @import tidyr
#' @author Rivka Lipkovitz
#' @export
fdid_prepare <- function(data,
                             Y_label,
                             X_labels = NULL,
                             G_label,
                             unit_label,
                             time_label,
                             cluster_label = NULL) {

  library(dplyr)
  library(tidyr)

  # 1. Pivot the outcome to wide
  wide_data <- data %>%
    tidyr::pivot_wider(
      id_cols     = tidyselect::all_of(unit_label),
      names_from  = tidyselect::all_of(time_label),
      values_from = tidyselect::all_of(Y_label),
      names_prefix = ""  # keep the time values as column names
    )

  # 2. Average the time-varying covariates using older dplyr syntax
  #    Summarize only the X_labels
  covar_data <- data %>%
    group_by(!!rlang::sym(unit_label)) %>%
    dplyr::summarise_at(
      .vars = X_labels,
      .funs = ~ mean(.x, na.rm = TRUE)
    ) %>%
    ungroup()

  # 3. Extract G_label and cluster_label (if present).
  #    We'll take the first occurrence for each unit
  #    (assuming they do not vary over time).
  #    Then we'll join these back to covar_data.
  if (!is.null(cluster_label)) {
    group_cluster_data <- data %>%
      distinct(
        !!rlang::sym(unit_label),
        !!rlang::sym(G_label),
        !!rlang::sym(cluster_label)
      )
    covar_data <- covar_data %>%
      left_join(group_cluster_data, by = unit_label)
  } else {
    group_data <- data %>%
      distinct(
        !!rlang::sym(unit_label),
        !!rlang::sym(G_label)
      )
    covar_data <- covar_data %>%
      left_join(group_data, by = unit_label)
  }

  # 4. Join these averaged covariates + group info back to wide outcome data
  wide_data <- wide_data %>%
    left_join(covar_data, by = unit_label)

  # 5. Rename columns:
  #    - First column => "unit"
  #    - Covariates => x1, x2, ...
  #    - G_label => "G"
  #    - cluster_label => "c"
  #    (using older rename_at / rename_with might need caution in older dplyr)
  colnames(wide_data)[1] <- "unit"

  # figure out which columns are the newly-averaged covariates
  # (the same names as X_labels), so we can rename them to x1, x2, ...
  x_cols <- X_labels
  # we also rename them in the same order they appear in X_labels
  for (i in seq_along(x_cols)) {
    old_col <- x_cols[i]
    new_col <- paste0("x", i)
    wide_data <- wide_data %>%
      rename(!!new_col := !!rlang::sym(old_col))
  }

  # rename G_label => "G"
  wide_data <- wide_data %>%
    rename(G = !!rlang::sym(G_label))

  # rename cluster_label => "c"
  if (!is.null(cluster_label)) {
    wide_data <- wide_data %>%
      rename(c = !!rlang::sym(cluster_label))
  }

  return(wide_data)
}
