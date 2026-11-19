#' Analyze Chemical Genomics (CGS) Data and Perform GSEA
#'
#' This function analyzes chemical genomics (CGS) data, filters by compound mechanism of action (MOA),
#' prepares the data for Gene Set Enrichment Analysis (GSEA), and performs GSEA for specified cell lines and metrics.
#'
#' @param dt_metrics A data.table containing screening data. Requires columns: `drug_moa`, `DrugName`, `CellLineName`,
#' and columns specified in the `metrics` argument.
#' @param metrics A character vector specifying the response metrics to analyze:
#' "x_mean", "x_AOC_range", "xc50", "ec50", "x_max".
#' @param cl_name An optional string specifying a single cell line to analyze. If NULL (default),
#' all cell lines in the data are analyzed. Should be NULL if \code{cellline1} and \code{cellline2} are provided.
#' @param cellline1 An optional string representing the first cell line name.
#' Should be specified alongside \code{cellline2}.
#' @param cellline2 An optional string representing the second cell line name.
#' Should be specified alongside \code{cellline1}.
#' @param normalization_type A string with normalization types to be selected, one of:
#' "GR" ("GRvalue") or "RV" ("RelativeViability").
#' Passed to \code{\link{prep_dt_response_metric_diff}}.
#'
#' @details Each cell line is analysed at each co-treatment concentration it was screened at, so cell lines
#' may use different concentrations. Drugs without a finite value of a metric are left out of the
#' GSEA for that metric.
#'
#' @return A list of results, where each element corresponds to a cell line or cell line difference.
#' A cell line screened at several co-treatment concentrations gets one element per concentration,
#' named \code{"<cell line> (co-treatment <concentration>_<dilution drug>)"}.
#' Each result contains:
#' \itemize{
#'   \item \code{fgsea}: GSEA results for the specified metrics,
#'   \item \code{metrics_diff}: The prepped data.table used for the GSEA analysis,
#'   \item \code{moa_list}: A list of DrugNames grouped by drug_moa.
#' }
#' @keywords cgs_plots
#'
#' @author Bartosz Czech \email{czech.bartosz@@external.gene.com}
#'
#' @examples
#' dt_metrics <- qs2::qs_read(system.file("testdata/cgs_data.qs2", package = "gDRplots"))
#' analyze_cgs(dt_metrics, metrics = "xc50", cl_name = "CellLineName_1")

#' cellline1 <- "CellLineName_1"
#' cellline2 <- "CellLineName_2"
#' analyze_cgs(dt_metrics, "xc50", cl_name = NULL, cellline1, cellline2)
#'
#' @export
analyze_cgs <- function(dt_metrics,
                        metrics,
                        cl_name = NULL,
                        cellline1 = NULL,
                        cellline2 = NULL,
                        normalization_type = "RV") {

  # identifiers
  drug_moa <- gDRutils::get_env_identifiers("drug_moa")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cellline <- gDRutils::get_env_identifiers("cellline_name")

  # asserts
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_character(metrics, any.missing = FALSE)
  checkmate::assert_subset(metrics, choices = c("x_mean", "x_AOC_range", "xc50", "ec50", "x_max"), empty.ok = FALSE)
  checkmate::assert_true(all(metrics %in% names(dt_metrics)))
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))

  if (!is.null(cellline1) || !is.null(cellline2)) {
    checkmate::assert_string(cellline1, null.ok = FALSE)
    checkmate::assert_string(cellline2, null.ok = FALSE)
    checkmate::assert_subset(cellline1, choices = as.character(unique(dt_metrics[[cellline]])), empty.ok = FALSE)
    checkmate::assert_subset(cellline2, choices = as.character(unique(dt_metrics[[cellline]])), empty.ok = FALSE)
  } else {
    checkmate::assert_string(cl_name, null.ok = TRUE)
    checkmate::assert_subset(cl_name, choices = as.character(unique(dt_metrics[[cellline]])), empty.ok = TRUE)
  }

  # filter out unwanted drug moa
  dt_metrics <- dt_metrics[!eval(drug_moa) %in% c("unknown", "Unknown"), ]

  # take care of Inf and NaN values in IC50 metrics
  if (any(is.infinite(dt_metrics$xc50))) {
    message("In `dt_metrics` some xc50 values are infinite.")
    dt_metrics <- dt_metrics[!is.infinite(get("xc50")), ]
  }

  # prepare the data with specified metric differences
  metrics_diff <- prep_dt_response_metric_diff(dt_metrics,
                                               metric = metrics,
                                               d_name = NULL,
                                               d_name2 = NULL,
                                               cellline1 = cellline1,
                                               cellline2 = cellline2,
                                               normalization_type = normalization_type,
                                               additional_cols = drug_moa)

  two_lines <- !is.null(cellline1) && !is.null(cellline2)
  if (two_lines) {
    original <- grep("_cellline_diff$", names(metrics_diff), value = TRUE)
    data.table::setnames(metrics_diff, original, gsub("_cellline_diff$", "", original))
  } else {
    metrics_diff <- .cgs_diff_by_cotrt(metrics_diff)
  }

  # create a list of DrugNames grouped by drug_moa, filtering out those with less than 4 unique drugs
  moa_list <- lapply(split(metrics_diff[[drug_name]], metrics_diff[[drug_moa]]), unique)
  moa_list <- moa_list[vapply(moa_list, length, FUN.VALUE = numeric(1)) > 3]

  metrics_diff <- metrics_diff[eval(drug_moa) %chin% names(moa_list)]

  if (two_lines) {
    # Perform FGSEA on the difference between the specified cell lines
    unit_name <- paste0(cellline1, " vs. ", cellline2)
    list_results <- lapply(metrics, function(metric) {
      .cgs_fgsea(metrics_diff, metric, moa_list, drug_name, drug_moa, unit_name)
    })
    names(list_results) <- metrics
    results <- list(list(fgsea = list_results,
                         metrics_diff = metrics_diff,
                         moa_list = moa_list))
    names(results) <- unit_name
    return(results)
  }

  # Analyze all or a single specified cell line, at each co-treatment concentration it has
  cell_lines <- if (!is.null(cl_name)) {
    cl_name
  } else {
    unique(metrics_diff[[cellline]])
  }
  units <- unique(metrics_diff[get(cellline) %chin% cell_lines, c(cellline, "cotrt"), with = FALSE])
  several <- duplicated(units[[cellline]]) | duplicated(units[[cellline]], fromLast = TRUE)
  unit_names <- units[[cellline]]
  unit_names[several] <- paste0(unit_names[several], " (co-treatment ", units$cotrt[several], ")")

  results <- lapply(seq_len(NROW(units)), function(i) {
    data_subset <- metrics_diff[get(cellline) == units[[cellline]][i] & get("cotrt") == units$cotrt[i]]
    data_subset$normalization_type <- normalization_type
    list_results <- lapply(metrics, function(metric) {
      .cgs_fgsea(data_subset, metric, moa_list, drug_name, drug_moa, unit_names[i])
    })
    names(list_results) <- metrics
    list(fgsea = list_results,
         metrics_diff = data_subset,
         moa_list = moa_list)
  })
  names(results) <- unit_names
  results
}

#' Compare two screening conditions of a chemical genomics screen and perform GSEA
#'
#' Compares the response to each library drug between two conditions, each a cell line at one
#' co-treatment concentration (\code{0} for the drug alone), and runs GSEA over the drug classes
#' (\code{drug_moa}). The conditions may use different cell lines and different co-treatment
#' concentrations, for example a resistant line kept in the anchor drug against its parental line
#' in standard medium. The difference is \code{condition1} minus \code{condition2}.
#'
#' @param dt_metrics A data.table of combination metrics, as for \code{\link{analyze_cgs}}, with
#'   \code{cotrt_value}.
#' @param condition1,condition2 Named lists with \code{cell_line} (a value of \code{CellLineName})
#'   and \code{cotrt_value} (a co-treatment concentration of that cell line, \code{0} for the drug
#'   alone).
#' @param metrics A character vector of metrics to compare: "x_mean", "x_AOC_range", "xc50",
#'   "ec50", "x_max". \code{xc50} is compared on the log10 scale.
#' @param normalization_type A string, "RV" or "GR".
#'
#' @return A list of one element, named \code{"<condition1> vs. <condition2>"}, with \code{fgsea},
#'   \code{metrics_diff} and \code{moa_list} as returned by \code{\link{analyze_cgs}}, so
#'   \code{\link{plot_cgs_ranking}} can draw it.
#' @keywords cgs_plots
#'
#' @examples
#' dt_metrics <- qs2::qs_read(system.file("testdata/cgs_data.qs2", package = "gDRplots"))
#' analyze_cgs_conditions(dt_metrics,
#'                        condition1 = list(cell_line = "CellLineName_1", cotrt_value = 1),
#'                        condition2 = list(cell_line = "CellLineName_2", cotrt_value = 0))
#'
#' @export
analyze_cgs_conditions <- function(dt_metrics,
                                   condition1,
                                   condition2,
                                   metrics = "x_mean",
                                   normalization_type = "RV") {

  drug_moa <- gDRutils::get_env_identifiers("drug_moa")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cellline <- gDRutils::get_env_identifiers("cellline_name")
  norm_col <- gDRutils::get_env_identifiers("normalization_type")

  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_character(metrics, any.missing = FALSE)
  checkmate::assert_subset(metrics, choices = c("x_mean", "x_AOC_range", "xc50", "ec50", "x_max"), empty.ok = FALSE)
  checkmate::assert_true(all(c(metrics, "cotrt_value", cellline, drug_name, drug_moa) %in% names(dt_metrics)))
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  for (cond in list(condition1, condition2)) {
    checkmate::assert_list(cond, names = "unique")
    checkmate::assert_names(names(cond), permutation.of = c("cell_line", "cotrt_value"))
    checkmate::assert_choice(cond$cell_line, choices = as.character(unique(dt_metrics[[cellline]])))
    checkmate::assert_number(cond$cotrt_value, lower = 0)
  }

  # computed outside `[` so that the column names do not shadow the arguments
  keep <- dt_metrics[[norm_col]] == normalization_type &
    !dt_metrics[[drug_moa]] %in% c("unknown", "Unknown")
  dt <- dt_metrics[keep]
  if ("fit_source" %in% names(dt)) {
    dt <- dt[dt$fit_source == "gDR"]
  }
  # non-finite xc50 values are left out per metric by .cgs_fgsea(), so other metrics keep those drugs
  if ("xc50" %in% metrics) {
    dt$xc50 <- log10(dt$xc50)
  }

  # mean response per drug in one condition
  id_cols <- c(drug_name, drug_moa)
  .condition_metrics <- function(cond) {
    rows <- dt[[cellline]] == cond$cell_line
    tol <- sqrt(.Machine$double.eps)
    sel <- dt[rows & abs(dt$cotrt_value - cond$cotrt_value) < tol]
    if (NROW(sel) == 0) {
      stop(sprintf("No %s data for '%s' at co-treatment %s; its co-treatment values are: %s.",
                   normalization_type, cond$cell_line, cond$cotrt_value,
                   toString(sort(unique(dt$cotrt_value[rows])))), call. = FALSE)
    }
    sel[, lapply(.SD, mean, na.rm = TRUE), by = id_cols, .SDcols = metrics]
  }
  .condition_label <- function(cond) {
    if (cond$cotrt_value == 0) {
      return(sprintf("%s (single agent)", cond$cell_line))
    }
    sprintf("%s (co-treatment %s)", cond$cell_line, cond$cotrt_value)
  }

  first <- .condition_metrics(condition1)
  second <- .condition_metrics(condition2)
  metrics_diff <- merge(first, second, by = id_cols, suffixes = c("_c1", "_c2"))
  for (m in metrics) {
    metrics_diff[[m]] <- metrics_diff[[paste0(m, "_c1")]] - metrics_diff[[paste0(m, "_c2")]]
  }
  metrics_diff[[norm_col]] <- rep(normalization_type, NROW(metrics_diff))

  moa_list <- lapply(split(metrics_diff[[drug_name]], metrics_diff[[drug_moa]]), unique)
  moa_list <- moa_list[vapply(moa_list, length, FUN.VALUE = numeric(1)) > 3]
  in_moa <- metrics_diff[[drug_moa]] %chin% names(moa_list)
  metrics_diff <- metrics_diff[in_moa]

  unit_name <- paste0(.condition_label(condition1), " vs. ", .condition_label(condition2))
  list_results <- lapply(metrics, function(metric) {
    .cgs_fgsea(metrics_diff, metric, moa_list, drug_name, drug_moa, unit_name)
  })
  names(list_results) <- metrics
  results <- list(list(fgsea = list_results, metrics_diff = metrics_diff, moa_list = moa_list))
  names(results) <- unit_name
  results
}

#' Reshape co-treatment differences to one row per drug and co-treatment concentration
#'
#' \code{prep_dt_response_metric_diff} returns one difference column per metric and
#' co-treatment concentration (e.g. \code{RV_gDR_x_mean_cotrt_diff_0.03_drug_1}). A cell line
#' screened at one concentration has \code{NA} in the columns of the others, so the columns are
#' stacked into one column per metric, with the concentration in column \code{cotrt}. Rows
#' without any value at a concentration are dropped.
#'
#' @param metrics_diff data.table returned by \code{prep_dt_response_metric_diff}.
#'
#' @return data.table with the identifier columns, \code{cotrt} and one column per metric.
#' @keywords internal
.cgs_diff_by_cotrt <- function(metrics_diff) {
  pattern <- "^.*gDR_(log10_)?(.*)_cotrt_diff_(.*)$"
  diff_cols <- grep("_cotrt_diff_", names(metrics_diff), value = TRUE)
  diff_cols <- diff_cols[!grepl("_fittings$", diff_cols)]
  metric_of <- sub(pattern, "\\2", diff_cols)
  cotrt_of <- sub(pattern, "\\3", diff_cols)
  meta_cols <- names(metrics_diff)[!grepl("cotrt|_fittings$", names(metrics_diff))]

  stacked <- lapply(unique(cotrt_of), function(ct) {
    cols <- diff_cols[cotrt_of == ct]
    out <- metrics_diff[, c(meta_cols, cols), with = FALSE]
    data.table::setnames(out, cols, metric_of[cotrt_of == ct])
    has_value <- rowSums(!is.na(out[, metric_of[cotrt_of == ct], with = FALSE])) > 0
    out <- out[has_value]
    out$cotrt <- rep(ct, NROW(out))
    out
  })
  data.table::rbindlist(stacked, use.names = TRUE, fill = TRUE)
}

#' Run GSEA on one metric of one cell line (or cell line comparison)
#'
#' Drugs without a finite value of the metric are left out, since \code{fgsea} accepts only
#' finite statistics. When no drug is left, an empty result is returned.
#'
#' @param dt data.table with the metric column, drug names and MOAs.
#' @param metric string, the metric column.
#' @param moa_list named list of drug names per MOA.
#' @param name_col,moa_col strings, the drug name and MOA columns.
#' @param unit_name string, the cell line (or comparison) used in messages.
#'
#' @return data.table with the \code{fgsea} result and a \code{median} column; \code{leadingEdge}
#'   holds the leading-edge drugs as a comma-separated string.
#' @keywords internal
.cgs_fgsea <- function(dt, metric, moa_list, name_col, moa_col, unit_name) {
  metric_values <- dt[[metric]]
  names(metric_values) <- dt[[name_col]]
  finite <- is.finite(metric_values)
  if (!all(finite)) {
    message(sprintf("%s, %s: %d drug(s) without a finite value left out of GSEA",
                    unit_name, metric, sum(!finite)))
  }
  metric_values <- metric_values[finite]

  if (length(metric_values) == 0) {
    message(sprintf("%s, %s: no drugs with a finite value, GSEA skipped", unit_name, metric))
    return(data.table::data.table(pathway = character(0), pval = numeric(0), padj = numeric(0),
                                  log2err = numeric(0), ES = numeric(0), NES = numeric(0),
                                  size = integer(0), leadingEdge = character(0), median = numeric(0)))
  }

  fgsea_result <- purrr::quietly(fgsea::fgsea)(pathways = moa_list,
                                               stats = metric_values,
                                               maxSize = 500,
                                               minSize = 4,
                                               nPermSimple = 1e5)$result

  medians <- dt[, list(median = stats::median(.SD[[1]], na.rm = TRUE)), by = moa_col, .SDcols = metric]
  median_values <- stats::setNames(medians$median, medians[[moa_col]])

  fgsea_result$median <- median_values[fgsea_result$pathway]
  # a comma-separated string instead of a list column, so the table can be written to xlsx
  fgsea_result$leadingEdge <- vapply(fgsea_result$leadingEdge, paste, character(1), collapse = ", ")
  data.table::setorder(fgsea_result, -NES)
  fgsea_result
}


#' Plot Chemical Genomics Screen GSEA Results
#'
#' Generates a ggplot2 visualization of chemical genomics screening data, highlighting GSEA results.
#'
#' @param results A list object returned from `analyze_cgs`.
#' @param cl_name A string specifying the cell line included in the \code{results} to prepare a visualization.
#' @param metric A string specifying the metric included in the \code{results} to prepare a visualization.
#' @param padj_threshold A numeric value specifying the threshold for filtering significant GSEA results
#' based on adjusted p-value.
#' @param top_results_no_sig A numeric value specifying the number of top results to plot
#' if there are no significant values.
#' @param max_results_with_sig A numeric value specifying the maximum number of results
#' to plot when there are more than this number of significant values.
#'
#' @return A ggplot2 object with cgs results
#'
#' @keywords cgs_plots
#'
#' @author Bartosz Czech \email{czech.bartosz@@external.gene.com}
#'
#' @examples
#' dt_metrics <- qs2::qs_read(system.file("testdata/cgs_data.qs2", package = "gDRplots"))
#' results <- analyze_cgs(dt_metrics, metrics = c("xc50"), cl_name = "CellLineName_1")
#' plot_cgs_ranking(results,
#'   cl_name = "CellLineName_1",
#'   metric = "xc50",
#'   padj_threshold = 0.1,
#'   top_results_no_sig = 5,
#'   max_results_with_sig = 15)
#' @export
plot_cgs_ranking <- function(results,
                             cl_name,
                             metric,
                             padj_threshold = 0.1,
                             top_results_no_sig = 5,
                             max_results_with_sig = 15) {

  # identifiers
  drug_moa <- gDRutils::get_env_identifiers("drug_moa")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  norm_type <- gDRutils::get_env_identifiers("normalization_type")

  # asserts
  checkmate::assert_list(results)
  checkmate::assert_string(cl_name, null.ok = FALSE)
  checkmate::assert_subset(cl_name, choices = names(results))
  checkmate::assert_string(metric, null.ok = FALSE)
  checkmate::assert_subset(metric, choices = names(results[[cl_name]]$fgsea), empty.ok = FALSE)


  # extract relevant data
  metrics_diff <- results[[cl_name]]$metrics_diff
  fgsea_results <- results[[cl_name]]$fgsea[[metric]]
  moa_groups_drugs <- results[[cl_name]]$moa_list

  # prepare data for plotting
  plot_data <- data.table::copy(metrics_diff)
  plot_data$x_pos <- NROW(plot_data) - rank(plot_data[[metric]]) + 1
  stats <- plot_data[[metric]]
  stats <- pmin(2, pmax(-2, stats))
  names(stats) <- plot_data[[drug_name]]
  stats <- stats[!is.na(plot_data[[drug_moa]])]

  norm_type <- unique(metrics_diff[[norm_type]])

  # filter significant GSEA results
  gsea_sign <- fgsea_results[padj < padj_threshold & !pathway %in% c("", "unknown")]
  if (NROW(gsea_sign) == 0) {
    # if there are no significant values, plot only top_results_no_sig top results
    gsea_sign <- fgsea_results[pval < sort(pval)[top_results_no_sig]]
  } else if (NROW(gsea_sign) > max_results_with_sig) {
    # if there are more than max_results_with_sig significant values, plot only max_results_with_sig top results
    gsea_sign <- utils::head(gsea_sign[order(padj)], max_results_with_sig)
  }
  gsea_sign[, y_pos := seq_len(.N)]
  gsea_sign[NES < 0, y_pos := -(seq_len(.N))]

  # pre-calculate statistics
  stats <- plot_data[[metric]]
  threshold_count <- sum(stats > mean(stats))
  yrange <- diff(range(stats))
  mean_effect <- mean(stats)

  # create the ggplot object
  plt <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_pos, y = !!rlang::sym(metric))) +
    ggplot2::geom_col(color = .get_setting("EDGE_COLOR")) +
    ggplot2::labs(title = cl_name,
                  y = bquote(~ Delta ~ .(metric) ~ "for" ~ .(norm_type)),
                  x = "Ranked drugs",
                  caption = sprintf("Top results with FDR < %.2f are shown. If no results meet this threshold,
                                    the top %d results by p-value are displayed.",
                                    padj_threshold, top_results_no_sig)
                  ) +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = 0,
                        color = .get_setting("HLINE_COLOR")) +
    ggplot2::geom_hline(yintercept = mean_effect, color = "black") +
    ggplot2::geom_segment(x = threshold_count, xend = threshold_count,
                          y = 0, yend = mean_effect + 0.2 * yrange,
                          color = "black") +
    ggplot2::annotate(geom = "text", x = threshold_count, y = mean_effect + 0.25 * yrange,
                      label = sprintf("Mean effect = %.2f", mean_effect),
                      hjust = 0, color = "black") +
    ggplot2::coord_cartesian(xlim = c(-2, NROW(plot_data) + 3),
                             ylim = c(-1.01 * yrange - 0.15 * yrange * NROW(gsea_sign), yrange + 0.01),
                             expand = FALSE, clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::unit(c(1, 16, 1, 1), "lines"))

  # define color palettes for the loop
  loop_colors <- gDRplots::get_qual_colors(NROW(gsea_sign))

  for (i in seq_len(NROW(gsea_sign))) {
    pathway <- gsea_sign$pathway[i]
    x <- plot_data[get(drug_name) %in% moa_groups_drugs[[pathway]]]$x_pos
    stats_moa <- plot_data[get(drug_name) %in% moa_groups_drugs[[pathway]]][[metric]]

    median_moa <- stats::median(stats_moa)
    count_above_median <- sum(stats > median_moa)

    current_color <- loop_colors[i]

    plt <- plt +
      ggplot2::geom_segment(
        data = data.frame(x = x),
        ggplot2::aes(x = x, xend = x),
        y = -yrange - (0.15 * yrange * (i - 1)),
        yend = -yrange - (0.15 * yrange * i),
        linewidth = 0.8, inherit.aes = FALSE, color = current_color) +
      ggplot2::annotate(
        geom = "text",
        x = max(plot_data$x_pos) + 10,
        y = -yrange - (0.15 * yrange * (i - 0.5)),
        label = sprintf("%s, NES=%.2f, FDR=%.2g", gsub("_", " ", pathway), gsea_sign$NES[i], gsea_sign$padj[i]),
        hjust = 0, color = current_color) +
      ggplot2::geom_segment(
        x = count_above_median - 0.5,
        xend = count_above_median - 0.5,
        y = median_moa, yend = -(gsea_sign$y_pos[i] + 0.5 * sign(gsea_sign$y_pos[i])) * 0.15 * yrange,
        color = current_color) +
      ggplot2::annotate(
        geom = "label",
        x = count_above_median,
        y = -(gsea_sign$y_pos[i] + 0.5 * sign(gsea_sign$y_pos[i])) * 0.185 * yrange,
        label = sprintf(" %s median = %.2f ", pathway, median_moa),
        hjust = 1 * (gsea_sign$NES[i] < 0),
        color = current_color,
        fill = "white", alpha = 0.65)
  }
  return(plt)
}
