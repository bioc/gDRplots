context("Test plot-chemical_genomics")

# Load test data

metrics_data <- qs2::qs_read(system.file("testdata/cgs_data.qs2", package = "gDRplots"))

test_that("analyze_cgs works correctly", {
  # test with a single cell line and single metric
  results <- analyze_cgs(metrics_data, metrics = "xc50", cl_name = "CellLineName_1")
  expect_is(results, "list")
  expect_equal(names(results), "CellLineName_1")
  expect_is(results$CellLineName_1$fgsea$xc50, "data.table")
  expect_gt(NROW(results$CellLineName_1$fgsea$xc50), 0) # Check for results
  expect_is(results$CellLineName_1$metrics_diff, "data.table")
  expect_gt(NROW(results$CellLineName_1$metrics_diff), 0)
  expect_is(results$CellLineName_1$moa_list, "list")
  expect_gt(length(results$CellLineName_1$moa_list), 0)


  # test with multiple metrics
  results <- analyze_cgs(metrics_data, metrics = c("xc50", "x_max"), cl_name = "CellLineName_1")
  expect_equal(names(results$CellLineName_1$fgsea), c("xc50", "x_max"))
  expect_is(results$CellLineName_1$fgsea$xc50, "data.table")
  expect_is(results$CellLineName_1$fgsea$x_max, "data.table")

  # test normalization type
  results1 <- analyze_cgs(metrics_data, metrics = "xc50", cl_name = "CellLineName_1", normalization_type = "RV")
  results2 <- analyze_cgs(metrics_data, metrics = "xc50", cl_name = "CellLineName_1", normalization_type = "GR")

  expect_false(identical(results1$CellLineName_1$metrics_diff$xc50, results2$CellLineName_1$metrics_diff$xc50))

  # no cell line selection
  results3 <- analyze_cgs(metrics_data, metrics = "xc50", cl_name = NULL)
  expect_is(results3, "list")
  expect_true(setequal(names(results3), unique(metrics_data$CellLineName)))

  # cell lines diff
  results4 <- analyze_cgs(metrics_data,
                          metrics = "xc50",
                          cl_name = NULL,
                          cellline1 = "CellLineName_1",
                          cellline2 = "CellLineName_2",
                          normalization_type = "RV")
  expect_is(results4, "list")
  expect_equal(names(results4), "CellLineName_1 vs. CellLineName_2")

  # scenario: CellLineName column is a factor
  metrics_data_factor <- data.table::copy(metrics_data)
  metrics_data_factor[["CellLineName"]] <- factor(metrics_data_factor[["CellLineName"]])

  results_factor <- analyze_cgs(metrics_data_factor, metrics = "xc50", cl_name = "CellLineName_1")
  expect_equal(names(results_factor), "CellLineName_1")
  expect_is(results_factor$CellLineName_1$fgsea$xc50, "data.table")
  expect_gt(NROW(results_factor$CellLineName_1$fgsea$xc50), 0)
  expect_is(results_factor$CellLineName_1$metrics_diff, "data.table")
  expect_gt(NROW(results_factor$CellLineName_1$metrics_diff), 0)

  results4_factor <- analyze_cgs(metrics_data_factor,
                                 metrics = "xc50",
                                 cl_name = NULL,
                                 cellline1 = "CellLineName_1",
                                 cellline2 = "CellLineName_2",
                                 normalization_type = "RV")
  expect_equal(names(results4_factor), "CellLineName_1 vs. CellLineName_2")
})

test_that("plot_cgs_ranking works correctly", {

  results <- analyze_cgs(metrics_data, metrics = c("xc50", "x_max"), cl_name = "CellLineName_1")
  plt <- plot_cgs_ranking(results, cl_name = "CellLineName_1", metric = "xc50")
  expect_is(plt, "ggplot")

  # rest with a different metric
  results <- analyze_cgs(metrics_data, metrics = c("xc50", "x_max"), cl_name = "CellLineName_1")
  plt <- plot_cgs_ranking(results, cl_name = "CellLineName_1", metric = "x_max")
  expect_is(plt, "ggplot")

  # test case where no significant GSEA results exist (using a dummy metric)
  results$CellLineName_1$fgsea$xc50$padj <- 1  # Set all padj values to 1 to simulate no significance
  plt <- plot_cgs_ranking(results, cl_name = "CellLineName_1", metric = "xc50") # Should still produce a plot
  expect_is(plt, "ggplot")
})

test_that("analyze_cgs handles cell lines screened at different co-treatment concentrations", {
  cellline <- gDRutils::get_env_identifiers("cellline_name")
  dt <- data.table::copy(metrics_data)
  # the second cell line gets another anchor concentration than the first
  dt[get(cellline) == "CellLineName_2" & cotrt_value == 1, cotrt_value := 2]

  results <- analyze_cgs(dt, metrics = c("xc50", "x_mean"))
  expect_setequal(names(results), c("CellLineName_1", "CellLineName_2"))
  for (cl in names(results)) {
    expect_true(all(is.finite(results[[cl]]$metrics_diff$x_mean)))
    expect_s3_class(results[[cl]]$fgsea$x_mean, "data.table")
  }
  expect_gt(NROW(results$CellLineName_1$fgsea$x_mean), 0)
  expect_equal(unique(results$CellLineName_2$metrics_diff$cotrt), "2_drug_1")
  expect_equal(unique(results$CellLineName_1$metrics_diff$cotrt), "1_drug_1")
})

test_that("analyze_cgs analyses a cell line at each of its co-treatment concentrations", {
  cellline <- gDRutils::get_env_identifiers("cellline_name")
  dt <- data.table::copy(metrics_data)[get(cellline) == "CellLineName_1"]
  second <- data.table::copy(dt[cotrt_value == 1])
  second[, `:=`(cotrt_value = 2, x_mean = x_mean * 0.9)]
  dt <- rbind(dt, second)

  results <- analyze_cgs(dt, metrics = "x_mean")
  expect_setequal(names(results), c("CellLineName_1 (co-treatment 1_drug_1)",
                                    "CellLineName_1 (co-treatment 2_drug_1)"))
  expect_false(isTRUE(all.equal(results[[1]]$metrics_diff$x_mean, results[[2]]$metrics_diff$x_mean)))
})

test_that(".cgs_fgsea leaves out drugs without a finite value", {
  dt <- data.table::data.table(DrugName = paste0("d", 1:12),
                               drug_moa = rep(c("A", "B", "C"), each = 4),
                               x_mean = c(seq(-1, 1, length.out = 10), NA, Inf))
  moa_list <- split(dt$DrugName, dt$drug_moa)
  expect_message(res <- .cgs_fgsea(dt, "x_mean", moa_list, "DrugName", "drug_moa", "cl"),
                 "2 drug\\(s\\) without a finite value")
  expect_s3_class(res, "data.table")
  expect_true("median" %in% names(res))
  expect_type(res$leadingEdge, "character")

  dt$x_mean <- NA_real_
  expect_message(empty <- .cgs_fgsea(dt, "x_mean", moa_list, "DrugName", "drug_moa", "cl"),
                 "GSEA skipped")
  expect_equal(NROW(empty), 0)
})

test_that("analyze_cgs_conditions compares two conditions of different cell lines and co-treatments", {
  cond1 <- list(cell_line = "CellLineName_1", cotrt_value = 1)
  cond2 <- list(cell_line = "CellLineName_2", cotrt_value = 0)
  res <- analyze_cgs_conditions(metrics_data, cond1, cond2, metrics = c("x_mean", "xc50"))

  expect_equal(names(res), "CellLineName_1 (co-treatment 1) vs. CellLineName_2 (single agent)")
  expect_named(res[[1]], c("fgsea", "metrics_diff", "moa_list"))
  expect_named(res[[1]]$fgsea, c("x_mean", "xc50"))
  expect_s3_class(res[[1]]$fgsea$x_mean, "data.table")
  expect_type(res[[1]]$fgsea$x_mean$leadingEdge, "character")

  # the difference is condition1 minus condition2, drug by drug
  d <- res[[1]]$metrics_diff
  expect_equal(d$x_mean, d$x_mean_c1 - d$x_mean_c2)
  drug <- d$DrugName[1]
  rv <- metrics_data[normalization_type == "RV" & fit_source == "gDR" & DrugName == drug]
  expected <- mean(rv[CellLineName == "CellLineName_1" & cotrt_value == 1, x_mean]) -
    mean(rv[CellLineName == "CellLineName_2" & cotrt_value == 0, x_mean])
  expect_equal(d$x_mean[1], expected)

  plt <- plot_cgs_ranking(res, names(res), "x_mean")
  expect_s3_class(plt, "ggplot")
})

test_that("analyze_cgs_conditions reports a condition without data", {
  cond1 <- list(cell_line = "CellLineName_1", cotrt_value = 5)
  cond2 <- list(cell_line = "CellLineName_2", cotrt_value = 0)
  expect_error(analyze_cgs_conditions(metrics_data, cond1, cond2), "its co-treatment values are: 0, 1")
  expect_error(analyze_cgs_conditions(metrics_data, list(cell_line = "nope", cotrt_value = 0), cond2),
               "Must be element")
  expect_error(analyze_cgs_conditions(metrics_data, list(cell_line = "CellLineName_1"), cond2), "Names must")
})
