# Every spatial metric must return the same columns, in the same order, differing
# only in the name of its "Observed" column. That is what lets downstream code --
# plotting, joining to clinical data, dplyr::bind_rows() across metrics -- treat
# them uniformly instead of special-casing each one.
#
# Before 2.0.0 there were effectively four different schemas: Theoretical CSR vs
# Theoretical G vs Theoretical g; Anchor/Counted vs From/To; Degree of Clustering
# vs Degree of Correlation vs Degree of Interaction; and
# Permuted_larger_than_Observed vs "Permutations Larger than Observed". This file
# is what stops that drift recurring, so add new metrics here rather than
# loosening it.

metric_fixture <- function() {
  mif <- create_mif(
    clinical_data = example_clinical %>%
      dplyr::mutate(deidentified_id = as.character(deidentified_id)),
    sample_data = example_summary %>%
      dplyr::mutate(deidentified_id = as.character(deidentified_id)),
    spatial_list = example_spatial["TMA3_[9,K].tif"],
    patient_id = "deidentified_id", sample_id = "deidentified_sample")
  list(mif = mif,
       mnames = c("CD3..Opal.570..Positive", "CD8..Opal.520..Positive"),
       r = seq(0, 60, 10))
}

# Every metric, the derived slot it writes, its Observed column, and whether it is
# bivariate. Adding a metric to the package means adding a row here.
metric_table <- function() {
  list(
    list(name = "ripleys_k",            slot = "univariate_Count",
         observed = "Observed K", bivariate = FALSE,
         run = function(f) ripleys_k(f$mif, mnames = f$mnames, r_range = f$r,
                                     permute = FALSE, workers = 1, overwrite = TRUE)),
    list(name = "bi_ripleys_k",         slot = "bivariate_Count",
         observed = "Observed K", bivariate = TRUE,
         run = function(f) bi_ripleys_k(f$mif, mnames = f$mnames, r_range = f$r,
                                        permute = FALSE, workers = 1, overwrite = TRUE)),
    list(name = "NN_G",                 slot = "univariate_NN",
         observed = "Observed G", bivariate = FALSE,
         run = function(f) NN_G(f$mif, mnames = f$mnames, r_range = f$r,
                                num_permutations = 2, workers = 1, overwrite = TRUE)),
    list(name = "bi_NN_G",              slot = "bivariate_NN",
         observed = "Observed G", bivariate = TRUE,
         run = function(f) bi_NN_G(f$mif, mnames = f$mnames, r_range = f$r,
                                   num_permutations = 2, workers = 1, overwrite = TRUE)),
    list(name = "pair_correlation",     slot = "univariate_pair_correlation",
         observed = "Observed g", bivariate = FALSE,
         run = function(f) pair_correlation(f$mif, mnames = f$mnames, r_range = f$r,
                                            num_permutations = 2, workers = 1, overwrite = TRUE)),
    list(name = "bi_pair_correlation",  slot = "bivariate_pair_correlation",
         observed = "Observed g", bivariate = TRUE,
         run = function(f) bi_pair_correlation(f$mif, mnames = f$mnames, r_range = f$r,
                                               num_permutations = 2, workers = 1, overwrite = TRUE)),
    list(name = "interaction_variable", slot = "interaction_variable",
         observed = "Observed Interaction", bivariate = TRUE,
         run = function(f) interaction_variable(f$mif, mnames = f$mnames, r_range = f$r,
                                                num_permutations = 2, workers = 1, overwrite = TRUE))
  )
}

test_that("every metric returns the canonical columns in canonical order", {
  f <- metric_fixture()
  for (m in metric_table()) {
    out <- suppressWarnings(suppressMessages(m$run(f)))$derived[[m$slot]]
    expect_true(is.data.frame(out), info = m$name)
    want <- c(standard_metric_cols("deidentified_sample", m$observed, m$bivariate), "Run")
    expect_identical(names(out), want, info = m$name)
  }
})

test_that("metric results from different metrics bind together cleanly", {
  # The practical payoff of one schema.
  f <- metric_fixture()
  tabs <- lapply(metric_table(), function(m) {
    out <- suppressWarnings(suppressMessages(m$run(f)))$derived[[m$slot]]
    as.data.frame(out)
  })
  combined <- dplyr::bind_rows(tabs)
  expect_gt(nrow(combined), 0)
  # Shared columns must genuinely be shared, not silently duplicated with variant
  # spellings.
  for (nm in c("r", "iter", "Theoretical CSR", "Permuted CSR", "Exact CSR",
               "Permutations Larger than Observed", "Degree of Clustering Theoretical",
               "Degree of Clustering Permutation", "Degree of Clustering Exact", "Run")) {
    expect_true(nm %in% names(combined), info = nm)
  }
  # No leftovers from the old per-metric naming.
  for (gone in c("Theoretical G", "Theoretical g", "Permuted G", "Permuted g",
                 "Permuted Interaction", "From", "To",
                 "Permuted_larger_than_Observed",
                 "Degree of Correlation Theoretical", "Degree of Correlation Permuted",
                 "Degree of Interaction Permuted")) {
    expect_false(gone %in% names(combined), info = gone)
  }
})

test_that("Exact CSR is populated only for Ripley's K", {
  f <- metric_fixture()
  for (m in metric_table()) {
    out <- suppressWarnings(suppressMessages(m$run(f)))$derived[[m$slot]]
    if (grepl("ripleys_k", m$name)) {
      expect_true(any(!is.na(out$`Exact CSR`)), info = m$name)
    } else {
      expect_true(all(is.na(out$`Exact CSR`)), info = m$name)
    }
  }
})

test_that("derived slot names are spelled correctly on the append path", {
  # Three functions used to write to a misspelled slot when overwrite = FALSE:
  # univaraite_pair_correlation, bivaraite_pair_correlation, and
  # mif$derived_intraction_variable. Appended runs therefore vanished.
  f <- metric_fixture()
  for (m in metric_table()) {
    first  <- suppressWarnings(suppressMessages(m$run(f)))
    f2 <- f; f2$mif <- first
    # Re-run with overwrite = FALSE by mutating the fixture's mif in place.
    second <- suppressWarnings(suppressMessages({
      g <- f2
      switch(m$name,
        ripleys_k            = ripleys_k(g$mif, mnames = g$mnames, r_range = g$r,
                                          permute = FALSE, workers = 1, overwrite = FALSE),
        bi_ripleys_k         = bi_ripleys_k(g$mif, mnames = g$mnames, r_range = g$r,
                                             permute = FALSE, workers = 1, overwrite = FALSE),
        NN_G                 = NN_G(g$mif, mnames = g$mnames, r_range = g$r,
                                     num_permutations = 2, workers = 1, overwrite = FALSE),
        bi_NN_G              = bi_NN_G(g$mif, mnames = g$mnames, r_range = g$r,
                                        num_permutations = 2, workers = 1, overwrite = FALSE),
        pair_correlation     = pair_correlation(g$mif, mnames = g$mnames, r_range = g$r,
                                                num_permutations = 2, workers = 1, overwrite = FALSE),
        bi_pair_correlation  = bi_pair_correlation(g$mif, mnames = g$mnames, r_range = g$r,
                                                   num_permutations = 2, workers = 1, overwrite = FALSE),
        interaction_variable = interaction_variable(g$mif, mnames = g$mnames, r_range = g$r,
                                                    num_permutations = 2, workers = 1, overwrite = FALSE))
    }))
    # The append must land in the SAME slot, with Run incremented, and must not
    # create any new slot.
    expect_identical(sort(names(second$derived)), sort(names(first$derived)), info = m$name)
    expect_identical(sort(unique(second$derived[[m$slot]]$Run)), c(1, 2), info = m$name)
    expect_equal(nrow(second$derived[[m$slot]]),
                 2 * nrow(first$derived[[m$slot]]), info = m$name)
  }
})
