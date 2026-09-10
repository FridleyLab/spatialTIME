# The golden-value test here is the single most important assertion added in 2.0.0.
# Before it, marker_freq_diff() returned a p-value from a contingency table whose
# second row was the compartment TOTAL rather than the count of marker-NEGATIVE
# cells, so the margin double-counted the positives and every p-value the function
# had ever produced was wrong. Nothing caught it because the function had no test.

test_that("p-values equal fisher.test on the positives-vs-negatives table", {
  spat <- example_spatial[["TMA3_[9,K].tif"]]
  mif  <- example_mif_full()
  markers <- c("CD3..Opal.570..Positive", "CD8..Opal.520..Positive",
               "FOXP3..Opal.620..Positive")

  out <- marker_freq_diff(mif, classifier = "Classifier.Label",
                          ref_level = "Tumor", diff_level = "Stroma",
                          mnames = markers, overwrite = TRUE)$derived$frequency_difference

  for (m in markers) {
    pos <- tapply(spat[[m]], spat$Classifier.Label, sum)[c("Tumor", "Stroma")]
    tot <- as.numeric(table(spat$Classifier.Label)[c("Tumor", "Stroma")])
    want <- stats::fisher.test(
      rbind(positive = as.numeric(pos), negative = tot - as.numeric(pos))
    )$p.value
    expect_equal(out[[paste0(m, "_p.value")]], want, info = m)
  }
})

test_that("the wrong table is no longer produced", {
  # Pin the specific regression: positives-vs-totals gives 7.949897e-07 on this
  # marker where positives-vs-negatives gives 6.001431e-07. Assert we return the
  # latter and NOT the former, so a re-introduction cannot pass.
  spat <- example_spatial[["TMA1_[3,B].tif"]]
  mif  <- example_mif_full("TMA1_[3,B].tif")
  m <- "CD3..Opal.570..Positive"
  got <- marker_freq_diff(mif, classifier = "Classifier.Label", ref_level = "Tumor",
                          diff_level = "Stroma", mnames = m,
                          overwrite = TRUE)$derived$frequency_difference[[paste0(m, "_p.value")]]
  # RELATIVE comparisons. all.equal() falls back to an ABSOLUTE difference once the
  # values are smaller than `tolerance`, so at tolerance = 1e-6 these two p-values --
  # which differ by 32% -- would compare as equal and the test would prove nothing.
  rel <- function(a, b) abs(a - b) / abs(b)
  expect_lt(rel(got, 6.001431e-07), 1e-4)   # the correct answer
  expect_gt(rel(got, 7.949897e-07), 0.1)    # the old, wrong one
})

test_that("a marker whose name is a prefix of another gets its own 2x2 table", {
  # dplyr::contains(marker) is a SUBSTRING match, so the shorter marker's table
  # used to absorb the longer marker's counts and fisher.test() was handed a 3x2
  # table whose r x c p-value was stored as that marker's 2x2 result. Both of these
  # are real columns of the shipped data.
  mif <- example_mif_full("TMA3_[9,K].tif")
  pair <- c("CD3..CD8.", "CD3..CD8..FOXP3.")
  together <- marker_freq_diff(mif, classifier = "Classifier.Label", ref_level = "Tumor",
                               diff_level = "Stroma", mnames = pair,
                               overwrite = TRUE)$derived$frequency_difference
  for (m in pair) {
    alone <- marker_freq_diff(mif, classifier = "Classifier.Label", ref_level = "Tumor",
                              diff_level = "Stroma", mnames = m,
                              overwrite = TRUE)$derived$frequency_difference
    expect_equal(together[[paste0(m, "_p.value")]],
                 alone[[paste0(m, "_p.value")]], info = m)
  }
})

test_that("the frequency difference column is ref minus diff, in percent", {
  spat <- example_spatial[["TMA3_[9,K].tif"]]
  mif  <- example_mif_full()
  m <- "CD3..Opal.570..Positive"
  out <- marker_freq_diff(mif, classifier = "Classifier.Label", ref_level = "Tumor",
                          diff_level = "Stroma", mnames = m,
                          overwrite = TRUE)$derived$frequency_difference
  pct <- tapply(spat[[m]], spat$Classifier.Label, mean) * 100
  expect_equal(out[[paste0("Tumor_M_Stroma ", m)]],
               unname(pct[["Tumor"]] - pct[["Stroma"]]))
  # And the per-level percent columns agree with that.
  expect_equal(out[[paste0("Tumor ", m, "%")]], unname(pct[["Tumor"]]))
})

test_that("missing classifier or level is rejected with a useful message", {
  mif <- example_mif()
  expect_error(
    marker_freq_diff(mif, classifier = "NotAColumn", ref_level = "Tumor",
                     diff_level = "Stroma", mnames = mnames_good()[1]),
    "classifier column not found"
  )
  expect_error(
    marker_freq_diff(mif, classifier = "Classifier.Label", ref_level = "Nope",
                     diff_level = "Stroma", mnames = mnames_good()[1]),
    "ref_level value not found"
  )
  expect_error(
    marker_freq_diff(mif, classifier = "Classifier.Label", ref_level = "Tumor",
                     diff_level = "Nope", mnames = mnames_good()[1]),
    "diff_level value not found"
  )
  expect_error(
    marker_freq_diff(mif, classifier = "Classifier.Label", ref_level = "Tumor",
                     diff_level = "Stroma", mnames = 1),
    "vector of marker names"
  )
  expect_error(
    marker_freq_diff(mif, classifier = c("a", "b"), ref_level = "Tumor",
                     diff_level = "Stroma", mnames = mnames_good()[1]),
    "classifier must be of length 1"
  )
})

test_that("one row per sample, carrying the sample id", {
  mif <- example_mif(which = 1:3, n_cells = 300)
  out <- marker_freq_diff(mif, classifier = "Classifier.Label", ref_level = "Tumor",
                          diff_level = "Stroma", mnames = mnames_good()[1],
                          overwrite = TRUE)$derived$frequency_difference
  expect_equal(nrow(out), 3)
  expect_setequal(out$deidentified_sample, names(mif$spatial))
})
