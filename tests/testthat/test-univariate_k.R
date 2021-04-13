context("univariate_k")

load("example_data.RData")

x <- create_mif(clinical_data = example_clinical,
                sample_data = example_summary,
                spatial_list = example_spatial,
                patient_id = "deidentified_id",
                sample_id = "deidentified_sample",
                clean_columns = TRUE)

test_that("returns proper number of elements ", {
  
  mnames <- c("foxp3_opal_620_positive", "cd3_opal_570_positive", "cd8_opal_520_positive",
              "pd1_opal_650_positive", "pdl1_opal_540_positive")
  y <- suppressWarnings(ripleys_k(x, id = "deidentified_sample", mnames = mnames, num_permutations = 5))

  expect_equal(nrow(y), length(mnames)*length(x$spatial)*3)
})
