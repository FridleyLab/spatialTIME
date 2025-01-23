library(spatialTIME)
library(dplyr)

mif <- create_mif(
  clinical_data = example_clinical,
  sample_data = example_summary,
  spatial_list = example_spatial,
  patient_id = "deidentified_id",
  sample_id = "deidentified_sample"
)

test_that("create_mif creates a MIF object", {
  expect_equal(length(mif), 6)
  expect_equal(mif$spatial, example_spatial)
})
