library(spatialTIME)
library(dplyr)

example_clinical_new = example_clinical %>% 
  mutate(deidentified_id = as.character(deidentified_id))
example_summary_new = example_summary %>% 
  mutate(deidentified_id = as.character(deidentified_id))
mif = create_mif(clinical_data = example_clinical_new,
                 sample_data = example_summary_new,
                 spatial_list = example_spatial,
                 patient_id = "deidentified_id", 
                 sample_id = "deidentified_sample")

test_that("create_mif creates a MIF object", {
  expect_equal(length(mif), 6)
  expect_equal(mif$clinical, example_clinical_new)
  expect_equal(mif$sample, example_summary_new)
  expect_equal(mif$spatial, example_spatial)
})
