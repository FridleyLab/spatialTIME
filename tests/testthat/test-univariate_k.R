context("univariate_k")

load("example_data.RData")

test_that("returns proper number of elements ", {
  # grabbing image tag names to make sample and clinical file -----
  spatial_names <- lapply(example_tma, function(x) {x$image.tag[[1]]})
  spatial_names <- unlist(spatial_names)
  # spatial_names <- gsub("\\.tif", "", spatial_names)
  
  example_tma <- lapply(example_tma, function(x){
    x <- x %>% janitor::clean_names()
  })
  
  names(example_tma) <- spatial_names
  
  set.seed(8675309)
  example_sample <- data.frame(image_tag = spatial_names,
                               patient_id = sample(c("patient_x", "patient_y", "patient_z"),
                                                   10, replace = TRUE))
  
  example_clinical <- data.frame(patient_id = c("patient_x", "patient_y", "patient_z"),
                                 covar_one = c("low", "high", "low"),
                                 covar_two = rnorm(3))
  
  x <- create_mif(spatial_list = example_tma, 
                  clinical_data = example_clinical, 
                  sample_data = example_sample,
                  clean_columns = TRUE)
  
  mnames <- c("foxp3_opal_620_positive", "cd3_opal_570_positive", "cd8_opal_520_positive",
              "pd1_opal_650_positive", "pdl1_opal_540_positive")
  
  y <- ripleys_k(x, id = "image_tag", mnames = mnames)
  
  expect_equal(nrow(y), length(mnames)*length(x$spatial)*3)
})
