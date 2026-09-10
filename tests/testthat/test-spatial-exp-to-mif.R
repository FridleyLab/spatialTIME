skip_if_not_installed("SpatialExperiment")
skip_if_not_installed("SummarizedExperiment")
skip_if_not_installed("S4Vectors")

# A SpatialExperiment built by hand. VectraPolarisData is deliberately NOT used:
# it is a ~GB Bioconductor data package and must never become a test dependency.

toy_spe <- function(n = 60, n_samples = 2) {
  ids <- rep(paste0("S", seq_len(n_samples)), length.out = n)
  set.seed(3)
  cd <- data.frame(
    sample_id      = ids,
    phenotype_cd3  = ifelse(rbinom(n, 1, 0.4) == 1, "CD3+", "CD3-"),
    phenotype_cd8  = ifelse(rbinom(n, 1, 0.3) == 1, "CD8+", "CD8-"),
    extra          = seq_len(n),
    stringsAsFactors = FALSE
  )
  coords <- cbind(cell_x_position = runif(n, 0, 100),
                  cell_y_position = runif(n, 0, 100))
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = matrix(0, nrow = 1, ncol = n)),
    colData = cd,
    spatialCoords = coords
  )
  S4Vectors::metadata(spe)$clinical_data <-
    data.frame(sample_id = unique(ids), stringsAsFactors = FALSE)
  spe
}

test_that("spatial_exp_to_mif converts a SpatialExperiment to a mif", {
  spe <- toy_spe()
  mif <- spatial_exp_to_mif(spe, markers = c("phenotype_cd3", "phenotype_cd8"),
                            patient_id = "sample_id", sample_id = "sample_id")
  expect_s3_class(mif, "mif")
  expect_length(mif$spatial, 2)
  # Marker columns are recoded to 0/1 by the +-regex.
  for (sp in mif$spatial) {
    for (m in c("phenotype_cd3", "phenotype_cd8")) {
      expect_true(all(sp[[m]] %in% c(0, 1)), info = m)
    }
  }
})

test_that("marker recoding follows marker_pos_regex", {
  spe <- toy_spe()
  cd <- SummarizedExperiment::colData(spe)
  expected <- sum(grepl("\\+", cd$phenotype_cd3))
  mif <- spatial_exp_to_mif(spe, markers = "phenotype_cd3",
                            patient_id = "sample_id", sample_id = "sample_id")
  expect_equal(sum(vapply(mif$spatial, function(x) sum(x$phenotype_cd3), numeric(1))),
               expected)
})

test_that("cols_to_keep is honoured", {
  spe <- toy_spe()
  mif <- spatial_exp_to_mif(spe, markers = "phenotype_cd3",
                            patient_id = "sample_id", sample_id = "sample_id",
                            cols_to_keep = "extra")
  expect_true("extra" %in% names(mif$spatial[[1]]))
})

test_that("a non-SpatialExperiment is rejected via inherits, not class ==", {
  # class(x) == "SpatialExperiment" errors or warns on a multi-element class vector;
  # inherits() is the correct test and also admits subclasses.
  expect_error(spatial_exp_to_mif(list(), markers = "m"), "class 'SpatialExperiment'")
  expect_error(spatial_exp_to_mif(data.frame(a = 1), markers = "m"),
               "class 'SpatialExperiment'")
})

test_that("bad column names are rejected with a useful message", {
  spe <- toy_spe()
  expect_error(
    spatial_exp_to_mif(spe, markers = "not_a_marker",
                       patient_id = "sample_id", sample_id = "sample_id"),
    "columns name in colData"
  )
  expect_error(
    spatial_exp_to_mif(spe, markers = "phenotype_cd3", x_coord = "nope",
                       patient_id = "sample_id", sample_id = "sample_id"),
    "x_coord"
  )
})
