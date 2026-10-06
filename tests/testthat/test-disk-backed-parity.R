# The contract for disk-backing: a disk-backed mif produces the SAME NUMBERS as the
# in-memory mif it was made from. Nothing else about the feature matters if this
# fails, so these run every metric both ways and compare the derived table.
#
# Two things make this the real safety net rather than a formality:
#
# 1. `mif_spatial()` ignores its `columns` argument for an in-memory mif, so the
#    in-memory path is byte-identical to what it was before any of this existed.
#    The consequence is that a projection which omits a column a metric actually
#    uses fails ONLY on the disk path -- these tests are where that shows up.
#
# 2. Seeds are drawn in the parent (see the note in ripleys_k()), so permuted
#    results must match too, and must match at more than one `workers`. `workers = 2`
#    is the highest CRAN allows.

# Deliberately more than one sample, so `mclapply()` has something to distribute and
# any per-worker state leaks between them.
two_sample_mif <- function() {
  example_mif(which = c("TMA3_[9,K].tif", "TMA1_[3,B].tif"), n_cells = 300)
}

expect_same_derived <- function(a, b, slot, info = NULL) {
  expect_named(a$derived, names(b$derived), info = info)
  expect_equal(a$derived[[slot]], b$derived[[slot]], info = info)
}


test_that("ripleys_k agrees between in-memory and disk-backed mifs", {
  m <- two_sample_mif()
  d <- disk_mif(m)
  for (w in c(1, 2)) for (p in c(FALSE, TRUE)) {
    run <- function(x) {
      set.seed(11)
      ripleys_k(x, mnames = mnames_good()[1:2], r_range = seq(0, 50, 10),
                permute = p, num_permutations = 5, workers = w, overwrite = TRUE)
    }
    expect_same_derived(run(m), run(d), "univariate_Count",
                        info = sprintf("workers=%d permute=%s", w, p))
  }
})

test_that("bi_ripleys_k agrees between in-memory and disk-backed mifs", {
  m <- two_sample_mif()
  d <- disk_mif(m)
  for (w in c(1, 2)) {
    run <- function(x) {
      set.seed(12)
      bi_ripleys_k(x, mnames = mnames_bivariate(), r_range = seq(0, 50, 10),
                   permute = FALSE, workers = w, overwrite = TRUE)
    }
    expect_same_derived(run(m), run(d), "bivariate_Count", info = paste("workers", w))
  }
})

test_that("NN_G and bi_NN_G agree between in-memory and disk-backed mifs", {
  m <- two_sample_mif()
  d <- disk_mif(m)
  run_uni <- function(x) {
    set.seed(13)
    NN_G(x, mnames = mnames_good()[1:2], r_range = seq(0, 50, 10),
         num_permutations = 5, workers = 2, overwrite = TRUE)
  }
  expect_same_derived(run_uni(m), run_uni(d), "univariate_NN")

  run_bi <- function(x) {
    set.seed(14)
    bi_NN_G(x, mnames = mnames_bivariate(), r_range = seq(0, 50, 10),
            num_permutations = 5, workers = 2, overwrite = TRUE)
  }
  expect_same_derived(run_bi(m), run_bi(d), "bivariate_NN")
})

test_that("pair_correlation and bi_pair_correlation agree", {
  m <- two_sample_mif()
  d <- disk_mif(m)
  run_uni <- function(x) {
    set.seed(15)
    pair_correlation(x, mnames = mnames_good()[1:2], r_range = seq(0, 50, 10),
                     num_permutations = 5, workers = 2, overwrite = TRUE)
  }
  expect_same_derived(run_uni(m), run_uni(d), "univariate_pair_correlation")

  run_bi <- function(x) {
    set.seed(16)
    bi_pair_correlation(x, mnames = mnames_bivariate(), r_range = seq(0, 50, 10),
                        num_permutations = 5, workers = 2, overwrite = TRUE)
  }
  expect_same_derived(run_bi(m), run_bi(d), "bivariate_pair_correlation")
})

test_that("interaction_variable agrees between in-memory and disk-backed mifs", {
  m <- two_sample_mif()
  d <- disk_mif(m)
  run <- function(x) {
    set.seed(17)
    interaction_variable(x, mnames = mnames_bivariate(), r_range = seq(0, 50, 10),
                         num_permutations = 5, workers = 2, overwrite = TRUE)
  }
  a <- run(m); b <- run(d)
  expect_equal(a$derived[[names(a$derived)[1]]], b$derived[[names(b$derived)[1]]])
})

test_that("the bivariate metrics accept a data-frame mnames on a disk-backed mif", {
  # `mnames` may be a two-column anchor/counted data frame (marker_combinations()).
  # A data frame is a list, so building the projection with `c(..., mnames, ...)`
  # produced a LIST, which dplyr::all_of() rejects -- the disk path errored while the
  # in-memory path worked, because mif_spatial() ignores `columns` in memory.
  m <- two_sample_mif()
  d <- disk_mif(m)
  pairs <- data.frame(anchor  = c("CD8..Opal.520..Positive", "CD3..Opal.570..Positive"),
                      counted = c("FOXP3..Opal.620..Positive", "CD8..Opal.520..Positive"),
                      stringsAsFactors = FALSE)
  run <- function(x) {
    set.seed(21)
    bi_ripleys_k(x, mnames = pairs, r_range = seq(0, 50, 10), permute = FALSE,
                 workers = 1, overwrite = TRUE)
  }
  expect_same_derived(run(m), run(d), "bivariate_Count")

  run_i <- function(x) {
    set.seed(22)
    interaction_variable(x, mnames = pairs, r_range = seq(0, 50, 10),
                         num_permutations = 5, workers = 1, overwrite = TRUE)
  }
  a <- run_i(m); b <- run_i(d)
  expect_equal(a$derived[[names(a$derived)[1]]], b$derived[[names(b$derived)[1]]])
})

test_that("marker_freq_diff agrees between in-memory and disk-backed mifs", {
  # This one needs no coordinates at all, so its projection asks for none.
  m <- two_sample_mif()
  d <- disk_mif(m)
  run <- function(x) {
    marker_freq_diff(x, classifier = "Classifier.Label", ref_level = "Tumor",
                     diff_level = "Stroma", mnames = mnames_good()[1:2],
                     overwrite = TRUE)
  }
  a <- run(m); b <- run(d)
  expect_equal(a$derived[[names(a$derived)[1]]], b$derived[[names(b$derived)[1]]])
})

test_that("plot_tissue_split works on a disk-backed mif after split_tissue", {
  # The compartment column plot_tissue_split() colours by lives in the OVERLAY, so
  # the projection has to union base + overlay columns or it comes back missing.
  m <- halfplane_mif()
  d <- disk_mif(m)
  args <- list(classifier = "Classifier.Label", class1 = "Tumor",
               class2 = "Stroma", sigma = 40, interface_width = 20, workers = 1)
  sm <- do.call(split_tissue, c(list(m), args))
  sd <- do.call(split_tissue, c(list(d), args))
  for (comp in c("refined_density_compartment", "density_compartment")) {
    pm <- plot_tissue_split(sm, compartment = comp, workers = 1)
    pd <- plot_tissue_split(sd, compartment = comp, workers = 1)
    expect_length(pd, length(pm))
    expect_s3_class(pd[[1]], "ggplot")
    expect_equal(pd[[1]]$data, pm[[1]]$data, info = comp)
  }
})

test_that("dixons_s agrees between in-memory and disk-backed mifs", {
  m <- toy_mif(list(S1 = toy_spatial("S1", n = 120, markers = c(A = 30, B = 30))))
  d <- disk_mif(m)
  run <- function(x) {
    set.seed(18)
    dixons_s(x, mnames = c("A", "B"), num_permutations = 20, workers = 1,
             overwrite = TRUE)
  }
  a <- run(m); b <- run(d)
  expect_equal(a$derived$Dixon_Z, b$derived$Dixon_Z)
  expect_equal(a$derived$Dixon_C, b$derived$Dixon_C)
})

test_that("split_tissue agrees, including the columns it writes back", {
  m <- halfplane_mif()
  d <- disk_mif(m)
  run <- function(x) {
    split_tissue(x, classifier = "Classifier.Label", class1 = "Tumor",
                 class2 = "Stroma", sigma = 40, interface_width = 20,
                 workers = 1, overwrite = TRUE)
  }
  a <- run(m); b <- run(d)

  expect_equal(a$sample, b$sample)
  expect_equal(a$derived$density_boundary, b$derived$density_boundary)
  # The three written columns come back through the overlay on the disk side.
  sa <- a$spatial[[1]]
  sb <- b$spatial[[1]]
  for (col in c("density_compartment", "refined_density_compartment",
                "density_score")) {
    expect_equal(sb[[col]], sa[[col]], info = col)
    # Factors, not characters -- levels and their order have to survive parquet.
    expect_equal(levels(sb[[col]]), levels(sa[[col]]), info = col)
  }
  # And the base columns are untouched by the overlay.
  expect_equal(sb[colnames(m$spatial[[1]])], sa[colnames(m$spatial[[1]])])
})

test_that("a disk-backed mif survives saveRDS and still computes", {
  m <- two_sample_mif()
  d <- disk_mif(m)
  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(d, f)
  d2 <- readRDS(f)

  expect_s3_class(d2$spatial, "mif_store")
  expect_equal(names(d2$spatial), names(d$spatial))
  set.seed(19)
  a <- ripleys_k(m, mnames = mnames_good()[1], r_range = seq(0, 30, 10),
                 permute = FALSE, workers = 1, overwrite = TRUE)
  set.seed(19)
  b <- ripleys_k(d2, mnames = mnames_good()[1], r_range = seq(0, 30, 10),
                 permute = FALSE, workers = 1, overwrite = TRUE)
  expect_equal(a$derived$univariate_Count, b$derived$univariate_Count)
})

test_that("the v1.4.0 numeric baseline holds through a disk-backed mif", {
  # Same assertion as test-regression-vs-baseline.R makes for the in-memory path,
  # so disk-backing is pinned to v1.4.0's numbers rather than merely to whatever the
  # current in-memory code happens to produce.
  path <- test_path("fixtures", "baseline-v1.4.0.rds")
  skip_if_not(file.exists(path), "v1.4.0 baseline fixture not present")
  b <- readRDS(path)
  mk <- attr(b, "markers"); rr <- attr(b, "r_range")

  mem <- create_mif(
    clinical_data = example_clinical %>%
      dplyr::mutate(deidentified_id = as.character(deidentified_id)),
    sample_data = example_summary %>%
      dplyr::mutate(deidentified_id = as.character(deidentified_id)),
    spatial_list = example_spatial[attr(b, "sample")],
    patient_id = "deidentified_id", sample_id = "deidentified_sample")

  for (ec in c("translation", "isotropic")) {
    key <- c(translation = "ripleys_k.exact", isotropic = "ripleys_k.iso")[[ec]]
    old <- b[[key]]$univariate_Count
    old <- old[order(old$Marker, old$r), ]
    new <- ripleys_k(disk_mif(mem), mnames = mk, r_range = rr, permute = FALSE,
                     workers = 1, edge_correction = ec,
                     overwrite = TRUE)$derived$univariate_Count
    new <- new[order(new$Marker, new$r), ]
    expect_equal(new$`Observed K`, old$`Observed K`, tolerance = 1e-9, info = ec)
    expect_equal(new$`Theoretical CSR`, old$`Theoretical CSR`, info = ec)
    expect_equal(new$`Exact CSR`, old$`Exact CSR`, tolerance = 1e-9, info = ec)
  }
})


test_that("every metric's derived slot is byte-identical after a store round trip", {
  # The parity tests above compare the mif a metric RETURNS. This one compares what
  # survives to disk, which until 2.0.0 was nothing: write_derived() mutated
  # mif$derived in memory only, so open_mif() after a metric run found an empty
  # `derived` and the work was silently gone.
  dir <- withr::local_tempdir()
  root <- file.path(dir, "s.mif")
  m <- two_sample_mif()
  d <- mif_to_disk(m, root)

  rr <- seq(0, 30, 10)
  runs <- list(
    univariate_Count = function(x) ripleys_k(x, mnames = mnames_good()[1], r_range = rr,
                                             permute = FALSE, workers = 1, overwrite = TRUE),
    univariate_NN    = function(x) NN_G(x, mnames = mnames_good()[1], r_range = rr,
                                       num_permutations = 3, workers = 1, overwrite = TRUE),
    bivariate_Count  = function(x) bi_ripleys_k(x, mnames = mnames_bivariate(), r_range = rr,
                                                permute = FALSE, workers = 1, overwrite = TRUE)
  )

  for (slot in names(runs)) {
    set.seed(31)
    d <- runs[[slot]](d)
    reopened <- open_mif(root)
    expect_identical(reopened$derived[[slot]], d$derived[[slot]], info = slot)
    # Slot order is the manifest's, not list.files()' alphabetical one.
    expect_identical(names(reopened$derived), names(d$derived), info = slot)
  }

  # And the numbers still match the in-memory mif they came from.
  set.seed(31)
  mem <- runs$univariate_Count(m)
  expect_equal(open_mif(root)$derived$univariate_Count,
               mem$derived$univariate_Count)
})


test_that("the store tracks whichever mif last ran a metric against it", {
  # A metric writes to the store as well as to the mif it returns, so two mifs
  # pointing at one directory are not independent: `a <- ripleys_k(xd, ...)` leaves
  # the in-memory xd alone, as copy-on-modify implies, but the store now matches `a`.
  # This is the one place a disk-backed mif departs from ordinary R semantics, and it
  # is documented in ?mif_to_disk and the vignette. Pinned here because it is a
  # consequence of write_derived() syncing, not an intended feature in itself -- if
  # the sync is ever made conditional this test says what changes.
  dir <- withr::local_tempdir()
  root <- file.path(dir, "s.mif")
  xd <- mif_to_disk(two_sample_mif(), root)
  mk <- mnames_good()[1:2]

  xd <- ripleys_k(xd, mnames = mk, r_range = seq(0, 40, 10), permute = FALSE,
                  workers = 1, overwrite = TRUE)
  # Assigned back to the same name: object and store agree.
  expect_identical(open_mif(root)$derived$univariate_Count,
                   xd$derived$univariate_Count)

  # Assigned elsewhere: the store follows the RUN, not the variable.
  set.seed(42)
  a <- ripleys_k(xd, mnames = mk[1], r_range = seq(0, 20, 10), permute = TRUE,
                 num_permutations = 5, workers = 1, overwrite = TRUE)
  stored <- open_mif(root)$derived$univariate_Count
  expect_identical(stored, a$derived$univariate_Count)
  expect_false(identical(stored, xd$derived$univariate_Count))
  # And the in-memory xd really is untouched, rather than merely different.
  expect_equal(nrow(xd$derived$univariate_Count),
               length(mk) * 5L * 2L)   # 2 markers x 5 radii x 2 samples

  # Only the slot being written is replaced; inherited slots carry forward.
  xd <- NN_G(xd, mnames = mk[1], r_range = seq(0, 20, 10), num_permutations = 2,
             workers = 1, overwrite = TRUE)
  reopened <- open_mif(root)
  expect_setequal(names(reopened$derived), c("univariate_Count", "univariate_NN"))
  expect_identical(reopened$derived$univariate_NN, xd$derived$univariate_NN)
})
