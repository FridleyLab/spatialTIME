# The store itself: round-trip fidelity, the row-order guarantee the whole design
# rests on, manifest validation, and the compatibility surface that keeps third-party
# code (mxfda, scSpatialSIM) and the package's own vignettes working.

store_of <- function(mif = example_mif()) {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  mif_to_disk(mif, file.path(dir, "s.mif"))
}

# `example_mif()` thins with `x[sort(sample(nrow(x), n)), ]`, which leaves row names
# like 19, 29, 86. A columnar store has nowhere to keep those and nothing in the
# package reads them (see the "What is and is not preserved" section of
# ?mif_to_disk), so comparisons of spatial data normalise them on both sides. The
# normalisation itself is asserted separately, below, so it cannot drift silently
# into hiding a real difference.
reset_rn <- function(x) {
  if (is.data.frame(x)) {
    rownames(x) <- NULL
    return(x)
  }
  lapply(x, reset_rn)
}


# ---------------------------------------------------------------------------
# Row order. This is the one that matters most.
# ---------------------------------------------------------------------------

test_that("read_parquet is order-stable across column subsets", {
  # Everything in this package is positional: k_from_pairs() indexes marker masks
  # into a pair list built in file order, and split_tissue() assigns columns back by
  # position. If coordinates came back in one order and a marker column in another,
  # markers would attach to the wrong cells and the numbers would be plausible and
  # wrong. Asserted rather than assumed so that an arrow upgrade which changes the
  # guarantee fails here instead of silently corrupting results.
  spat <- toy_spatial("S1", n = 5000, markers = c(A = 900, B = 700))
  f <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(spat, f)

  full <- as.data.frame(arrow::read_parquet(f))
  expect_identical(full$XMin, spat$XMin)

  one <- as.data.frame(arrow::read_parquet(f, col_select = dplyr::all_of("XMin")))
  two <- as.data.frame(arrow::read_parquet(f,
    col_select = dplyr::all_of(c("XMin", "A"))))
  expect_identical(one$XMin, full$XMin)
  expect_identical(two$XMin, full$XMin)
  expect_identical(two$A, full$A)

  # Repeated reads agree with each other too.
  expect_identical(as.data.frame(arrow::read_parquet(f,
    col_select = dplyr::all_of("A")))$A, full$A)
})

test_that("open_dataset does NOT preserve row order, which is why we never use it", {
  # Documents the hazard this design exists to avoid. If a future arrow makes
  # open_dataset order-stable this test fails, and that is the right moment to
  # reconsider the reader rule in R/utils-mif-store.R -- deliberately, not by
  # accident.
  # Reordering happens at multiples of arrow's 2^15 read-batch size and needs a file
  # big enough to be read in several batches -- 600k x 2 columns is 7 MB and ~0.03 s
  # to write, and reproduces. It is scheduling-dependent, not guaranteed (1M x 4
  # columns did not reproduce in the same session), so a negative result skips rather
  # than fails: a false alarm here would be worse than a missed canary, and the rule
  # itself is enforced by the test above.
  n <- 600000L
  d <- data.frame(seq = seq_len(n), val = seq_len(n) / 2)
  f <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(d, f)

  full <- as.data.frame(arrow::read_parquet(f))
  expect_identical(full$seq, d$seq)
  ds <- as.data.frame(dplyr::collect(dplyr::select(arrow::open_dataset(f), "seq")))
  # Same cells either way -- only the order is in question, which is what makes it
  # silent corruption rather than a visible error.
  expect_setequal(ds$seq, full$seq)
  if (identical(ds$seq, full$seq)) {
    skip("this arrow build returned Dataset rows in file order; reader rule unchanged")
  }
  expect_false(identical(ds$seq, full$seq))
  # The divergence lands on a read-batch boundary, which is the mechanism.
  expect_equal(which(ds$seq != full$seq)[1] - 1L, 32768L * 18L, tolerance = 32768)
})

test_that("a projected read lines up with the full read row for row", {
  d <- store_of()
  full <- d$spatial[[1]]
  proj <- mif_spatial(d, 1, c("XMin", "YMin", mnames_good()[1]))
  expect_identical(proj$XMin, full$XMin)
  expect_identical(proj[[mnames_good()[1]]], full[[mnames_good()[1]]])
})


test_that("no package code resizes arrow's thread pool", {
  # Resizing arrow's thread pool from inside an mclapply child DEADLOCKS on arrow
  # 25.0.1 / R 4.6.1: the run hangs at 0% CPU with no error. An earlier version of
  # pq_read() did exactly that, and it cost nothing to lose -- pinning measured
  # slower than arrow's default anyway.
  #
  # Checked by reading the sources rather than by forking-and-hoping, because a test
  # that reproduces the bug would HANG the suite (and a CRAN check) rather than fail
  # it. If you need to limit arrow's threads, do it in the parent before calling a
  # metric; children inherit that safely.
  r_dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "sources not available (installed package)")
  src <- unlist(lapply(list.files(r_dir, pattern = "[.]R$", full.names = TRUE),
                       readLines, warn = FALSE))
  code <- grep("^\\s*#", src, value = TRUE, invert = TRUE)
  expect_length(grep("set_cpu_count|set_io_thread_count", code), 0)
})

test_that("reading a store inside an mclapply fork works and leaves arrow alone", {
  # Reading parquet in a forked child is safe -- verified on arrow 23.0.1.2 and
  # 25.0.1 -- so this one does fork. It is the thread-pool resize above, not the
  # read, that was the hazard.
  d <- store_of(example_mif(which = c("TMA3_[9,K].tif", "TMA1_[3,B].tif"),
                            n_cells = 120))
  before <- c(arrow::cpu_count(), arrow::io_thread_count())
  res <- parallel::mclapply(seq_len(length(d$spatial)), function(i) {
    list(n = nrow(mif_spatial(d, i, "XMin")), cpu = arrow::cpu_count())
  }, mc.cores = 2, mc.preschedule = FALSE)

  expect_false(any(vapply(res, inherits, logical(1), "try-error")))
  expect_equal(vapply(res, function(x) x$n, numeric(1)),
               unname(attr(d$spatial, "nrow")))
  # Children inherited the parent's setting, and the parent's is unchanged.
  expect_equal(vapply(res, function(x) x$cpu, numeric(1)),
               rep(before[1], length(res)))
  expect_equal(c(arrow::cpu_count(), arrow::io_thread_count()), before)
})


# ---------------------------------------------------------------------------
# Round-trip
# ---------------------------------------------------------------------------

test_that("collect_mif inverts mif_to_disk exactly", {
  m <- example_mif(which = c("TMA3_[9,K].tif", "TMA1_[3,B].tif"), n_cells = 200)
  back <- collect_mif(store_of(m))
  expect_identical(reset_rn(back$spatial), reset_rn(m$spatial))
  expect_identical(back$clinical, m$clinical)
  expect_identical(back$sample, m$sample)
  expect_identical(back$patient_id, m$patient_id)
  expect_identical(back$sample_id, m$sample_id)
})

test_that("row names are normalised to 1:nrow, and that is the only difference", {
  # Documents the one documented fidelity gap (?mif_to_disk, "What is and is not
  # preserved"), so reset_rn() above cannot quietly be hiding anything else.
  m <- example_mif(which = "TMA3_[9,K].tif", n_cells = 200)
  expect_false(identical(rownames(m$spatial[[1]]),
                         as.character(seq_len(nrow(m$spatial[[1]])))))
  back <- collect_mif(store_of(m))$spatial[[1]]
  expect_identical(rownames(back), as.character(seq_len(nrow(back))))
  # Everything except row names is untouched.
  expect_equal(back, m$spatial[[1]], ignore_attr = "row.names")
  expect_identical(colnames(back), colnames(m$spatial[[1]]))
})

test_that("factor columns and their levels survive the round-trip", {
  spat <- toy_spatial("S1", n = 100, markers = c(A = 20))
  # An unused level and an NA, the two cases that silently degrade.
  spat$f <- factor(rep(c("x", "y", NA, "x"), length.out = 100),
                   levels = c("y", "x", "unused"))
  m <- toy_mif(list(S1 = spat))
  back <- collect_mif(store_of(m))
  expect_identical(back$spatial$S1$f, spat$f)
  expect_identical(levels(back$spatial$S1$f), c("y", "x", "unused"))
})

test_that("collect_mif can take a subset, by name or position", {
  m <- example_mif(which = c("TMA3_[9,K].tif", "TMA1_[3,B].tif"), n_cells = 100)
  d <- store_of(m)
  expect_named(collect_mif(d, samples = "TMA1_[3,B].tif")$spatial, "TMA1_[3,B].tif")
  expect_length(collect_mif(d, samples = 1)$spatial, 1)
  expect_error(collect_mif(d, samples = "nope"), "No sample named")
})

test_that("open_mif reopens a store faithfully, including derived slots", {
  m <- example_mif(n_cells = 150)
  set.seed(4)
  m <- ripleys_k(m, mnames = mnames_good()[1], r_range = seq(0, 20, 10),
                 permute = FALSE, workers = 1, overwrite = TRUE)
  dir <- withr::local_tempdir()
  root <- file.path(dir, "s.mif")
  mif_to_disk(m, root)
  o <- open_mif(root)
  expect_s3_class(o, "mif")
  expect_identical(o$derived$univariate_Count, m$derived$univariate_Count)
  expect_identical(reset_rn(collect_mif(o)$spatial), reset_rn(m$spatial))
})

test_that("a moved store still opens, because paths are relative", {
  dir <- withr::local_tempdir()
  m <- example_mif(n_cells = 100)
  mif_to_disk(m, file.path(dir, "a.mif"))
  file.rename(file.path(dir, "a.mif"), file.path(dir, "b.mif"))
  o <- open_mif(file.path(dir, "b.mif"))
  expect_identical(reset_rn(collect_mif(o)$spatial), reset_rn(m$spatial))
})


# ---------------------------------------------------------------------------
# Validation: every failure names the sample and the discrepancy
# ---------------------------------------------------------------------------

test_that("mif_to_disk requires an explicit path and refuses to clobber", {
  m <- example_mif(n_cells = 50)
  expect_error(mif_to_disk(m), "must name the directory")
  dir <- withr::local_tempdir()
  root <- file.path(dir, "s.mif")
  mif_to_disk(m, root)
  expect_error(mif_to_disk(m, root), "already exists")
  expect_no_error(mif_to_disk(m, root, overwrite = TRUE))
})

test_that("open_mif rejects a directory that is not a store", {
  dir <- withr::local_tempdir()
  expect_error(open_mif(file.path(dir, "nope")), "No such directory")
  expect_error(open_mif(dir), "not a mif store")
})

test_that("open_mif detects a missing spatial file", {
  dir <- withr::local_tempdir()
  root <- file.path(dir, "s.mif")
  d <- mif_to_disk(example_mif(n_cells = 50), root)
  unlink(store_paths(d$spatial, 1))
  expect_error(open_mif(root), "Store is incomplete")
})

test_that("open_mif detects a file whose size no longer matches", {
  dir <- withr::local_tempdir()
  root <- file.path(dir, "s.mif")
  d <- mif_to_disk(example_mif(n_cells = 50), root)
  # Rewrite with fewer rows: a plausible truncation, still valid parquet.
  arrow::write_parquet(d$spatial[[1]][1:10, ], store_paths(d$spatial, 1))
  expect_error(open_mif(root), "has changed on disk")
})

test_that("open_mif detects a row count that disagrees with the manifest", {
  dir <- withr::local_tempdir()
  root <- file.path(dir, "s.mif")
  d <- mif_to_disk(example_mif(n_cells = 50), root)
  mp <- file.path(root, "manifest.json")
  man <- jsonlite::fromJSON(mp, simplifyVector = FALSE)
  man$samples[[1]]$nrow <- 999999
  writeLines(jsonlite::toJSON(man, auto_unbox = TRUE, pretty = TRUE, null = "null"), mp)
  expect_error(open_mif(root), "rows, file has")
})

test_that("open_mif refuses a future format version", {
  dir <- withr::local_tempdir()
  root <- file.path(dir, "s.mif")
  mif_to_disk(example_mif(n_cells = 50), root)
  mp <- file.path(root, "manifest.json")
  man <- jsonlite::fromJSON(mp, simplifyVector = FALSE)
  man$format <- 99L
  writeLines(jsonlite::toJSON(man, auto_unbox = TRUE, pretty = TRUE, null = "null"), mp)
  expect_error(open_mif(root), "format version 99")
})

test_that("asking for a column that does not exist says so, and names the sample", {
  d <- store_of()
  expect_error(mif_spatial(d, 1, "not_a_column"), "not found in sample")
})

test_that("create_mif validates a character spatial_list", {
  m <- example_mif(n_cells = 50)
  d <- store_of(m)
  p <- store_paths(d$spatial, 1)
  clin <- data.frame(pid = "p1")
  samp <- data.frame(pid = "p1", sid = names(d$spatial)[1],
                     not_there = names(d$spatial)[1])

  expect_error(
    create_mif(clin, samp, spatial_list = c(a = file.path(tempdir(), "no.parquet")),
               patient_id = "pid", sample_id = "sid"),
    "not found")
  expect_error(
    create_mif(clin, samp, spatial_list = character(0),
               patient_id = "pid", sample_id = "sid"),
    "empty character vector")
  csv <- withr::local_tempfile(fileext = ".csv")
  writeLines("a,b", csv)
  expect_error(
    create_mif(clin, samp, spatial_list = c(a = csv),
               patient_id = "pid", sample_id = "sid"),
    "must name parquet files")
  # sample_id must actually be a column of the file.
  expect_error(
    create_mif(clin, samp, spatial_list = stats::setNames(p, "x"),
               patient_id = "pid", sample_id = "not_there"),
    "has no column named")
})

test_that("reference mode points at existing files without copying them", {
  m <- example_mif(which = c("TMA3_[9,K].tif", "TMA1_[3,B].tif"), n_cells = 120)
  d <- store_of(m)
  paths <- store_paths(d$spatial)
  ref <- create_mif(
    clinical_data = data.frame(deidentified_id = "p1"),
    sample_data = data.frame(deidentified_id = "p1",
                             deidentified_sample = names(d$spatial)),
    spatial_list = stats::setNames(paths, names(d$spatial)),
    patient_id = "deidentified_id", sample_id = "deidentified_sample")

  expect_s3_class(ref$spatial, "mif_store")
  expect_true(is.na(attr(ref$spatial, "root")))
  expect_identical(reset_rn(collect_mif(ref)$spatial), reset_rn(m$spatial))
  set.seed(5)
  a <- ripleys_k(m, mnames = mnames_good()[1], r_range = seq(0, 20, 10),
                 permute = FALSE, workers = 1, overwrite = TRUE)
  set.seed(5)
  b <- ripleys_k(ref, mnames = mnames_good()[1], r_range = seq(0, 20, 10),
                 permute = FALSE, workers = 1, overwrite = TRUE)
  expect_equal(b$derived$univariate_Count, a$derived$univariate_Count)
})

test_that("an unnamed reference vector takes names from the sample_id column", {
  m <- example_mif(which = "TMA3_[9,K].tif", n_cells = 120)
  d <- store_of(m)
  ref <- create_mif(
    clinical_data = data.frame(deidentified_id = "p1"),
    sample_data = data.frame(deidentified_id = "p1",
                             deidentified_sample = names(d$spatial)),
    spatial_list = unname(store_paths(d$spatial)),
    patient_id = "deidentified_id", sample_id = "deidentified_sample")
  expect_identical(names(ref$spatial), "TMA3_[9,K].tif")
})

test_that("a file holding two samples is refused when names must be derived", {
  # Every metric takes the window to be the convex hull of all of a file's cells, so
  # two samples in one file is silently wrong rather than merely untidy.
  two <- rbind(toy_spatial("S1", n = 60), toy_spatial("S2", n = 60))
  f <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(two, f)
  expect_error(
    create_mif(data.frame(deidentified_id = "p1"),
               data.frame(deidentified_id = "p1", deidentified_sample = "S1"),
               spatial_list = f,
               patient_id = "deidentified_id", sample_id = "deidentified_sample"),
    "contains 2 distinct values")
})


# ---------------------------------------------------------------------------
# Compatibility: the surface CRAN users and reverse dependencies already use
# ---------------------------------------------------------------------------

test_that("a disk-backed mif keeps the same six slots", {
  d <- store_of()
  expect_named(d, c("clinical", "sample", "spatial", "derived",
                    "patient_id", "sample_id"))
  expect_s3_class(d, "mif")
})

test_that("length, names and [[ work on the spatial slot as they always did", {
  m <- example_mif(which = c("TMA3_[9,K].tif", "TMA1_[3,B].tif"), n_cells = 100)
  d <- store_of(m)
  # The vignettes do exactly this (deriving_functions.Rmd:31, intro.Rmd:361), and
  # mxfda/scSpatialSIM index the slot directly.
  expect_equal(length(d$spatial), length(m$spatial))
  expect_identical(names(d$spatial), names(m$spatial))
  expect_identical(colnames(d$spatial[[1]]), colnames(m$spatial[[1]]))
  expect_identical(reset_rn(d$spatial[["TMA1_[3,B].tif"]]),
                   reset_rn(m$spatial[["TMA1_[3,B].tif"]]))
  expect_identical(table(d$spatial[[1]]$Classifier.Label),
                   table(m$spatial[[1]]$Classifier.Label))
  expect_error(d$spatial[["nope"]], "No sample named")
})

test_that("print methods work without reading spatial data", {
  d <- store_of()
  expect_output(print(d$spatial), "mif_store")
  expect_output(print(d$spatial), "cells")
  expect_output(print(d), "spatial data frames were found")
})

test_that("subsetting a store keeps it a store", {
  m <- example_mif(which = c("TMA3_[9,K].tif", "TMA1_[3,B].tif"), n_cells = 100)
  d <- store_of(m)
  s <- d$spatial[1]
  expect_s3_class(s, "mif_store")
  expect_length(s, 1)
  expect_identical(reset_rn(s[[1]]), reset_rn(m$spatial[[1]]))
})

test_that("subscripting a store accepts names, positions, negatives and logicals", {
  # The logical case is why store_subscript() exists: `as.integer(c(FALSE, TRUE))` is
  # `c(0L, 1L)`, and R drops the 0, so a naive conversion returned sample 1 for
  # `store[c(FALSE, TRUE)]` -- wrong data with no error.
  m <- example_mif(which = c("TMA3_[9,K].tif", "TMA1_[3,B].tif"), n_cells = 80)
  s <- store_of(m)$spatial
  nms <- names(s)
  for (idx in list(1, 2, c(1, 2), -1, -2, c(TRUE, FALSE), c(FALSE, TRUE),
                   c(TRUE, TRUE))) {
    expect_identical(names(s[idx]), nms[idx],
                     info = paste(class(idx), paste(idx, collapse = ",")))
  }
  # Names select themselves.
  expect_identical(names(s[nms[2]]), nms[2])
  expect_identical(names(s[nms]), nms)
  expect_identical(names(s[rev(nms)]), rev(nms))
  expect_error(s[5], "out of range")
  expect_error(s["nope"], "No sample named")
})

test_that("collect_mif takes the same subscript kinds as the store", {
  m <- example_mif(which = c("TMA3_[9,K].tif", "TMA1_[3,B].tif"), n_cells = 80)
  d <- store_of(m)
  expect_named(collect_mif(d, samples = c(FALSE, TRUE))$spatial, names(m$spatial)[2])
  expect_named(collect_mif(d, samples = -1)$spatial, names(m$spatial)[2])
  expect_named(collect_mif(d, samples = 2)$spatial, names(m$spatial)[2])
})

test_that("sample names that sanitise to the same stem get distinct files", {
  # "A/1" and "A 1" both sanitise to "A_1", so without disambiguation one sample's
  # parquet would overwrite the other's while the manifest listed both.
  sp <- list(toy_spatial("A/1", n = 40), toy_spatial("A 1", n = 50))
  names(sp) <- c("A/1", "A 1")
  m <- toy_mif(sp)
  d <- store_of(m)
  expect_length(unique(store_paths(d$spatial)), 2)
  expect_identical(nrow(d$spatial[["A/1"]]), 40L)
  expect_identical(nrow(d$spatial[["A 1"]]), 50L)

  # And the same holds for the subset store, and for overlay files.
  dir <- withr::local_tempdir()
  out <- subset_mif(d, classifier = "Classifier.Label", level = "Tumor",
                    markers = "A", path = file.path(dir, "sub.mif"))
  expect_length(unique(store_paths(out$spatial)), length(out$spatial))
  expect_no_error(open_mif(file.path(dir, "sub.mif")))
})

test_that("overlay columns are visible to a projected read", {
  m <- halfplane_mif()
  dir <- withr::local_tempdir()
  d <- mif_to_disk(m, file.path(dir, "s.mif"))
  sd <- split_tissue(d, classifier = "Classifier.Label", class1 = "Tumor",
                     class2 = "Stroma", sigma = 40, interface_width = 20,
                     workers = 1)
  # Asking for an overlay column by name must return it, not drop it.
  got <- mif_spatial(sd, 1, c("XMin", "refined_density_compartment"))
  expect_named(got, c("XMin", "refined_density_compartment"))
  expect_true(is.factor(got$refined_density_compartment))
  # And spatial_columns() must keep it in the projection it builds.
  expect_true("refined_density_compartment" %in%
    spatial_columns(sd, extra = "refined_density_compartment", i = 1))
  expect_true("refined_density_compartment" %in% mif_spatial_colnames(sd, 1))
})

test_that("two disk-backed mifs merge, and mixing modes is refused", {
  m1 <- example_mif(which = "TMA3_[9,K].tif", n_cells = 100)
  m2 <- example_mif(which = "TMA1_[3,B].tif", n_cells = 100)
  d1 <- store_of(m1); d2 <- store_of(m2)

  merged <- merge_mifs(list(d1, d2), check.names = TRUE)
  expect_s3_class(merged$spatial, "mif_store")
  expect_length(merged$spatial, 2)
  expect_identical(reset_rn(merged$spatial[[2]]), reset_rn(m2$spatial[[1]]))

  expect_error(merge_mifs(list(d1, m2), check.names = TRUE),
               "disk-backed and in-memory")
})

test_that("subset_mif on a disk-backed mif needs a path and writes a store", {
  m <- example_mif(which = "TMA3_[9,K].tif", n_cells = 300)
  d <- store_of(m)
  expect_error(
    subset_mif(d, classifier = "Classifier.Label", level = "Tumor",
               markers = mnames_good()[1]),
    "needs a `path`")

  dir <- withr::local_tempdir()
  out <- subset_mif(d, classifier = "Classifier.Label", level = "Tumor",
                    markers = mnames_good()[1], path = file.path(dir, "sub.mif"))
  expect_s3_class(out$spatial, "mif_store")
  ref <- subset_mif(m, classifier = "Classifier.Label", level = "Tumor",
                    markers = mnames_good()[1])
  expect_identical(names(out$spatial), names(ref$spatial))
  expect_identical(reset_rn(collect_mif(out)$spatial), reset_rn(ref$spatial))
  expect_equal(out$sample, ref$sample)
  # And it is a real store, reopenable on its own.
  expect_no_error(open_mif(file.path(dir, "sub.mif")))
})

test_that("the store index stays small and does not grow with cohort size", {
  # The point of the exercise: the parent process holds metadata, not cells.
  m <- example_mif(which = names(example_spatial), n_cells = NULL)
  d <- store_of(m)
  expect_lt(as.numeric(object.size(d$spatial)),
            as.numeric(object.size(m$spatial)) / 100)
  # A homogeneous cohort stores one schema, not one per sample.
  expect_false(is.list(attr(d$spatial, "columns")))
})

test_that("split_tissue on a disk-backed mif writes an overlay, not the base file", {
  m <- halfplane_mif()
  dir <- withr::local_tempdir()
  root <- file.path(dir, "s.mif")
  d <- mif_to_disk(m, root)
  base_before <- file.size(store_paths(d$spatial, 1))

  out <- split_tissue(d, classifier = "Classifier.Label", class1 = "Tumor",
                      class2 = "Stroma", sigma = 40, interface_width = 20,
                      workers = 1)

  expect_true(file.exists(file.path(root, "overlay")))
  expect_identical(file.size(store_paths(d$spatial, 1)), base_before)
  expect_true(all(c("density_compartment", "refined_density_compartment",
                    "density_score") %in% colnames(out$spatial[[1]])))
  # The overlay is recorded, so a reopened store still has the columns.
  o <- open_mif(root)
  expect_true("density_score" %in% colnames(o$spatial[[1]]))
  expect_identical(o$spatial[[1]]$density_compartment,
                   out$spatial[[1]]$density_compartment)
  # Re-running replaces the overlay columns rather than duplicating them.
  again <- split_tissue(o, classifier = "Classifier.Label", class1 = "Tumor",
                        class2 = "Stroma", sigma = 40, interface_width = 20,
                        workers = 1, overwrite = TRUE)
  expect_equal(sum(colnames(again$spatial[[1]]) == "density_score"), 1)
})

test_that("split_tissue refuses to write derived columns in reference mode", {
  # Reference mode points at the user's own files; adding a covariate must not
  # modify them, and there is no store to put an overlay in.
  m <- halfplane_mif()
  d <- store_of(m)
  ref <- create_mif(
    clinical_data = data.frame(deidentified_id = "p1"),
    sample_data = data.frame(deidentified_id = "p1", deidentified_sample = "S1"),
    spatial_list = stats::setNames(store_paths(d$spatial), "S1"),
    patient_id = "deidentified_id", sample_id = "deidentified_sample")
  expect_error(
    split_tissue(ref, classifier = "Classifier.Label", class1 = "Tumor",
                 class2 = "Stroma", sigma = 40, interface_width = 20, workers = 1),
    "nowhere to write")
})
