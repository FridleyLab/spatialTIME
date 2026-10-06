# Writing, opening and collecting a disk-backed mif.
#
# What goes in which format, and why it is not all parquet
# -------------------------------------------------------
# Spatial frames are parquet. They are the only slot that is large, the only one
# read per-worker, and the only one that benefits from column projection -- which is
# the whole point of the exercise. Factor columns round-trip exactly through
# parquet, including unused levels and NAs (verified), which matters because
# `split_tissue()` writes two factor columns per cell.
#
# `clinical`, `sample` and `derived` are RDS. They are small -- one row per patient,
# one per sample, and a few thousand rows per metric table -- so projection buys
# nothing, while RDS preserves R types and, crucially, ATTRIBUTES exactly.
# `split_tissue()` hangs a `call_info` attribute on `derived$density_boundary`
# (R/split_tissue.R:323-335) and `merge_mifs()` reads it back
# (R/merge_mifs.R:132-141); parquet would silently drop it. `derived` also holds
# list-valued slots (`density_boundary`, `spatial_plots`) that are not tabular at
# all. One format for that whole group is simpler than a per-slot rule, and costs
# nothing measurable.
#
# The integrity stamp
# -------------------
# The manifest records each spatial file's byte size, row count and column schema,
# and `open_mif()` re-probes every file and compares. That catches truncation, a
# half-finished write, a replaced file and a schema change, at the cost of one
# footer read per sample (about a millisecond). It deliberately does NOT hash file
# contents: md5 over 2.19 GB takes seconds to minutes and would make opening a
# cohort feel broken, and `digest` is not a dependency. It also deliberately does
# not record mtime, which changes on a perfectly valid copy or restore and would
# produce false alarms on exactly the operation the relative-path layout exists to
# support.


#' Write a mif's spatial data to disk
#'
#' @description
#' Materialises an in-memory `mif` into a self-describing directory and returns a
#' disk-backed `mif` that reads from it. Metric functions then load only the
#' columns they need, for the one sample a worker is handling, instead of the
#' parent process holding every sample resident for the whole run.
#'
#' @details
#' Use this when the spatial data no longer comfortably fits in memory, which in
#' practice means whole slide images. On a 283-slide cohort of 342 million cells the
#' in-memory spatial slot is 34.2 GB against 2.19 GB of parquet, and since the
#' parent process must hold all of it while `mclapply()` forks, a multi-core run
#' cannot start at all. Disk-backing makes peak memory scale with `workers` rather
#' than with the number of samples.
#'
#' It is not a general speed-up. For [ripleys_k()] the per-worker pair list
#' overtakes the spatial frame at roughly `max(r_range) = 30`, so at the default
#' `r_range` the frame is a minor part of a worker's footprint; what disk-backing
#' removes is the parent's copy of the whole cohort.
#'
#' The store is a plain directory you can inspect, copy and archive:
#' \preformatted{
#' cohort.mif/
#'   manifest.json      format, ids, and per-sample rows/columns/size
#'   clinical.rds  sample.rds
#'   spatial/  <sample>.parquet    one file per sample, never rewritten
#'   overlay/  <sample>.parquet    columns added later by split_tissue()
#'   derived/  <slot>.rds
#' }
#' Paths inside `manifest.json` are relative, so the whole directory can be moved
#' or copied and still open.
#'
#' @section What is and is not preserved:
#' Column names, column order, types and factor levels -- including unused levels
#' and `NA`s -- round-trip exactly, so `collect_mif(mif_to_disk(x, p))` returns the
#' spatial data `x` had.
#'
#' The spatial slot is iterable, so `sapply(mif$spatial, nrow)` and
#' `lapply(mif$spatial, ...)` work the same way they do on an in-memory `mif`. **But
#' iterating reads every sample**, because `lapply()` materialises its input before
#' applying anything, so peak memory becomes the whole cohort — the cost
#' disk-backing exists to avoid. On a cohort that does not fit, iterate over names
#' instead and only one sample is resident at a time:
#'
#' \preformatted{
#' # reads everything at once
#' sapply(mif$spatial, nrow)
#'
#' # reads one sample at a time
#' vapply(names(mif$spatial), function(s) nrow(mif$spatial[[s]]), numeric(1))
#' }
#'
#' Writing to the slot is refused outright — `mif$spatial[[i]] <- x` and
#' `length(mif$spatial) <- n` error rather than destroying the index. Use
#' [collect_mif()] to bring the data into memory if you need to modify it.
#'
#' @section The store is shared state:
#' A metric writes its results into the store as well as into the `mif` it returns,
#' so they survive a reopen. One consequence is worth knowing, because it is the one
#' place a disk-backed `mif` does not behave like an ordinary R object:
#'
#' \preformatted{
#' a <- ripleys_k(xd, ...)                      # result goes to `a`
#' identical(open_mif(store)$derived, xd$derived)   # FALSE
#' identical(open_mif(store)$derived, a$derived)    # TRUE
#' }
#'
#' The in-memory `xd` is untouched, exactly as copy-on-modify implies — but `xd` and
#' `a` point at the same directory, and the metric wrote there. The store reflects
#' whichever `mif` most recently ran a metric against it.
#'
#' Nothing is silently dropped: each run carries forward the derived slots it
#' inherited, so only the slot being written is replaced. Assigning the result back
#' to the same name (`xd <- ripleys_k(xd, ...)`) keeps object and store in step. To
#' hold two parameterisations at once, give them separate stores rather than separate
#' variables.
#'
#' **Row names are not preserved**; they come back as `1:nrow`. A columnar file has
#' nowhere to put them, and storing them would mean an extra column per sample for
#' something no function in this package reads -- every metric indexes cells by
#' position. Non-default row names only arise from having subset a data frame
#' (`spat[keep, ]`), and `dplyr` resets them anyway, so this is visible only if you
#' compare with [identical()] after such a subset. Use
#' `all.equal(check.attributes = FALSE)`, or reset them yourself, if you need that
#' comparison to hold.
#'
#' @param mif object of class `mif` created with [create_mif()].
#' @param path directory to write the store to. **Required, with no default** --
#'   a store can be tens of gigabytes, so it is always written where you say and
#'   never somewhere chosen for you. `~` is expanded.
#' @param overwrite whether to replace an existing store at `path`.
#'
#' @return object of class `mif` whose `spatial` slot reads from `path`.
#' @seealso [open_mif()] to reopen one, [collect_mif()] to load it back into memory.
#' @export
#' @examples
#' x <- create_mif(clinical_data = spatialTIME::example_clinical,
#'   sample_data = spatialTIME::example_summary,
#'   spatial_list = spatialTIME::example_spatial[1],
#'   patient_id = "deidentified_id", sample_id = "deidentified_sample")
#' store <- file.path(tempdir(), "example.mif")
#' xd <- mif_to_disk(x, store, overwrite = TRUE)
#' xd
#' unlink(store, recursive = TRUE)
mif_to_disk <- function(mif, path, overwrite = FALSE) {
  if (!inherits(mif, "mif")) {
    stop("mIF should be of class `mif` created with function `create_mif()`",
         call. = FALSE)
  }
  # No default for `path`, following `sigma` in split_tissue(): a silent default
  # here would write tens of GB somewhere the caller never named.
  if (missing(path) || is.null(path)) {
    stop("`path` must name the directory to write the store to, for example\n",
         "  mif_to_disk(mif, \"~/projects/cohort.mif\")", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a single non-empty directory path.", call. = FALSE)
  }
  n <- length(mif$spatial)
  if (!n || (n == 1L && !is_mif_store(mif$spatial) && identical(mif$spatial[[1]], NA))) {
    stop("This mif has no spatial data to write.", call. = FALSE)
  }

  root <- path.expand(path)
  if (file.exists(root)) {
    if (!overwrite) {
      stop("A store already exists at \"", root, "\".\n",
           "  Pass `overwrite = TRUE` to replace it, or choose another `path`.",
           call. = FALSE)
    }
    unlink(root, recursive = TRUE)
  }
  dir.create(file.path(root, "spatial"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "derived"), recursive = TRUE, showWarnings = FALSE)

  samples <- names(mif$spatial)
  if (is.null(samples) || anyNA(samples) || !all(nzchar(samples))) {
    stop("Every element of the spatial slot must be named to be written to disk.",
         call. = FALSE)
  }
  stems <- sample_file_stems(samples)

  rel <- character(n); nrows <- numeric(n)
  cols <- vector("list", n); types <- vector("list", n); bytes <- numeric(n)
  for (k in seq_len(n)) {
    spat <- mif_spatial(mif, k)
    if (!is.data.frame(spat)) {
      stop("Sample \"", samples[k], "\" is not a data frame.", call. = FALSE)
    }
    check_storable_columns(spat, samples[k])
    rel[k] <- file.path("spatial", paste0(stems[k], ".parquet"))
    pq_write(spat, file.path(root, rel[k]))
    p <- pq_probe(file.path(root, rel[k]))
    nrows[k] <- p$nrow; cols[[k]] <- p$columns; types[[k]] <- p$types
    bytes[k] <- file.size(file.path(root, rel[k]))
  }

  store <- new_mif_store(
    files = rel, samples = samples, nrow = nrows,
    columns = collapse_schema(cols), types = collapse_schema(types),
    overlay = NA_character_, root = root
  )
  attr(store, "overlay_columns") <- rep(list(character(0)), n)

  mif$spatial <- store
  write_store_tables(root, mif)
  write_manifest(root, mif, store, bytes)
  mif
}


#' Refuse column types that parquet cannot carry back unchanged
#'
#' Checked before the first byte is written, so a cohort does not get half-written
#' and then rejected.
#'
#' `raw` is the dangerous one: `arrow::write_parquet()` accepts it and
#' `collect_mif()` hands back an `integer`, silently, with no warning anywhere. A
#' list column round-trips as an `AsIs` column whose contents are not
#' [identical()] to the input. `complex` is refused by arrow itself, but with
#' "Cannot infer type from vector", which names neither the column nor the sample.
#'
#' Factors, `logical`, `integer`, `double`, `character`, `Date` and `POSIXct`
#' (including `tzone`) are all exact, so they are not mentioned here. `integer64`
#' is allowed through but comes back as `integer` when every value fits in 32 bits
#' -- documented in `?mif_to_disk` rather than refused, because it is lossless for
#' the values concerned.
#'
#' @param spat one sample's spatial data frame. @param sample its name.
#' @keywords internal
#' @noRd
check_storable_columns <- function(spat, sample) {
  is_bad <- function(col) {
    (is.list(col) && !is.data.frame(col)) || is.raw(col) || is.complex(col)
  }
  bad <- names(spat)[vapply(spat, is_bad, logical(1))]
  if (!length(bad)) return(invisible(TRUE))
  kind <- vapply(spat[bad], function(col) {
    if (is.raw(col)) "raw" else if (is.complex(col)) "complex" else "list"
  }, character(1))
  stop("Sample \"", sample, "\" has ", if (length(bad) > 1) "columns" else "a column",
       " that cannot be stored on disk: ",
       paste0("\"", bad, "\" (", kind, ")", collapse = ", "), ".\n",
       "  Parquet would not return ", if (length(bad) > 1) "them" else "it",
       " unchanged. Drop the column", if (length(bad) > 1) "s" else "",
       ", or keep this mif in memory.", call. = FALSE)
}


#' Write a store's non-spatial slots
#'
#' `clinical`, `sample` and every `derived` slot, as RDS. Shared by `mif_to_disk()`,
#' `subset_mif()` and `sync_manifest()` -- `subset_mif()` writing a store without
#' these is exactly the bug that made its output unopenable.
#' @keywords internal
#' @noRd
write_store_tables <- function(root, mif) {
  saveRDS(mif$clinical, file.path(root, "clinical.rds"))
  saveRDS(mif$sample,   file.path(root, "sample.rds"))
  dir.create(file.path(root, "derived"), recursive = TRUE, showWarnings = FALSE)
  for (slot in names(mif$derived)) {
    saveRDS(mif$derived[[slot]],
            file.path(root, "derived", paste0(sample_file_stem(slot), ".rds")))
  }
  invisible(NULL)
}


#' Reopen a mif store written by `mif_to_disk()`
#'
#' @description
#' Reads and validates `manifest.json`, re-probes every spatial file, and returns
#' the `mif`. Validation happens up front so that a moved, truncated or
#' half-written store fails immediately with a message naming the sample, rather
#' than as a confusing error partway through a long computation.
#'
#' @param path directory containing the store.
#' @return object of class `mif`.
#' @seealso [mif_to_disk()], [collect_mif()]
#' @export
open_mif <- function(path) {
  if (missing(path) || !is.character(path) || length(path) != 1L) {
    stop("`path` must be a single directory path.", call. = FALSE)
  }
  root <- path.expand(path)
  mp <- file.path(root, "manifest.json")
  if (!dir.exists(root)) {
    stop("No such directory: \"", root, "\".", call. = FALSE)
  }
  if (!file.exists(mp)) {
    stop("\"", root, "\" is not a mif store: no manifest.json.\n",
         "  Stores are created with `mif_to_disk(mif, path)`.", call. = FALSE)
  }
  m <- jsonlite::fromJSON(mp, simplifyVector = FALSE)

  if (!identical(as.integer(m$format), MIF_STORE_FORMAT)) {
    stop("This store uses mif format version ", m$format, ", but this version of ",
         "spatialTIME reads version ", MIF_STORE_FORMAT, ".\n",
         "  It was written by spatialTIME ", m$spatialTIME_version %||% "(unknown)",
         "; install a matching version, or rebuild the store.", call. = FALSE)
  }

  n <- length(m$samples)
  if (!n) stop("This store's manifest lists no samples.", call. = FALSE)
  samples <- vapply(m$samples, function(s) s$sample, character(1))
  rel     <- vapply(m$samples, function(s) s$file, character(1))
  nrows   <- vapply(m$samples, function(s) as.numeric(s$nrow), numeric(1))
  bytes   <- vapply(m$samples, function(s) as.numeric(s$bytes), numeric(1))
  cols    <- lapply(m$samples, function(s) as.character(unlist(s$columns)))
  types   <- lapply(m$samples, function(s) as.character(unlist(s$types)))

  ov  <- vapply(m$samples, function(s) s$overlay %||% NA_character_, character(1))
  ovc <- lapply(m$samples, function(s) as.character(unlist(s$overlay_columns)))

  # Re-probe rather than trust. This is the point of the manifest.
  for (k in seq_len(n)) {
    validate_store_file(file.path(root, rel[k]), rel = rel[k], sample = samples[k],
                        what = "spatial file", nrow_expect = nrows[k],
                        cols_expect = cols[[k]], bytes_expect = bytes[k])
    # A sample id column that is absent or misnamed fails much later and far from
    # the cause -- inside add_cell_centres() or a dplyr::select() in whichever
    # metric was called. Catch it here, where the message can say so.
    if (!m$sample_id %in% cols[[k]]) {
      stop("Sample \"", samples[k], "\" has no column named \"", m$sample_id,
           "\", which this store's manifest gives as its `sample_id`.\n",
           "  Available: ", abbreviate_names(cols[[k]]), call. = FALSE)
    }
  }

  # Overlays used to get file.exists() and nothing else, so a truncated or
  # wrong-schema overlay opened cleanly and then failed deep inside a read. There
  # is no recorded byte size for them (that would need a format bump), but the row
  # count must match the sample and the columns must match what the manifest says
  # the overlay added -- both already in format 1.
  for (k in which(!is.na(ov))) {
    validate_store_file(file.path(root, ov[k]), rel = ov[k], sample = samples[k],
                        what = "overlay", nrow_expect = nrows[k],
                        cols_expect = ovc[[k]], bytes_expect = NULL)
  }

  store <- new_mif_store(
    files = rel, samples = samples, nrow = nrows,
    columns = collapse_schema(cols), types = collapse_schema(types),
    overlay = ov, root = root
  )
  attr(store, "overlay_columns") <- ovc

  # Manifest order is authoritative: list.files() sorts alphabetically, which would
  # silently reorder derived slots on every reopen. And the manifest is the list of
  # what SHOULD be there, so read from it rather than from the directory -- an
  # intersect() with whatever happens to be on disk turns a half-copied store into
  # a mif that quietly lost a metric table.
  derived <- list()
  known <- as.character(unlist(m$derived))
  for (slot in known) {
    derived[[slot]] <- read_store_rds(root, file.path("derived",
                                                      paste0(sample_file_stem(slot), ".rds")))
  }

  structure(list(clinical = read_store_rds(root, "clinical.rds"),
                 sample   = read_store_rds(root, "sample.rds"),
                 spatial  = store,
                 derived  = derived,
                 patient_id = m$patient_id,
                 sample_id  = m$sample_id),
            class = "mif")
}


#' Probe one store file and compare it to what the manifest claims
#'
#' Shared by the base-file and overlay loops of [open_mif()] and by [verify_mif()],
#' so a file can never be validated one way in one place and another way elsewhere
#' -- which is how overlays ended up with only an existence check.
#'
#' @param f absolute path. @param rel path as recorded, for the message.
#' @param sample sample name, for the message. @param what "spatial file" or "overlay".
#' @param nrow_expect,cols_expect expected row count and column names.
#' @param bytes_expect expected size, or `NULL` where none is recorded.
#' @keywords internal
#' @noRd
validate_store_file <- function(f, rel, sample, what, nrow_expect, cols_expect,
                                bytes_expect = NULL) {
  if (!file.exists(f)) {
    stop("Store is incomplete: ", what, " for sample \"", sample,
         "\" should be at \"", rel, "\", which does not exist.", call. = FALSE)
  }
  if (!is.null(bytes_expect)) {
    got <- file.size(f)
    if (!identical(as.numeric(got), as.numeric(bytes_expect))) {
      stop("Sample \"", sample, "\" has changed on disk: manifest records ",
           format(bytes_expect, big.mark = ","), " bytes, ", what, " is ",
           format(got, big.mark = ","), " bytes.\n",
           "  The store was modified outside spatialTIME; rebuild it with ",
           "`mif_to_disk()`.", call. = FALSE)
    }
  }
  p <- pq_probe(f)
  if (!identical(as.numeric(p$nrow), as.numeric(nrow_expect))) {
    stop("Sample \"", sample, "\" (", what, "): manifest says ",
         format(nrow_expect, big.mark = ","), " rows, file has ",
         format(p$nrow, big.mark = ","), ".", call. = FALSE)
  }
  if (!identical(p$columns, as.character(cols_expect))) {
    stop("Sample \"", sample, "\" (", what, "): columns on disk do not match the ",
         "manifest.\n  manifest: ", paste(cols_expect, collapse = ", "),
         "\n  file:     ", paste(p$columns, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}


#' Read one of a store's RDS slots, reporting a missing one as a store problem
#'
#' A bare `readRDS()` on an absent file gives a raw `gzfile` warning followed by
#' "cannot open the connection", which names neither the store nor what was
#' missing. A half-copied store is exactly the case this has to explain.
#' @keywords internal
#' @noRd
read_store_rds <- function(root, rel) {
  f <- file.path(root, rel)
  if (!file.exists(f)) {
    stop("Store is incomplete: \"", rel, "\" is missing from \"", root, "\".\n",
         "  The store was not finished writing, or was only partially copied.",
         call. = FALSE)
  }
  readRDS(f)
}


#' Re-probe every spatial file a mif points at
#'
#' @description
#' Checks that each file still has the row count and columns the mif expects, and
#' reports the first that does not.
#'
#' [open_mif()] does this automatically, but two cases are not covered by it. A mif
#' built with `create_mif(spatial_list = <parquet paths>)` references your files
#' directly and has no manifest, so there is nothing to validate at construction
#' and nothing stops those files changing afterwards. And a long-running session
#' holds a store open across hours of computation, during which the files can be
#' replaced underneath it. In both cases the recorded row count goes stale, and
#' because the metrics index cells by position a stale count means markers attached
#' to the wrong cells rather than an error.
#'
#' @param mif object of class `mif`. In-memory mifs have nothing to verify and pass
#'   trivially, so this is safe to call unconditionally.
#' @return invisibly `TRUE`; errors naming the first sample that fails.
#' @seealso [mif_to_disk()], [open_mif()]
#' @export
#' @examples
#' x <- create_mif(clinical_data = spatialTIME::example_clinical,
#'   sample_data = spatialTIME::example_summary,
#'   spatial_list = spatialTIME::example_spatial[1],
#'   patient_id = "deidentified_id", sample_id = "deidentified_sample")
#' store <- file.path(tempdir(), "verify.mif")
#' xd <- mif_to_disk(x, store, overwrite = TRUE)
#' verify_mif(xd)
#' unlink(store, recursive = TRUE)
verify_mif <- function(mif) {
  if (!inherits(mif, "mif")) {
    stop("mIF should be of class `mif` created with function `create_mif()`",
         call. = FALSE)
  }
  if (!is_disk_mif(mif)) return(invisible(TRUE))
  x <- mif$spatial
  for (k in seq_along(x)) {
    validate_store_file(
      store_paths(x, k), rel = unname(unclass(x)[k]), sample = names(x)[k],
      what = "spatial file", nrow_expect = attr(x, "nrow")[[k]],
      cols_expect = store_columns(x, k), bytes_expect = NULL)
    ov <- store_overlay_path(x, k)
    if (!is.na(ov)) {
      validate_store_file(
        ov, rel = attr(x, "overlay")[[k]], sample = names(x)[k], what = "overlay",
        nrow_expect = attr(x, "nrow")[[k]],
        cols_expect = attr(x, "overlay_columns")[[k]], bytes_expect = NULL)
    }
  }
  invisible(TRUE)
}


#' Load a disk-backed mif's spatial data back into memory
#'
#' @description
#' The inverse of [mif_to_disk()], named after `dplyr::collect()` for the same
#' reason: it is the point where lazily-referenced data is materialised.
#'
#' @param mif object of class `mif`.
#' @param samples optionally a subset of sample names or positions to load, so that
#'   a few samples can be inspected interactively without materialising a cohort
#'   that does not fit in memory.
#' @return object of class `mif` with an ordinary named list of data frames in its
#'   `spatial` slot.
#' @seealso [mif_to_disk()], [open_mif()]
#' @export
collect_mif <- function(mif, samples = NULL) {
  if (!inherits(mif, "mif")) {
    stop("mIF should be of class `mif` created with function `create_mif()`",
         call. = FALSE)
  }
  if (!is_disk_mif(mif)) {
    if (is.null(samples)) return(mif)
    # Validated the same way as the disk path, and with the same message. Plain
    # `list[samples]` returns a NULL element named NA for a name that is not there,
    # so `collect_mif(m, samples = "typo")` used to hand back a one-sample mif whose
    # sample was NULL -- on the disk path the identical call errors.
    mif$spatial <- mif$spatial[list_subscript(mif$spatial, samples)]
    return(mif)
  }
  store <- mif$spatial
  idx <- if (is.null(samples)) seq_along(store) else store_subscript(store, samples)
  out <- lapply(idx, function(k) store_read_sample(store, k))
  names(out) <- names(store)[idx]
  mif$spatial <- out
  mif
}


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

#' Write manifest.json
#' @keywords internal
#' @noRd
write_manifest <- function(root, mif, store, bytes) {
  n <- length(store)
  ovc <- attr(store, "overlay_columns")
  m <- list(
    format = MIF_STORE_FORMAT,
    backend = "parquet",
    spatialTIME_version = as.character(utils::packageVersion("spatialTIME")),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    patient_id = mif$patient_id,
    sample_id = mif$sample_id,
    derived = as.list(names(mif$derived)),
    samples = lapply(seq_len(n), function(k) list(
      sample  = names(store)[k],
      file    = unname(unclass(store)[k]),
      nrow    = attr(store, "nrow")[[k]],
      bytes   = bytes[[k]],
      columns = as.list(store_columns(store, k)),
      types   = as.list(store_types(store, k)),
      overlay = attr(store, "overlay")[[k]],
      overlay_columns = as.list(if (is.null(ovc)) character(0) else ovc[[k]])
    ))
  )
  writeLines(jsonlite::toJSON(m, auto_unbox = TRUE, pretty = TRUE, null = "null"),
             file.path(root, "manifest.json"))
  invisible(file.path(root, "manifest.json"))
}


#' Rewrite manifest.json to match a store whose overlays have changed
#'
#' `mif_spatial_set()` updates the in-memory store as it writes overlay files; this
#' persists that so the next `open_mif()` finds them. Called once at the end of an
#' operation rather than per sample, because rewriting a 283-sample manifest inside
#' a per-sample loop is pointless I/O.
#' @keywords internal
#' @noRd
sync_manifest <- function(mif) {
  if (!is_disk_mif(mif)) return(invisible(NULL))
  root <- attr(mif$spatial, "root")
  if (is.na(root)) return(invisible(NULL))
  bytes <- vapply(seq_along(mif$spatial),
                  function(k) file.size(store_paths(mif$spatial, k)), numeric(1))
  write_store_tables(root, mif)
  write_manifest(root, mif, mif$spatial, bytes)
}


#' Store one schema instead of N identical ones
#'
#' A homogeneous cohort -- which is the normal case, since samples come off one
#' imaging platform -- would otherwise carry 283 copies of the same 21 column
#' names, which is most of the weight of an index that is supposed to be
#' negligible.
#' @keywords internal
#' @noRd
collapse_schema <- function(x) {
  if (!length(x)) return(character(0))
  first <- x[[1]]
  if (all(vapply(x, identical, logical(1), first))) first else x
}


#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
