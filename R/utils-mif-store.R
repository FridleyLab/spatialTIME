# Disk-backed spatial storage for mif objects.
#
# Why this exists
# ---------------
# `mif$spatial` is normally a named list of in-memory data frames, and every metric
# forks over it with `parallel::mclapply(seq_along(mif$spatial), ...)`. The closure
# captures the whole mif, so the parent process holds every sample resident for the
# entire run. On a 283-slide whole-slide cohort (342,267,952 cells, 21 columns) that
# is 34.2 GB as R data frames against 2.19 GB of parquet on disk -- 100.0 bytes per
# row, measured -- and it is spent to no purpose: a metric reads only the coordinate
# columns, the sample-id column, and the one or two marker columns of the moment.
# 17 of those 21 columns are never touched on a given pass.
#
# A `mif_store` replaces that list with a manifest of per-sample files, so each
# worker reads only the columns it needs for the sample it was given. Peak memory
# becomes O(workers x one sample) instead of O(cohort).
#
# What this does NOT fix: for Ripley's K the per-worker pair list from `k_pairs()`
# overtakes the spatial frame at about max(r_range) = 30, and at the default
# r_range = seq(0, 100, 1) it is roughly 12x the frame (1.23 GB against 0.10 GB on a
# 1,000,977-cell slide). Disk-backing removes the O(cohort) parent term, not the
# O(n x lambda x r^2) per-worker term. Measured: parent peak RSS -64%, per-worker
# peak unchanged or slightly higher. Bounding the pair list is what `big` and
# `chunk_index()` are for; see R/utils-k-engine.R.
#
# THE READER RULE: read_parquet() only, never open_dataset()
# ---------------------------------------------------------
# `arrow::open_dataset() |> select() |> collect()` DOES NOT PRESERVE ROW ORDER. On
# a 1,000,977-row single-row-group file it diverges from the full read at index
# 32,769 -- the 2^15 record-batch boundary -- with 935,427 of 1,000,977 positions
# differing and an identical sorted multiset. Across all 283 files of the cohort it
# matched file order in only 72.
#
# That is disqualifying here, because this package is positional throughout:
# `k_from_pairs()` indexes marker masks into a pair list built in file order
# (R/utils-k-engine.R:185), and `split_tissue()` assigns its three columns back by
# position. A reordered marker column would attach markers to the wrong cells and
# return plausible, wrong numbers.
#
# `arrow::read_parquet(file, col_select = ...)` IS order-stable: verified 283/283
# files (1 to 4 row groups), across repeated reads, across differing column
# subsets, and inside `mclapply(mc.preschedule = FALSE)` children. So every read in
# this file goes through `pq_read()` below and nothing else. Do not "modernise" it
# to a Dataset or to lazy dplyr -- tests/testthat/test-mif-store.R asserts the
# order property precisely so that an arrow upgrade which changes it is caught.


#' On-disk format version
#'
#' Bumped when the store layout changes incompatibly. `open_mif()` refuses a store
#' whose manifest records a different value rather than guessing.
#' @keywords internal
#' @noRd
MIF_STORE_FORMAT <- 1L


# ---------------------------------------------------------------------------
# Backend shim
#
# Parquet is the only spatial backend (arrow is in Imports). These three functions
# exist anyway so that the reader rule above lives in exactly one place instead of
# being repeated at every call site, and so that an `rds` backend -- which would
# lose column projection but keep the O(cohort) win -- is a small, local addition
# if arrow ever proves unavailable somewhere.
# ---------------------------------------------------------------------------

#' Read a parquet file, optionally projecting to some columns
#'
#' @param file path to a parquet file.
#' @param columns character vector of columns, or `NULL` for all of them.
#' @return a base `data.frame`, rows in file order.
#' @keywords internal
#' @noRd
pq_read <- function(file, columns = NULL) {
  if (!file.exists(file)) {
    stop("Spatial file not found: ", file, call. = FALSE)
  }
  # col_select is tidyselect, so a character vector has to go through all_of();
  # bare `columns` would be interpreted as a column *name* in the data mask.
  out <- if (is.null(columns)) {
    arrow::read_parquet(file)
  } else {
    arrow::read_parquet(file, col_select = dplyr::all_of(columns))
  }
  # Metrics index with `[[`, `$<-` and `[` and hand coordinates to spatstat; a
  # tibble's stricter `[` and its list-column tolerance only get in the way here.
  as.data.frame(out, stringsAsFactors = FALSE)
}


# NEVER resize arrow's thread pool from inside an mclapply child.
#
# An earlier version of this file called `arrow::set_cpu_count(1)` /
# `set_io_thread_count(1)` in each forked worker, as insurance against `workers x
# cores` threads competing for `cores` on a large-core HPC. It DEADLOCKS on arrow
# 25.0.1 / R 4.6.1: the child inherits a copy of the thread pool's mutex state from
# the moment of the fork, and resizing the pool then waits forever on a lock no
# thread in the child holds. The symptom is the whole run hanging at 0% CPU with no
# error -- the worst possible failure mode. Reproduced in isolation: reading parquet
# in a forked child is fine (verified on arrow 23.0.1.2 and 25.0.1), and adding only
# the `set_cpu_count(1)` call is what hangs.
#
# It was never worth it anyway: pinning measured SLOWER than leaving arrow's default
# alone (0.31 s against 0.26 s for 24 files on 8 workers). If a huge-core machine
# ever does thrash, set the thread count in the PARENT before calling a metric --
# `arrow::set_cpu_count(2)` -- which children inherit safely across the fork.


#' Write one spatial frame to parquet
#' @keywords internal
#' @noRd
pq_write <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(x, file)
  invisible(file)
}


#' Row count and schema of a parquet file, without reading any data
#'
#' Reads the footer only -- about a millisecond even for a 20 MB file -- which is
#' what makes it affordable to probe every file when a store is opened.
#'
#' @return list with `nrow`, `columns` and `types`.
#' @keywords internal
#' @noRd
pq_probe <- function(file) {
  if (!file.exists(file)) {
    stop("Spatial file not found: ", file, call. = FALSE)
  }
  reader <- arrow::ParquetFileReader$create(file)
  schema <- reader$GetSchema()
  n <- schema$num_fields
  list(
    nrow    = as.numeric(reader$num_rows),
    columns = schema$names,
    # field() is 0-indexed.
    types   = vapply(seq_len(n), function(i) schema$field(i - 1L)$type$ToString(),
                     character(1))
  )
}


# ---------------------------------------------------------------------------
# The mif_store class
#
# Underneath, a mif_store IS a named character vector of file paths. That is
# deliberate rather than incidental: it makes `length()` and `names()` correct with
# no methods at all, and those are exactly the two things the rest of the package
# and its reverse dependencies already call on `mif$spatial` (R/print.R:21,
# R/merge_mifs.R:74-80). Only `[[`, `[`, `c` and `print` need overriding.
#
# Paths are relative to `root` for a store written by mif_to_disk(), and absolute
# in reference mode where the user pointed create_mif() at files they already had.
# Relative-plus-root is what lets a store be moved on disk and still open.
# ---------------------------------------------------------------------------

#' Construct a mif_store
#'
#' @param files character vector of paths; relative to `root` when `root` is set.
#' @param samples sample names; become `names()` of the result.
#' @param nrow integer row count per sample.
#' @param columns either one character vector shared by every sample, or a list of
#'   per-sample character vectors. The shared form is not just an optimisation: a
#'   homogeneous 283-sample cohort would otherwise carry 283 copies of the same 21
#'   names, which is most of what a "the index is tiny" claim is spending.
#' @param types same shape as `columns`.
#' @param overlay per-sample overlay path (relative to `root`), `NA` where none.
#' @param root store root, or `NA_character_` for reference mode.
#' @keywords internal
#' @noRd
new_mif_store <- function(files, samples, nrow, columns, types,
                          overlay = NA_character_, root = NA_character_,
                          format = MIF_STORE_FORMAT) {
  structure(
    stats::setNames(as.character(files), as.character(samples)),
    nrow    = nrow,
    columns = columns,
    types   = types,
    overlay = rep_len(overlay, length(files)),
    root    = root,
    format  = format,
    class   = "mif_store"
  )
}


#' Is this a disk-backed store / a disk-backed mif?
#' @keywords internal
#' @noRd
is_mif_store <- function(x) inherits(x, "mif_store")

#' @keywords internal
#' @noRd
is_disk_mif <- function(mif) is_mif_store(mif$spatial)


#' Resolve a store entry to an absolute path
#'
#' @param x a mif_store.
#' @param i index, or `NULL` for all of them.
#' @keywords internal
#' @noRd
store_paths <- function(x, i = NULL) {
  p <- unclass(x)
  if (!is.null(i)) p <- p[i]
  root <- attr(x, "root")
  if (is.na(root)) return(unname(p))
  unname(file.path(root, p))
}


#' Resolve a store entry's overlay file, or NA
#' @keywords internal
#' @noRd
store_overlay_path <- function(x, i) {
  ov <- attr(x, "overlay")[[i]]
  if (is.na(ov)) return(NA_character_)
  root <- attr(x, "root")
  if (is.na(root)) return(ov)
  file.path(root, ov)
}


#' Per-sample schema, handling the shared and per-sample forms
#' @keywords internal
#' @noRd
store_columns <- function(x, i) {
  cols <- attr(x, "columns")
  if (is.list(cols)) cols[[i]] else cols
}

#' @keywords internal
#' @noRd
store_types <- function(x, i) {
  ty <- attr(x, "types")
  if (is.list(ty)) ty[[i]] else ty
}


#' Normalise a sample index to an integer position
#'
#' Accepts a position or a sample name, because `mif$spatial[["S1"]]` has to keep
#' working for the same reason `mif$spatial[[1]]` does.
#' @keywords internal
#' @noRd
store_index <- function(x, i) {
  if (length(i) != 1L) {
    stop("A mif store must be indexed one sample at a time.", call. = FALSE)
  }
  if (is.character(i)) {
    pos <- match(i, names(x))
    if (is.na(pos)) {
      stop("No sample named \"", i, "\" in this mif. Available: ",
           paste0("\"", names(x), "\"", collapse = ", "), ".", call. = FALSE)
    }
    return(pos)
  }
  pos <- as.integer(i)
  if (is.na(pos) || pos < 1L || pos > length(x)) {
    stop("Sample index ", i, " is out of range; this mif has ", length(x),
         " sample", if (length(x) == 1L) "" else "s", ".", call. = FALSE)
  }
  pos
}


#' Normalise a multi-sample subscript to positive positions
#'
#' Accepts names, positive positions, negative positions and a logical mask, the way
#' `[` on a vector does.
#'
#' `seq_along(x)[i]` rather than `as.integer(i)`: the latter turns a logical mask
#' into 0/1 (`c(FALSE, TRUE)` becomes `c(0L, 1L)`), and because R drops 0 from a
#' subscript that silently returns sample 1 instead of sample 2 -- wrong data, no
#' error. Indexing `seq_along()` handles all four subscript kinds in one line.
#'
#' @keywords internal
#' @noRd
subscript_samples <- function(nms, n, i) {
  if (is.character(i)) {
    pos <- match(i, nms)
    if (anyNA(pos)) {
      stop("No sample named ", paste0("\"", i[is.na(pos)], "\"", collapse = ", "),
           " in this mif.", call. = FALSE)
    }
    return(pos)
  }
  pos <- seq_len(n)[i]
  if (anyNA(pos)) {
    stop("Sample subscript is out of range; this mif has ", n,
         " sample", if (n == 1L) "" else "s", ".", call. = FALSE)
  }
  pos
}

#' @keywords internal
#' @noRd
store_subscript <- function(x, i) subscript_samples(names(x), length(x), i)

#' The same subscript rules for an in-memory spatial slot
#'
#' So that `collect_mif(samples = )` validates identically whichever representation
#' it is handed. Sharing one implementation is the point: when the in-memory branch
#' used plain `list[i]`, an unknown name produced a `NULL` element named `NA` while
#' the disk branch errored, and nothing made the two converge.
#' @keywords internal
#' @noRd
list_subscript <- function(x, i) subscript_samples(names(x), length(x), i)


#' Read one sample, applying the overlay
#'
#' Overlay columns shadow base columns of the same name. That is what lets
#' `split_tissue(overwrite = TRUE)` replace its three columns without rewriting the
#' base file, and it is why the overlay is read second.
#'
#' @keywords internal
#' @noRd
store_read_sample <- function(x, i, columns = NULL) {
  i <- store_index(x, i)

  base_cols <- store_columns(x, i)
  ov_path   <- store_overlay_path(x, i)
  ov_cols   <- if (is.na(ov_path)) character(0) else attr(x, "overlay_columns")[[i]]

  # Default to everything, base order first then whatever the overlay added, which
  # is the column order an in-memory mif would have after the same operations.
  want <- if (is.null(columns)) union(base_cols, ov_cols) else unique(columns)

  unknown <- setdiff(want, c(base_cols, ov_cols))
  if (length(unknown)) {
    stop("Column", if (length(unknown) > 1) "s" else "", " not found in sample \"",
         names(x)[i], "\": ", paste0("\"", unknown, "\"", collapse = ", "), ".\n",
         "  Available: ", abbreviate_names(c(base_cols, ov_cols)), call. = FALSE)
  }

  from_ov   <- intersect(want, ov_cols)
  from_base <- setdiff(want, from_ov)

  out <- if (length(from_base)) {
    pq_read(store_paths(x, i), from_base)
  } else {
    # Every requested column lives in the overlay. Still need the right row count.
    data.frame(row.names = seq_len(attr(x, "nrow")[[i]]))
  }
  if (length(from_ov)) {
    out <- cbind(out, pq_read(ov_path, from_ov), stringsAsFactors = FALSE)
  }
  out[, want, drop = FALSE]
}


#' @export
`[[.mif_store` <- function(x, i, ...) store_read_sample(x, i)


# Subset a store, keeping it a store, so that collect_mif(samples = ) and any future
# per-sample subsetting do not have to reach into the attributes by hand. Documented
# as a comment rather than with roxygen: an Rd page for `[.mif_store` would be noise
# in the reference index, and @export here exists only to register the S3 method.
#' @export
`[.mif_store` <- function(x, i, ...) {
  i <- store_subscript(x, i)
  cols <- attr(x, "columns")
  ty   <- attr(x, "types")
  ovc  <- attr(x, "overlay_columns")
  out <- new_mif_store(
    files   = unclass(x)[i],
    samples = names(x)[i],
    nrow    = attr(x, "nrow")[i],
    columns = if (is.list(cols)) cols[i] else cols,
    types   = if (is.list(ty)) ty[i] else ty,
    overlay = attr(x, "overlay")[i],
    root    = attr(x, "root"),
    format  = attr(x, "format")
  )
  if (!is.null(ovc)) attr(out, "overlay_columns") <- ovc[i]
  out
}


# Make the slot iterable, which it looks like it already is.
#
# `lapply()`, `sapply()`, `vapply()`, `Map()` and `purrr::map()` all begin by
# calling as.list() on anything that is not already a list. Without a method here
# that returned the underlying character vector of FILE PATHS, so
# `sapply(mif$spatial, nrow)` gave a list of NULLs on a disk-backed mif and the row
# counts on an in-memory one -- silently, with no error, which is the worst
# available outcome and the one failure class this store is otherwise careful to
# avoid. The same trap is recorded for `colnames()` in the comment on
# mif_spatial_colnames() below; the package works around it internally by iterating
# `seq_along(mif$spatial)` and going through mif_spatial(), but user code reasonably
# expects a named list of data frames to behave like one.
#
# The hook is the `is.object()` clause in base lapply():
#
#   if (!is.vector(X) || is.object(X)) X <- as.list(X)
#
# A classed object always takes that branch, so defining as.list() here is enough to
# make every *apply function work. The method has to be REGISTERED, not merely
# defined -- an `@export` tag that has not been through document() leaves no
# S3method() line in NAMESPACE, dispatch silently falls through to as.list.default,
# and the paths come back again looking exactly like the original bug.
#
# THIS READS EVERY SAMPLE. There is no way around that: lapply() materialises
# as.list()'s result before applying anything, so the elements cannot be read one at
# a time and released. Promises do not help -- they cache on forcing, and a list
# cannot hold them anyway. Peak memory is therefore the whole cohort, which is the
# thing disk-backing exists to avoid. For a cohort that does not fit, iterate over
# names and stay bounded to one sample:
#
#   for (s in names(mif$spatial)) {
#     spat <- mif$spatial[[s]]      # one sample, released next iteration
#     ...
#   }
#
# or use vapply() over names() rather than over the slot. Both are noted in
# ?mif_to_disk and in the disk-backed vignette.
#' @export
as.list.mif_store <- function(x, ...) {
  out <- lapply(seq_along(x), function(i) store_read_sample(x, i))
  names(out) <- names(x)
  out
}


# Concatenate stores. merge_mifs() merges spatial slots with do.call(c, .), so two
# disk-backed mifs need this. Mixing a disk-backed and an in-memory mif is refused
# there rather than here, where there is no room for a message naming which mif was
# the odd one out.
#' @export
c.mif_store <- function(...) {
  parts <- list(...)
  if (!all(vapply(parts, is_mif_store, logical(1)))) {
    stop("Cannot combine a disk-backed spatial slot with an in-memory one.\n",
         "  Use `collect_mif()` on the disk-backed mif first, or `mif_to_disk()` ",
         "on the in-memory one.", call. = FALSE)
  }
  # Two stores written by versions that disagree about the on-disk layout cannot be
  # read by one set of accessors. Previously the result silently took the first
  # part's format and claimed to be that.
  fmts <- unique(vapply(parts, function(p) as.integer(attr(p, "format")), integer(1)))
  if (length(fmts) > 1L) {
    stop("Cannot combine mif stores written in different on-disk formats (",
         paste(fmts, collapse = ", "), ").\n",
         "  Rebuild the older one with `mif_to_disk()`.", call. = FALSE)
  }
  roots <- unique(vapply(parts, function(p) attr(p, "root"), character(1)))
  # Different roots mean the relative paths are no longer comparable, so absolutise.
  absolute <- length(roots) > 1L
  files <- unlist(lapply(parts, function(p) if (absolute) store_paths(p) else unclass(p)),
                  use.names = FALSE)
  cols <- lapply(seq_along(parts), function(k)
    lapply(seq_along(parts[[k]]), function(i) store_columns(parts[[k]], i)))
  tys <- lapply(seq_along(parts), function(k)
    lapply(seq_along(parts[[k]]), function(i) store_types(parts[[k]], i)))
  out <- new_mif_store(
    files   = files,
    samples = unlist(lapply(parts, names), use.names = FALSE),
    nrow    = unlist(lapply(parts, function(p) attr(p, "nrow")), use.names = FALSE),
    columns = unlist(cols, recursive = FALSE),
    types   = unlist(tys, recursive = FALSE),
    overlay = unlist(lapply(parts, function(p) attr(p, "overlay")), use.names = FALSE),
    root    = if (absolute) NA_character_ else roots,
    format  = attr(parts[[1]], "format")
  )
  ovc <- unlist(lapply(parts, function(p) {
    o <- attr(p, "overlay_columns")
    if (is.null(o)) rep(list(character(0)), length(p)) else o
  }), recursive = FALSE)
  attr(out, "overlay_columns") <- ovc
  out
}


# A store is a named character vector underneath, which is what makes length() and
# names() correct for free -- but it also means the ordinary ways of writing to a
# named list land on that vector instead of on the data. `mif$spatial[1] <- list(df)`
# silently turned the store into a plain list, and `length(mif$spatial) <- 0` into a
# plain character vector; both lose every sample with no error. The assignment forms
# are therefore refused outright, and `$` is supported because the vignettes and
# reverse dependencies read `mif$spatial$SampleName`.
#' @export
`$.mif_store` <- function(x, name) store_read_sample(x, name)

# Tab completion for `mif$spatial$`. Without this it offers nothing at all:
# utils:::.DollarNames.default short-circuits on `if (is.atomic(x) || is.symbol(x))
# character()`, and a store is an atomic character vector underneath. Costs no file
# access -- sample names come from names(), which the index already holds -- so
# completing on a 283-slide store is as cheap as on one core.
#
# Prefix matching with startsWith() rather than the grep() the default would use,
# because these names are image tags: `TMA3_[9,K].tif`,
# `Peres_P3_110212 A3_Scan1.unmixed.qptiff`. A partially typed one is not a valid
# regular expression -- grep("TMA3_[", ...) is an error, not a non-match -- so regex
# matching would turn completion into a thrown condition partway through typing a
# perfectly ordinary sample name.
#
# The @importFrom is load-bearing, and `utils::` prefixing is not a substitute for
# it. Registering an S3 method needs the GENERIC present in this namespace --
# S3method(.DollarNames, mif_store) is resolved at load time -- so without the
# import the package does not load at all: "object '.DollarNames' not found whilst
# loading namespace 'spatialTIME'". devtools::load_all() does not catch this,
# because it registers methods more permissively than a real namespace load; only
# R CMD check does.
#' @importFrom utils .DollarNames
#' @export
.DollarNames.mif_store <- function(x, pattern = "") {
  nms <- names(x)
  if (!length(pattern) || !nzchar(pattern)) return(nms)
  nms[startsWith(nms, pattern)]
}

store_readonly_msg <- function() {
  paste0("A disk-backed spatial slot is read-only.\n",
         "  Use `collect_mif(mif)` to bring it into memory first, or `split_tissue()` ",
         "to add per-cell columns through the overlay.")
}

#' @export
`[[<-.mif_store` <- function(x, i, value) stop(store_readonly_msg(), call. = FALSE)

#' @export
`[<-.mif_store` <- function(x, i, value) stop(store_readonly_msg(), call. = FALSE)

#' @export
`$<-.mif_store` <- function(x, name, value) stop(store_readonly_msg(), call. = FALSE)

#' @export
`length<-.mif_store` <- function(x, value) stop(store_readonly_msg(), call. = FALSE)


#' @export
print.mif_store <- function(x, ...) {
  root <- attr(x, "root")
  cat("<mif_store>", length(x), "samples,",
      format(sum(attr(x, "nrow")), big.mark = ","), "cells\n")
  cat("  ", if (is.na(root)) "referenced files (no store root)" else paste0("root: ", root),
      "\n", sep = "")
  n <- min(3L, length(x))
  if (n) {
    cat("  ", paste0(names(x)[seq_len(n)], " (",
                     format(attr(x, "nrow")[seq_len(n)], big.mark = ","), ")",
                     collapse = ", "),
        if (length(x) > n) paste0(", ... and ", length(x) - n, " more") else "",
        "\n", sep = "")
  }
  invisible(x)
}


# ---------------------------------------------------------------------------
# Accessors used by the metric functions
# ---------------------------------------------------------------------------

#' One sample's spatial data
#'
#' The single read path for every metric. For an in-memory mif this returns the
#' stored data frame and **ignores `columns` entirely**, which is deliberate: the
#' in-memory path stays byte-for-byte what it was before disk-backing existed, so
#' no already-published result can move. The consequence is that a `columns` list
#' missing something a metric actually uses fails only on the disk path, which is
#' what tests/testthat/test-disk-backed-parity.R is there to catch.
#'
#' @param mif a mif.
#' @param i sample position or name.
#' @param columns character vector to project to, or `NULL` for everything.
#' @keywords internal
#' @noRd
mif_spatial <- function(mif, i, columns = NULL) {
  sp <- mif$spatial
  if (!is_mif_store(sp)) return(sp[[i]])
  store_read_sample(sp, i, columns)
}


#' Sample names and count, without reading any spatial data
#'
#' `length()` and `names()` already do the right thing for both representations;
#' these exist so that call sites read as intentional rather than as an accident of
#' the store being a character vector underneath.
#' @keywords internal
#' @noRd
mif_samples <- function(mif) names(mif$spatial)

#' @keywords internal
#' @noRd
n_mif_samples <- function(mif) length(mif$spatial)


#' Column names of one sample, without reading any data
#'
#' Needed wherever the package inspects the shape of the spatial data rather than
#' its contents -- `split_tissue()`'s already-run check is the main one. Going
#' through `colnames(mif$spatial[[i]])` there would read the whole sample on a
#' disk-backed mif, and iterating `vapply(mif$spatial, colnames, ...)` would be
#' worse than slow: a store is a character vector underneath, so `colnames()` of an
#' element is NULL and the check would silently pass for every sample.
#'
#' @keywords internal
#' @noRd
mif_spatial_colnames <- function(mif, i) {
  sp <- mif$spatial
  if (!is_mif_store(sp)) return(colnames(sp[[i]]))
  pos <- store_index(sp, i)
  ovc <- attr(sp, "overlay_columns")
  union(store_columns(sp, pos),
        if (is.null(ovc)) character(0) else ovc[[pos]])
}


#' The columns a metric needs from the spatial data
#'
#' Derives the projection from what each metric already knows, so the list is
#' computed in one place rather than hand-maintained at nine call sites. Mirrors
#' the contract in `add_cell_centres()`: either `xloc`/`yloc` name columns, or the
#' XMin/XMax/YMin/YMax bounding box is required.
#'
#' `mnames` is intersected with what the sample actually has rather than demanded,
#' because `ripleys_k()` selects markers with `dplyr::any_of()` and several metrics
#' tolerate a marker being absent from a given sample. Genuinely missing markers
#' are diagnosed by the metrics themselves, which have the context to say so well.
#'
#' @param mif a mif.
#' @param mnames marker columns, or `NULL`.
#' @param xloc,yloc coordinate columns, or both `NULL` for the bounding box.
#' @param extra anything else the caller needs (a classifier, a cell-type column).
#' @param i sample index, used only to intersect against that sample's schema.
#' @param coords `FALSE` for a metric that needs no cell positions at all
#'   (`marker_freq_diff()` is the only one), so that it does not demand
#'   `XMin`/`XMax`/`YMin`/`YMax` from data that legitimately has no bounding box.
#' @keywords internal
#' @noRd
spatial_columns <- function(mif, mnames = NULL, xloc = NULL, yloc = NULL,
                            extra = NULL, i = NULL, coords = TRUE) {
  if (is.null(xloc) != is.null(yloc)) {
    stop("`xloc` and `yloc` must either both be NULL or both name a column.",
         call. = FALSE)
  }
  coords <- if (!isTRUE(coords)) character(0)
            else if (is.null(xloc)) c("XMin", "XMax", "YMin", "YMax")
            else c(xloc, yloc)
  # `mnames` is a two-column anchor/counted DATA FRAME for the bivariate metrics --
  # see marker_combinations() and the `inherits(mnames, "data.frame")` branch in
  # bi_ripleys_k(). A data frame IS a list, so `c(..., mnames, ...)` would build a
  # list, `%in%` would deparse each column to `c("A", "B")` and match nothing, and
  # dplyr::all_of() would reject the result outright.
  flatten <- function(x) if (is.null(x)) NULL else as.character(unlist(x, use.names = FALSE))
  want <- unique(c(mif$sample_id, coords, flatten(mnames), flatten(extra)))
  want <- want[!is.na(want) & nzchar(want)]

  if (is.null(i) || !is_disk_mif(mif)) return(want)
  # Drop only the optional parts (markers/extras) that this sample does not have;
  # keep coords and sample_id so their absence is still an error, raised by
  # add_cell_centres() with its much better message.
  #
  # mif_spatial_colnames(), not store_columns(): the latter is the BASE parquet
  # schema, so anything split_tissue() wrote to the overlay would be dropped here and
  # then be missing from the frame. That breaks plot_tissue_split(), which asks for
  # `refined_density_compartment` by name.
  have <- mif_spatial_colnames(mif, i)
  keep <- want %in% have | want %in% c(mif$sample_id, coords)
  want[keep]
}


#' Write new columns for one sample into its overlay
#'
#' The base spatial file is never rewritten. Appending `split_tissue()`'s three
#' columns to a 283-slide cohort by rewriting would touch 2.19 GB to add about 12
#' MB, and in reference mode those files are the user's own primary data -- a
#' function that adds a covariate has no business mutating them. So new columns go
#' to a sidecar keyed by sample, and `store_read_sample()` cbinds them back.
#'
#' @param mif a mif.
#' @param i sample position or name.
#' @param new_columns a data frame of the columns to add, with as many rows as the
#'   sample has.
#' @return the mif, with its store updated.
#' @keywords internal
#' @noRd
mif_spatial_set <- function(mif, i, new_columns) {
  if (!is_disk_mif(mif)) {
    for (nm in names(new_columns)) mif$spatial[[i]][[nm]] <- new_columns[[nm]]
    return(mif)
  }
  x   <- mif$spatial
  pos <- store_index(x, i)
  n   <- attr(x, "nrow")[[pos]]
  if (nrow(new_columns) != n) {
    stop("Cannot write ", nrow(new_columns), " values into sample \"",
         names(x)[pos], "\", which has ", n, " cells.", call. = FALSE)
  }
  root <- attr(x, "root")
  if (is.na(root)) {
    stop("This mif references spatial files directly, so there is nowhere to write ",
         "derived columns without modifying your source data.\n",
         "  Use `mif_to_disk(mif, path)` to make a store first.", call. = FALSE)
  }

  # Merge with whatever the overlay already holds, so two functions can each add
  # columns and a re-run can replace its own.
  # Mirror the base file's name rather than re-deriving a stem from the sample name:
  # the base name already went through sample_file_stems() and so is unique, whereas
  # sample_file_stem() alone can map two distinct sample names onto one file and make
  # two samples share an overlay.
  ov_rel  <- file.path("overlay", basename(unclass(x)[[pos]]))
  ov_path <- file.path(root, ov_rel)
  merged <- if (file.exists(ov_path)) {
    prev <- pq_read(ov_path)
    prev[names(new_columns)] <- NULL
    if (ncol(prev)) cbind(prev, new_columns, stringsAsFactors = FALSE) else new_columns
  } else {
    new_columns
  }
  pq_write(merged, ov_path)

  attr(x, "overlay")[[pos]] <- ov_rel
  ovc <- attr(x, "overlay_columns")
  if (is.null(ovc)) ovc <- rep(list(character(0)), length(x))
  ovc[[pos]] <- names(merged)
  attr(x, "overlay_columns") <- ovc
  mif$spatial <- x
  mif
}


#' Join names for an error message without printing a wall of text
#'
#' A HALO export carries 51 columns and a real one carries more, so listing them all
#' buries the name that was actually wrong.
#' @keywords internal
#' @noRd
abbreviate_names <- function(x, n = 12L) {
  if (length(x) <= n) return(paste0(paste(x, collapse = ", "), "."))
  paste0(paste(x[seq_len(n)], collapse = ", "),
         ", ... and ", length(x) - n, " more.")
}


#' A filesystem-safe file stem for a sample name
#'
#' Sample names here are image tags like `TMA3_[9,K].tif` and
#' `Peres_P3_110212 A3_Scan1.unmixed.qptiff`, which carry spaces, brackets, commas
#' and dots. The manifest keeps the real name; only the file name is sanitised, and
#' a positional suffix is added by the caller when two names would collide.
#' @keywords internal
#' @noRd
sample_file_stem <- function(name) {
  stem <- gsub("[^A-Za-z0-9._-]+", "_", name)
  stem <- gsub("^[._]+|[._]+$", "", stem)
  ifelse(nzchar(stem), stem, "sample")
}


#' Unique, filesystem-safe file stems for a set of sample names
#'
#' Sanitising can map two distinct sample names onto one stem (`A/1` and `A 1` both
#' become `A_1`), which would silently make one sample's file overwrite another's.
#' Disambiguate positionally instead.
#' @keywords internal
#' @noRd
sample_file_stems <- function(names_) {
  stems <- sample_file_stem(names_)
  dup <- duplicated(stems) | duplicated(stems, fromLast = TRUE)
  if (any(dup)) stems[dup] <- paste0(stems[dup], "__", seq_along(stems)[dup])
  stems
}
