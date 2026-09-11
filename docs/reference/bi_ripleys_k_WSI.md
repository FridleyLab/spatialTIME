# Bivariate Ripley's K for Whole Slide Images (removed)

**Removed in spatialTIME 2.0.0.** Use
[`bi_ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/bi_ripleys_k.md)
instead, which now handles whole slide images directly. Calling this
function raises an error.

## Usage

``` r
bi_ripleys_k_WSI(
  mif,
  mnames,
  r_range = 0:100,
  edge_correction = "translation",
  num_permutations = 50,
  permute = FALSE,
  keep_permutation_distribution = FALSE,
  overwrite = FALSE,
  workers = 1,
  big = 1000,
  nlarge = 1000,
  xloc = NULL,
  yloc = NULL,
  ...
)
```

## Arguments

- mif, mnames, r_range, edge_correction, num_permutations, permute:

  Ignored.

- keep_permutation_distribution, overwrite, workers, big, nlarge, xloc,
  yloc:

  Ignored.

- ...:

  Ignored.

## Value

Nothing. Always throws an error directing you to
[`bi_ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/bi_ripleys_k.md).

## Details

`bi_ripleys_k_WSI()` existed because
[`bi_ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/bi_ripleys_k.md)
used to build an n-by-n distance matrix, which is impossible at
whole-slide scale, and so a separate tiled implementation was maintained
alongside it. The two functions were otherwise near-duplicates, and the
tiled one traded accuracy for memory: above a cell-count threshold it
replaced the requested edge correction with `"none"`, its tiled branch
hard-coded translation regardless of the `edge_correction` argument, and
it silently returned garbage for `edge_correction = "isotropic"` because
that branch of its inner loop returned `NULL`.

As of 2.0.0 that trade-off is gone.
[`bi_ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/bi_ripleys_k.md)
only ever materialises the cell pairs closer than `max(r_range)`, so
whole-slide images are handled directly – at 200,000 cells roughly 75 MB
rather than 319 GB – with the requested edge correction always honoured
and results equal to
[`spatstat.explore::Kcross()`](https://rdrr.io/pkg/spatstat.explore/man/Kcross.html)
to floating-point precision.

Migration is a rename. `big` and `nlarge` have no equivalent because
they no longer affect the statistic;
[`bi_ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/bi_ripleys_k.md)
takes a `big` argument that bounds peak memory only.


    # before
    bi_ripleys_k_WSI(mif, mnames, r_range = 0:100, big = 1000, nlarge = 1000)

    # after
    bi_ripleys_k(mif, mnames, r_range = 0:100)

## See also

[`bi_ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/bi_ripleys_k.md)
