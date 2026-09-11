# Univariate Pair Correlation Function

The pair correlation function g(r) is the derivative of Ripley's K, so
it measures clustering *at* a radius rather than cumulatively up to it.
It is correspondingly slower to compute and noisier at small r.

`xloc` and `yloc`, if `NULL`, are taken as the midpoints of
`XMin`/`XMax` and `YMin`/`YMax`.

## Usage

``` r
pair_correlation(
  mif,
  mnames,
  r_range = NULL,
  num_permutations = 100,
  edge_correction = "translation",
  keep_permutation_distribution = FALSE,
  workers = 1,
  overwrite = FALSE,
  xloc = NULL,
  yloc = NULL,
  ...
)
```

## Arguments

- mif:

  object of class `mif`

- mnames:

  character vector of marker names

- r_range:

  numeric vector including 0. If `NULL`, `spatstat` chooses the range.

- num_permutations:

  integer number of permutations used to estimate CSR

- edge_correction:

  edge correction passed to
  [`spatstat.explore::pcf()`](https://rdrr.io/pkg/spatstat.explore/man/pcf.html)

- keep_permutation_distribution:

  boolean; keep each permutation's result or average them to one row per
  marker and radius

- workers:

  integer number of CPU cores used to process samples in parallel

- overwrite:

  boolean; replace an existing `univariate_pair_correlation` slot rather
  than appending it as a new `Run`

- xloc, yloc:

  the x and y columns giving cell centres. If left `NULL`, `XMin`,
  `XMax`, `YMin` and `YMax` must be present.

- ...:

  other parameters passed to
  [`spatstat.explore::pcf()`](https://rdrr.io/pkg/spatstat.explore/man/pcf.html),
  plus support for deprecated argument names (see Details).

## Value

`mif` object with the `univariate_pair_correlation` derived slot filled
or appended to

## Details

`keep_perm_dis` is accepted as a deprecated alias for
`keep_permutation_distribution`.

## Output columns

As of 2.0.0 this function returns the same columns as
[`ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/ripleys_k.md),
with `Observed g` in place of `Observed K`. `Theoretical g` and
`Permuted g` are now `Theoretical CSR` and `Permuted CSR`,
`Degree of Correlation *` is now `Degree of Clustering *`, and
`Permuted_larger_than_Observed` is now
`Permutations Larger than Observed`. `Exact CSR` is present but always
`NA`: the closed-form CSR shortcut available to
[`ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/ripleys_k.md)
has not been implemented for the pair correlation function yet.
