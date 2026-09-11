# Bivariate Interaction Variable

For each anchor cell, find the distance to its nearest counted cell;
then report the cumulative percentage of cells whose nearest counted
neighbour falls within each radius. Cells positive for both markers of a
pair are excluded from that pair.

## Usage

``` r
interaction_variable(
  mif,
  mnames,
  r_range = NULL,
  num_permutations = 100,
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

  a character vector, or a two-column data frame of anchor/counted
  markers to assess

- r_range:

  numeric vector of radii at which to evaluate the interaction variable

- num_permutations:

  integer number of permutations used to derive the interaction estimate
  under CSR

- keep_permutation_distribution:

  boolean; keep each permutation's result or average them to one row per
  marker pair and radius

- workers:

  integer number of CPU cores used to process samples in parallel

- overwrite:

  boolean; replace an existing `interaction_variable` slot rather than
  appending it as a new `Run`

- xloc, yloc:

  the x and y columns giving cell centres. If left `NULL`, `XMin`,
  `XMax`, `YMin` and `YMax` must be present.

- ...:

  support for deprecated argument names (see Details).

## Value

object of class `mif` with the `interaction_variable` derived slot
filled

## Details

Single-cell spatial-protein metric introduced by Steinhart et al.,
[doi:10.1158/1541-7786.MCR-21-0411](https://doi.org/10.1158/1541-7786.MCR-21-0411)
.

`keep_perm_dis` is accepted as a deprecated alias for
`keep_permutation_distribution`.

## How the percentage is normalised

The numerator counts anchor cells (one nearest-neighbour distance per
anchor cell), but the denominator is the **total** number of anchor plus
counted cells in the pair. `Observed Interaction` therefore approaches
`100 * n_anchor / (n_anchor + n_counted)` rather than 100 as r grows,
and its ceiling depends on the relative abundance of the two markers.

This is preserved exactly as implemented before 2.0.0 so that existing
results remain comparable – it is documented here rather than changed,
because redefining a published metric is not a refactoring decision. If
you want a quantity that reaches 100%, divide by the anchor count:
`Observed Interaction * (n_anchor + n_counted) / n_anchor`.

## Output columns

As of 2.0.0 this function returns the same columns as
[`bi_ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/bi_ripleys_k.md),
with `Observed Interaction` in place of `Observed K`. `From`/`To` are
now `Anchor`/`Counted`, `Permuted Interaction` is now `Permuted CSR`,
and `Degree of Interaction Permuted` is now
`Degree of Clustering Permutation`. `Theoretical CSR` and `Exact CSR`
are present but always `NA`: this metric has no closed-form null, which
is why it is permutation-only.
