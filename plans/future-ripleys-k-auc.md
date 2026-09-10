# Deferred capability: `ripleys_k_auc()`

Recorded during the 2.0.0 cleanup so it is not lost. Source was `R_dev/ripleys_k_auc.R`,
deleted in commit "chore: remove R_dev/ and doc/ build artifacts". Recover the original with:

```
git show <sha-of-that-commit>^:R_dev/ripleys_k_auc.R
```

## What it did

Took a `mif` that had already been through `ripleys_k()`, collapsed
`mif$derived$univariate_Count` to one curve per `Run`/`iter`/`Label`/`Marker`, and integrated
every value column across `r` to produce a single scalar per curve, written to a new
`mif$derived$univariate_AUC` slot with `"<column> AUC"` names.

This is genuinely useful: it reduces a K curve to one number per sample/marker, which is what
you actually feed into a downstream survival or regression model.

## Why it was parked

It called `flux::auc()`, and `flux` is not in `DESCRIPTION` — so the function could never have
run from an installed copy of the package.

## Why it is not simply revived in 2.0.0

Three reasons, all fixable, none appropriate for a cleanup release:

1. **The dependency is unnecessary.** `flux::auc()` is a trapezoid rule. Base R does it in one
   line, so reviving this should not add a dependency:
   ```r
   auc_trapz <- function(x, y) {
     keep <- !is.na(x) & !is.na(y)
     x <- x[keep]; y <- y[keep]
     if (length(x) < 2L) return(NA_real_)
     o <- order(x); x <- x[o]; y <- y[o]
     sum(diff(x) * (utils::head(y, -1L) + utils::tail(y, -1L)) / 2)
   }
   ```
2. **It referenced columns that do not exist.** `across(\`Theoretical K\`:\`Degree of Clustering
   Theoretical\`)` — there has never been a `Theoretical K` column; it is `Theoretical CSR`. The
   column-range selection would error. It also matched `iter == "Estimate"` while `ripleys_k()`
   writes `"Estimater"` (sic) / `"Estimator"` depending on the branch, and read a `Label` column
   that `ripleys_k()` renames to the user's `sample_id` before returning.
3. **It should generalise.** Post-2.0.0 every metric shares one output schema, so this belongs as
   a single `metric_auc(mif, slot)` that works for K, G, pcf and interaction alike — not as a
   K-only function. That is a natural follow-on to the unified schema, not part of it.

## Recommended shape when revived

```
metric_auc(mif, slot = "univariate_Count", ...)
```

- Integrate every numeric value column over `r` within each
  `Run` / `iter` / `<sample_id>` / `Marker`-or-`Anchor`+`Counted` group.
- Use the base trapezoid helper above; add no dependency.
- Write to `mif$derived$<slot>_AUC`.
- Because the 2.0.0 schema is identical across metrics, one implementation covers
  `univariate_Count`, `bivariate_Count`, `univariate_NN`, `bivariate_NN`,
  `univariate_pair_correlation`, `bivariate_pair_correlation`, `interaction_variable`.
