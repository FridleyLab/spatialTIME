# Symbols that R CMD check would otherwise report as undefined globals.
#
# Three kinds of thing end up here:
#   1. Column names referenced through dplyr/tidyr non-standard evaluation,
#      including backtick-quoted multi-word output columns like `Theoretical CSR`.
#   2. Spatial-file column names that come from the user's data (XMin/XMax/...).
#   3. A few bare function names (across, bind_rows, mutate, ...) reached via
#      `import(dplyr)` rather than being namespace-qualified at the call site.
#
# Entries removed in 2.0.0 because the only code referencing them was deleted:
#   areapp                       -> get_exactK() (also the reason it could never run)
#   Na, Naa, Nab, Nb, Nba, Nbb   -> dix_s_z()
#   "    Obs.Count"              -> dix_s_z() (note the four leading spaces)
#   Permuted K, Range            -> compute_metrics() and the gen-1 drivers
#
# Keep this list honest: when you delete code, re-check whether an entry here went
# with it. A stale entry silently hides a genuine undefined-variable bug -- that is
# exactly how `areapp` in get_exactK() and the undefined `W` in bi_NN_G() survived.
utils::globalVariables(c("r", "label", "Marker", "theo", "Positive", "anchor",
                         "counted", "xloc", "yloc", "XMin", "XMax",
                         "YMin", "YMax", "Var1", "Var2", "Theoretical CSR",
                         "iter", "Permuted G", "Permuted CSR",
                         "Observed G", "Observed K", "Observed", ":=", ".", "is",
                         'Exact CSR', 'Label', 'Var3', 'across',
                         'bind_rows', 'dist', 'distinct', 'mutate',
                         'relocate', 'rename', 'select', 'un', '.data',
                         'Anchor', 'Counted', 'Direction', 'From',
                         'Observed Interaction', 'Observed g',
                         'Permuted Interaction', 'Permuted g',
                         'Permuted_larger_than_Observed', 'Run',
                         'Theoretical G', 'To', 'W', 'cell',
                         'marker', 'marks', 'indicator'))
