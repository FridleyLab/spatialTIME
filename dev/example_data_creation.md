
## Generating example data

Randomly generated example data for demonstrating the functions
available in this package.

``` r
library(magrittr)

set.seed(8675309)

tma_one <- tibble::tibble(
  sample = rep("example_one", 150),
  XMin = sample(100:1450, 150),
  XMax = XMin + sample(10:30,150, replace = TRUE),
  YMin = sample(100:1450, 150),
  YMax = YMin + sample(10:30,150, replace = TRUE),
  markers = sample(c("CD3", "CD8", "FOXP3", "PD1", "PCK", "DAPI"),
                   150, replace = TRUE),
  indicator = rep(1, 150),
  cell_type = sample(c("Stroma", "Tumor"), 150, replace = TRUE)
) %>% tidyr::pivot_wider(names_from = markers, values_from = indicator, 
                         values_fill = 0)

tma_two <- tibble::tibble(
  sample = rep("example_two", 150),
  XMin = sample(100:1450, 150),
  XMax = XMin + sample(10:30,150, replace = TRUE),
  YMin = sample(100:1450, 150),
  YMax = YMin + sample(10:30,150, replace = TRUE),
  markers = sample(c("CD3", "CD8", "FOXP3", "PD1", "PCK", "DAPI"),
                   150, replace = TRUE),
  indicator = rep(1, 150),
  cell_type = sample(c("Stroma", "Tumor"), 150, replace = TRUE)
) %>% tidyr::pivot_wider(names_from = markers, values_from = indicator, 
                         values_fill = 0)

tma_three <- tibble::tibble(
  sample = rep("example_three", 150),
  XMin = sample(100:1450, 150),
  XMax = XMin + sample(10:30,150, replace = TRUE),
  YMin = sample(100:1450, 150),
  YMax = YMin + sample(10:30,150, replace = TRUE),
  markers = sample(c("CD3", "CD8", "FOXP3", "PD1", "PCK", "DAPI"),
                   150, replace = TRUE),
  indicator = rep(1, 150),
  cell_type = sample(c("Stroma", "Tumor"), 150, replace = TRUE)
) %>% tidyr::pivot_wider(names_from = markers, values_from = indicator, 
                         values_fill = 0)

tma_four <- tibble::tibble(
  sample = rep("example_four", 150),
  XMin = sample(100:1450, 150),
  XMax = XMin + sample(10:30,150, replace = TRUE),
  YMin = sample(100:1450, 150),
  YMax = YMin + sample(10:30,150, replace = TRUE),
  markers = sample(c("CD3", "CD8", "FOXP3", "PD1", "PCK", "DAPI"),
                   150, replace = TRUE),
  indicator = rep(1, 150),
  cell_type = sample(c("Stroma", "Tumor"), 150, replace = TRUE)
) %>% tidyr::pivot_wider(names_from = markers, values_from = indicator, 
                         values_fill = 0)

if_data <- list(tma_one, tma_two, tma_three, tma_four)

# store 
usethis::use_data(if_data, overwrite = TRUE)
```

    ## ✓ Setting active project to '/Users/creedja/Documents/github/spatialIHC'

    ## ✓ Saving 'if_data' to 'data/if_data.rda'

    ## ● Document your data (see 'https://r-pkgs.org/data.html')
