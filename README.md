# Okavango AdDSWE

Code accompanying the manuscript *"A 41-year satellite record reveals abrupt
redistribution of the Okavango Delta's flood pulse since 2012."*

This repository provides the Google Earth Engine and Python code used to generate the
Adapted Dynamic Surface Water Extent (AdDSWE) inundation record for the Okavango Delta
and to reproduce the analyses and figures in the manuscript.

> **Status: under construction.** The repository is being populated. Full documentation,
> the script-to-figure map, and reproducibility instructions are forthcoming.

## Repository layout

```
addswe/        Core AdDSWE algorithm (classification + Landsat compositing)
generation/    GEE notebooks that produce the inundation and ancillary products
analysis/      Notebooks that reproduce the manuscript figures/tables from the archived data
```

## Data

The AdDSWE products and derived datasets are archived on Dryad:
https://doi.org/10.5061/dryad.msbcc2gct

## License

Code is released under the MIT License (see `LICENSE`). The archived data on Dryad are
released under CC0 1.0.
