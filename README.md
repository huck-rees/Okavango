# Okavango AdDSWE

Code accompanying the manuscript *"A 41-year satellite record reveals abrupt
redistribution of the Okavango Delta's flood pulse since 2012."*

This repository provides the Google Earth Engine and Python code used to (1) generate the
**Adapted Dynamic Surface Water Extent (AdDSWE)** inundation record for the Okavango
Delta, Botswana, and (2) reproduce the analyses and figures in the manuscript.

AdDSWE is a 41-year (June 1984 – December 2025), 30-m, monthly, multiclass surface-water
record derived from the Landsat Collection 2 archive. It augments the standard Dynamic
Surface Water Extent (DSWE) confidence framework (Jones 2019) with a dynamic
shortwave-infrared (SWIR2) threshold test ("Test 6") that detects water beneath
vegetation, classifying every pixel into five confidence classes:

| Class | Meaning |
|-------|---------|
| 0 | Not water / dry land |
| 1 | Low-confidence inundation (vegetated / partial) |
| 2 | Partial / water under canopy |
| 3 | Moderate-confidence open water |
| 4 | High-confidence open water |

Open water = classes 3–4; total inundation = classes 1–4.

## Data availability

The AdDSWE products and all derived datasets are archived on Dryad:

- **Dryad:** https://doi.org/10.5061/dryad.msbcc2gct
- **Interactive viewer:** https://okavango-water-mask-viewer.streamlit.app/

The `analysis/` notebooks are designed to run against the Dryad file layout — download the
archive and point each notebook's input paths at your local copy. See the Dryad `README`
for file-level descriptions, variable definitions, units, and coordinate reference systems.

## Repository layout

```
addswe/          Core AdDSWE algorithm (importable package)
  classification.py   spectral indices, 5 DSWE tests, dynamic SWIR2 Test 6, decision tree, morphological filter
  composites.py       Landsat monthly compositing + iterative temporal gap-filling
generation/      GEE notebooks that produce the inundation and ancillary products
analysis/        Notebooks that reproduce the manuscript figures/tables from the archived data
environment.yml  Conda environment specification
```

## Script → figure / table map

| Notebook | Manuscript output |
|----------|-------------------|
| `generation/landsat_addswe_monthly_30m.ipynb` | AdDSWE + original DSWE products (inputs to Figs. 2, 4, 7) |
| `generation/chirps_precipitation.ipynb` | Precipitation series (Fig. 7, Supp. Fig. 4) |
| `generation/pekel_watersurface.ipynb` | JRC/Pekel comparison series (Supp. Fig. 4a) |
| `analysis/monthly_record_fig2.ipynb` | Figure 2 + `LS_AdDSWE_monthly_summary` |
| `analysis/distributary_total_inundation_fig3.ipynb` | Figure 3, Table 1 (total inundation) |
| `analysis/distributary_open_water_suppfig1.ipynb` | Supplementary Fig. 1, Table 1 (open water) |
| `analysis/mean_difference_map_fig4.ipynb` | Figure 4 (7-band difference raster) |
| `analysis/centroid_migration_fig5.ipynb` | Figure 5 + `centroid_weighted_all.csv` |
| `analysis/open_water_validation_supp.ipynb` | Fig. 7a–d, Supp. Fig. 5, Supp. Table 1 |
| `analysis/qc_mask_stats_suppnote2.ipynb` | Supplementary Note 2 (QC statistics) |

## Reproducibility & requirements

The two stages differ in how they can be re-run:

- **Generation** (`generation/`, and the product-generation sections of the Fig. 4 and
  Fig. 5 notebooks) requires a Google Earth Engine account and project, authentication,
  and long-running export tasks. These notebooks document exactly how the products were
  produced; they are **not** turnkey.
- **Analysis** (`analysis/`) reproduces the manuscript figures and tables from the
  **archived Dryad data** and is runnable by anyone once the archive is downloaded.

The core AdDSWE algorithm is in the `addswe` package so it can be read, cited, and reused
independently of the notebooks. Call `ee.Initialize()` before using the Earth Engine
functions.

### Installation

```bash
conda env create -f environment.yml
conda activate okavango
```

Then, for Earth Engine functionality, authenticate once with `earthengine authenticate`.

## Note on `Study_Compare_MaxInundation.csv`

The cross-product comparison table (Dryad `inundation/Study_Compare_MaxInundation.csv`;
used in Fig. 7 and Supp. Fig. 4) is **compiled from multiple sources** rather than produced
by a single script. Its columns are drawn from:

| Column | Source |
|--------|--------|
| AdDSWE Total Inundation (JASO / MJJASO) | this study (`generation/` + `analysis/`) |
| AdDSWE Open Water (MJJASO) | this study |
| Landsat-SWIR threshold (JASO) | Inman & Lyons (2020) |
| MODIS | Thito et al. (2016) |
| NOAA AVHRR | Gumbricht et al. (2004); McCarthy et al. (2003) |
| MODIS-SWIR threshold | Wolski et al. (2017); Hofmann et al. (2024) |
| JRC Global Surface Water | Pekel et al. (2016) — via `generation/pekel_watersurface.ipynb` |
| Total Precipitation (mm) | CHIRPS — via `generation/chirps_precipitation.ipynb` |
| Annual Peak Discharge at Mohembo | Global Runoff Data Centre (GRDC) |

## Citation

If you use this code or the AdDSWE products, please cite:

- **Manuscript:** Rees, H., Crompton, O., Lourenço, M., WinklerPrins, L., Larsen, L.,
  Caylor, K., Ganti, V. *A 41-year satellite record reveals abrupt redistribution of the
  Okavango Delta's flood pulse since 2012.* (in preparation; DOI to be added upon publication).
- **Data:** Dryad, https://doi.org/10.5061/dryad.msbcc2gct
- **Code:** this repository (archived release DOI to be added via Zenodo).

## License

Code is released under the MIT License (see `LICENSE`). The archived data on Dryad are
released under CC0 1.0.
