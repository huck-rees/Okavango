"""AdDSWE: Adapted Dynamic Surface Water Extent for the Okavango Delta.

Core algorithm (spectral indices, the five DSWE tests, the dynamic SWIR2 Test 6, the
decision tree, and morphological filtering) and Landsat monthly compositing/gap-filling,
extracted from the AdDSWE monthly generator notebook.

Confidence classes: 0 = not water, 1 = low, 2 = partial / water-under-canopy,
3 = moderate open water, 4 = high-confidence open water. Open water = classes 3-4;
total inundation = classes 1-4.

Call ``ee.Initialize()`` before using these functions.
"""

from .composites import (
    maskL8sr,
    getLandsatCollection,
    rescale,
    load_roi,
    get_filled_composite_before_dswe,
)
from .classification import (
    Mndwi,
    Mbsrv,
    Mbsrn,
    Ndvi,
    Awesh,
    find_bimodal_trough,
    calculate_dynamic_swir2_threshold,
    Dswe,
    morphological_filter,
    Dswe_with_Test6,
)

__all__ = [
    # composites
    "maskL8sr",
    "getLandsatCollection",
    "rescale",
    "load_roi",
    "get_filled_composite_before_dswe",
    # classification
    "Mndwi",
    "Mbsrv",
    "Mbsrn",
    "Ndvi",
    "Awesh",
    "find_bimodal_trough",
    "calculate_dynamic_swir2_threshold",
    "Dswe",
    "morphological_filter",
    "Dswe_with_Test6",
]
