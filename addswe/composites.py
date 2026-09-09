"""Landsat monthly compositing and temporal gap-filling for AdDSWE.

Builds cloud/shadow-masked, reflectance-rescaled Landsat Collection 2 monthly median
composites over a region of interest, with iterative temporal gap-filling that widens
the compositing window by one month at a time until a target spatial coverage is met.
The per-pixel number of months of expansion is returned as a QC ("expansion") mask.

Functions are transcribed verbatim from the AdDSWE monthly generator notebook; only the
module-level imports have been added. Call ``ee.Initialize()`` before using them.
"""

import ee
import geopandas as gpd
from datetime import datetime
from dateutil.relativedelta import relativedelta


def maskL8sr(image):
    """
    Masks out clouds and cloud shadows in Landsat 8/9 imagery using the BQA band.
    Uses bits 8–9 for cloud confidence and bits 10–11 for shadow confidence.

    Args:
        image: ee.Image, the input Landsat image.

    Returns:
        ee.Image: The masked image with cloud and shadow pixels removed.
    """
    qa = image.select('BQA')

    # Cloud confidence (bits 8–9): mask if medium or high
    cloud_conf = qa.rightShift(8).bitwiseAnd(3)
    cloud_ok = cloud_conf.lte(1)

    # Shadow confidence (bits 10–11): mask if medium or high
    shadow_conf = qa.rightShift(10).bitwiseAnd(3)
    shadow_ok = shadow_conf.lte(1)

    mask = cloud_ok.And(shadow_ok)
    return image.updateMask(mask)


def getLandsatCollection():
    """
    Merges Landsat 5, 7, 8, and 9 collections (Tier 1, Collection 2 SR)
    and standardizes the band names for consistent analysis.

    Returns:
        ee.ImageCollection: A merged collection of standardized Landsat images.
    """
    # Define the band mappings for each Landsat version
    bn9 = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL']
    bn8 = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL']
    bn7 = ['SR_B1', 'SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL']
    bn5 = ['SR_B1', 'SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL']
    # Standardized names for all bands
    standard_bands = ['uBlue', 'Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2', 'BQA']

    # Fetch and rename bands in the Landsat collections
    ls5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").select(bn5, standard_bands)
    ls7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterDate('1999-04-15', '2003-05-30').select(bn7, standard_bands)
    ls8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").select(bn8, standard_bands)
    ls9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").select(bn9, standard_bands)

    # Merge all collections
    merged_collection = ls5.merge(ls7).merge(ls8).merge(ls9)

    return merged_collection


def rescale(image):
    """
    Rescale the reflectance values of Landsat imagery to allow for use of Landsat Collection 2 in DSWE.

    Parameters:
    image (ee.Image): The input image with bands to be rescaled.

    Returns:
    ee.Image: The image with rescaled bands added.
    """
    bns = ['uBlue', 'Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2']
    optical_bands = image.select(bns).multiply(0.0000275).add(-0.2)
    return image.addBands(optical_bands, None, True)


def load_roi(shapefile_path):
    """Load ROI from a shapefile and return as an EE Geometry."""
    gdf = gpd.read_file(shapefile_path)
    return ee.Geometry.Polygon(gdf.unary_union.__geo_interface__["coordinates"])


def get_filled_composite_before_dswe(start_date, end_date, roi, max_months=6, fill_threshold=0.95):
    def get_composite(s, e):
        collection = (getLandsatCollection()
                      .map(maskL8sr)
                      .map(rescale)
                      .filterDate(s, e)
                      .filterBounds(roi))

        ids = collection.aggregate_array('system:index').distinct().sort()
        id_string = ids.join(",")  # Convert to comma-separated string
        return collection.median().clip(roi).set({'source_ids': id_string})

    def get_coverage(image):
        valid = image.select("Green").mask()
        stats = valid.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=roi,
            scale=30,
            maxPixels=1e13
        )
        return ee.Number(stats.get("Green"))

    base = get_composite(start_date, end_date)
    last_used_composite = base

    filled = base
    expansion_mask = base.select("Green").mask().Not().multiply(0).rename("expansion_mask")
    coverage = get_coverage(filled)

    for offset in range(1, max_months + 1):
        if coverage.gte(fill_threshold).getInfo():
            break

        delta = relativedelta(months=offset)
        expanded_start = (datetime.strptime(start_date, "%Y-%m-%d") - delta).strftime("%Y-%m-%d")
        expanded_end = (datetime.strptime(end_date, "%Y-%m-%d") + delta).strftime("%Y-%m-%d")
        extra = get_composite(expanded_start, expanded_end)
        last_used_composite = extra

        update_mask = base.select("Green").mask().Not()
        new_expansion = extra.select("Green").mask().And(update_mask).multiply(offset).rename("expansion_mask")
        expansion_mask = expansion_mask.where(expansion_mask.eq(0).And(new_expansion.neq(0)), new_expansion)

        filled = extra.blend(filled)
        coverage = get_coverage(filled)

    # Attach source scene metadata
    filled = filled.set({'source_ids': last_used_composite.get('source_ids')})
    return filled, expansion_mask
