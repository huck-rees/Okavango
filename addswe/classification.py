"""Adapted Dynamic Surface Water Extent (AdDSWE) classification.

Implements the AdDSWE algorithm: the standard five DSWE spectral tests and decision
tree (Jones 2019 thresholds, unmodified), the dynamic SWIR2 "Test 6" that upgrades
confidence for water beneath vegetation, and a morphological filter that removes small,
isolated low-confidence blobs while preserving any blob containing open water.

Confidence classes: 0 = not water, 1 = low, 2 = partial/water-under-canopy,
3 = moderate open water, 4 = high-confidence open water. Open water = classes 3–4;
total inundation = classes 1–4.

Functions are transcribed verbatim from the AdDSWE monthly generator notebook; only the
module-level imports have been added. Call ``ee.Initialize()`` before using them.
"""

import os
import logging

import ee
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, argrelextrema
import matplotlib.pyplot as plt


# Normalized Difference Water Index (MNDWI)
def Mndwi(image):
    """
    Calculate the Modified Normalized Difference Water Index (MNDWI) for a given image.

    Parameters:
    image (ee.Image): The input image.

    Returns:
    ee.Image: The resulting image with the MNDWI band named 'mndwi'.
    """
    return image.normalizedDifference(['Green', 'Swir1']).rename('mndwi')


# Modified Bare Soil Reflectance Variables
def Mbsrv(image):
    """
    Calculate the Modified Bare Soil Reflectance Variable (MBSRV) for a given image.

    Parameters:
    image (ee.Image): The input image.

    Returns:
    ee.Image: The resulting image with the MBSRV band named 'mbsrv'.
    """
    return image.select(['Green']).add(image.select(['Red'])).rename('mbsrv')


def Mbsrn(image):
    """
    Calculate the Modified Bare Soil Reflectance Normalized (MBSRN) for a given image.

    Parameters:
    image (ee.Image): The input image.

    Returns:
    ee.Image: The resulting image with the MBSRN band named 'mbsrn'.
    """
    return image.select(['Nir']).add(image.select(['Swir1'])).rename('mbsrn')


# Normalized Difference Vegetation Index (NDVI)
def Ndvi(image):
    """
    Calculate the Normalized Difference Vegetation Index (NDVI) for a given image.

    Parameters:
    image (ee.Image): The input image.

    Returns:
    ee.Image: The resulting image with the NDVI band named 'ndvi'.
    """
    return image.normalizedDifference(['Nir', 'Red']).rename('ndvi')


# Automated Water Extraction Index (AWESH)
def Awesh(image):
    """
    Calculate the Automated Water Extraction Index (AWEsh) for a given image.

    Parameters:
    image (ee.Image): The input image with the necessary bands for MBSRN calculation.

    Returns:
    ee.Image: The resulting image with the AWEsh band named 'awesh'.
    """
    return image.expression(
        'Blue + 2.5 * Green + (-1.5) * mbsrn + (-0.25) * Swir2',
        {
            'Blue': image.select(['Blue']),
            'Green': image.select(['Green']),
            'mbsrn': Mbsrn(image).select(['mbsrn']),
            'Swir2': image.select(['Swir2'])
        }
    ).rename('awesh')


def find_bimodal_trough(histogram_data, band_name='SWIR1', smoothing_sigma=2):
    """
    Find the trough (local minimum) between two peaks in a bimodal distribution.

    This implements the concept from Inman & Lyons (2020) of finding the natural
    boundary between wet and dry pixels in the SWIR reflectance histogram. When
    you plot SWIR values for a wetland image, you typically see two "humps"
    (peaks): one for water/wet areas (low reflectance) and one for dry land
    (high reflectance). The valley (trough) between these peaks represents the
    natural separation point.

    Parameters:
    -----------
    histogram_data : dict
        GEE histogram output with structure:
        {'BandName': {'bucketMeans': [...], 'histogram': [...]}}
    band_name : str
        Name of the band ('SWIR1' or 'SWIR2')
    smoothing_sigma : float
        Gaussian smoothing parameter to reduce noise in the histogram
        Higher values = smoother curve but may miss subtle features

    Returns:
    --------
    float : The reflectance value at the trough (threshold)
    dict : Diagnostic information about the distribution
    """

    # Extract histogram data
    band_data = histogram_data[band_name]
    means = np.array(band_data['bucketMeans'])  # Reflectance values (x-axis)
    counts = np.array(band_data['histogram'])    # Pixel counts (y-axis)

    # Smooth the histogram to reduce noise
    # Think of this like drawing a smooth curve through scattered points
    counts_smooth = gaussian_filter1d(counts, sigma=smoothing_sigma)

    # Find peaks (the two "humps" in the histogram)
    # prominence ensures we only find significant peaks, not small bumps
    peaks, peak_properties = find_peaks(
        counts_smooth,
        prominence=np.max(counts_smooth) * 0.1  # Peak must be at least 10% of tallest peak
    )

    # If we don't find two clear peaks, fall back to percentile method
    if len(peaks) < 2:
        # Use the 30th percentile as a conservative wet/dry boundary
        cumsum = np.cumsum(counts)
        total = cumsum[-1]
        threshold_idx = np.where(cumsum >= total * 0.30)[0][0]
        threshold = means[threshold_idx]

        diagnostics = {
            'method': 'percentile_fallback',
            'threshold': threshold,
            'peaks_found': len(peaks),
            'wet_mode': None,
            'dry_mode': None,
            'reason': 'Bimodal structure not clear - using 30th percentile'
        }

        return threshold, diagnostics

    # Sort peaks by reflectance value (left to right on histogram)
    peak_indices = peaks[np.argsort(means[peaks])]

    # The first peak (leftmost) = wet mode (low reflectance)
    # The second peak (rightmost) = dry mode (high reflectance)
    wet_peak_idx = peak_indices[0]
    dry_peak_idx = peak_indices[1] if len(peak_indices) > 1 else peak_indices[0]

    # Find the lowest point (trough) between the two peaks
    search_range = counts_smooth[wet_peak_idx:dry_peak_idx+1]
    local_minima = argrelextrema(search_range, np.less)[0]

    if len(local_minima) > 0:
        # Take the deepest minimum (lowest point in the valley)
        trough_local_idx = local_minima[np.argmin(search_range[local_minima])]
        trough_idx = wet_peak_idx + trough_local_idx
        threshold = means[trough_idx]
    else:
        # Fallback: use midpoint between the two peaks
        threshold = (means[wet_peak_idx] + means[dry_peak_idx]) / 2

    # Package diagnostic information
    diagnostics = {
        'method': 'bimodal_trough',
        'threshold': threshold,
        'wet_mode': means[wet_peak_idx],           # Reflectance of wet peak
        'wet_mode_count': int(counts[wet_peak_idx]),  # Height of wet peak
        'dry_mode': means[dry_peak_idx],           # Reflectance of dry peak
        'dry_mode_count': int(counts[dry_peak_idx]),  # Height of dry peak
        'peaks_found': len(peaks),
        'trough_position': (threshold - means[wet_peak_idx]) / (means[dry_peak_idx] - means[wet_peak_idx])  # 0-1 scale
    }

    return threshold, diagnostics


def calculate_dynamic_swir2_threshold(image, roi, min_swir2=0.04, max_swir2=0.15,
                                       save_plot=True, output_dir=None,
                                       year=None, month=None):
    """
    Calculate dynamic SWIR2 threshold for a given image based on its histogram.

    This function analyzes the distribution of SWIR2 reflectance values across
    your study area and finds the natural separation between wet and dry pixels.
    The threshold is scene-specific and adapts to seasonal flooding conditions.

    This simplified version returns only SWIR2 threshold for use in Test 6
    (vegetated inundation enhancement).

    Parameters:
    -----------
    image : ee.Image
        The Landsat composite image
    roi : ee.Geometry
        Region of interest (your study area boundary)
    min_swir2, max_swir2 : float
        Safety constraints on SWIR2 threshold
        Default range: 0.04 to 0.15 (4% to 15% reflectance)
    save_plot : bool, optional (default=True)
        Whether to save histogram plot with threshold
    output_dir : str, optional (default=None)
        Directory to save plots. If None, saves to current directory
    year : int, optional
        Year for plot filename
    month : int, optional
        Month for plot filename

    Returns:
    --------
    float : The calculated SWIR2 threshold value
    """

    # Extract SWIR2 band
    swir2 = image.select(['Swir2'])

    # Get histogram from Google Earth Engine
    hist_dict = swir2.reduceRegion(
        reducer=ee.Reducer.histogram(maxBuckets=100),
        geometry=roi,
        scale=30,
        maxPixels=1e13
    ).getInfo()

    # Prepare histogram data for analysis
    swir2_hist = {'Swir2': hist_dict['Swir2']}

    # Find the trough (natural boundary) in histogram
    swir2_threshold, swir2_diag = find_bimodal_trough(swir2_hist, 'Swir2')

    # Apply safety constraints to prevent unreasonable values
    swir2_threshold_clipped = np.clip(swir2_threshold, min_swir2, max_swir2)

    # Create plot if requested
    if save_plot:
        import matplotlib.pyplot as plt
        import os

        # Extract histogram data
        means = np.array(hist_dict['Swir2']['bucketMeans'])
        counts = np.array(hist_dict['Swir2']['histogram'])

        # Create figure
        plt.figure(figsize=(10, 6))
        plt.bar(means, counts, width=(means[1] - means[0]) * 0.8,
                color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.5)

        # Add threshold line
        plt.axvline(swir2_threshold_clipped, color='red', linestyle='--',
                   linewidth=2, label=f'Threshold = {swir2_threshold_clipped:.4f}')

        # Add wet and dry mode lines if available
        if swir2_diag.get('wet_mode') is not None:
            plt.axvline(swir2_diag['wet_mode'], color='blue', linestyle=':',
                       linewidth=1.5, alpha=0.7, label=f'Wet Mode = {swir2_diag["wet_mode"]:.4f}')
        if swir2_diag.get('dry_mode') is not None:
            plt.axvline(swir2_diag['dry_mode'], color='orange', linestyle=':',
                       linewidth=1.5, alpha=0.7, label=f'Dry Mode = {swir2_diag["dry_mode"]:.4f}')

        # Labels and title
        plt.xlabel('SWIR2 Reflectance', fontsize=12, fontweight='bold')
        plt.ylabel('Pixel Count', fontsize=12, fontweight='bold')

        if year and month:
            plt.title(f'SWIR2 Histogram with Dynamic Threshold\n{year}-{month:02d}',
                     fontsize=14, fontweight='bold')
        else:
            plt.title('SWIR2 Histogram with Dynamic Threshold',
                     fontsize=14, fontweight='bold')

        plt.legend(loc='upper right', fontsize=10)
        plt.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()

        # Determine output directory
        if output_dir is None:
            output_dir = os.getcwd()
        else:
            os.makedirs(output_dir, exist_ok=True)

        # Create filename
        if year and month:
            base_filename = f'SWIR2_threshold_{year}_{month:02d}'
        else:
            base_filename = 'SWIR2_threshold'

        # Save as PNG and JPEG
        png_path = os.path.join(output_dir, f'{base_filename}.png')
        jpeg_path = os.path.join(output_dir, f'{base_filename}.jpeg')

        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        plt.savefig(jpeg_path, dpi=300, bbox_inches='tight', format='jpeg')
        plt.close()

        print(f"Plots saved:")
        print(f"  PNG:  {png_path}")
        print(f"  JPEG: {jpeg_path}")

    return float(swir2_threshold_clipped)


# Decision Tree for Surface Water Extent (DSWE)
def Dswe(image):
    """
    Calculate the Decision Tree for Surface Water Extent (DSWE) for a given image.

    Parameters:
    image (ee.Image): The input image with bands required for the DSWE calculation.

    Returns:
    ee.Image: The resulting image with the DSWE classification named 'dswe'.
    """
    mndwi = Mndwi(image)
    mbsrv = Mbsrv(image)
    mbsrn = Mbsrn(image)
    awesh = Awesh(image)
    swir1 = image.select(['Swir1'])
    nir = image.select(['Nir'])
    ndvi = Ndvi(image)
    blue = image.select(['Blue'])
    swir2 = image.select(['Swir2'])

    # Decision tree thresholds
    t1 = mndwi.gt(0.124)
    t2 = mbsrv.gt(mbsrn)
    t3 = awesh.gt(0)
    t4 = (mndwi.gt(-0.44)
          .And(swir1.lt(0.09))
          .And(nir.lt(0.15))
          .And(ndvi.lt(0.7)))
    t5 = (mndwi.gt(-0.5)
          .And(blue.lt(0.1))
          .And(swir1.lt(0.3))
          .And(swir2.lt(0.1))
          .And(nir.lt(0.25)))

    # Combine results using weights to create unique classes
    t = t1.add(t2.multiply(10)).add(t3.multiply(100)).add(t4.multiply(1000)).add(t5.multiply(10000))

    # Define DSWE classification levels
    noWater = t.eq(0).Or(t.eq(1)).Or(t.eq(10)).Or(t.eq(100)).Or(t.eq(1000))
    hWater = t.eq(1111).Or(t.eq(10111)).Or(t.eq(11011)).Or(t.eq(11101)).Or(t.eq(11110)).Or(t.eq(11111))
    mWater = (t.eq(111).Or(t.eq(1011)).Or(t.eq(1101)).Or(t.eq(1110))
              .Or(t.eq(10011)).Or(t.eq(10101)).Or(t.eq(10110))
              .Or(t.eq(11001)).Or(t.eq(11010)).Or(t.eq(11100)))
    pWetland = t.eq(11000)
    lWater = (t.eq(11).Or(t.eq(101)).Or(t.eq(110)).Or(t.eq(1001))
              .Or(t.eq(1010)).Or(t.eq(1100)).Or(t.eq(10000))
              .Or(t.eq(10001)).Or(t.eq(10010)).Or(t.eq(10100)))

    # Assign classification levels to DSWE
    iDswe = (noWater.multiply(0)
             .add(hWater.multiply(4))
             .add(mWater.multiply(3))
             .add(pWetland.multiply(2))
             .add(lWater.multiply(1)))

    return iDswe.rename(['dswe'])


def morphological_filter(dswe_image, size_threshold=50, max_class_threshold=2,
                         roi=None, return_diagnostics=True):
    """
    Remove isolated blobs of low-confidence water classifications that are completely
    surrounded by dry pixels. Preserves any blob containing high-confidence water pixels
    (class > 2) regardless of size, ensuring the main floodplain "megablob" is retained.

    A blob is removed if:
    1. It is smaller than size_threshold (in pixels), AND
    2. All pixels in the blob are class 1 or 2 (max value <= max_class_threshold)

    Parameters:
    -----------
    dswe_image : ee.Image
        DSWE classification image (0=no water, 1=low, 2=partial, 3=moderate, 4=high)
    size_threshold : int, optional (default=50)
        Maximum blob size (in pixels) eligible for removal
        At 30m resolution: 50 pixels = 4.5 hectares
    max_class_threshold : int, optional (default=2)
        Maximum DSWE class value - blobs with ANY pixel > this are always preserved
        Default of 2 means blobs containing class 3 or 4 are kept regardless of size
    roi : ee.Geometry, optional
        Region of interest for calculating diagnostics (if return_diagnostics=True)
    return_diagnostics : bool, optional (default=True)
        Whether to calculate and return diagnostic statistics

    Returns:
    --------
    ee.Image : Filtered DSWE classification
    dict (optional) : Diagnostic statistics if return_diagnostics=True

    Notes:
    ------
    - Uses 8-connectivity (diagonal neighbors count as connected)
    - Blobs touching the image boundary are treated the same as interior blobs
    - The main Okavango floodplain is preserved because it contains class 3-4 pixels
    """

    # Step 1: Create binary mask of all wet pixels (any class > 0)
    wet_mask = dswe_image.gt(0)

    # Step 2: Label connected components
    # Use 8-connectivity (diagonal neighbors connect) to avoid fragmenting natural wetlands
    # maxSize is the tile size for processing - must be <= 1024
    labeled = wet_mask.connectedComponents(
        connectedness=ee.Kernel.square(1),  # 8-connectivity
        maxSize=256  # Tile size for processing (not max blob size!)
    )

    # Step 3: Add labels band to DSWE image for connected components reduction
    dswe_with_labels = dswe_image.addBands(labeled.select('labels'))

    # Step 4: Calculate statistics for each blob
    # reduceConnectedComponents maps the blob-level statistic back to every pixel in that blob
    # So blob_max will be an image where each pixel has its blob's maximum DSWE value
    blob_max = dswe_with_labels.reduceConnectedComponents(
        reducer=ee.Reducer.max(),
        labelBand='labels'
    )

    blob_count = dswe_with_labels.reduceConnectedComponents(
        reducer=ee.Reducer.count(),
        labelBand='labels'
    )

    # Step 5: Identify pixels belonging to blobs that should be removed
    # Each pixel now knows its blob's size and max class value
    # A pixel should be removed if its blob is small AND low-confidence
    small_blobs = blob_count.select('dswe').lt(size_threshold)
    low_confidence_only = blob_max.select('dswe').lte(max_class_threshold)
    removal_mask = small_blobs.And(low_confidence_only)

    # Step 6: Apply filter
    # Set pixels in removable blobs to 0 (no water)
    # All other pixels remain unchanged
    filtered_dswe = dswe_image.where(removal_mask, 0)

    # Step 7: Calculate diagnostics if requested
    diagnostics = None
    if return_diagnostics and roi is not None:
        try:
            # Count total pixels changed
            changed_pixels = dswe_image.neq(filtered_dswe).And(dswe_image.mask())

            # Calculate statistics
            original_stats = dswe_image.gt(0).reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=roi,
                scale=30,
                maxPixels=1e13
            ).getInfo()

            filtered_stats = filtered_dswe.gt(0).reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=roi,
                scale=30,
                maxPixels=1e13
            ).getInfo()

            original_water_pixels = original_stats.get('dswe', 0)
            filtered_water_pixels = filtered_stats.get('dswe', 0)
            pixels_removed = original_water_pixels - filtered_water_pixels
            area_removed_km2 = pixels_removed * 0.0009  # 30m pixels = 0.0009 km²

            # Count pixels by original class that were removed
            class_1_removed = dswe_image.eq(1).And(changed_pixels).reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=roi,
                scale=30,
                maxPixels=1e13
            ).getInfo().get('dswe', 0)

            class_2_removed = dswe_image.eq(2).And(changed_pixels).reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=roi,
                scale=30,
                maxPixels=1e13
            ).getInfo().get('dswe', 0)

            # These should always be zero if filter works correctly
            class_3_removed = dswe_image.eq(3).And(changed_pixels).reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=roi,
                scale=30,
                maxPixels=1e13
            ).getInfo().get('dswe', 0)

            class_4_removed = dswe_image.eq(4).And(changed_pixels).reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=roi,
                scale=30,
                maxPixels=1e13
            ).getInfo().get('dswe', 0)

            diagnostics = {
                'pixels_removed': int(pixels_removed) if pixels_removed else 0,
                'area_removed_km2': round(area_removed_km2, 2) if area_removed_km2 else 0.0,
                'class_1_pixels_removed': int(class_1_removed) if class_1_removed else 0,
                'class_2_pixels_removed': int(class_2_removed) if class_2_removed else 0,
                'class_3_pixels_removed': int(class_3_removed) if class_3_removed else 0,
                'class_4_pixels_removed': int(class_4_removed) if class_4_removed else 0,
                'size_threshold_used': size_threshold,
                'max_class_threshold_used': max_class_threshold,
                'percent_water_removed': round(100 * pixels_removed / original_water_pixels, 2) if original_water_pixels > 0 else 0.0
            }

        except Exception as e:
            logging.warning(f"Could not calculate diagnostics: {e}")
            diagnostics = {
                'pixels_removed': 0,
                'area_removed_km2': 0.0,
                'class_1_pixels_removed': 0,
                'class_2_pixels_removed': 0,
                'class_3_pixels_removed': 0,
                'class_4_pixels_removed': 0,
                'size_threshold_used': size_threshold,
                'max_class_threshold_used': max_class_threshold,
                'percent_water_removed': 0.0,
                'error': str(e)
            }

    # Add metadata to filtered image
    filtered_dswe = filtered_dswe.set({
        'morphological_filter_applied': True,
        'blob_size_threshold': size_threshold,
        'blob_max_class_threshold': max_class_threshold
    })

    if return_diagnostics:
        return filtered_dswe, diagnostics
    else:
        return filtered_dswe


def Dswe_with_Test6(image, roi, min_swir2=0.04, max_swir2=0.15,
                     save_plot=True, output_dir=None, year=None, month=None):
    """
    Calculate DSWE classification with Test 6 enhancement for vegetated inundation.

    This function applies the standard DSWE algorithm, then upgrades class confidence
    for pixels that also pass Test 6 (SWIR2 < dynamic threshold). This enhancement
    is designed to better capture water beneath dense vegetation (e.g., papyrus swamps)
    where traditional spectral indices may fail but SWIR2 still indicates moisture.

    Upgrade logic:
    - Class 0 (No Water) + Test 6 pass → Class 1 (Low Water)
    - Class 1 (Low Water) + Test 6 pass → Class 2 (Partial Wetland)
    - Class 2 (Partial Wetland) + Test 6 pass → Class 3 (Moderate Water)
    - Class 3 (Moderate Water) + Test 6 pass → Class 4 (High Water)
    - Class 4 (High Water) → Remains Class 4 (no change)

    Parameters:
    -----------
    image : ee.Image
        Landsat composite with standard bands
    roi : ee.Geometry
        Region of interest for threshold calculation
    min_swir2, max_swir2 : float
        Safety constraints on SWIR2 threshold
        Default range: 0.04 to 0.15 (4% to 15% reflectance)
    save_plot : bool, optional (default=True)
        Whether to save SWIR2 histogram plot
    output_dir : str, optional (default=None)
        Directory to save plots
    year : int, optional
        Year for metadata and plot filename
    month : int, optional
        Month for metadata and plot filename

    Returns:
    --------
    tuple : (upgraded_classification, original_classification, swir2_threshold)
        - upgraded_classification: ee.Image with Test 6 upgrades applied
        - original_classification: ee.Image with standard DSWE (for comparison)
        - swir2_threshold: float, the calculated SWIR2 threshold value
    """

    # Step 1: Run standard DSWE algorithm
    original_dswe = Dswe(image)

    # Step 2: Calculate dynamic SWIR2 threshold
    swir2_threshold = calculate_dynamic_swir2_threshold(
        image, roi,
        min_swir2=min_swir2,
        max_swir2=max_swir2,
        save_plot=save_plot,
        output_dir=output_dir,
        year=year,
        month=month
    )

    # Step 3: Create Test 6 (SWIR2 < threshold)
    swir2 = image.select(['Swir2'])
    test6 = swir2.lt(swir2_threshold)

    # Step 4: Apply upgrade logic
    # Start with original classification
    upgraded_dswe = original_dswe

    # Upgrade class 0 → 1 if Test 6 passes
    upgraded_dswe = upgraded_dswe.where(
        original_dswe.eq(0).And(test6),
        1
    )

    # Upgrade class 1 → 2 if Test 6 passes
    upgraded_dswe = upgraded_dswe.where(
        original_dswe.eq(1).And(test6),
        2
    )

    # Upgrade class 2 → 3 if Test 6 passes
    upgraded_dswe = upgraded_dswe.where(
        original_dswe.eq(2).And(test6),
        3
    )

    # Upgrade class 3 → 4 if Test 6 passes
    upgraded_dswe = upgraded_dswe.where(
        original_dswe.eq(3).And(test6),
        4
    )

    # Class 4 remains unchanged (no .where() operation needed)

    # Step 5: Add metadata to both images
    metadata = {
        'swir2_threshold': swir2_threshold,
        'test6_applied': True,
        'algorithm': 'DSWE_with_Test6_vegetated_enhancement'
    }

    if year is not None:
        metadata['year'] = year
    if month is not None:
        metadata['month'] = month

    upgraded_dswe = upgraded_dswe.set(metadata).rename(['dswe'])
    original_dswe = original_dswe.set({
        'swir2_threshold': swir2_threshold,
        'test6_applied': False,
        'algorithm': 'DSWE_standard'
    }).rename(['dswe'])

    return upgraded_dswe, original_dswe, swir2_threshold
