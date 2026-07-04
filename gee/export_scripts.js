/**
 * export_scripts.js
 * Exports Hansen Global Forest Watch tree cover loss/gain data.
 * Combines with NDVI to add a structural deforestation/recovery dimension.
 */

var AOI = ee.Geometry.Rectangle([138.0, -37.0, 146.0, -30.0]);

// Hansen v1.11 — annual tree cover loss year (0 = no loss, 1-23 = year of loss)
var hansen = ee.Image('UMD/hansen/global_forest_change_2023_v1_11');

var treeCover2000 = hansen.select('treecover2000');  // % canopy cover in 2000
var lossYear      = hansen.select('lossyear');        // Year of loss (2001-2023)
var gain          = hansen.select('gain');            // Gain 2000-2012 (binary)

// Mask to study region
var cover = treeCover2000.clip(AOI);
var loss  = lossYear.clip(AOI);

// Export tree cover (for baseline)
Export.image.toDrive({
  image:          cover,
  description:    'hansen_treecover_2000_baseline',
  folder:         'GEE_NDVI_Exports',
  fileNamePrefix: 'hansen_cover_2000',
  region:         AOI,
  scale:          500,
  crs:            'EPSG:4326',
  maxPixels:      1e10
});

// Export loss year raster
Export.image.toDrive({
  image:          loss,
  description:    'hansen_loss_year_2001_2023',
  folder:         'GEE_NDVI_Exports',
  fileNamePrefix: 'hansen_lossyear',
  region:         AOI,
  scale:          500,
  crs:            'EPSG:4326',
  maxPixels:      1e10
});

print('Hansen exports queued. Check Tasks tab.');