/**
 * site_analysis.js
 * Generates annual NDVI composites for spatial visualisation.
 * Exports GeoTIFFs to Drive for use with rasterio/geopandas in notebooks.
 */

var YEARS    = [2010, 2015, 2019, 2022, 2024];
var AOI      = ee.Geometry.Rectangle([138.0, -37.0, 146.0, -30.0]); // SA + NSW

var modis = ee.ImageCollection('MODIS/061/MOD13Q1')
  .select('NDVI')
  .map(function(img) {
    // Quality mask: keep pixels with VI quality = 00 (good data)
    var qa = img.select('SummaryQA');
    var mask = qa.eq(0).or(qa.eq(1)); // 0=good, 1=marginal
    return img.updateMask(mask).multiply(0.0001);
  });

// Annual max-value composite (reduces cloud/atmospheric noise)
YEARS.forEach(function(year) {
  var annual = modis
    .filterDate(year + '-01-01', year + '-12-31')
    .max()  // Max-value composite
    .clip(AOI);
  
  Export.image.toDrive({
    image:       annual,
    description: 'ndvi_annual_' + year,
    folder:      'GEE_NDVI_Exports',
    fileNamePrefix: 'ndvi_' + year,
    region:      AOI,
    scale:       500,        // 500m for manageable file sizes
    crs:         'EPSG:4326',
    maxPixels:   1e10
  });
});

// Visualise in GEE map panel (for sanity-checking)
var viz = {min: 0, max: 0.8, palette: ['#d73027','#fee08b','#1a9850']};
Map.addLayer(modis.filterDate('2024-01-01','2024-12-31').max().clip(AOI), viz, 'NDVI 2024');
Map.centerObject(AOI, 6);