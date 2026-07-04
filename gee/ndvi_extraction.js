/**
 * ndvi_extraction.js
 * ==================
 * Project: Reforestation NDVI Recovery Analysis
 * Author:  Mambwe Nsanje
 * 
 * PURPOSE
 * -------
 * Extracts mean NDVI time series from NASA MODIS MOD13Q1 for each
 * Australian study site (treated vs control). Exports a single merged
 * CSV to Google Drive for download into data/processed/.
 *
 * PRODUCT DETAILS
 * ---------------
 * Collection : MODIS/061/MOD13Q1
 * Band       : NDVI (raw integer, scale factor = 0.0001 → range -1 to 1)
 * Resolution : 250 m spatial, 16-day temporal composite
 * Period     : 2010-01-01 → 2024-12-31
 * QA filter  : SummaryQA ∈ {0, 1} (good + marginal data only)
 *
 * HOW TO USE
 * ----------
 * 1. Paste this script into https://code.earthengine.google.com
 * 2. Click  Run  — check the Console for record counts
 * 3. Click  Tasks  tab → Run each queued export
 * 4. When complete, download from Google Drive →
 *    place CSV at:  data/processed/ndvi_modis_all_sites.csv
 *
 * UPDATE BOUNDING BOXES
 * ---------------------
 * Replace the placeholder bbox values below with real coordinates
 * from Coexist Australia project maps or the CAPAD protected areas
 * database: https://www.dcceew.gov.au/environment/land/nrs/science/capad
 */


// ─────────────────────────────────────────────────────────────────────────────
// 1. STUDY SITE DEFINITIONS
// ─────────────────────────────────────────────────────────────────────────────
// Format: ee.Geometry.Rectangle([minLon, minLat, maxLon, maxLat])
// Southern hemisphere → latitudes are negative

var SITES = {
  // ── Treated sites (known reforestation / rewilding intervention) ───────────
  coexist_mallee_treated: {
    geom:              ee.Geometry.Rectangle([140.5, -35.8, 141.5, -34.8]),
    type:              'treated',
    intervention_year: 2017,
    biome:             'semi_arid_mallee',
    notes:             'SA Mallee — native mallee scrub restoration'
  },
  mungo_rewild_treated: {
    geom:              ee.Geometry.Rectangle([143.0, -33.2, 144.0, -32.2]),
    type:              'treated',
    intervention_year: 2019,
    biome:             'semi_arid_mulga',
    notes:             'NSW Mungo region — rewilding project'
  },

  // ── Control sites (same biome, no known intervention) ─────────────────────
  coexist_mallee_control: {
    geom:              ee.Geometry.Rectangle([142.5, -35.8, 143.5, -34.8]),
    type:              'control',
    intervention_year: null,
    biome:             'semi_arid_mallee',
    notes:             'SA Mallee control — no intervention recorded'
  },
  mungo_rewild_control: {
    geom:              ee.Geometry.Rectangle([143.0, -31.5, 144.0, -30.5]),
    type:              'control',
    intervention_year: null,
    biome:             'semi_arid_mulga',
    notes:             'NSW semi-arid control — no intervention recorded'
  }
};


// ─────────────────────────────────────────────────────────────────────────────
// 2. QUALITY MASK FUNCTION
// ─────────────────────────────────────────────────────────────────────────────
/**
 * applyQAMask()
 * Retains only pixels where SummaryQA = 0 (good) or 1 (marginal).
 * Removes cloud-contaminated, snow/ice, and fill pixels.
 * Also applies the valid NDVI range mask (-2000 to 10000 raw integer).
 *
 * SummaryQA values:
 *   0 → Good data, use with confidence
 *   1 → Marginal data, useful but look at other QA info
 *   2 → Snow/ice — DO NOT USE for vegetation analysis
 *   3 → Cloudy   — DO NOT USE
 */
function applyQAMask(image) {
  var qa        = image.select('SummaryQA');
  var qaMask    = qa.lte(1);                          // Keep 0 and 1 only

  var ndvi      = image.select('NDVI');
  var rangeMask = ndvi.gte(-2000).and(ndvi.lte(10000)); // Valid raw range

  return image
    .updateMask(qaMask)
    .updateMask(rangeMask);
}


// ─────────────────────────────────────────────────────────────────────────────
// 3. LOAD AND PREPARE MODIS COLLECTION
// ─────────────────────────────────────────────────────────────────────────────
var modis = ee.ImageCollection('MODIS/061/MOD13Q1')
  .filterDate('2010-01-01', '2024-12-31')
  .select(['NDVI', 'SummaryQA'])
  .map(applyQAMask)
  .map(function(image) {
    // Apply scale factor: raw integer → true NDVI float (-1 to 1)
    return image
      .select('NDVI')
      .multiply(0.0001)
      .copyProperties(image, image.propertyNames());
  });

print('MODIS collection loaded. Image count:', modis.size());


// ─────────────────────────────────────────────────────────────────────────────
// 4. EXTRACT NDVI STATISTICS PER SITE
// ─────────────────────────────────────────────────────────────────────────────
/**
 * extractSiteStats()
 * For each MODIS 16-day image, computes spatial statistics over the
 * site bounding box and returns a FeatureCollection of one row per date.
 *
 * Outputs per row:
 *   date, year, month, doy        — temporal identifiers
 *   ndvi_mean, ndvi_median        — central tendency
 *   ndvi_std, ndvi_p25, ndvi_p75  — spread / IQR
 *   pixel_count                   — valid (non-masked) pixels
 *   site_id, site_type            — site metadata
 *   intervention_year, biome      — analysis metadata
 */
function extractSiteStats(siteName, siteMeta) {
  var geom = siteMeta.geom;

  var siteCollection = modis.map(function(image) {
    var stats = image.reduceRegion({
      reducer: ee.Reducer.mean()
                 .combine({ reducer2: ee.Reducer.median(),  sharedInputs: true })
                 .combine({ reducer2: ee.Reducer.stdDev(),  sharedInputs: true })
                 .combine({ reducer2: ee.Reducer.percentile([25, 75]), sharedInputs: true })
                 .combine({ reducer2: ee.Reducer.count(),   sharedInputs: true }),
      geometry:  geom,
      scale:     250,         // Native MODIS resolution — do not upsample
      maxPixels: 1e9,
      bestEffort: true        // Slightly loosens scale if needed to avoid timeout
    });

    // Pixel count: how many valid pixels contributed to this observation?
    // Low count = most pixels were masked (cloud/qa) → treat with caution
    var pixelCount = image.reduceRegion({
      reducer:  ee.Reducer.count(),
      geometry: geom,
      scale:    250,
      maxPixels: 1e9
    }).get('NDVI');

    return ee.Feature(null, {
      // ── Temporal ──────────────────────────────────────────────────────────
      'date':              image.date().format('YYYY-MM-dd'),
      'year':              image.date().get('year'),
      'month':             image.date().get('month'),
      'doy':               image.date().getRelative('day', 'year'),

      // ── NDVI statistics ───────────────────────────────────────────────────
      'ndvi_mean':         stats.get('NDVI_mean'),
      'ndvi_median':       stats.get('NDVI_median'),
      'ndvi_std':          stats.get('NDVI_stdDev'),
      'ndvi_p25':          stats.get('NDVI_p25'),
      'ndvi_p75':          stats.get('NDVI_p75'),
      'pixel_count':       pixelCount,

      // ── Site metadata ─────────────────────────────────────────────────────
      'site_id':           siteName,
      'site_type':         siteMeta.type,
      'intervention_year': siteMeta.intervention_year,
      'biome':             siteMeta.biome,
      'notes':             siteMeta.notes
    });
  });

  return siteCollection;
}


// ─────────────────────────────────────────────────────────────────────────────
// 5. LOOP OVER ALL SITES AND MERGE
// ─────────────────────────────────────────────────────────────────────────────
var siteNames = Object.keys(SITES);
var allFeatures = ee.FeatureCollection([]);

siteNames.forEach(function(name) {
  var fc = extractSiteStats(name, SITES[name]);
  allFeatures = allFeatures.merge(fc);
  print('Queued site:', name);
});

print('Total feature count (all sites × all dates):', allFeatures.size());
print('Sample record (first feature):', allFeatures.first());


// ─────────────────────────────────────────────────────────────────────────────
// 6. SANITY CHECK — PLOT NDVI FOR ONE SITE IN GEE CONSOLE
// ─────────────────────────────────────────────────────────────────────────────
// Visualise the Mallee treated site time series directly in GEE
// before exporting — confirms extraction is working as expected.

var mallee_ts = allFeatures
  .filter(ee.Filter.eq('site_id', 'coexist_mallee_treated'))
  .sort('date');

print('Mallee treated — first 5 rows:', mallee_ts.limit(5));

// Map layer: 2024 NDVI mosaic for visual inspection
var ndvi2024 = modis
  .filterDate('2024-01-01', '2024-12-31')
  .mean()
  .clip(ee.Geometry.Rectangle([138.0, -37.0, 146.0, -30.0]));

Map.addLayer(
  ndvi2024,
  { min: 0.05, max: 0.65, palette: ['#d73027', '#fee08b', '#91cf60', '#1a9850'] },
  'NDVI 2024 mean'
);

// Show site bounding boxes on map
siteNames.forEach(function(name) {
  var color = SITES[name].type === 'treated' ? '1a9850' : 'd73027';
  Map.addLayer(
    ee.FeatureCollection([ee.Feature(SITES[name].geom, { label: name })]),
    { color: color },
    name
  );
});

Map.centerObject(SITES['coexist_mallee_treated'].geom, 6);


// ─────────────────────────────────────────────────────────────────────────────
// 7. EXPORT TO GOOGLE DRIVE
// ─────────────────────────────────────────────────────────────────────────────
/**
 * IMPORTANT: After clicking Run above, go to the Tasks tab (top right)
 * and click Run next to the export task.
 *
 * The export lands in Google Drive → GEE_NDVI_Exports/
 * Download the CSV and place it at:
 *   data/processed/ndvi_modis_all_sites.csv
 *
 * Expected file size: ~2 MB for 4 sites × ~350 16-day periods
 */

Export.table.toDrive({
  collection:     allFeatures,
  description:    'ndvi_modis_all_sites_2010_2024',
  folder:         'GEE_NDVI_Exports',
  fileNamePrefix: 'ndvi_modis_all_sites',
  fileFormat:     'CSV',
  selectors: [
    // Column order matches what notebooks/01_data_exploration.ipynb expects
    'date', 'year', 'month', 'doy',
    'ndvi_mean', 'ndvi_median', 'ndvi_std', 'ndvi_p25', 'ndvi_p75',
    'pixel_count',
    'site_id', 'site_type', 'intervention_year', 'biome', 'notes'
  ]
});

print('──────────────────────────────────────────────────────');
print('Export task queued: ndvi_modis_all_sites_2010_2024');
print('→ Go to Tasks tab and click Run to start the export.');
print('──────────────────────────────────────────────────────');