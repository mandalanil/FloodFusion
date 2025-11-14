/*******************************************************************************
 *
 * Title: Advanced Flood Mapping and Analysis Application (Final Polished Version)
 *
 * Description: An interactive GEE application for flood mapping using a
 * fused Sentinel-1 and Sentinel-2 approach with a Random Forest classifier.
 * This version extracts and displays only the flooded pixels.
 *
 * NOTE: This code assumes you have imported the 'Nepal' asset in your Imports tab.
 *
 ******************************************************************************/

//================================================================================
// === CONFIGURATION & CONSTANTS ===
//================================================================================

var CONFIG = {
  // Data sources
  S1_COLLECTION: 'COPERNICUS/S1_GRD',
  S2_COLLECTION: 'COPERNICUS/S2_SR',
  DEM: 'USGS/SRTMGL1_003',

  // Default analysis parameters
  DEFAULT_START_DATE: '2021-06-01',
  DEFAULT_END_DATE: '2021-07-31',
  DEFAULT_S1_ORBIT: 'DESCENDING',
  DEFAULT_RF_TREES: 500,
  DEFAULT_TRAINING_SPLIT: 0.7, // 70% for training, 30% for validation
  DEFAULT_SLOPE_THRESHOLD: 5, // degrees
  DEFAULT_CONNECTIVITY_THRESHOLD: 8, // pixels
  
  // --- USER-DEFINED DEFAULTS ---
  DEFAULT_TRAINING_ASSET: 'users/srijal2023/Melamchi_points_water', // <-- Set your default asset
  DEFAULT_CLASS_COLUMN: 'Planet_flo', // <-- Set your default class property name
  // -----------------------------

  // Visualization parameters
  VIS_S1_VV: {
    min: -25,
    max: 0,
    bands: ['VV_Filtered']
  },
  VIS_S1_FALSE_COLOR: {
    min: [-20, -25, 0.5],
    max: [0, -5, 5],
    bands: ['VV_Filtered', 'VH_Filtered', 'Ratio_Filtered'] // VV, VH, VV/VH Ratio
  },
  VIS_S2_RGB: {
    min: 0.0,
    max: 0.3,
    bands: ['B4', 'B3', 'B2']
  },
  VIS_CLASSIFICATION: {
    palette: ['#0000FF'] // Blue for Flood/Water
  },
  LEGEND_INFO: {
    'Flood/Water': '#0000FF'
  },
  AOI_STYLE: {
    color: 'red',
    fillColor: '00000000'
  }
};

//================================================================================
// === PROCESSING LIBRARY ===
//================================================================================

/**
 * Applies a Refined Lee speckle filter to a Sentinel-1 image.
 */
function refinedLee(img) {
  var bandNames = img.bandNames();
  var weights = ee.List.repeat(ee.List.repeat(1, 7), 7);
  var kernel = ee.Kernel.fixed(7, 7, weights, 3, 3, false);
  var mean = img.reduceNeighborhood(ee.Reducer.mean(), kernel);
  var variance = img.reduceNeighborhood(ee.Reducer.variance(), kernel);
  var center_pixel = img.select(bandNames);
  var b = variance.divide(mean.multiply(mean));
  var T1 = ee.Image(1.0).subtract(b).max(0);
  var T2 = T1.divide(ee.Image(1.0).add(T1));
  var new_pixel = center_pixel.multiply(T2).add(mean.multiply(ee.Image(1.0).subtract(T2)));
  return new_pixel.rename(bandNames);
}


/**
 * Masks clouds in a Sentinel-2 SR image.
 */
function maskS2srClouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000)
    .select('B.*')
    .copyProperties(image, ['system:time_start']);
}

/**
 * Creates an analysis-ready Sentinel-1 composite.
 */
function processS1_ARD(startDate, endDate, aoi) {
  var s1col = ee.ImageCollection(CONFIG.S1_COLLECTION)
    .filterDate(startDate, endDate)
    .filterBounds(aoi)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('orbitProperties_pass', CONFIG.DEFAULT_S1_ORBIT))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));

  var s1_composite = ee.Image(ee.Algorithms.If(
    s1col.size().gt(0),
    s1col.median().clip(aoi),
    ee.Image().rename(['VV', 'VH'])
  ));

  var vv_filtered = ee.Image(ee.Algorithms.If(
    s1_composite.bandNames().contains('VV'),
    refinedLee(s1_composite.select('VV')).rename('VV_Filtered'),
    ee.Image().rename('VV_Filtered')
  ));
  
  var vh_filtered = ee.Image(ee.Algorithms.If(
    s1_composite.bandNames().contains('VH'),
    refinedLee(s1_composite.select('VH')).rename('VH_Filtered'),
    ee.Image().rename('VH_Filtered')
  ));

  var ratio = ee.Image(ee.Algorithms.If(
    vv_filtered.bandNames().size().gt(0).and(vh_filtered.bandNames().size().gt(0)),
    vv_filtered.divide(vh_filtered).rename('Ratio_Filtered'),
    ee.Image().rename('Ratio_Filtered')
  ));

  return ee.Image.cat([vv_filtered, vh_filtered, ratio]);
}

/**
 * Creates an analysis-ready Sentinel-2 composite.
 */
function processS2_ARD(startDate, endDate, aoi) {
  var s2col = ee.ImageCollection(CONFIG.S2_COLLECTION)
    .filterDate(startDate, endDate)
    .filterBounds(aoi)
    .map(maskS2srClouds);

  return ee.Image(ee.Algorithms.If(
      s2col.size().gt(0),
      s2col.median().clip(aoi),
      ee.Image()
  ));
}


//================================================================================
// === UI SETUP ===
//================================================================================

// --- Main Panels ---
var mainPanel = ui.Panel({
  style: {
    width: '350px',
    padding: '8px'
  }
});
ui.root.insert(0, mainPanel);

var map = ui.root.widgets().get(1);
map.style().set('cursor', 'crosshair');


// --- Title and Description ---
mainPanel.add(ui.Label({
  value: 'FloodFusion',
  style: {
    fontSize: '20px',
    fontWeight: 'bold'
  }
}));
mainPanel.add(ui.Label('This app performs a fused Sentinel-1 and Sentinel-2 analysis to map flood extent using a Random Forest classifier.'));

// --- Section 1: Analysis Scope ---
mainPanel.add(ui.Label({
  value: '1. Analysis Scope',
  style: {
    fontWeight: 'bold',
    fontSize: '16px',
    margin: '10px 0 4px 0'
  }
}));
mainPanel.add(ui.Label('Select Start & End Dates:'));
var startDateBox = ui.Textbox({
  value: CONFIG.DEFAULT_START_DATE,
  style: {
    width: '150px'
  }
});
var endDateBox = ui.Textbox({
  value: CONFIG.DEFAULT_END_DATE,
  style: {
    width: '150px'
  }
});
var datePanel = ui.Panel([startDateBox, endDateBox], ui.Panel.Layout.flow('horizontal'));
mainPanel.add(datePanel);


mainPanel.add(ui.Label('Define Area of Interest (AOI):'));
var drawingTools = map.drawingTools();
drawingTools.setShown(false);
while (drawingTools.layers().length() > 0) {
  drawingTools.layers().remove(drawingTools.layers().get(0));
}
var dummyGeometry = ui.Map.GeometryLayer({
  geometries: null,
  name: 'geometry',
  color: 'red'
});
drawingTools.layers().add(dummyGeometry);

drawingTools.onDraw(function() {
  drawingTools.stop();
});

function clearAndDraw(shape) {
  resetApp(false);  
  drawingTools.setShape(shape);
  drawingTools.draw();
}

// --- MODIFIED DRAW BUTTONS PANEL to include Default AOI ---
var drawButtons = ui.Panel({
  widgets: [
    ui.Button({
      label: '⬛ Rectangle',
      onClick: function() {
        clearAndDraw('rectangle');
      },
      style: {
        stretch: 'horizontal'
      }
    }),
    ui.Button({
      label: '🔺 Polygon',
      onClick: function() {
        clearAndDraw('polygon');
      },
      style: {
        stretch: 'horizontal'
      }
    }),
    ui.Button({ 
      label: '📍 Default AOI',
      onClick: loadDefaultAoi,
      style: {
        stretch: 'horizontal',
        backgroundColor: '#D6EAF8' 
      }
    })
  ],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {
    stretch: 'horizontal'
  }
});
mainPanel.add(drawButtons);
// --- END MODIFIED PANEL ---

// --- Section 2: Classification Parameters ---
mainPanel.add(ui.Label({
  value: '2. Classification Parameters',
  style: {
    fontWeight: 'bold',
    fontSize: '16px',
    margin: '10px 0 4px 0'
  }
}));
mainPanel.add(ui.Label('Training Data (Shapefile GEE Asset ID):'));
var trainingAssetBox = ui.Textbox({
  value: CONFIG.DEFAULT_TRAINING_ASSET, // Use the default value
  style: {
    width: '95%'
  }
});
mainPanel.add(trainingAssetBox);

var fetchColumnsButton = ui.Button({
  label: 'Fetch Columns',
  onClick: fetchAndPopulateColumns,
  style: {
    stretch: 'horizontal',
    margin: '4px 0 0 0'
  }
});
mainPanel.add(fetchColumnsButton);

mainPanel.add(ui.Label('Select Class Property Column:'));
var columnSelectDropdown = ui.Select({
  items: [],
  placeholder: 'Fetch columns first',
  disabled: true,
  style: {
    width: '95%'
  }
});
mainPanel.add(columnSelectDropdown);


mainPanel.add(ui.Label('Random Forest Trees:'));
var rfTreesBox = ui.Textbox({
  value: CONFIG.DEFAULT_RF_TREES,
  style: {
    width: '95%'
  }
});
mainPanel.add(rfTreesBox);

// --- Section 3: Post-Processing ---
mainPanel.add(ui.Label({
  value: '3. Post-Processing',
  style: {
    fontWeight: 'bold',
    fontSize: '16px',
    margin: '10px 0 4px 0'
  }
}));
mainPanel.add(ui.Label('Slope Threshold (degrees, 0-30):'));
var slopeSlider = ui.Slider({
  min: 0,
  max: 30,
  value: CONFIG.DEFAULT_SLOPE_THRESHOLD,
  step: 1,
  style: { stretch: 'horizontal' }
});
mainPanel.add(slopeSlider);

mainPanel.add(ui.Label('Min. Flood Patch Size (pixels, 0-50):'));
var connectivitySlider = ui.Slider({
  min: 0,
  max: 50,
  value: CONFIG.DEFAULT_CONNECTIVITY_THRESHOLD,
  step: 1,
  style: { stretch: 'horizontal' }
});
mainPanel.add(connectivitySlider);

// --- Execution and Status ---
var runButton = ui.Button({
  label: 'Run Analysis',
  onClick: runAnalysis,
  style: {
    stretch: 'horizontal',
    fontWeight: 'bold',
    backgroundColor: '#E74C3C',
    color: 'black',
    margin: '10px 0 0 0'
  }
});
var statusLabel = ui.Label({
  value: 'Status: Ready.',
  style: {
    margin: '8px 0',
    color: 'gray'
  }
});
mainPanel.add(runButton);
mainPanel.add(statusLabel);

// --- Results Panel ---
var resultsPanel = ui.Panel({
  style: {
    padding: '0 8px'
  }
});
mainPanel.add(resultsPanel);

// --- Legend Panel ---
var legendPanel = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px',
    backgroundColor: 'rgba(255, 255, 255, 0.8)'
  }
});

// --- Area Display Panel ---
var areaPanelLabel = ui.Label('Flooded Area: N/A', {
  fontWeight: 'bold',
  fontSize: '14px',
  margin: '0',
});
var areaPanel = ui.Panel({
  widgets: [areaPanelLabel],
  style: {
    position: 'bottom-right',
    padding: '8px 15px',
    backgroundColor: 'rgba(255, 255, 255, 0.8)'
  }
});
map.add(areaPanel);


function buildLegend(title, legendInfo) {
  legendPanel.clear();
  var legendTitle = ui.Label({
    value: title,
    style: {
      fontWeight: 'bold',
      fontSize: '18px',
      margin: '0 0 4px 0',
      padding: '0'
    }
  });
  legendPanel.add(legendTitle);
  for (var key in legendInfo) {
    var name = key;
    var color = legendInfo[key];
    var colorBox = ui.Label({
      style: {
        backgroundColor: color,
        padding: '8px',
        margin: '0 0 4px 0'
      }
    });
    var description = ui.Label({
      value: name,
      style: {
        margin: '0 0 4px 6px'
      }
    });
    var legendItemPanel = ui.Panel([colorBox, description], ui.Panel.Layout.flow('horizontal'));
    legendPanel.add(legendItemPanel);
  }
}
map.add(legendPanel);

//================================================================================
// === APP LOGIC & EVENT HANDLERS ===
//================================================================================

/**
 * Resets the application state.
 */
function resetApp(clearAoi) {
  map.layers().reset();
  if (clearAoi) {
      drawingTools.layers().get(0).geometries().reset();
  }
  resultsPanel.clear();
  legendPanel.clear();
  areaPanelLabel.setValue('Flooded Area: N/A');
  statusLabel.setValue('Status: Ready. Please draw an Area of Interest (AOI).').style().set('color', 'gray');
  runButton.setDisabled(false);
}

/**
 * Loads and displays the default AOI from the Nepal asset.
 * This function uses .evaluate() to handle the computed server-side geometry.
 */
function loadDefaultAoi() {
  // Check if the global 'Nepal' asset is available
  if (typeof Nepal === 'undefined') {
    handleError('The "Nepal" asset is not imported. Please add it to your Imports section.');
    return;
  }
  
  resetApp(true); // Reset the map and clear any existing AOI geometry
  drawingTools.stop();
  statusLabel.setValue('Status: Loading default AOI (Melamchi)...').style().set('color', 'orange');

  // 1. Define the computed Earth Engine Geometry object (server-side)
  // FIX: Using the corrected column name 'GaPa_NaPa'
  var filteredFeatures = Nepal.filter(ee.Filter.eq('GaPa_NaPa', 'Melamchi'));
  
  // Use ee.Algorithms.If to check if the filtered collection is empty, preventing a 'null' error
  var defaultAoiEE = ee.Algorithms.If(
    filteredFeatures.size().gt(0),
    filteredFeatures.geometry().bounds(),
    ee.Geometry.Point(85.58, 27.83).buffer(1) // Fallback to a small point near Melamchi
  );
  defaultAoiEE = ee.Geometry(defaultAoiEE); // Cast the result

  // 2. Compute the geometry on the server and use the result client-side.
  defaultAoiEE.evaluate(function(geojson, error) {
    // Note: The check for 'geojson.coordinates' ensures the default AOI was actually found,
    // otherwise, the fallback geometry is used.
    if (error || !geojson.coordinates) {
      handleError('Could not find Melamchi in the asset. Default point used. Check "GaPa_NaPa" value.');
    }
    
    // Get the geometry layer from the drawing tools.
    var geometryLayer = drawingTools.layers().get(0);
    
    // Reset the geometries array in the layer.
    geometryLayer.geometries().reset();
    
    // Add the new client-side GeoJSON object to the layer.
    geometryLayer.geometries().add(geojson);
    
    // Center the map on the new AOI (using the server-side object for centering).
    map.centerObject(defaultAoiEE, 12);

    statusLabel.setValue('Status: Default AOI (Melamchi) loaded. Ready to run.').style().set('color', 'blue');
  });
}


/**
 * Fetches property names from the user-provided asset.
 */
function fetchAndPopulateColumns() {
  fetchColumnsButton.setDisabled(true);
  statusLabel.setValue('Status: Fetching columns...').style().set('color', 'orange');
  columnSelectDropdown.setPlaceholder('Fetching...');

  var assetId = trainingAssetBox.getValue();
  if (!assetId) {
    handleError('Please enter a training data Asset ID first.');
    fetchColumnsButton.setDisabled(false);
    return;
  }

  var fc;
  try {
    fc = ee.FeatureCollection(assetId);
  } catch (e) {
    handleError('Invalid Asset ID. Could not load FeatureCollection.');
    fetchColumnsButton.setDisabled(false);
    return;
  }

  var first = fc.first();
  var propNames = first.propertyNames();

  propNames.evaluate(function(names, error) {
    if (error) {
      handleError('Could not fetch column names: ' + error);
      columnSelectDropdown.setPlaceholder('Fetch failed');
    } else {
      var userNames = names.filter(function(name) {
        return name !== 'system:index';
      });
      
      columnSelectDropdown.items().reset(userNames);
      columnSelectDropdown.setDisabled(false);
      
      // Check if the default column exists in the list and set it
      var defaultColumnExists = userNames.indexOf(CONFIG.DEFAULT_CLASS_COLUMN) !== -1;
      
      if (defaultColumnExists) {
        columnSelectDropdown.setValue(CONFIG.DEFAULT_CLASS_COLUMN);
        statusLabel.setValue('Status: Columns fetched. Default class selected.').style().set('color', 'blue');
      } else {
        columnSelectDropdown.setPlaceholder('Select a column');
        statusLabel.setValue('Status: Columns fetched. Select class property.').style().set('color', 'blue');
      }
    }
    fetchColumnsButton.setDisabled(false);
  });
}


/**
 * Main function to orchestrate the analysis workflow.
 */
function runAnalysis() {
  var aoi = drawingTools.layers().get(0).getEeObject();
  
  if (!aoi) {
    handleError('Please draw an Area of Interest (AOI) first.');
    return;
  }

  runButton.setDisabled(true);
  statusLabel.setValue('Status: Processing...').style().set('color', 'orange');
  resultsPanel.clear();
  legendPanel.clear();
  areaPanelLabel.setValue('Flooded Area: Processing...');
  
  var startDate = ee.Date(startDateBox.getValue());
  var endDate = ee.Date(endDateBox.getValue());
  var trainingAssetId = trainingAssetBox.getValue();
  var classColumn = columnSelectDropdown.getValue();
  var rfTrees = ee.Number.parse(String(rfTreesBox.getValue()));
  var slopeThreshold = slopeSlider.getValue();
  var connectivityThreshold = connectivitySlider.getValue();

  if (!trainingAssetId || !classColumn) {
    handleError('Please provide a Training Asset and select a Class Column.');
    return;
  }
  
  map.centerObject(aoi, 11);
  
  var empty = ee.Image().byte();
  var aoiOutline = empty.paint({
    featureCollection: ee.FeatureCollection(aoi),
    color: 1,
    width: 3
  });
  map.addLayer(aoiOutline, {palette: 'FF0000'}, 'Area of Interest');
  
  drawingTools.layers().get(0).geometries().reset();

  statusLabel.setValue('Status: Processing satellite data...');
  var s1_image = processS1_ARD(startDate, endDate, aoi);
  var s2_image = processS2_ARD(startDate, endDate, aoi);

  var bandCounts = ee.Dictionary({
    s1: s1_image.bandNames().size(),
    s2: s2_image.bandNames().size()
  });

  bandCounts.evaluate(function(counts, error) {
    if (error) {
      handleError('Could not verify input data: ' + error);
      return;
    }
    if (counts.s1 === 0 || counts.s2 === 0) {
      handleError('No Sentinel-1 or Sentinel-2 images found for the criteria.');
      return;
    }

    var s2_bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'];
    var stackedImage = s2_image.select(s2_bands).addBands(s1_image);

    statusLabel.setValue('Status: Loading training data...');
    var trainingDataRaw;
    try {
      trainingDataRaw = ee.FeatureCollection(trainingAssetId);
    } catch (e) {
      handleError('Could not load training data. Check Asset ID.');
      return;
    }

    statusLabel.setValue('Status: Sampling training data...');
    var allSampledPoints = stackedImage.sampleRegions({
      collection: trainingDataRaw,
      properties: [classColumn],
      scale: 10,
      tileScale: 8
    });

    allSampledPoints.size().evaluate(function(numSampledPoints, error) {
      if (error || numSampledPoints === 0) {
        handleError('No valid training data found. Points may be in cloudy areas or outside image extent.');
        return;
      }
      
      var sampledWithRandom = allSampledPoints.randomColumn('random');
      var trainingSet = sampledWithRandom.filter(ee.Filter.lt('random', CONFIG.DEFAULT_TRAINING_SPLIT));
      var validationSet = sampledWithRandom.filter(ee.Filter.gte('random', CONFIG.DEFAULT_TRAINING_SPLIT));

      statusLabel.setValue('Status: Training classifier...');
      var classifier = ee.Classifier.smileRandomForest(rfTrees)
        .train({
          features: trainingSet,
          classProperty: classColumn,
          inputProperties: stackedImage.bandNames()
        });

      statusLabel.setValue('Status: Classifying image...');
      var classified = stackedImage.classify(classifier);

      var dem = ee.Image(CONFIG.DEM);
      var slope = ee.Terrain.slope(dem);
      var slopeMask = slope.lte(slopeThreshold);

      var floodPixels = classified.eq(1);

      var finalFloodPixels = ee.Image(ee.Algorithms.If(
          ee.Number(connectivityThreshold).gt(0),
          floodPixels.connectedPixelCount({maxSize: 100, eightConnected: true})
                       .gte(connectivityThreshold),
          floodPixels
      ));

      var finalClassification = classified.where(floodPixels, finalFloodPixels)
                                         .updateMask(slopeMask);
      
      var floodLayer = finalClassification.eq(1).selfMask();
      
      // Add result layers to map
      map.addLayer(s2_image, CONFIG.VIS_S2_RGB, 'Sentinel-2 RGB', false);
      map.addLayer(s1_image.select(['VV_Filtered', 'VH_Filtered', 'Ratio_Filtered']), CONFIG.VIS_S1_FALSE_COLOR, 'Sentinel-1 False Color', false);
      
      map.addLayer(floodLayer, CONFIG.VIS_CLASSIFICATION, 'Flooded Area');
      buildLegend('Legend', CONFIG.LEGEND_INFO);

      statusLabel.setValue('Status: Assessing accuracy...');
      var validation = validationSet.classify(classifier);
      var confusionMatrix = validation.errorMatrix(classColumn, 'classification');
      
      displayResults(confusionMatrix, aoi, finalClassification);
    });
  });
}

/**
 * Displays final results and accuracy metrics.
 */
function displayResults(confusionMatrix, aoi, finalClassification) {
  
  resultsPanel.clear();
  resultsPanel.add(ui.Label({
    value: 'Analysis Results',
    style: { fontWeight: 'bold', fontSize: '16px' }
  }));
  
  statusLabel.setValue('Status: Calculating area...').style().set('color', 'orange');
  
  var floodAreaImage = finalClassification.eq(1).multiply(ee.Image.pixelArea());
  var areaMetrics = ee.Dictionary({
      floodArea: floodAreaImage.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: aoi,
      scale: 10,
      maxPixels: 1e13,
      tileScale: 4
    }).get('classification'),
    aoiArea: aoi.area({'maxError': 1})
  });

  areaMetrics.evaluate(function(areaResults, error){
    if(error){
      handleError('Could not calculate area. ' + error);
      return;
    }

    var floodHa = ee.Number(areaResults.floodArea).divide(10000).getInfo();
    var aoiHa = ee.Number(areaResults.aoiArea).divide(10000).getInfo();

    resultsPanel.add(ui.Label('AOI Area: ' + aoiHa.toFixed(2) + ' ha'));
    resultsPanel.add(ui.Label('Mapped Flood Area: ' + floodHa.toFixed(2) + ' ha'));
    areaPanelLabel.setValue('Flooded Area: ' + floodHa.toFixed(2) + ' ha');

    statusLabel.setValue('Status: Finalizing results...').style().set('color', 'orange');
    
    confusionMatrix.accuracy().evaluate(function(accuracy, error) {
      if (error) {
        handleError('Could not compute accuracy metrics.');
        return;
      }
      
      resultsPanel.add(ui.Label('Accuracy Assessment', { fontWeight: 'bold', margin: '8px 0 4px 0' }));
      resultsPanel.add(ui.Label('Overall Accuracy: ' + (accuracy * 100).toFixed(2) + '%'));
      
      // Calculate and display Kappa Coefficient
      confusionMatrix.kappa().evaluate(function(kappa, kappaError) {
        if (!kappaError) {
          resultsPanel.add(ui.Label('Kappa Coefficient: ' + kappa.toFixed(3)));
        } else {
          print('Kappa Error:', kappaError);
        }
      });
      
      var floodDownloadLayer = finalClassification.eq(1).selfMask();
      floodDownloadLayer.getDownloadURL({
        name: 'flood_area_extraction',
        region: aoi,
        scale: 10,
        format: 'GEO_TIFF'
      }, function(url, failure) {
          if(failure){
            resultsPanel.add(ui.Label('Download Error: ' + failure, {color: 'red'}));
          } else {
            var downloadLink = ui.Label({
              value: 'Download Flood Area Map (GeoTIFF)',
              style: { color: 'blue', textDecoration: 'underline', margin: '8px 0' },
              targetUrl: url
            });
            resultsPanel.add(downloadLink);
          }
          
          statusLabel.setValue('Status: Complete.').style().set('color', 'green');
          runButton.setDisabled(false);
          
          var layers = map.layers();
          var layersToRemove = [];
          layers.forEach(function(layer) {
            if (layer.getName() === 'Area of Interest') {
              layersToRemove.push(layer);
            }
          });
          layersToRemove.forEach(function(layer) {
            map.layers().remove(layer);
          });
          
          var finalAoiOutline = ee.Image().byte().paint({
            featureCollection: ee.FeatureCollection(aoi),
            color: 1,
            width: 3
          });
          map.addLayer(finalAoiOutline, {palette: 'red'}, 'Final AOI');
      });
    });
  });
}


/**
 * Handles application errors by updating the UI.
 */
function handleError(message) {
  statusLabel.setValue('Error: ' + message).style().set('color', 'red');
  runButton.setDisabled(false);
}

//================================================================================
// === INITIALIZATION ===
//================================================================================

// Center the map on Melamchi, Nepal by default.
var melamchiCoords = {lon: 85.58, lat: 27.83};
var defaultZoom = 12;
map.setCenter(melamchiCoords.lon, melamchiCoords.lat, defaultZoom);
statusLabel.setValue('Status: Ready. Please draw an Area of Interest (AOI).');