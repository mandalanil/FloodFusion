# FloodFusion: Advanced Flood Mapping and Analysis 🌊

**FloodFusion** is an interactive Google Earth Engine (GEE) application for mapping flood extent.  
It uses a data fusion approach that combines Sentinel-1 (SAR) and Sentinel-2 (optical) satellite imagery, classifying the data with a Random Forest algorithm to identify and display flooded pixels.

- **GEE App Link:** [**Click here to open the FloodFusion Application**](https://your-gee-app-url-here.com)  
- **Demo Training Data:** [**Sample Training Points (GEE Asset)**](https://earthengine.google.com/asset/users/your_username/your_sample_asset_name)

---

## Key Features 🚀

- **Data Fusion:** Combines the all-weather capabilities of **Sentinel-1** radar with the multispectral detail of **Sentinel-2** optical data for robust flood detection.
- **Machine Learning:** Uses a **Random Forest classifier** for accurate land and water discrimination.
- **Bring Your Own Data:** Accepts your own **GEE Feature Collection asset** for training, offering full flexibility.
- **Advanced Filtering:** Removes pixels on steep slopes (unlikely flood zones) and isolated noisy pixels, improving final accuracy.
- **Interactive UI:** A control panel lets you define your area, time frame, and analysis settings without coding.
- **Results & Export:** Calculates total flooded area (hectares), reports model accuracy, and allows **GeoTIFF download**.

---

## How to Use the App 🗺️

### Step 1: Define Your Scope
1. **Draw an AOI:** Use the `⬛ Rectangle` or `🔺 Polygon` tools to draw your analysis area on the map.
2. **Select Dates:** Enter the **Start Date** and **End Date** in `YYYY-MM-DD` format.

### Step 2: Configure the Classifier
1. **Provide Training Data:** Paste your **GEE Asset ID** for training points.  
   - Must be a `FeatureCollection` with integer labels (**1 = Flood/Water**, **0 = Non-Flood**).  
   - Asset must have public read permissions.
2. **Fetch & Select Label Column:** Click **Fetch Columns**, then choose the column containing the class labels.
3. **Set RF Trees:** Adjust the number of Random Forest trees (default `500`).

### Step 3: Apply Post-Processing Filters
1. **Slope Threshold:** Exclude slopes above a set degree (default: `5°`).
2. **Minimum Patch Size:** Remove small, isolated areas (default: `8` connected pixels).

### Step 4: Run and Get Results
1. **Run Analysis:** Click the red **Run Analysis** button.
2. **Review Results:**  
   - Flooded areas appear in blue; toggle other layers as needed.  
   - Side panel shows `Mapped Flood Area (ha)` and `Overall Accuracy`.  
   - A **GeoTIFF download** link will be available.

---

## Technical Details ⚙️

**Primary Data Sources:**
- **Sentinel-1:** `COPERNICUS/S1_GRD` radar data, filtered with **Refined Lee**.
- **Sentinel-2:** `COPERNICUS/S2_SR` optical data, with cloud masking applied.
- **Topography:** `USGS/SRTMGL1_003` DEM for slope masking.

**Workflow:**
1. Collect S1 & S2 images within the AOI and date range.
2. Stack bands from both sensors into one composite image.
3. Sample the stack with the provided training dataset.
4. Split into **70% training** / **30% validation**.
5. Train `smileRandomForest` classifier and classify the stack.
6. Apply slope and patch size filters.
7. Add flood layer to the map and calculate accuracy metrics.
