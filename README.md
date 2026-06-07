# GeoSPEI-Framework
# ======================
Monsoon SPEI Computation and Visualization Framework
Developed by: Swapnajit Chowdhury
Description:
This script computes seasonal SPEI (Standardized Precipitation-Evapotranspiration Index) using CHIRPS precipitation and ERA5-Land evapotranspiration datasets in Google Earth Engine and generates multi-temporal visualization outputs.
Users may adapt this workflow to their own study areas
and analysis periods.
# 1. RUN THIS FIRST IF IN COLAB OR JUPYTER: 
# !pip install cartopy rasterio geemap 
import ee 
import geemap 
import os 
import glob 
import rasterio 
import numpy as np 
import matplotlib.pyplot as plt 
!pip install cartopy 
import cartopy 
import cartopy.crs as ccrs 
import cartopy.feature as cfeature 
from matplotlib import colors 
# 2. Initialize Earth Engine 
try: 
ee.Initialize() 
except Exception: 
ee.Authenticate() 
ee.Initialize() 
# --- Parameters & Study Area --- 
wb = ee.FeatureCollection("FAO/GAUL/2015/level1") \ 
.filter(ee.Filter.eq('ADM1_NAME', 'West Bengal')) \ 
.geometry() 
start_year, end_year = 2015, 2025 
baseline_start, baseline_end = 2010, 2019 
start_month, end_month = 6, 9 
scale = 10000 
# --- Extraction Functions --- 
def get_precip(year): 
return ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY").filterBounds(wb) \ 
.filter(ee.Filter.calendarRange(year, year, 'year')) \ 
.filter(ee.Filter.calendarRange(start_month, end_month, 'month')) \ 
.sum().rename('P') 
def get_et(year): 
return ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR").filterBounds(wb) \ 
.filter(ee.Filter.calendarRange(year, year, 'year')) \ 
.filter(ee.Filter.calendarRange(start_month, end_month, 'month')) \ 
.select('total_evaporation_sum') \ 
.sum().multiply(-1000).rename('ET') 
def balance(year): 
return get_precip(year).subtract(get_et(year)).rename('D') 
# --- SPEI Logic --- 
baseline_list = ee.List.sequence(baseline_start, baseline_end).map(lambda y: balance(y)) 
baseline_col = ee.ImageCollection.fromImages(baseline_list) 
mean = baseline_col.mean() 
std = baseline_col.reduce(ee.Reducer.stdDev()).where(ee.Image(0), 1) 
# --- Process & Export --- 
output_folder = "SPEI_WB" 
os.makedirs(output_folder, exist_ok=True) 
for y in range(start_year, end_year + 1): 
spei_img = balance(y).subtract(mean).divide(std).rename('SPEI') 
filename = f"{output_folder}/SPEI_{y}.tif" 
if not os.path.exists(filename): 
print(f"Downloading data for {y}...") 
geemap.ee_export_image(spei_img, filename=filename, scale=scale, region=wb) 
# --- Visualization --- 
files = sorted(glob.glob(f"{output_folder}/SPEI_*.tif")) 
ncols = 4 
nrows = int(np.ceil(len(files) / ncols)) 
fig = plt.figure(figsize=(18, 5 * nrows)) 
cmap = plt.get_cmap('RdYlBu') 
norm = colors.TwoSlopeNorm(vcenter=0, vmin=-2.5, vmax=2.5) 
for i, file in enumerate(files): 
ax = plt.subplot(nrows, ncols, i + 1, projection=ccrs.PlateCarree()) 
with rasterio.open(file) as src: 
data = src.read(1) 
data = np.ma.masked_where(data == src.nodata, data) 
extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top] 
im = ax.imshow(data, extent=extent, cmap=cmap, norm=norm, 
transform=ccrs.PlateCarree(), origin='upper') 
ax.add_feature(cfeature.BORDERS, linewidth=0.8, edgecolor='black') 
ax.add_feature(cfeature.COASTLINE, linewidth=0.8) 
# Accurate state boundaries for West Bengal context 
states = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', 
scale='10m', facecolor='none') 
ax.add_feature(states, edgecolor='gray', linewidth=0.5) 
gl = ax.gridlines(draw_labels=True, linestyle=':', alpha=0.6) 
gl.top_labels = gl.right_labels = False 
gl.xlabel_style = {'size': 8} 
gl.ylabel_style = {'size': 8} 
ax.set_title(f"Monsoon SPEI {2015 + i}", fontsize=12, fontweight='bold', pad=10) 
# Adjust spacing 
plt.subplots_adjust(wspace=0.15, hspace=0.3, top=0.92, bottom=0.15) 
# Colorbar 
cbar_ax = fig.add_axes([0.3, 0.08, 0.4, 0.015]) 
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal') 
cbar.set_label("Standardized Precipitation-Evapotranspiration Index (SPEI)", fontsize=12) 
# Global Title 
plt.suptitle("SPEI Analysis", fontsize=22, fontweight='bold') 
# Attribution Text 
attribution_text = "Data: UCSB CHIRPS & ECMWF ERA5-Land\nVisualization: Swapnajit Chowdhury" 
fig.text(0.98, 0.02, attribution_text, fontsize=10, color='black', 
ha='right', va='bottom', style='italic', fontweight='medium') 
# Save and Show 
plt.savefig("SPEI.png", dpi=300, bbox_inches='tight') 
print("SPEI_.png") 
plt.show()

