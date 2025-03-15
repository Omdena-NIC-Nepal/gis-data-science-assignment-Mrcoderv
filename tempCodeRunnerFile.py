import geopandas as gpd
import rasterio
import numpy as np
import matplotlib.pyplot as plt
from rasterstats import zonal_stats

# Load administrative regions with corrected path
admin_regions = gpd.read_file("nepal_climate_data/nepal_admin_regions.gpkg")

# Function to process raster data
def process_raster(raster_path):
    with rasterio.open(raster_path) as src:
        data = src.read()
        annual_avg = np.mean(data, axis=0)
        return annual_avg, src.transform

# Process temperature data
temp_2020, transform = process_raster("nepal_climate_data/nepal_temperature_2020.tif")
temp_2050, _ = process_raster("nepal_climate_data/nepal_temperature_2050.tif")
temp_diff = temp_2050 - temp_2020

# Zonal statistics function
def calculate_zonal_stats(geodf, raster_array, affine_transform):
    stats = zonal_stats(
        geodf,
        raster_array,
        affine=affine_transform,
        stats=['mean', 'median', 'min', 'max'],
        geojson_out=True
    )
    return gpd.GeoDataFrame.from_features(stats)

# Calculate statistics
temp_stats_2020 = calculate_zonal_stats(admin_regions, temp_2020, transform)
temp_stats_2050 = calculate_zonal_stats(admin_regions, temp_2050, transform)

# Merge results
admin_regions['temp_2020'] = temp_stats_2020['mean']
admin_regions['temp_2050'] = temp_stats_2050['mean']
admin_regions['temp_change'] = admin_regions['temp_2050'] - admin_regions['temp_2020']

# Visualization 1: Temperature Change Map
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
admin_regions.plot(column='temp_change', legend=True, ax=ax,
                   legend_kwds={'label': "Temperature Change (°C)"})
plt.title("Projected Temperature Change 2020-2050 by Administrative Region")
plt.savefig("temperature_change_map.png")
plt.close()

# Visualization 2: Climate Vulnerability vs Temperature Change
fig, ax = plt.subplots(figsize=(10, 6))
admin_regions.plot.scatter(x='climate_vulnerability_2050', y='temp_change', ax=ax)
plt.title("Climate Vulnerability vs Temperature Change")
plt.xlabel("Climate Vulnerability Index (2050)")
plt.ylabel("Temperature Change (°C)")
plt.savefig("vulnerability_vs_temp.png")
plt.close()

# River Flow Analysis
rivers = gpd.read_file("nepal_climate_data/nepal_rivers.gpkg")
plt.figure(figsize=(10, 6))
rivers['flow_reduction_pct'].plot(kind='bar')
plt.xticks(range(len(rivers)), rivers['name'])
plt.title("Projected River Flow Reductions by 2050")
plt.ylabel("Flow Reduction (%)")
plt.savefig("river_flow_reductions.png")
plt.close()

# EDA Statistics
print("\nTemperature Statistics (2020-2050):")
print(f"Average Temperature Change: {admin_regions['temp_change'].mean():.2f}°C")
print(f"Max Regional Change: {admin_regions['temp_change'].max():.2f}°C")
print(f"Min Regional Change: {admin_regions['temp_change'].min():.2f}°C")

print("\nClimate Vulnerability Statistics (2050):")
print(admin_regions['climate_vulnerability_2050'].describe())

# Glacier Analysis
glaciers = gpd.read_file("nepal_climate_data/nepal_glaciers.gpkg")
print("\nGlacier Retreat Statistics:")
print(f"Average Retreat Rate 2020: {glaciers['retreat_2020'].mean():.1f} m/yr")
print(f"Projected 2050 Rate: {glaciers['retreat_2050'].mean():.1f} m/yr")