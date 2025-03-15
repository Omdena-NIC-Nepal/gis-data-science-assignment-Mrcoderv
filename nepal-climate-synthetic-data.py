import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import os
from scipy.ndimage import gaussian_filter

# Create output directory
output_dir = "nepal_climate_data"
os.makedirs(output_dir, exist_ok=True)

# Nepal's approximate bounding box
nepal_bounds = {
    "min_lon": 80.0,
    "max_lon": 88.2,
    "min_lat": 26.3,
    "max_lat": 30.5
}

# Parameters for raster data
width = 400  # number of pixels in x direction
height = 200  # number of pixels in y direction
x_res = (nepal_bounds["max_lon"] - nepal_bounds["min_lon"]) / width
y_res = (nepal_bounds["max_lat"] - nepal_bounds["min_lat"]) / height

# Create transform for raster data
transform = from_origin(nepal_bounds["min_lon"], nepal_bounds["max_lat"], x_res, y_res)

# Function to create synthetic temperature data with warming trend
def create_temperature_data(baseline_year=2020, projection_year=2050):
    """
    Create synthetic temperature data with a warming trend.
    Returns both historical and projected temperature rasters.
    """
    # Nepal's elevation affects temperature, so create an elevation-like base
    # Higher elevations in the north, lower in the south
    y_coords = np.linspace(0, 1, height)
    elevation_factor = np.repeat(y_coords[:, np.newaxis], width, axis=1)
    
    # Add some random variation for mountains
    mountains = gaussian_filter(np.random.rand(height, width), sigma=5) * 2
    elevation = elevation_factor + mountains
    
    # Baseline temperature decreases with elevation
    baseline_temp = 30 - elevation * 25  # Higher elevation = lower temperature
    
    # Add seasonal pattern (warmer in summer, cooler in winter)
    seasonal_pattern = np.sin(np.linspace(0, 2 * np.pi, 12))
    
    # Create monthly data for baseline year
    baseline_data = []
    for month in range(12):
        # Add seasonal variation
        month_temp = baseline_temp + seasonal_pattern[month] * 8
        
        # Add some random noise
        noise = np.random.normal(0, 1, (height, width)) * 0.5
        month_temp += noise
        
        baseline_data.append(month_temp)
    
    # Create projected data with warming trend
    # Assume 0.04°C warming per year on average (varies spatially)
    years_diff = projection_year - baseline_year
    warming_trend = np.random.normal(0.04 * years_diff, 0.01 * years_diff, (height, width))
    
    # Higher warming in mountains
    warming_trend = warming_trend * (1 + elevation * 0.5)
    
    projected_data = []
    for month in range(12):
        # Base temperature with warming trend
        month_temp = baseline_data[month] + warming_trend
        
        # Add some random noise
        noise = np.random.normal(0, 1, (height, width)) * 0.5
        month_temp += noise
        
        projected_data.append(month_temp)
    
    return baseline_data, projected_data

# Function to create synthetic precipitation data with changing patterns
def create_precipitation_data(baseline_year=2020, projection_year=2050):
    """
    Create synthetic precipitation data with changing patterns.
    Returns both historical and projected precipitation rasters.
    """
    # Create elevation-like base (for orographic effect)
    y_coords = np.linspace(0, 1, height)
    elevation_factor = np.repeat(y_coords[:, np.newaxis], width, axis=1)
    mountains = gaussian_filter(np.random.rand(height, width), sigma=5) * 2
    elevation = elevation_factor + mountains
    
    # Monsoon season affects precipitation (June-September has higher precipitation)
    monthly_factors = np.array([0.2, 0.3, 0.4, 0.6, 0.8, 2.0, 2.5, 2.2, 1.5, 0.7, 0.4, 0.2])
    
    # Baseline precipitation (affected by elevation - orographic effect)
    baseline_precip_annual = 1500 + elevation * 1000  # mm per year
    
    # Create monthly data for baseline year
    baseline_data = []
    for month in range(12):
        # Monthly precipitation based on annual and monthly factor
        month_precip = baseline_precip_annual * monthly_factors[month] / np.sum(monthly_factors)
        
        # Add some random spatial variation
        spatial_var = np.random.normal(1, 0.2, (height, width))
        month_precip = month_precip * spatial_var
        
        baseline_data.append(month_precip)
    
    # Create projected data with changing patterns
    # Climate change effects: more intense monsoon, drier dry season
    years_diff = projection_year - baseline_year
    
    # Monsoon intensification factor (more rain in wet months, less in dry)
    monsoon_intensification = 1 + (monthly_factors > 1) * 0.1 * years_diff / 30
    monsoon_intensification -= (monthly_factors < 0.5) * 0.05 * years_diff / 30
    
    projected_data = []
    for month in range(12):
        # Adjust precipitation based on climate change patterns
        month_precip = baseline_data[month] * monsoon_intensification[month]
        
        # Add some random spatial variation specific to future projections
        spatial_var = np.random.normal(1, 0.25, (height, width))
        month_precip = month_precip * spatial_var
        
        projected_data.append(month_precip)
    
    return baseline_data, projected_data

# Function to save raster data
def save_raster(data, output_path, description):
    """Save data as a GeoTIFF raster file."""
    with rasterio.open(
        output_path,
        'w',
        driver='GTiff',
        height=height,
        width=width,
        count=len(data),
        dtype=data[0].dtype,
        crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
        transform=transform,
    ) as dst:
        for i, layer in enumerate(data):
            dst.write(layer, i + 1)
        dst.descriptions = [f"Month {i+1}" for i in range(len(data))]
        dst.update_tags(description=description)
    
    print(f"Saved: {output_path}")

# Function to create vector data: administrative regions (simplified)
def create_admin_regions():
    """Create simplified vector dataset of Nepal's administrative regions."""
    # Simplified provinces of Nepal
    provinces = [
        {"name": "Province 1", "geometry": Polygon([
            (87.0, 26.5), (88.0, 27.0), (87.8, 28.0), (87.0, 28.5), (86.5, 27.5), (87.0, 26.5)
        ])},
        {"name": "Province 2", "geometry": Polygon([
            (85.0, 26.5), (87.0, 26.5), (86.5, 27.5), (85.0, 27.2), (85.0, 26.5)
        ])},
        {"name": "Bagmati", "geometry": Polygon([
            (85.0, 27.2), (86.5, 27.5), (86.2, 28.5), (85.0, 28.2), (85.0, 27.2)
        ])},
        {"name": "Gandaki", "geometry": Polygon([
            (83.5, 27.5), (85.0, 27.2), (85.0, 28.2), (84.0, 29.0), (83.5, 28.0), (83.5, 27.5)
        ])},
        {"name": "Lumbini", "geometry": Polygon([
            (82.0, 27.0), (83.5, 27.5), (83.5, 28.0), (82.5, 28.2), (82.0, 27.5), (82.0, 27.0)
        ])},
        {"name": "Karnali", "geometry": Polygon([
            (81.0, 28.0), (82.5, 28.2), (83.5, 28.0), (82.5, 29.5), (81.5, 29.0), (81.0, 28.0)
        ])},
        {"name": "Sudurpashchim", "geometry": Polygon([
            (80.2, 28.5), (81.0, 28.0), (81.5, 29.0), (80.5, 30.0), (80.2, 28.5)
        ])}
    ]
    
    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(provinces, crs="EPSG:4326")
    
    # Add random climate vulnerability index (for demonstration)
    np.random.seed(42)
    gdf["climate_vulnerability_2020"] = np.random.uniform(0.3, 0.7, len(gdf))
    
    # Project higher vulnerability for 2050
    gdf["climate_vulnerability_2050"] = gdf["climate_vulnerability_2020"] * np.random.uniform(1.2, 1.5, len(gdf))
    
    # Ensure values are between 0 and 1
    gdf["climate_vulnerability_2050"] = np.minimum(gdf["climate_vulnerability_2050"], 1.0)
    
    return gdf

# Function to create vector data: river network (simplified)
def create_river_network():
    """Create simplified vector dataset of Nepal's major rivers."""
    # Simplified major rivers
    rivers = [
        {"name": "Koshi", "flow_2020": 2200, "flow_2050": 2000, "geometry": 
         Polygon([
             (86.0, 26.5), (86.1, 26.5), (87.0, 28.0), (86.9, 28.0), (86.0, 26.5)
         ])},
        {"name": "Gandaki", "flow_2020": 1800, "flow_2050": 1650, "geometry": 
         Polygon([
             (84.0, 26.5), (84.1, 26.5), (85.0, 28.5), (84.9, 28.5), (84.0, 26.5)
         ])},
        {"name": "Karnali", "flow_2020": 1600, "flow_2050": 1400, "geometry": 
         Polygon([
             (81.0, 26.5), (81.1, 26.5), (82.0, 29.0), (81.9, 29.0), (81.0, 26.5)
         ])},
        {"name": "Mahakali", "flow_2020": 1200, "flow_2050": 1050, "geometry": 
         Polygon([
             (80.0, 26.5), (80.1, 26.5), (80.5, 29.0), (80.4, 29.0), (80.0, 26.5)
         ])}
    ]
    
    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(rivers, crs="EPSG:4326")
    
    # Add flow reduction percentage
    gdf["flow_reduction_pct"] = ((gdf["flow_2020"] - gdf["flow_2050"]) / gdf["flow_2020"]) * 100
    
    return gdf

# Function to create vector data: glacier retreat points (simplified)
def create_glacier_data():
    """Create simplified vector dataset of Nepal's glacier retreat points."""
    # Simulated glacier monitoring points in Nepal's high mountains
    np.random.seed(42)
    n_points = 20
    
    # Create points in the northern region (higher altitude)
    lons = np.random.uniform(nepal_bounds["min_lon"], nepal_bounds["max_lon"], n_points)
    lats = np.random.uniform(nepal_bounds["max_lat"] - 1.5, nepal_bounds["max_lat"], n_points)
    
    # Create data for the points
    glacier_data = []
    for i in range(n_points):
        # Glacier retreat in meters per year (varying between 10-25m/year for 2020)
        retreat_2020 = np.random.uniform(10, 25)
        
        # Projected retreat for 2050 (increase by 30-80%)
        increase_factor = np.random.uniform(1.3, 1.8)
        retreat_2050 = retreat_2020 * increase_factor
        
        glacier_data.append({
            "id": f"GL{i+1:02d}",
            "geometry": Point(lons[i], lats[i]),
            "retreat_2020": retreat_2020,
            "retreat_2050": retreat_2050,
            "increase_pct": (increase_factor - 1) * 100
        })
    
    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(glacier_data, crs="EPSG:4326")
    
    return gdf

# --- Generate all data ---

# 1. Generate raster data
print("Generating raster data...")
baseline_temp, projected_temp = create_temperature_data(2020, 2050)
baseline_precip, projected_precip = create_precipitation_data(2020, 2050)

# Save temperature data
save_raster(
    baseline_temp, 
    os.path.join(output_dir, "nepal_temperature_2020.tif"),
    "Nepal monthly temperature (°C) for 2020"
)
save_raster(
    projected_temp, 
    os.path.join(output_dir, "nepal_temperature_2050.tif"),
    "Nepal projected monthly temperature (°C) for 2050"
)

# Save precipitation data
save_raster(
    baseline_precip, 
    os.path.join(output_dir, "nepal_precipitation_2020.tif"),
    "Nepal monthly precipitation (mm) for 2020"
)
save_raster(
    projected_precip, 
    os.path.join(output_dir, "nepal_precipitation_2050.tif"),
    "Nepal projected monthly precipitation (mm) for 2050"
)

# 2. Generate vector data
print("Generating vector data...")

# Administrative regions
admin_regions = create_admin_regions()
admin_regions.to_file(os.path.join(output_dir, "nepal_admin_regions.gpkg"), driver="GPKG")
print(f"Saved: {os.path.join(output_dir, 'nepal_admin_regions.gpkg')}")

# River network
rivers = create_river_network()
rivers.to_file(os.path.join(output_dir, "nepal_rivers.gpkg"), driver="GPKG")
print(f"Saved: {os.path.join(output_dir, 'nepal_rivers.gpkg')}")

# Glacier retreat points
glaciers = create_glacier_data()
glaciers.to_file(os.path.join(output_dir, "nepal_glaciers.gpkg"), driver="GPKG")
print(f"Saved: {os.path.join(output_dir, 'nepal_glaciers.gpkg')}")

# 3. Create metadata
metadata = {
    "title": "Nepal Climate Change Synthetic Dataset",
    "description": "Synthetic climate data for Nepal showing baseline (2020) and projected (2050) conditions",
    "created_date": datetime.now().strftime("%Y-%m-%d"),
    "spatial_coverage": f"Nepal ({nepal_bounds['min_lon']}, {nepal_bounds['min_lat']}) to ({nepal_bounds['max_lon']}, {nepal_bounds['max_lat']})",
    "temporal_coverage": "2020 (baseline) and 2050 (projection)",
    "raster_resolution": f"{x_res:.4f} degrees (~{x_res * 111}km at equator)",
    "projection": "EPSG:4326 (WGS84)",
    "scenario": "Synthetic data approximating RCP4.5-like scenario",
    "variables": {
        "temperature": "Monthly average temperature in degrees Celsius",
        "precipitation": "Monthly total precipitation in millimeters",
        "admin_regions": "Administrative regions with climate vulnerability index",
        "rivers": "Major river systems with projected flow changes",
        "glaciers": "Glacier monitoring points with retreat rates"
    },
    "notes": "This is synthetic data for demonstration purposes. It should not be used for actual climate research or planning."
}

# Save metadata as JSON
import json
with open(os.path.join(output_dir, "metadata.json"), "w") as f:
    json.dump(metadata, f, indent=2)
print(f"Saved: {os.path.join(output_dir, 'metadata.json')}")

print("\nData generation complete. Files saved in the 'nepal_climate_data' directory.")