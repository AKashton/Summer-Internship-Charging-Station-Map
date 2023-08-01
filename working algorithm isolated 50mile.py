import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from decimal import Decimal, ROUND_HALF_UP
# Load the GeoJSON file with charging station locations
#EV_charge_gdf

Future_Charge_gdf = gpd.read_file('/content/drive/MyDrive/Data for web mapping project/Under_Construction_EV_Charge.geojson')

#merge future stations with current stations
charging_stations_gdf = gpd.GeoDataFrame(pd.concat([EV_charge_gdf, Future_Charge_gdf], ignore_index=True))

# List to store marker points
marker_points = []

for road_coordinates in rounded_coordinates:
    total_distance = 0.0
    prev_lon, prev_lat, _ = road_coordinates[0]
    segment_distance = 0.0
    segment_interval = 1.5 # Check for charging stations every 1.5 miles
    charging_station_found = False

    for index in range(0, len(road_coordinates)):
        lon, lat, _ = road_coordinates[index]
        distance = haversine(prev_lon, prev_lat, lon, lat)

        segment_distance += distance
        total_distance += distance
        prev_lon, prev_lat = lon, lat

        # Check for charging stations every segment_interval (e.g., 1.5 miles)
        if segment_distance >= segment_interval:
            # Calculate the bounding box for the current segment with a radius of 1.5 miles
            lon_min, lat_min, lon_max, lat_max = create_bounding_box(lon, lat)

            # Filter charging stations within the bounding box
            charging_stations_within_bbox = charging_stations_gdf.cx[lon_min:lon_max, lat_min:lat_max]

            # Check if any charging stations are within the bounding box
            if not charging_stations_within_bbox.empty:
                 charging_station_found = True  # Set the flag if charging stations are found

            # Reset the segment distance for the next segment
            segment_distance = 0.0

            if charging_station_found:
                # Reset total_distance if charging station is found within 1.5 miles and go to the next segment
                total_distance = 0.0
                charging_station_found = False

        if total_distance >= 50:  # If accumulated distance is 50 miles or more
            # Create a marker point at the current coordinates
            marker_points.append(Point(lon, lat))

            # Adjust total_distance to avoid skipping points
            total_distance = total_distance % 50

# Create a GeoDataFrame from the marker points
gdf_markers = gpd.GeoDataFrame(geometry=marker_points)