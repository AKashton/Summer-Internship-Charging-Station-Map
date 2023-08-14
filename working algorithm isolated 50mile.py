!pip install geopandas pyogrio fiona folium matplotlib mapclassify

import geopandas as gpd
import pandas as pd
import math
from decimal import Decimal
import numpy as np
from shapely.geometry import mapping

roads_gdf = gpd.read_file('/content/drive/MyDrive/Data for web mapping project/AlaskaHW.geojson')
EV_charge_gdf = gpd.read_file('/content/drive/MyDrive/Data for web mapping project/AK_charge_fast.geojson')


# List to store all the individual coordinates
all_coordinates = []

# Iterate over all the MultiLineString geometries in the GeoDataFrame
for geom in roads_gdf['geometry']:
    # Convert the MultiLineString geometry to a dictionary representation
    geom_dict = mapping(geom)

    # Extract individual coordinates from the dictionary
    coordinates = geom_dict['coordinates']

    # Append the coordinates to the 'all_coordinates' list
    all_coordinates.extend(coordinates)



def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the distance between two points on the Earth's surface using the haversine formula.

    Parameters:
        lat1, lon1: Latitude and longitude of the first point in degrees.
        lat2, lon2: Latitude and longitude of the second point in degrees.
        radius: Radius of the Earth in the desired units (default is kilometers).

    Returns:
        The distance between the two points in the same units as the radius. converted to miles for US (Alaska calculations)
    """
    R =6371.0 #Earth's radius in km

    # Convert latitude and longitude from degrees to radians
    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2)
    lon2_rad = np.radians(lon2)

    # Differences in coordinates
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad

    # Haversine formula
    a = np.sin(dlat / 2)**2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    # Calculate the distance in kilometers
    distance = R * c
    #convert to miles
    distance_mile = distance * 0.62137

    return distance_mile

"""
the create_bounding_box function creates a box around the road coordinates to identify other charging stations
within 1.5 miles of the road system. (the NEVI plan specifies that charging stations must be within 1 mile of fuel corridors or major road systems.)

I used 1.5 miles instead of 1 mile to account for accuray errors and rounding errors.

Paramateres:
  lon, lat, alt, radius: coordinates of the current section of road that the algorithim is at, units are in degrees

Returns:
  a bounding box with an radius of 1.5 miles from the coordinates of the road.
"""



def create_bounding_box(lon, lat, _=0, radius=1.5):

    # Conversion factor from miles to degrees (approximate at the equator)
    MILES_TO_DEGREES = Decimal(radius / 69.17)
    # Ensure lon and lat are Decimal objects for consistency
    lon = Decimal(str(lon))
    lat = Decimal(str(lat))

    # Calculate the distance covered by one degree of latitude at the given latitude
    distance_per_degree_latitude = MILES_TO_DEGREES * Decimal(math.cos(math.radians(lat)))

    # Calculate the bounding box coordinates around the point (lon, lat) with a 1.5-mile radius
    lon_min = lon - MILES_TO_DEGREES
    lat_min = lat - distance_per_degree_latitude
    lon_max = lon + MILES_TO_DEGREES
    lat_max = lat + distance_per_degree_latitude

    return lon_min, lat_min, lon_max, lat_max

# Function to round a single coordinate list
def round_coordinate(coord, decimal_places = 6):
    return [round(coord[0], decimal_places), round(coord[1], decimal_places), round(coord[2], decimal_places)]

# Function to round the coordinates in a nested list
def round_coordinates(coord_list):
    return [round_coordinate(coord) for coord in coord_list]

# Apply the rounding function to all the coordinates in the all_coordinates array
rounded_coordinates = [round_coordinates(road) for road in all_coordinates]


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
