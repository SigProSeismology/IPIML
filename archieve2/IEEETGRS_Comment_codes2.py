import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, Point

# Read the points data from CSV file
points = pd.read_csv('/home/shazam/Desktop/catalog/GOG_NEW.csv')

# Define the polygon
polygon1 = Polygon([(-1.208, 6.013), (-1.015, 6.569), (-0.925, 6.536), (-1.112, 5.977)])

gdf_points = gpd.GeoDataFrame(points, geometry=gpd.points_from_xy(points.longitude, points.latitude))
points_in_polygon = gdf_points[gdf_points.geometry.within(polygon1)]
points_in_polygon.to_csv('/home/shazam/Desktop/catalog/cluster1.csv', index=False)


polygon2 = Polygon([(-0.0773,6.0556), (-0.1325,6.233), (-0.0713,6.25), (-0.0162,6.0781)])
gdf_points = gpd.GeoDataFrame(points, geometry=gpd.points_from_xy(points.longitude, points.latitude))
points_in_polygon = gdf_points[gdf_points.geometry.within(polygon2)]
points_in_polygon.to_csv('/home/shazam/Desktop/catalog/cluster2.csv', index=False)





polygon3 = Polygon([(-2.11, 5.38), (-1.97, 5.27), (-2.123, 5.057), (-2.267, 5.16)])
gdf_points = gpd.GeoDataFrame(points, geometry=gpd.points_from_xy(points.longitude, points.latitude))
points_in_polygon = gdf_points[gdf_points.geometry.within(polygon3)]
points_in_polygon.to_csv('/home/shazam/Desktop/catalog/cluster3.csv', index=False)










polygon4 = Polygon([(-0.587, 5.82), (-0.209, 5.877), (-0.19, 5.767), (-0.57, 5.701)])
gdf_points = gpd.GeoDataFrame(points, geometry=gpd.points_from_xy(points.longitude, points.latitude))
points_in_polygon = gdf_points[gdf_points.geometry.within(polygon4)]
points_in_polygon.to_csv('/home/shazam/Desktop/catalog/cluster4.csv', index=False)




