# import csv
# from datetime import datetime, timedelta
# from shapely.geometry import Point, Polygon
#
# # Function to convert string datetime to datetime object for the second catalog file
# def convert_to_datetime_format2(datetime_str):
#     datetime_str = datetime_str.split('.')[0]  # Remove the fractional part
#     return datetime.strptime(datetime_str, '%Y%m%d%H%M%S')
#
# def convert_to_datetime(datetime_str):
#     return datetime.strptime(datetime_str, '%Y/%m/%d %H:%M:%S')
#
# # Function to check if the time difference is less than 2 seconds
# def is_time_difference_less(starttime1, starttime2):
#     time1 = convert_to_datetime_format2(starttime1)
#     time2 = convert_to_datetime(starttime2)
#     time_diff = abs(time1 - time2)
#     return time_diff <= timedelta(seconds=200)
#
# # Path to the first catalog CSV file
# catalog_file1 = "/home/shazam/Desktop/Catalog_output_SEQT/IPIML.csv"
# # Path to the second catalog CSV file
# catalog_file2 = "/home/shazam/Desktop/Catalog_output_SEQT/SEQT.csv"
# # Path to the output CSV file
# output_file = "/home/shazam/Desktop/Catalog_output_SEQT/matching_events.csv"
#
# # Define the polygon points
# polygon_points = [
#     [5.451, -2.135],
#     [5.298, -1.881],
#     [4.96, -2.068],
#     [5.103, -2.318]
# ]
#
# # Create a polygon object from the polygon points
# polygon = Polygon(polygon_points)
#
# # Lists to store matching events
# matching_events = []
#
# # Read and process the first catalog file
# with open(catalog_file1, 'r') as file:
#     reader = csv.DictReader(file)
#     for row in reader:
#         latitude = float(row['latitude'])
#         longitude = float(row['longitude'])
#
#         # Check if the event is within the polygon
#         if polygon.contains(Point(latitude, longitude)):
#             starttime = row['starttime']
#             matching_events.append(starttime)
#
# # Read and process the second catalog file
# matching_points = []
# with open(catalog_file2, 'r') as file:
#     reader = csv.DictReader(file)
#     for row in reader:
#         starttime = row['starttime']
#
#         # Check if the event's starttime is within 2 seconds of any matching event from the first catalog
#         for matching_event in matching_events:
#             if is_time_difference_less(starttime, matching_event):
#                 matching_points.append(row)
#
# # Write the matching event points to the output CSV file
# fieldnames = ['starttime', 'Latitude', 'Longitude', 'depth', 'mag']
# with open(output_file, 'w', newline='') as file:
#     writer = csv.DictWriter(file, fieldnames=fieldnames)
#     writer.writeheader()
#     writer.writerows(matching_points)
#
#
# # Path to the first catalog CSV file
# catalog_file1 = "/home/shazam/Desktop/Catalog_output_SEQT/IPIML.csv"
# # Path to the output CSV file
# output_file = "/home/shazam/Desktop/Catalog_output_SEQT/filtered_events.csv"
#
# # Read and process the first catalog file
# filtered_events = []
#
# with open(catalog_file1, 'r') as file:
#     reader = csv.DictReader(file)
#     for row in reader:
#         latitude = float(row['latitude'])
#         longitude = float(row['longitude'])
#
#         # Check if the event is within the polygon
#         if polygon.contains(Point(latitude, longitude)):
#             filtered_events.append(row)
import csv
from shapely.geometry import Point, Polygon

# Function to convert string datetime to datetime object for the second catalog file
def convert_to_datetime_format2(datetime_str):
    datetime_str = datetime_str.split('.')[0]  # Remove the fractional part
    return datetime.strptime(datetime_str, '%Y%m%d%H%M%S')

def convert_to_datetime(datetime_str):
    return datetime.strptime(datetime_str, '%Y/%m/%d %H:%M:%S')

# Function to check if the time difference is less than 2 seconds
def is_time_difference_less(starttime1, starttime2):
    time1 = convert_to_datetime_format2(starttime1)
    time2 = convert_to_datetime(starttime2)
    time_diff = abs(time1 - time2)
    return time_diff <= timedelta(seconds=200)

# Path to the first catalog CSV file
catalog_file1 = "/home/shazam/Desktop/Catalog_output_SEQT/IPIML.csv"
# Path to the output CSV file
output_file = "/home/shazam/Desktop/Catalog_output_SEQT/filtered.csv"

# Define the polygon points
polygon_points = [
    [5.451, -2.135],
    [5.298, -1.881],
    [4.96, -2.068],
    [5.103, -2.318]
]

# Create a polygon object from the polygon points
polygon = Polygon(polygon_points)

# Read and process the first catalog file
filtered_events = []

with open(catalog_file1, 'r') as file:
    reader = csv.DictReader(file)
    fieldnames = reader.fieldnames  # Read the header
    for row in reader:
        latitude = float(row['latitude'])
        longitude = float(row['longitude'])

        # Check if the event is within the polygon
        if polygon.contains(Point(latitude, longitude)):
            filtered_events.append(row)

# Write the filtered events to the output CSV file
with open(output_file, 'w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(filtered_events)
