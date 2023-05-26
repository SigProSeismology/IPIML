import csv
from datetime import datetime, timedelta

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
# Path to the second catalog CSV file
catalog_file2 = "/home/shazam/Desktop/Catalog_output_SEQT/SEQT.csv"
# Path to the output CSV file
output_file = "/home/shazam/Desktop/Catalog_output_SEQT/matching_events.csv"


# Lists to store matching events
matching_events = []

# Read and process the first catalog file
with open(catalog_file1, 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        latitude = float(row['latitude'])
        longitude = float(row['longitude'])

        # Check if the event is within the specified latitude and longitude range
        if -2.3 <= longitude <= -1.9 and 5.0 <= latitude <= 5.5:
            starttime = row['starttime']
            matching_events.append(starttime)


# Read and process the second catalog file
matching_points = []
with open(catalog_file2, 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        starttime = row['starttime']

        # Check if the event's starttime is within 2 seconds of any matching event from the first catalog
        for matching_event in matching_events:
            if is_time_difference_less(starttime, matching_event):
                matching_points.append(row)

# Write the matching event points to the output CSV file
fieldnames = ['starttime', 'Latitude', 'Longitude', 'depth', 'mag']
with open(output_file, 'w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(matching_points)




import csv

# Path to the first catalog CSV file
catalog_file1 = "/home/shazam/Desktop/Catalog_output_SEQT/IPIML.csv"
# Path to the output CSV file
output_file = "/home/shazam/Desktop/Catalog_output_SEQT/filtered_events.csv"

# Read and process the first catalog file
filtered_events = []

with open(catalog_file1, 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        latitude = float(row['latitude'])
        longitude = float(row['longitude'])

        # Check if the event is within the specified latitude and longitude range
        if -2.3 <= longitude <= -1.8 and 5.0 <= latitude <= 6.0:
            filtered_events.append(row)

# Write the filtered events to the output CSV file
fieldnames = list(filtered_events[0].keys())

with open(output_file, 'w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(filtered_events)










#
# import pandas as pd
#
# # Read the CSV file
# csv_path = '/home/shazam/Desktop/Catalog_output_SEQT/SEQT.csv'
# df = pd.read_csv(csv_path)
#
# # Define the polygon coordinates
# polygon_coordinates = [(-1.046, 5.491), (-1.046, 6.072), (-1.604, 6.024), (-1.632, 5.489)]
#
# # Apply depth and location filters
# depth_filtered_events = df[
#     (df['Longitude'] >= polygon_coordinates[0][0]) &
#     (df['Longitude'] <= polygon_coordinates[2][0]) &
#     (df['Latitude'] >= polygon_coordinates[0][1]) &
#     (df['Latitude'] <= polygon_coordinates[1][1]) &
#     (df['depth'] >= 0) &
#     (df['depth'] <= 50)
# ]
#
# # Write filtered events to a new CSV file
# output_path = '/home/shazam/Desktop/Catalog_output_SEQT/SEQT_extra.csv'
# depth_filtered_events.to_csv(output_path, index=False)
#
# print("Filtered events have been written to:", output_path)

#
# import csv
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
#         # Check if the event is within the specified latitude and longitude range
#         if -2.3 <= longitude <= -1.9 and 5.0 <= latitude <= 5.5:
#             filtered_events.append(row)
#
# # Write the filtered events to the output CSV file
# fieldnames = list(filtered_events[0].keys())
#
# with open(output_file, 'w', newline='') as file:
#     writer = csv.DictWriter(file, fieldnames=fieldnames)
#     writer.writeheader()
#     writer.writerows(filtered_events)
