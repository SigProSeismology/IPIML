import csv
from datetime import datetime

with open('/home/shazam/Desktop/catalog/GOG.csv', 'r') as f_input, open('/home/shazam/Desktop/catalog/GOG_NEW.csv', 'w', newline='') as f_output:
    csv_input = csv.DictReader(f_input)
    fieldnames = ['longitude', 'latitude', 'year', 'month', 'day', 'mag', 'depth_km', 'hour', 'minute']
    csv_output = csv.DictWriter(f_output, fieldnames=fieldnames)
    csv_output.writeheader()

    for row in csv_input:
        if row['time'] != '':
            time_str = row['time']
            time_obj = datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%S.%fZ')
            csv_output.writerow({
                'longitude': row['longitude'],
                'latitude': row['latitude'],
                'year': time_obj.year,
                'month': time_obj.month,
                'day': time_obj.day,
                'mag': row['mag'],
                'depth_km': row['depth_km'],
                'hour': time_obj.hour,
                'minute': time_obj.minute
            })












