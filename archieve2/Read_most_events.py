import pandas as pd

# Read the Excel file
file_path = '/home/shazam/Desktop/IPIML_catalog.xlsx'
column_name = 'time'
df = pd.read_excel(file_path)

# Convert the 'date' column to datetime format
df[column_name] = pd.to_datetime(df[column_name])

# Extract the day from the datetime column
df['day'] = df[column_name].dt.date

# Count the number of events per day
event_counts = df['day'].value_counts()

# Find the top 10 days with the most events
top_10_days = event_counts.head(6)

print("Top 10 days with the most events:")
for day, num_events in top_10_days.items():
    print(day, "- Number of events:", num_events)
