from obspy import read_events
import numpy


a= numpy.load('/home/shazam/PycharmProjects/IPIMLTEST/Final/SEqTCatalogs/seqt_real_e_dict.npy')

# Path to the EVT file
file_path = "/path/to/file.EVT"

# Read events from file
catalog = read_events(file_path, format="EVT")

# Loop through events in catalog
for event in catalog.events:
    print("Event:", event)

    # Loop through picks for the event
    for pick in event.picks:
        print("Pick:", pick)
