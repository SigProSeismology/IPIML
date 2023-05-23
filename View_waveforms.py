import os
from obspy import read, Stream
import matplotlib.pyplot as plt

# set the directory path to search for mseed files
dir_path = '/home/shazam/PycharmProjects/IPIMLTEST/AAAAA/MIL_OUTPUTS/data_ML/primary_events/mseeds/2013-03-14-09-39-07'

# initialize an empty Stream object to merge the data
st = Stream()

# search for all .mseed files in the directory and its subdirectories
for root, dirs, files in os.walk(dir_path):
    for file in files:
        if file.endswith('.mseed'):
            file_path = os.path.join(root, file)
            # read the mseed file and add it to the Stream object
            st += read(file_path)

# merge the data
st.merge()
st.detrend('demean')
st.detrend('simple')
# st.filter('bandpass', freqmin=0, freqmax=50, corners=2,zerophase=True)
st.taper(max_percentage=0.001, type='cosine', max_length=1)  # to avoid anormaly at bounday

# plot the data
st.plot()