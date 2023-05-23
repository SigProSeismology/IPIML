from obspy import  read
import obspy

parent_dir='/home/shazam/PycharmProjects/IPIMLTEST/A/ori_data'

stream = obspy.Stream()

# st=read("")
# st.trim(dt, dt + 5)
#
import os

from obspy.core.utcdatetime import UTCDateTime

# dt = UTCDateTime("2013-04-09T13:07:00.183642Z")


for subdir, dirs, files in os.walk(parent_dir):
    for file in files:
        if file.endswith('.mseed'):
            # Construct the full file path
            file_path = os.path.join(subdir, file)
            # print(file_path)
            # Use ObsPy to read the file
            stream +=read(file_path)
       # read(file)
# stream.detrend('demean')
# stream.detrend('simple')
stream.detrend("spline", order=3, dspline=1000)
# stream.filter('bandpass', freqmin=0.8, freqmax=50)
# stream.trim(dt-20,dt+120)


stream.plot()