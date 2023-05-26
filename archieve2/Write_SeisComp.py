from obspy import read
import os
stream = read('/home/shazam/PycharmProjects/malmi/data/seismic_data/SDS/*.mseed')
# stream.plot()

SDS1 = '{network}.{station}..{channel}'
SDS2 ='{year}.{julday}.mseed'
# Parent_dir='./SDS_converted'


# Set the parent directory to search in
parent_dir = '/home/shazam/PycharmProjects/IPIMLTEST/TEST_server/ori_data'
writefolder= '/home/shazam/PycharmProjects/IPIMLTEST/OnDiskData/GhanaOnDiskSeiscomp'


# if os.path.exists(Parent_dir):
#     pass
# else:
#     os.makedirs(Parent_dir)

import os
from obspy import read

def WriteSeiscomp(stream, writefolder):

    for tr in stream:
        fname = '{}/{}/{}/{}/{}/{}.{}..{}.{}.{}.{}'.format(writefolder, str(stream.traces[0].stats.starttime.year),str(stream.traces[0].stats.network)
                                          ,str(stream.traces[0].stats.station),str(stream.traces[0].stats.channel),
                                          str(stream.traces[0].stats.network),str(stream.traces[0].stats.station),str(stream.traces[0].stats.channel),str(stream.traces[0].stats.starttime.year),str(stream.traces[0].stats.starttime.julday),'mseed')
        path = '{}/{}/{}/{}/{}/'.format(writefolder, str(stream.traces[0].stats.starttime.year),str(stream.traces[0].stats.network)
                                          ,str(stream.traces[0].stats.station),str(stream.traces[0].stats.channel))
        os.makedirs(path, exist_ok=True)
        stream.write(fname, 'MSEED')


# Loop through all subdirectories and files in the parent directory
for subdir, dirs, files in os.walk(parent_dir):
    for file in files:
        # Check if the file has a .mseed extension
        if file.endswith('.mseed'):
            # Construct the full file path
            file_path = os.path.join(subdir, file)
            # print(file_path)
            # Use ObsPy to read the file
            st = read(file_path)

            # Use ObsPy's WriteSeiscomp function to write the data to a file
            WriteSeiscomp(st, writefolder)

    # GH.AKOS..HHE.2012.279.mseed
    # / Year / NET / STA / CHAN.TYPE / NET.STA.LOC.CHAN.TYPE.YEAR.DAY
    # str(tr.stats.starttime.julday)




