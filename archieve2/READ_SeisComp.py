from obspy.clients.filesystem.sds import Client
client = Client("/my/SDS/archive/root")


from obspy import read
stream = read()

SDS = 'test_{channel}.mseed'
for tr in stream:
    fname = SDS.format(**tr.stats)
    tr.write(fname, 'MSEED')
