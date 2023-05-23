import os
from obspy.core import read
from datetime import datetime


class Waveforms:

    def __init__(self, event_path, extension='*', comps=['E','N','Z'], freq=None, station_idmode='network.station'):
        if not os.path.isdir(event_path):
            raise ValueError('Error: data path does not exist')
        try:
            self.load_waveforms(event_path, extension, comps, freq, station_idmode)
        except:
            raise WaveformLoadingError('Error: data not read for the event: %s' %(event_path))
        self.station_list()

    def station_list(self):
        data_stalist=[]
        for comp in (self.stream).keys():
            for sta in (self.stream[comp]).keys():
                if sta not in data_stalist:
                    data_stalist.append(sta)
        self.data_stations=set(data_stalist)

    def load_waveforms(self, event_path, extension, comps, freq, station_idmode):
        files=os.path.join(event_path,extension)
        traces=read(files)
        
        if freq:
            traces.detrend('demean')
            traces.detrend('linear')
            if len(freq) == 1:
                traces.filter("highpass", freq=freq[0])
            elif len(freq) == 2:
                traces.filter("bandpass", freqmin=freq[0], freqmax=freq[1])
        
        self.stream={}
        for comp in comps:
            self.stream[comp]={}
            for tr in traces:
                if tr.stats.channel[-1]==comp:
                    dtime=datetime.strptime(str(tr.stats.starttime),"%Y-%m-%dT%H:%M:%S.%fZ")
                    if station_idmode == 'station':
                        station_id = "{}".format(tr.stats.station)
                    elif station_idmode == 'network.station':
                        station_id = "{}.{}".format(tr.stats.network, tr.stats.station)
                    elif station_idmode == 'network.station.location':
                        station_id = "{}.{}.{}".format(tr.stats.network, tr.stats.station, tr.stats.location)
                    elif station_idmode == 'network.station.location.instrument':
                        station_id = "{}.{}.{}.{}".format(tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel[:-1])
                    else:
                        raise ValueError
                    self.stream[comp][station_id]=[dtime, tr.stats.delta, tr.data]


class WaveformLoadingError(Exception):
    pass
