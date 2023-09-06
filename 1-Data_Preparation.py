#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 2023

@authors: Xiao Zhuowei, Hamzeh Mohammadigheymasi, Peidong Shi
"""

"""
Download data from IRIS using obspy and EqT libraries
Convert the OnDisk data format from the SeisComp format to the EQT format
"""
import sys
from obspy.io.stationxml.core import _read_stationxml, _write_stationxml
import xml.etree.ElementTree as ET
import glob
from pandas import to_datetime
import fnmatch
import obspy
sys.path.append('./src/S_EqT_codes/src/EqT_libs')
sys.path.append('./srcmalmi')
from downloader2 import makeStationList, downloadMseeds
import argparse
from datetime import datetime
import os
import json
import pandas as pd
from pathlib import Path
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.inventory import Inventory, Network, Station, Channel




def stream2EQTinput_raw(stream, dir_output, instrument_code=None, component_code=None, freqband=None, station_code=None):
    """
    This function is used to format the input obspy stream into the EQ-Transformer
    acceptable seismic data inputs.

    The three component seismic data of a station should be downloaded at the same time range.
    The output filename contains the time range of the data. For a perticular station,
    the starttime and endtime with a wider range of the three component is used as the
    unified time range in the output filename.
    So don't split different component data of the same station to differnt stream,
    they must be kept in the same stream and be checked for outputting. You can simply
    merge different streams to a final complete stream which contrain all stations or at
    least all components of the same station, and then pass the steam to this function.

    In general the seismic data (stream) span a day (data are usually downloaded daily).
    However, this function also accept data time range longer or smaller than a day.
    But using daily data segment is highly recommended, because by default the EQ-Transformer
    are set to process this kind of data (daily date segment). It now also works for longer or
    short time range. But there is no guarantee that the future updated version will also
    support this feature.

    Parameters
    ----------
    stream : obspy stream
        input seismic data.
    dir_output : str
        directory for outputting.
    instrument_code : list of str
        instrument_code for outputting, such as: ["HH", "BH", "EH", "SH", "HG", "HN"];
        only output data which have the listed instrument codes;
        if None or [], then searching for all avaliable instrument code in the input stream.
    component_code : list of str
        component_code for outputting, such as: ['Z','N','E','1','2','3'];
        only output data of the listed components;
        Note complete data should have at least three component data, such as ['Z','N','E'];
        if None or [], then searching for all avaliable component code in the input stream.
    freqband : list of float
        frequency range in Hz for filtering seismic data,
        e.g. [3, 45] meaning filter seismic data to 3-45 Hz.
        default is None, means no filtering.
    station_code : list of str, default is None
        specify the stations for output.
        If None or [], will output all avaliable stations.

    Returns
    -------
    None.

    Example
    -------
    dir_output = '/Users/human/eqt/examples/mseeds'
    stream2EQTinput(stream, dir_output)
    """

    timeformat = "%Y%m%dT%H%M%SZ"  # NOTE here output until second
    if not instrument_code:  # for None or [] or Flase will reset
        # no input instrument codes
        # search for all available instrument codes in the input stream data
        instrument_code = []
        for tr in stream:
            if tr.stats.channel[:-1] not in instrument_code:
                instrument_code.append(tr.stats.channel[:-1])
        del tr
    if not component_code:
        # no input component codes
        # search for all available component codes in the input stream data
        component_code = []
        for tr in stream:
            if tr.stats.channel[-1] not in component_code:
                component_code.append(tr.stats.channel[-1])
        del tr
    if not station_code:
        # no input station codes
        # scan all traces to get the station names
        station_code = []
        for tr in stream:
            sname = tr.stats.station
            if sname not in station_code:
                station_code.append(sname)
        del tr
    # for a particular station, first check starttime and endtime, then output data
    for ista in station_code:
        # select and output data for a perticular station
        stdata_ista = stream.select(station=ista)  # data for this station
        for iinstru in instrument_code:
            # select and output data for a perticular instrument code
            ista_save = False  # flag to indicate whether data of this station have been saved
            # scan different channels for getting a unified time range (choose the wider one) at a perticular station
            stdata = stdata_ista.select(channel=iinstru + '*')  # stream data of an instrument code
            if stdata.count() > 0:
                dcount = 0
                for tr in stdata:
                    if dcount == 0:
                        starttime = tr.stats.starttime
                        endtime = tr.stats.endtime
                    else:
                        starttime = min(starttime, tr.stats.starttime)
                        endtime = max(endtime, tr.stats.endtime)
                    dcount += 1
                # round datetime to the nearest second, and convert to the setted string format
                starttime_str = to_datetime(starttime.datetime).round('1s').strftime(timeformat)
                endtime_str = to_datetime(endtime.datetime).round('1s').strftime(timeformat)
                # Output data for each station and each channel
                # For a particular station, the three channel (if there are) share
                # the same time range in the final output filename.
                for icomp in component_code:  # not all component exist
                    trdata = stdata.select(component=icomp)  # stream data of a component
                    if trdata.count() > 0:
                        # if freqband is not None:
                            # filter data in specified frequency range
                            # note need to process in this way to avoide glitch after filtering
                            # trdata.detrend('demean')
                            # trdata.detrend('simple')
                            # trdata.filter('bandpass', freqmin=freqband[0], freqmax=freqband[1], corners=2,
                            #               zerophase=True)
                            # trdata.taper(max_percentage=0.001, type='cosine',
                            #              max_length=1)  # to avoid anormaly at bounday
                        # creat a folder for each station and output data in the folder
                        # the data from the same station are output to the same folder
                        dir_output_sta = os.path.join(dir_output, ista)
                        if not os.path.exists(dir_output_sta):
                            os.makedirs(dir_output_sta)
                        OfileName = trdata[0].id + '__' + starttime_str + '__' + endtime_str + '.mseed'
                        trdata.write(os.path.join(dir_output_sta, OfileName), format="MSEED")
                        ista_save = True
            if ista_save:
                break  # already save data for this station, no need to look for the next instrument code
    return
def makeStationListOnDisk(inventory, file_name,channel_list,
                    filter_network=[], filter_station=[]):
    """

   Uses fdsn to find available stations in a specific geographical location and time period.

   Parameters
   ----------
   client_list: list
       List of client names e.g. ["IRIS", "SCEDC", "USGGS"].

   min_lat: float
       Min latitude of the region.

   max_lat: float
       Max latitude of the region.

   min_lon: float
       Min longitude of the region.

   max_lon: float
       Max longitude of the region.

   start_time: str
       Start DateTime for the beginning of the period in "YYYY-MM-DDThh:mm:ss.f" format.

   end_time: str
       End DateTime for the beginning of the period in "YYYY-MM-DDThh:mm:ss.f" format.

   channel_list: str, default=[]
       A list containing the desired channel codes. Downloads will be limited to these channels based on priority. Defaults to [] --> all channels

   filter_network: str, default=[]
       A list containing the network codes that need to be avoided.

   filter_station: str, default=[]
       A list containing the station names that need to be avoided.

   Returns
   ----------
   stations_list.json: A dictionary containing information for the available stations.

    """
    station_list = {}
    for ev in inventory:
        net = ev.code
        for st in ev:
            station = st.code
            print(str(net) + "--" + str(station))
            elv = st.elevation
            lat = st.latitude
            lon = st.longitude
            new_chan = channel_list
            if len(new_chan) > 0 and (station not in station_list):
                station_list[str(station)] = {"network": net,
                                              "channels": list(set(new_chan)),
                                              "coords": [lat, lon, elv]
                                              }
    with open(file_name, 'w') as fp:
        json.dump(station_list, fp)
def format_SDS(seisdate, stainv, dir_seismic, dir_output, instrument_code=["HH", "BH", "EH", "SH", "HG", "HN"],
               location_code=['', '00', 'R1', 'BT', 'SF', '*'], freqband=None, split=False):
    """
    Format seismic data organized in SDS data structure so that the ouput data
    can be feed to various ML models.
    Seismic data sets are formated per station.
    Suitable for formatting large or long-duration data set.
    SDS fromat of data archiving structure:
        year/network_code/station_code/channel_code.D/network_code.station_code.location_code.channel_code.D.year.day_of_year
        for example: 2020/CH/VDR/HHZ.D/CH.VDR..HHZ.D.2020.234
    Instrument code has a higher priority than location code.
    Both instrument code list and location code list are priority code list, the
    program will try load only one instrument code and one location code, the code
    listed in front has higher priority.

    Parameters
    ----------
    seisdate : datetime.date
        the date of seismic data to be formated.
    stainv : obspy station inventory object.
        obspy station inventory containing the station information.
    dir_seismic : str
        path to the SDS archive directory.
    instrument_code : list of str
        the perfered list of instrument code of the input seismic data.
        We need the instrument code to look for data in SDS dirctory, this program will
        loop over this list until it can find data. Code listed first has higher priority.
        such as ["HH", "BH", "EH", "SH", "HG", "HN"] (in this example, 'HH' has the highest priority).
    dir_output : str
        directory for outputting seismic data,
        NOTE do not add '/' at the last.
    location_code : list of str, optional, default is ['','00','R1', 'BT', 'SF', '*'] (Note the last code '*' will match any location code it can find).
        the prefered list of location cods; specifying the perference order to load the data;
        For example: ['','00','R1'], in this situation '' will have the highest priority.
        If you only want load a specific location, just specify the perferred one, such as ['00'].
        If you don't want to spcify location code, use None which will use the first location_code where it can load data.
    freqband : list of float
        frequency range in Hz for filtering seismic data,
        e.g. [3, 45] meaning filter seismic data to 3-45 Hz.
        default is None, means no filtering.
    split: boolen or dict, default is False.
        whether to split the input continous data into unmasked traces without gaps.
        split['mask_value']: float, int or None
            input continous seismic data of the specified value will be recognized as gap,
            and will be masked and used to split the traces.
            This is good for filtering, because filter the contious data with
            0 (for example) filled gap will produce glitches. It is recommand
            to filter the data before merge the seismic data.
        split['minimal_continous_points'] : int
            this specifies that at least certain continuous points having the mask_value
            will be recognized as gap.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """

    if instrument_code is None:
        instrument_code = ['*']

    tdate = seisdate  # the date to be processed for SDS data files
    tyear = tdate.year  # year
    tday = tdate.timetuple().tm_yday  # day of the year

    for network in stainv:
        for station in network:
            # loop over each station for formatting input date set
            dir_stalevel = os.path.join(dir_seismic, str(tyear), network.code, station.code)  # station level

            if os.path.exists(dir_stalevel):
                # station folder exist

                for iinstru in instrument_code:
                    # loop over instrument code list to check and load data
                    dir_chalevel_want = os.path.join(dir_stalevel, iinstru + '*')
                    dir_chalevel = glob.glob(dir_chalevel_want)  # channel level
                    if len(dir_chalevel) == 0:
                        # folder of current instrument code does not exist
                        print("No data found for path: {}! Pass!".format(dir_chalevel_want))
                    elif len(dir_chalevel) <= 3:
                        # folder of current instrument code exists

                        # determine the location code
                        ilocation = None
                        if isinstance(location_code, list) and (len(location_code) == 1) and (location_code[0] != '*'):
                            # have a specific location code; only load data of that location
                            ilocation = location_code[0]
                        elif (location_code is None) or ((len(location_code) == 1) and (location_code[0] == '*')):
                            # no specifying location code list, use the first location code it can find
                            for dir_icha in dir_chalevel:
                                dir_datelevel = os.path.join(dir_icha, '*.{:03d}*'.format(tday))
                                sdatafile = glob.glob(dir_datelevel)
                                if len(sdatafile) > 0:
                                    ilocation = sdatafile[0].split(os.sep)[-1].split('.')[2]
                                    break
                        else:
                            # search avaliable location codes from the input location code preferece list
                            data_location_codes = []
                            for dir_icha in dir_chalevel:
                                dir_datelevel = os.path.join(dir_icha, '*.{:03d}*'.format(tday))
                                sdatafile = glob.glob(dir_datelevel)
                                for ifile in sdatafile:
                                    data_location_codes.append(ifile.split(os.sep)[-1].split('.')[2])
                            data_location_codes = list(set(data_location_codes))
                            for iicd in location_code:
                                location_code_filtered = fnmatch.filter(data_location_codes, iicd.upper())
                                if len(location_code_filtered) == 1:
                                    ilocation = location_code_filtered[0]
                                    print('Find data at the prefered station location code: {}.'.format(ilocation))
                                    break
                                elif len(location_code_filtered) > 1:
                                    ilocation = location_code_filtered[0]
                                    warnings.warn(
                                        'Find multiple location codes ({}) matching the current tested code {}. Choose the first one as the prefered station location code: {}.'
                                        .format(location_code_filtered, iicd, ilocation))
                                    break

                        stream = obspy.Stream()  # initilize an empty obspy stream
                        if ilocation is not None:
                            for dir_icha in dir_chalevel:
                                # loop over each channel folder to load data of the current station
                                dir_datelevel = os.path.join(dir_icha, '*.{}.*.{:03d}*'.format(ilocation,
                                                                                              tday))  # date and location level, the final filename, use day of the year to identify data
                                sdatafile = glob.glob(
                                    dir_datelevel)  # final seismic data filename for the specified station, component and date

                                if len(sdatafile) == 0:
                                    print("No data found for {}! Pass!".format(dir_datelevel))
                                elif len(sdatafile) == 1:
                                    print('Load data: {}.'.format(sdatafile[0]))
                                    stream += obspy.read(sdatafile[0])
                                else:
                                    raise ValueError(
                                        "More than one file exist: {}! This should not happen.".format(sdatafile))
                        else:
                            warnings.warn('Cannot find data from the input preferred location code list: {}.'.format(
                                location_code))

                        # output data for the current station
                        if stream.count() > 0:
                            # have at least one component data

                            if isinstance(split, dict):
                                stream = stream_split_gaps(stream, mask_value=split['mask_value'],
                                                           minimal_continous_points=split['minimal_continous_points'])
                            stream2EQTinput_raw(stream=stream, dir_output=dir_output, instrument_code=None,
                                            freqband=freqband)
                            break  # already find and output data for this instrument code, no need to look at the rest instrument codes
                            del stream
                        else:
                            warnings.warn(
                                'No data found at station {} for the specified instrument codes {}, date {} and location code {}!'.format(
                                    station.code, instrument_code, seisdate, location_code))
                            del stream
                    else:
                        warnings.warn(
                            'More than 3 folders ({}) found for the instrument code {}! Pass!'.format(dir_chalevel,
                                                                                                      iinstru))
            else:
                # station folder does not exist, no data
                warnings.warn('No data found for: {}! Pass!'.format(dir_stalevel))
    return
def read_stainv_csv(file_stainv):
    """
    Read the csv format station inventory file, and format it to obspy station inventory object.
    Each peculiar station is identified by network.station.location.instrument (such as: TA.N59A..BH)
    The input CSV file using ',' as the delimiter in which the first row
    is column name and must contain: 'network', 'station', 'latitude',
    'longitude', 'elevation'. Latitude and longitude are in decimal degree
    and elevation in meters relative to the sea-level (positive for up).
    Other optional colume names are:
    location: location code of station, such as "00", "", "01". Default: "".
    depth: the local depth or overburden of the instruments location in meter. Default: 0.
    instrument: instrument code, such as "SH", "HH", "BH", "FP";
    component: component code, shch as "ZNE", "Z12";
    Parameters
    ----------
    file_stainv : str
        filename of the input station inventory of csv format.
    Returns
    -------
    stainv : obspy station inventory object.
    """
    stainv = Inventory(networks=[])
    stadf = pd.read_csv(file_stainv, delimiter=',', encoding='utf-8',
                        header="infer", skipinitialspace=True)
    if ('instrument' in stadf.columns) and ('component' in stadf.columns):
        have_channels = True
    elif ('instrument' not in stadf.columns) and ('component' not in stadf.columns):
        have_channels = False
    else:
        raise ValueError(
            "Instrument code and component code must exist at the same time! You cannot present only one of them!")

    net_dict = {}
    for rid, row in stadf.iterrows():
        if row['network'] not in list(net_dict.keys()):
            # network not include in net_dict
            net = Network(code=row['network'], stations=[])
            net_dict[row['network']] = net

        sta = Station(code=row['station'], latitude=row['latitude'],
                      longitude=row['longitude'], elevation=row['elevation'])
        if have_channels:
            # add channel information
            if len(row['component']) != 3:
                raise ValueError("Must input three components! Current is {}!".format(row['component']))

            if ('location' in row) and (row['location'] is not None) and (not np.isnan(row['location'])):
                jlocation = row['location']
            else:
                jlocation = ""  # default location code

            if 'depth' in row:
                jdepth = row['depth']
            else:
                jdepth = 0  # default depth

            for icomp in row['component']:
                jcha = Channel(code=row['instrument'] + icomp, location_code=jlocation,
                               latitude=row['latitude'], longitude=row['longitude'],
                               elevation=row['elevation'], depth=jdepth)
                sta.channels.append(jcha)

        net_dict[row['network']].stations.append(sta)

    for inet in net_dict.keys():
        stainv.networks.append(net_dict[inet])

    return stainv
def format_ML_inputs(cfg):
    from ioseisdata import seisdata_format_4ML
    task_dir = './' + cfgs['Project'] + '/'
    if os.path.exists(task_dir):
        pass
    else:
        os.makedirs(task_dir)
    STAJSONPATH = task_dir + cfgs['InputData']['sta_json_name']
    DATASAVEPATH = task_dir + cfgs['InputData']['data_save_name']
    print('MALMI starts to format input data set for ML models:')
    DFMT = {}
    DFMT['seisdatastru_input']=cfg['InputData']['DataFormat']['seisdatastru_input']
    DFMT['dir_seismic_input'] = cfg['InputData']['DataFormat']['dir_seismic_input']
    DFMT['dir_seismic_output'] = DATASAVEPATH
    DFMT['seismic_date'] = date
    DFMT['stainv'] = station_inv
    DFMT['instrument_code'] = cfg['InputData']['DataFormat']['instrument_code']
    DFMT['freqband'] = cfg['InputData']['DataFormat']['freqband']
    DFMT['split'] = cfg['InputData']['DataFormat']['split']
    seisdata_format_4ML(DFMT=DFMT)
    gc.collect()
    print('MALMI_format_ML_inputs complete!')
    return

def import_xml_files_to_a_function(folder_path,station_inv_output):
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".xml"):
                file_path = os.path.join(root, file)
                invvvv=_read_stationxml(file_path)
                folder_to_save=station_inv_output

                fname = './{}/{}/{}/{}.{}.{}'.format(folder_to_save,'ori_dataxml' ,invvvv.networks[0].stations[0].code, invvvv.networks[0].code,invvvv.networks[0].stations[0].code,'xml')
                path = './{}/{}/{}/'.format(folder_to_save,'ori_dataxml', invvvv.networks[0].stations[0].code)
                os.makedirs(path, exist_ok=True)
                roott = ET.Element(' ')
                tree = ET.ElementTree(roott)
                tree.write(fname)
                _write_stationxml(invvvv, fname, validate=False,
                                  nsmap=None, level="response")
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='1-Data_Preparation')
    parser.add_argument('--config-file', dest='config_file', type=str, help='Configuration file path',
                        default='./Configuration_Parameters.json')
    args = parser.parse_args()
    with open(args.config_file, 'r') as f:
        cfgs = json.load(f)
    task_dir = './' + cfgs['Project'] + '/'
    if os.path.exists(task_dir):
        pass
    else:
        os.makedirs(task_dir)
    STIME = cfgs['InputData']['start_time']
    ETIME = cfgs['InputData']['end_time']
    STAJSONPATH = task_dir + cfgs['InputData']['sta_json_name']
    DATASAVEPATH = task_dir + cfgs['InputData']['data_save_name']
    if cfgs['InputData']['Datatype']== 'OnServer':
        MINLON = cfgs['InputData']['DataDownload']['minlon']
        MAXLON = cfgs['InputData']['DataDownload']['maxlon']
        MINLAT = cfgs['InputData']['DataDownload']['minlat']
        MAXLAT = cfgs['InputData']['DataDownload']['maxlat']
        CHANLIST = cfgs['InputData']['DataDownload']['channel_list']
        FILTERNETWORK = cfgs['InputData']['DataDownload']['exclude_network']
        CLIENTLIST = cfgs['InputData']['DataDownload']['client_list']
        makeStationList(client_list=["IRIS"],
                        min_lat=MINLAT,
                        max_lat=MAXLAT,
                        min_lon=MINLON,
                        max_lon=MAXLON,
                        start_time=STIME,
                        end_time=ETIME,
                        channel_list=CHANLIST,
                        filter_network=FILTERNETWORK,
                        filter_station=[])
        if os.path.exists(STAJSONPATH):
            os.remove(STAJSONPATH)
        os.rename('./station_list.json', STAJSONPATH)
        downloadMseeds(client_list=CLIENTLIST,
                       stations_json=STAJSONPATH,
                       output_dir=DATASAVEPATH,
                       start_time=STIME,
                       end_time=ETIME,
                       min_lat=MINLAT,
                       max_lat=MAXLAT,
                       min_lon=MINLON,
                       max_lon=MAXLON,
                       chunk_size=1,
                       channel_list=CHANLIST,
                       n_processor=8)
        # remove empty folders
        mseed_path = Path(DATASAVEPATH)
        for sub_path in mseed_path.glob('*'):
            if len(list(sub_path.glob('*'))) == 0:
                print('Remove Empty Folder: {}'.format(str(sub_path)))
                os.rmdir(str(sub_path))
    elif cfgs['InputData']['Datatype']== 'OnDisk':
        stanvv= read_stainv_csv(cfgs['InputData']['DataFormat']['station_inv'])
        parent_directory= cfgs['InputData']['DataFormat']['dir_seismic_input']
        ch_list = cfgs['InputData']['DataFormat']['instrument_code']
        fd=STIME;
        floatdate=UTCDateTime(fd)
        import_xml_files_to_a_function(cfgs['InputData']['DataFormat']['station_inv_input'],cfgs['Project'])
        while UTCDateTime(floatdate)<UTCDateTime(ETIME):
            makeStationListOnDisk(stanvv,task_dir+'station_list.json', ch_list)
            # assuming you have a UTCDateTime object called utc_dt
            dt = datetime.utcfromtimestamp(floatdate.timestamp)
            # stream = sds()
            search_dir=parent_directory
            format_SDS(seisdate=dt, stainv=stanvv,
                       dir_seismic=cfgs['InputData']['DataFormat']['dir_seismic_input'], dir_output=DATASAVEPATH,
                       instrument_code=['HH', 'BH', 'EH', 'SH', 'HG'], freqband=cfgs['InputData']['DataFormat']['freqband'], split=False)
            floatdate=floatdate+86400