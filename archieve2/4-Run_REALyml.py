import obspy
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
from obspy import read, UTCDateTime
from datetime import datetime
import os
import yaml
from pathlib import Path
import re
import sys
import argparse


def find_files_and_merge(directory_path, start_time_str, end_time_str, output_directory_path=None,freqband={}):
    start_time = UTCDateTime(start_time_str)
    end_time = UTCDateTime(end_time_str)

    current_day = UTCDateTime(start_time.year, start_time.month, start_time.day)
    end_day = UTCDateTime(end_time.year, end_time.month, end_time.day)

    if output_directory_path is not None:
        if not os.path.exists(output_directory_path):
            os.makedirs(output_directory_path)
    num_files_merged = 0
    while current_day <= end_day:
        current_day_str = current_day.strftime("%Y%m%d")
        file_pattern = f"*__{current_day_str}T*.mseed"
        file_path = os.path.join(directory_path, file_pattern)
        files = glob.glob(file_path)
        components = set([os.path.basename(file).split("..")[1][-45:-42] for file in files])
        for component in components:
            component_files = [file for file in files if os.path.basename(file).split("..")[1][-45:-42] == component]
            st = read(component_files[0])
            for file in component_files[1:]:
                st += read(file)
            st.merge(method=1, fill_value='latest')
            st.trim(start_time, end_time)
            num_files_merged += 1
            if output_directory_path is not None:
                output_file_name = f"{st[0].stats.network}.{st[0].stats.station}..{component}__{start_time_str}_{end_time_str}.mseed"
                output_file_path = os.path.join(output_directory_path, output_file_name)
                if freqband is not None:
                    st.detrend('demean')
                    st.detrend('simple')
                    st.filter('bandpass', freqmin=freqband[0], freqmax=freqband[1], corners=2,
                                  zerophase=True)
                    st.taper(max_percentage=0.001, type='cosine',
                                 max_length=1)  # to avoid anormaly at bounday
                st.write(output_file_path, format="MSEED")
        current_day += 86400
    return num_files_merged
"""
Workaround codes. Need Cleaning.
"""
def Primary_events(event_file,freqband):
    # read the file into a pandas dataframe
    df = pd.read_csv(event_file, delimiter='\s+', header=None)
    events=df.values
    for event in events:
        directory_path= './'+cfgs['InputData']['data_save_name']
        evedate='{}-{}-{}-{}'.format(str(event[1]), str(event[2]),str(event[3]), str(event[4]))
        dt = datetime.strptime(evedate, '%Y-%m-%d-%H:%M:%S.%f')
        evedate=dt.strftime('%Y-%m-%d-%H-%M-%S')
        output_directory_path= './'+cfgs['MIL']['control']['dir_output']+'/data_ML/primary_events/mseeds/'+evedate


            # 'primary_events'+ '/'+evedate
        # # convert datetime object to output format
        start_time_str = dt.strftime('%Y%m%dT%H%M%SZ')
        stt = UTCDateTime(start_time_str) - 30
        start_time_str= stt.strftime('%Y%m%dT%H%M%SZ')
        endtime = UTCDateTime(start_time_str)+150
        end_time_str = endtime.strftime('%Y%m%dT%H%M%SZ')
        subdirectories = [os.path.join(directory_path, d) for d in os.listdir(directory_path) if
                          os.path.isdir(os.path.join(directory_path, d))]
        for subdirectory in subdirectories:
            subdirectory_name = os.path.basename(subdirectory)
            subdirectory_output_path = os.path.join(output_directory_path, subdirectory_name)
            if not os.path.exists(subdirectory_output_path):
                os.makedirs(subdirectory_output_path)
            find_files_and_merge(subdirectory, start_time_str, end_time_str, subdirectory_output_path,freqband=freqband)

def convert2sec(t, t_ref):
    """
    convert UTCDatetime object to seconds
    Params:
    t       UTCDateTime     Time to be converted
    t_ref   UTCDateTime     Reference time
    """
    t_utc = UTCDateTime(t)
    t_ref_utc = UTCDateTime(t_ref)
    return t_utc - t_ref_utc
def runREAL(cfgs):
    """
    Run REAL Scripts
    """
    freqband=cfgs['MIL']['seismic']['freqband']
    if os.path.exists(cfgs['REAL']['seqt_catalog_dir']):
        pass
    else:
        os.makedirs(cfgs['REAL']['seqt_catalog_dir'])
    for idx in range(len(cfgs['REAL']['year'])):
        # copy temp perl file
        f_perl = open('../REAL_scripts/runREAL.pl', 'r')
        f_perl_source = f_perl.read()
        f_perl.close()
        f_perl_source = f_perl_source.replace('YEAR_KEY', cfgs['REAL']['year'][idx])
        f_perl_source = f_perl_source.replace('MON_KEY', cfgs['REAL']['mon'][idx])
        f_perl_source = f_perl_source.replace('DAY_KEY', cfgs['REAL']['day'][idx])
        f_perl_source = f_perl_source.replace('DIR_KEY', '\"' + cfgs['REAL']['seqt_dir'] + '\"')
        f_perl_source = f_perl_source.replace('STATION_KEY', cfgs['REAL']['station'])
        f_perl_source = f_perl_source.replace('TTIME_KEY', cfgs['REAL']['ttime'])
        f_perl_source = f_perl_source.replace('R_KEY', cfgs['REAL']['R'])
        f_perl_source = f_perl_source.replace('G_KEY', cfgs['REAL']['G'])
        f_perl_source = f_perl_source.replace('V_KEY', cfgs['REAL']['V'])
        f_perl_source = f_perl_source.replace('S_KEY', cfgs['REAL']['S'])
        f_perl_temp = open('../REAL_scripts/runREAL_temp.pl', 'w')
        f_perl_temp.write(f_perl_source)
        f_perl_temp.close()
        real_output = os.system('../REAL_scripts/runREAL_temp.pl')
        print('STATUS: {}'.format(real_output))
        os.rename('./catalog_sel.txt', '{}seqt_real_catalog_sel.txt'.format(cfgs['REAL']['seqt_catalog_dir']))
        os.rename('./phase_sel.txt', '{}seqt_real_phase_sel.txt'.format(cfgs['REAL']['seqt_catalog_dir']))
    if os.path.exists('{}seqt_real_catalog_sel.txt'.format(cfgs['REAL']['seqt_catalog_dir'])) and os.path.getsize('{}seqt_real_catalog_sel.txt'.format(cfgs['REAL']['seqt_catalog_dir'])) > 0:
        Primary_events('{}seqt_real_catalog_sel.txt'.format(cfgs['REAL']['seqt_catalog_dir']),freqband)
    else:
        print("No primary events to be imported to the Migration location step!")
    return

def merge_phasesel(cfgs):
    """
    Merge phase sel files
    """
    e_dict = dict()
    base_time = obspy.UTCDateTime(cfgs['REAL']['ref_time'])
    e_ID = None
    f_sel = open('{}/seqt_real_phase_sel.txt'.format(cfgs['REAL']['seqt_catalog_dir']), 'r')
    for line in f_sel.readlines():
        line_split = re.sub('\s{2,}', ' ', line).split(' ')
        if len(line_split) > 11:
            e_ID = '{}'.format(int(line_split[1]))
            e_dict[e_ID] = dict()
            real_time = base_time + float(line_split[6])
            e_dict[e_ID]['REAL_TIME'] = real_time
            e_dict[e_ID]['REAL_LAT'] = float(line_split[8])
            e_dict[e_ID]['REAL_LON'] = float(line_split[9])
            e_dict[e_ID]['REAL_DEP'] = float(line_split[10])
            e_dict[e_ID]['Picks'] = list()
        else:
            sta_name = line_split[1] + '.' + line_split[2]
            pick_type = line_split[3]
            pick_time = base_time + float(line_split[4])
            if pick_time - e_dict[e_ID]['REAL_TIME'] < 0.01:
                continue
            e_dict[e_ID]['Picks'].append([sta_name, pick_type, pick_time])
    f_sel.close()
    np.save('{}/seqt_real_e_dict.npy'.format(cfgs['REAL']['seqt_catalog_dir']), e_dict)
    return


def print_dict(e_dict):
    for key in e_dict.keys():
        print('E_ID: {} VELTime: {} LAT: {} LON: {} DEP: {}'.format(key,
                                                                    e_dict[key]['VELEST_TIME'],
                                                                    e_dict[key]['VELEST_LAT'],
                                                                    e_dict[key]['VELEST_LON'],
                                                                    e_dict[key]['VELEST_DEP']))

        print('REALTime: {} LAT: {} LON: {} DEP: {}'.format(e_dict[key]['REAL_TIME'],
                                                            e_dict[key]['REAL_LAT'],
                                                            e_dict[key]['REAL_LON'],
                                                            e_dict[key]['REAL_DEP']))
        for pick in e_dict[key]['Picks']:
            print(pick)
    return


def pad_empty_sta(cfgs):
    f = open(cfgs['REAL']['save_sta'], 'r')
    lines = f.readlines()
    f.close()
    save_folder = cfgs['REAL']['seqt_dir']
    for line in lines:
        splits = line.split(' ')
        sta_name = splits[3]
        net_name = splits[2]
        t_P_name = save_folder + net_name + '.' + sta_name + '.P.txt'
        t_S_name = save_folder + net_name + '.' + sta_name + '.S.txt'
        if os.path.exists(t_P_name):
            pass
        else:
            t_f = open(t_P_name, 'w')
            t_f.close()
        if os.path.exists(t_S_name):
            pass
        else:
            t_f = open(t_S_name, 'w')
            t_f.close()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='4-run_REAL')
    parser.add_argument('--config-file', dest='config_file',
                        type=str, help='Configuration file path', default='./Configuration_Parameters.yaml')
    args = parser.parse_args()
    cfgs = yaml.load(open(args.config_file, 'r'), Loader=yaml.SafeLoader)
    task_dir = './' + cfgs['Project'] + '/'
    os.chdir(task_dir)
    pad_empty_sta(cfgs)
    runREAL(cfgs)
    merge_phasesel(cfgs)




