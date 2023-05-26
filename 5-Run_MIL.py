#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 2023

@authors: Peidong Shi, Hamzeh Mohammadigheymasi
"""

"""
Runs Migration location method on primary list of events associated by REAL method
"""




import sys
sys.path.append('./src')  # user need to change the path here!
from mainMALMI import MALMI
import datetime
import json
import os
import argparse
from pathlib import Path
from catalog_plot import catalog_plot_depth, catalog_plot_otime



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='1_run_EqT')
    parser.add_argument('--config-file', dest='config_file', type=str, help='Configuration file path',
                        default='./Configuration_Parameters.json')
    args = parser.parse_args()
    with open(args.config_file, 'r') as f:
        cfgs = json.load(f)
    task_dir = './' + cfgs['Project'] + '/'
    # %% Initialize MALMI
    coseismiq = MALMI(cfgs=cfgs)
     # %% Run ML models to get continuous phase probabilities
    ML = {}
    ML['engine'] = cfgs['MIL']['ML']['engine']
    ML['model'] = cfgs['MIL']['ML']['model']
    ML['seqtpicks'] = './' + cfgs['Project'] + '/SEqTCatalogs/seqt_real_e_dict.npy'
    # ML['model'] =  '/home/shipe/codes/EQTransformer/ModelsAndSampleData/EqT_model.h5'  # path to a trained EQT model
    ML[
        'overlap'] = 0.9  # overlap rate of time window for generating probabilities. e.g. 0.6 means 60% of time window are overlapped
    coseismiq.generate_prob(ML)
    # %% Detect locatable events from continuous phase probabilities
    coseismiq.event_detect_ouput()
    # %% Migration location for each event
    coseismiq.migration()
    # %% Generate waveform plots for each event
    coseismiq.rsprocess_view()
    # # %% Delete the input continuous seismic data for ML models for saving disk space
    # CL = {}
    # CL['hdf5_prob'] = True
    # coseismiq.clear_interm(CL)
    # # extract catalog and plot
    CAT = {}
    CAT['dir_output'] = task_dir + cfgs['MIL']['CAT']['dir_output']
    CAT['fname'] = cfgs['MIL']['CAT']['fname']
    CAT['fformat'] = cfgs['MIL']['CAT']['fformat']
    catalog = coseismiq.get_catalog(CAT)
    # catalog_plot_depth(coseismiq.grid['pltregion'], catalog, depthrg=None, cmap="hot", sta_inv=coseismiq.stainv, mkregion=coseismiq.grid['mgregion'],
    #                     fname="../figure/catalog_depth.png", plot_stationname=False, eq_size=0.09, markers=None)
    # catalog_plot_otime(coseismiq.grid['pltregion'], catalog, time_ref=None, cmap="hot", sta_inv=coseismiq.stainv, mkregion=coseismiq.grid['mgregion'],
    #                     fname="../figure/catalog_time.png", plot_stationname=False, eq_size=0.09, markers=None)







