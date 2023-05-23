"""
Coordinate convention: X-axis -> East; Y-axis -> North; Z-axis -> vertical down.

"""


import os
import numpy as num
import matplotlib
matplotlib.use('Agg')  # set the backend before importing pyplot
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import datetime
import copy
import gc
from scipy import stats
from loki import ioformatting
#
from loki import traveltimes
from loki import waveforms
from loki import stacktraces
from loki import LatLongUTMconversion
import tt_processing                       # C
import location_t0                         # C  for multiplying the P- and S-stacking values using this
import location_4D                         # C  for obtaining 4D stacking matrix
#import location_t0_tmax                   # C  for also output max coherency over time samples
#import location_t0_plus                   # C  for adding the P- and S-stacking values using this


class Loki:
    """docstring for Loki"""

    def __init__(self, data_path, output_path, db_path, hdr_filename, mode='locator'):
        self.data_path = data_path
        self.output_path = output_path
        self.db_path = db_path
        self.hdr_filename = hdr_filename
        if mode == 'locator':
            self.data_tree, self.events = self.location_data_struct(self.data_path, self.output_path)
        elif mode == 'detector':
            self.data_tree, self.events = self.detection_data_struct(self.data_path, self.output_path)
        else:
            raise ValueError('mode must be "detector" or "locator"')

    def location_data_struct(self, data_path, output_path):
        events=[]
        data_tree=[]
        for root, dirs, files in os.walk(data_path):
            if (not dirs) and (files) and (files[0][-4:] != '.txt'):
                data_tree.append(root)
        data_tree.sort()
        events = [idtree.split('/')[-1] for idtree in data_tree]
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        return data_tree, events

    def detection_data_struct(self, data_path, output_path):
        events = []
        data_tree = []
        for root, dirs, files in os.walk(data_path):
            if not dirs:
                data_tree.append(root)
        data_tree.sort()
        events = [idtree.split('/')[-1] for idtree in data_tree]
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        return data_tree, events

    def location(self, extension='*', comp=['E', 'N', 'Z'], precision='single', **inputs):
        if 'tshortp_min' in inputs:
            # need to calculate STA/LTA for stacking
            STALTA = True
            tshortp_min = inputs['tshortp_min']
            tshortp_max = inputs['tshortp_max']
            tshorts_min = inputs['tshorts_min']
            tshorts_max = inputs['tshorts_max']
            slrat = inputs['slrat']
            ntrial = inputs['ntrial']
        
            tshortp = num.linspace(tshortp_min, tshortp_max, ntrial)
            tshorts = num.linspace(tshorts_min, tshorts_max, ntrial)
        else:
            # no need to calculate STA/LTA ratio for stacking
            STALTA = False
            ntrial = 1
            
        npr = inputs['npr']
        model = inputs['model']
        
        # whether filtering input data
        if 'freq' in inputs:
            freq = inputs['freq']
        else:
            freq = None  # no filtering
            
        # whether ouput the computed CFs
        if 'opsf' in inputs:
            opsf = inputs['opsf']
        else:
            opsf = False
            
        # whether to compute element wise power over the phase probabilities before stacking
        if 'ppower'  not in inputs:
            inputs['ppower'] = None
        
        # whether to output the migration volume
        if 'output_migv' not in inputs:
            inputs['output_migv'] = True
        
        # to calculate 3D (xyz) or 4D (xyzt) stacking matrix, default is 3D
        if 'migv_4D' not in inputs:
            inputs['migv_4D'] = False
        
        # set how to identify each unique station for waveforms.Waveforms
        # can be 'station', 'network.station', 'network.station.location', 'network.station.location.instrument'
        if 'station_idmode' in inputs:
            station_idmode = inputs['station_idmode']
        else:
            station_idmode = 'station'

        # create the file for outputting catalog
        ff = open(self.output_path+'/'+'catalogue', 'a')
        ff.close()
        
        # load traveltime data set
        tobj = traveltimes.Traveltimes(self.db_path, self.hdr_filename)
        tp = tobj.load_traveltimes('P', model, precision)
        ts = tobj.load_traveltimes('S', model, precision)

        for event_path in self.data_tree:
            wobj = waveforms.Waveforms(event_path=event_path, extension=extension, comps=comp, freq=freq, station_idmode=station_idmode)
            sobj = stacktraces.Stacktraces(tobj, wobj, **inputs)
            # event = sobj.evid
            event = event_path.split('/')[-1]

            print('Processing to the event folder: ', event_path, event)
            if os.path.isdir(self.output_path+'/'+event):
                continue
            else:
                os.mkdir(self.output_path+'/'+event)

            tp_modse, ts_modse = sobj.time_extractor(tp, ts)  # traveltime table in second
            tp_mod, ts_mod = tt_processing.tt_f2i(sobj.deltat, tp_modse, ts_modse, npr)  # traveltime table in time sample, for each imaging point traveltimes have substracted the minimal P traveltime

            cmax_pre = -1.0
            for i in range(ntrial):
                if STALTA:
                    # need to calculate STA/LTA from the characteristic funtion
                    # then stack the STA/LTA for imaging
                    nshort_p = round(tshortp[i]/sobj.deltat)
                    nshort_s = round(tshorts[i]/sobj.deltat)
                    obs_dataP, obs_dataS = sobj.loc_stalta(nshort_p, nshort_s, slrat, norm=1)
                else:
                    # no need to calculate STA/LTA 
                    # directly stack the characteristic function for imaging
                    obs_dataP = sobj.obs_dataV  # vertical -> P
                    obs_dataS = sobj.obs_dataH  # horizontal -> S

                if opsf:
                    # output the characteristic functions for stacking
                    datainfo = {}
                    datainfo['dt'] = sobj.deltat
                    datainfo['starttime'] = sobj.dtime_max
                    for ista, sta in enumerate(sobj.stations):
                        datainfo['station_name'] = sta
                        datainfo['channel_name'] = 'CFP'  # note maximum three characters, the last one must be 'P'
                        ioformatting.vector2trace(datainfo, obs_dataP[ista,:], self.output_path+'/'+event+'/cf/trial{}'.format(i))
                        datainfo['channel_name'] = 'CFS'  # note maximum three characters, the last one must be 'S'
                        ioformatting.vector2trace(datainfo, obs_dataS[ista,:], self.output_path+'/'+event+'/cf/trial{}'.format(i))

                if inputs['migv_4D']:
                    corrmatrix = location_4D.stacking(tp_mod, ts_mod, obs_dataP, obs_dataS, npr)
                    corrmatrix = num.reshape(corrmatrix, (sobj.ns, tobj.nx, tobj.ny, tobj.nz))
                    cmax_indx = num.unravel_index(num.argmax(corrmatrix, axis=None), corrmatrix.shape)
                    indx_grid = num.ravel_multi_index(cmax_indx[1:], corrmatrix.shape[1:])
                    indx_time = cmax_indx[0]
                else:
                    # corrmatrix is the 3D stacking matrix, in 1D format but can be 
                    # reformat to 3D format, each point saves the maximum stacking 
                    # value during this calculation time period
                    iloctime, corrmatrix = location_t0.stacking(tp_mod, ts_mod, obs_dataP, obs_dataS, npr)  # iloctime[0]: the grid index of the maximum stacking point; iloctime[1]: the time index at the maximum stacking point
                    #iloctime, corrmatrix, cohmaxtt = location_t0_tmax.stacking(tp_mod, ts_mod, obs_dataP, obs_dataS, npr)  # iloctime[0]: the grid index of the maximum stacking point; iloctime[1]: the time index at the maximum stacking point
                    indx_grid = iloctime[0]
                    indx_time = iloctime[1]
                    corrmatrix = num.reshape(corrmatrix,(tobj.nx,tobj.ny,tobj.nz))
                    
                evtpmin = num.amin(tp_modse[indx_grid,:])
                event_t0 = sobj.dtime_max + datetime.timedelta(seconds=indx_time*sobj.deltat) - datetime.timedelta(seconds=evtpmin)  # event origin time
                event_t0s = (event_t0).isoformat()
                cmax = num.max(corrmatrix)
                (ixloc, iyloc, izloc) = num.unravel_index(indx_grid,(tobj.nx,tobj.ny,tobj.nz))
                xloc = tobj.x[ixloc]
                yloc = tobj.y[iyloc]
                zloc = tobj.z[izloc]
                
                # output the current location result
                if ntrial > 1:
                    cmfilename = self.output_path+'/'+event+'/'+event
                else:
                    cmfilename = self.output_path+'/'+event+'/'+event_t0s
                out_file = open(cmfilename+'.loc', 'a')
                if STALTA:
                    out_file.write(str(i)+' '+str(xloc)+' '+str(yloc)+' '+str(zloc)+' '+str(cmax)+' '+str(nshort_p)+' '+str(nshort_s)+' '+str(slrat)+'\n')
                else:
                    out_file.write(str(i)+' '+str(xloc)+' '+str(yloc)+' '+str(zloc)+' '+str(cmax)+'\n')
                out_file.close()
                
                # save the stacked coherence matrix
                if inputs['output_migv']:
                    num.save(self.output_path+'/'+event+'/'+'corrmatrix_trial_'+str(i),corrmatrix)
                
                # # save the maximum stacked coherency over time
                # num.save(self.output_path+'/'+event+'/'+'cohmaxtt_trial_'+str(i),cohmaxtt)
                
                # plot migration profiles
                self.coherence_plot(self.output_path+'/'+event, corrmatrix, tobj.x, tobj.y, tobj.z, i)
            
                # output theoretical P- and S-wave arrivaltimes
                fname = cmfilename + '_trial{}.phs'.format(i)
                self.write_phasetime(sobj.stations, event_t0, tp_modse, ts_modse, indx_grid, fname)
                
                if cmax > cmax_pre:
                    event_t0s_final = copy.deepcopy(event_t0s)
                    cmax_pre = copy.deepcopy(cmax)
            
            self.catalogue_creation(event, event_t0s_final, tobj.lat0, tobj.lon0, ntrial, corrmatrix)
        print('Location process completed!!!')
        gc.collect()

    def catalogue_creation(self, event, event_t0s, lat0, lon0, ntrial, corrmatrix, refell=23):
        (zorig, eorig, norig) = LatLongUTMconversion.LLtoUTM(refell, lat0, lon0) #da adeguare in python 3
        if (ntrial > 1):
            ev_file = self.output_path+'/'+event+'/'+event+'.loc'
            data = num.loadtxt(ev_file)
            w = num.sum(data[:, 4])
            xb = ((num.dot(data[:, 1], data[:, 4])/w)*1000)+eorig
            yb = ((num.dot(data[:, 2], data[:, 4])/w)*1000)+norig
            late, lone = LatLongUTMconversion.UTMtoLL(refell, yb, xb, zorig)
            zb = num.dot(data[:, 3], data[:, 4])/w  # depth in km
            cb = num.mean(data[:, 4])  # the mean coherence over the ntrial realizations
            cmax = num.max(data[:, 4])  # the maximum coherence over the ntrial realizations
            merr = num.vstack((data[:, 1], data[:, 2], data[:, 3]))
            err = num.cov(merr)
            errmax = num.sqrt(num.max(num.linalg.eigvals(err)))
            
            f = open(self.output_path+'/'+'catalogue', 'a')
            f.write(event_t0s+'    '+str(late)+'   '+str(lone)+'   '+str(zb)+'   '+str(errmax)+'   '+str(cb)+'   '+str(cmax)+'\n')
            f.close()    
        else:
            ev_file = self.output_path+'/'+event+'/'+event_t0s+'.loc'
            data = num.loadtxt(ev_file)
            late, lone = LatLongUTMconversion.UTMtoLL(refell, (data[2]*1000)+norig, (data[1]*1000)+eorig, zorig)  # latitude, longitude
            evdepth_km = data[3]  # depth in km
            migv_max = data[4]  # the maximum coherence over the corrmatrix
            
            migv_std = num.std(corrmatrix, axis=None)  # the coherence standard deviation of the corrmatrix
            migv_median = num.median(corrmatrix, axis=None)  # the median coherence of the corrmatrix
            migv_mean = num.mean(corrmatrix, axis=None)  # the mean value of the migration data volume
            migv_min = num.amin(corrmatrix, axis=None)  # the minimal value of the migration data volume
            migv_MAD = stats.median_abs_deviation(corrmatrix, axis=None, scale=1, nan_policy='omit')  # median absolute deviation of the migration data volume
            migv_kurtosis = stats.kurtosis(corrmatrix, axis=None, nan_policy='omit')  # kurtosis of the migration data volume
            migv_skewness = stats.skew(corrmatrix, axis=None, nan_policy='omit')  # skewness of the migration data volume
            
            # nomalize corrmatrix, let minimal->0, maximum->1
            n1 = 0.0  # minimal limit
            n2 = 1.0  # maximum limit
            dmax = num.amax(corrmatrix, axis=None, keepdims=True)
            dmin = num.amin(corrmatrix, axis=None, keepdims=True)
            k = (n2-n1)/(dmax-dmin)
            b = (dmax*n1-dmin*n2)/(dmax-dmin)
            corrmatrix = k*corrmatrix + b
            
            migv_normstd = num.std(corrmatrix, axis=None)  # the coherence standard deviation of the corrmatrix
            migv_normMAD = stats.median_abs_deviation(corrmatrix, axis=None, scale=1, nan_policy='omit')  # median absolute deviation of the migration data volume
            migv_normkurtosis = stats.kurtosis(corrmatrix, axis=None, nan_policy='omit')  # kurtosis of the migration data volume
            migv_normskewness = stats.skew(corrmatrix, axis=None, nan_policy='omit')  # skewness of the migration data volume
            
            f = open(self.output_path+'/'+'catalogue', 'a')
            f.write(event_t0s+'    '+str(late)+'    '+str(lone)+'    '+str(evdepth_km)+'    '+str(migv_std)+'    '+str(migv_median)+'    '+str(migv_max)
                    +'    '+str(migv_mean)+'    '+str(migv_min)+'    '+str(migv_MAD)+'    '+str(migv_kurtosis)+'    '+str(migv_skewness)+'    '
                    +str(migv_normstd)+'    '+str(migv_normMAD)+'    '+str(migv_normkurtosis)+'    '+str(migv_normskewness)+'\n')
            f.close()

    def coherence_plot(self, event_path, corrmatrix, xax, yax, zax, itrial, normalization=False, figfmt='.png'):
        
        if corrmatrix.ndim == 4:
            # stacking matrix in form of (nt, nx, ny, nz)
            nt, nx, ny, nz = num.shape(corrmatrix)
            CXY = num.zeros([ny, nx])
            for i in range(ny):
                for j in range(nx):
                    CXY[i,j]=num.max(corrmatrix[:,j,i,:])
    
            CXZ = num.zeros([nz, nx])
            for i in range(nz):
                for j in range(nx):
                    CXZ[i, j] = num.max(corrmatrix[:,j,:,i])
    
            CYZ = num.zeros([nz, ny])
            for i in range(nz):
                for j in range(ny):
                    CYZ[i, j] = num.max(corrmatrix[:,:, j, i])
        elif corrmatrix.ndim == 3:
            # stacking matrix in form of (nx, ny, nz)
            nx, ny, nz = num.shape(corrmatrix)
            CXY = num.zeros([ny, nx])
            for i in range(ny):
                for j in range(nx):
                    CXY[i,j]=num.max(corrmatrix[j,i,:])
    
            CXZ = num.zeros([nz, nx])
            for i in range(nz):
                for j in range(nx):
                    CXZ[i, j] = num.max(corrmatrix[j,:,i])
    
            CYZ = num.zeros([nz, ny])
            for i in range(nz):
                for j in range(ny):
                    CYZ[i, j] = num.max(corrmatrix[:, j, i])

        if normalization:
            nrm = Normalize(vmin=0., vmax=1.)
        else:
            nrm = None

        fig = plt.figure(dpi=600)
        fig.suptitle('Coherence matrix X-Y', fontsize=14, fontweight='bold')
        ax = fig.gca()
        cmap = plt.cm.get_cmap('jet', 100)
        cs = plt.contourf(xax, yax, CXY, 20, cmap=cmap, interpolation='bilinear', norm=nrm)
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        cbar = plt.colorbar(cs)
        ax.set_aspect('equal')
        figname = os.path.join(event_path, 'Coherence_matrix_xy'+str(itrial)+figfmt)
        plt.savefig(figname, dpi=600)

        fig = plt.figure(dpi=600)
        fig.suptitle('Coherence matrix X-Z', fontsize=14, fontweight='bold')
        ax = fig.gca()
        cmap = plt.cm.get_cmap('jet', 100)
        cs = plt.contourf(xax, zax, CXZ, 20, cmap=cmap, interpolation='bilinear', norm=nrm)
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Z (km)')
        cbar = plt.colorbar(cs)
        ax.invert_yaxis()
        ax.set_aspect('equal')
        figname = os.path.join(event_path, 'Coherence_matrix_xz'+str(itrial)+figfmt)
        plt.savefig(figname, dpi=600)

        fig = plt.figure(dpi=600)
        fig.suptitle('Coherence matrix Y-Z', fontsize=14, fontweight='bold')
        ax = fig.gca()
        cmap = plt.cm.get_cmap('jet', 100)
        cs = plt.contourf(yax, zax, CYZ, 20, cmap=cmap, interpolation='bilinear', norm=nrm)
        ax.set_xlabel('Y (km)')
        ax.set_ylabel('Z (km)')
        cbar = plt.colorbar(cs)
        ax.invert_yaxis()
        ax.set_aspect('equal')
        figname = os.path.join(event_path, 'Coherence_matrix_yz'+str(itrial)+figfmt)
        plt.savefig(figname, dpi=600)
        plt.close("all")
        
    
    def write_phasetime(self, stations, event_t0, tp_modse, ts_modse, grididx, fname):
        """
        Calculate the theoretical arrival-times of P- and S-phases for the located
        event and output to a text file.

        Parameters
        ----------
        stations : list of str
            station names.
        event_t0 : datetime
            event origin time.
        tp_modse : numpy array, shape: n_stations*n_grids
            P-wave traveltime table in second.
        ts_modse : numpy array, shape: n_stations*n_grids
            S-wave traveltime table in second.
        grididx : int
            grid index where the seismic event is located.
        fname : str
            output filename including path.

        Returns
        -------
        None.

        """
        
        ofile = open(fname, 'a')
        ofile.write('# station    P_arrivaltime        S_arrivaltime \n')
        
        for ii, sta in enumerate(stations):
            # loop over each station to output the theoretical arrival-times for the P- and S-phases
            tp_tavt = event_t0 + datetime.timedelta(seconds=tp_modse[grididx,ii])  # P_arrival-time = event_origin_time + P_traveltime
            ts_tavt = event_t0 + datetime.timedelta(seconds=ts_modse[grididx,ii])  # S_arrival-time = event_origin_time + S_traveltime
            ofile.write(sta+' '+tp_tavt.isoformat()+' '+ts_tavt.isoformat()+'\n')
            ofile.flush()
            
        ofile.close()        
        
        return

