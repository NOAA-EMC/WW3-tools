import matplotlib
# matplotlib.use('Agg') # for backend plots, not for rendering in a window
import time
from time import strptime
from calendar import timegm
import pandas as pd
import xarray as xr
import netCDF4 as nc
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import sys
import pandas as pd
from matplotlib import ticker
# import pickle
import sys
import warnings; warnings.filterwarnings("ignore")
import tarfile
import math

def spec_ww3(*args):
    '''
    WAVEWATCH III, wave spectrum, netcdf (.nc) or text (.spec) format
    Input: file names (list of file names), and station names (list of station names)
    Output: list of dictionaries containing:
      time(seconds since 1970),time(datetime64),lat,lon; Arrays: freq,dfreq,pwst,d1sp,dire,dspec,wnds,wndd
    '''

    if len(args) < 2:
        sys.exit(' Two inputs are required: list of file names and list of station names')

    fnames = args[0]
    stnames = args[1]
    sk = 1
    if len(args) > 2:
        sk = int(args[2])
    if len(args) > 3:
        sys.exit(' Too many inputs')

    file_names = args[0]
    station_names = args[1]
    results = []

    for fname in fnames:
        for stname in stnames:
            try:
                # Text format (only one point allowed here, same as WW3/NOAA operational)
                fp = open(fname)
                nt = fp.read().count(stname)
                fp.close()
                del fp
                if nt >= 1:
                    # Open file and read the first parameters
                    fp = open(fname)
                    cabc = fp.readline()
                    cabc = cabc.strip().split()
                    nf = int(cabc[3])  # number of frequencies
                    nd = int(cabc[4])  # number of directions
                    npo = int(cabc[5])  # number of point outputs

                    freq = zeros(nf, 'f')
                    dire = zeros(nd, 'f')
                    dspec = zeros((nt, nf, nd), 'f')
                    adire = zeros(dire.shape)
                    adspec = zeros(dspec.shape)
                    mtime = np.zeros((nt), 'd')

                    # Frequencies --------------------
                    ncf = int(np.floor(nf/8))
                    rncf = int(np.round(8*((float(nf)/8)-ncf)))
                    k = 0
                    for i in range(0, ncf):
                        line = fp.readline()
                        line = line.strip().split()
                        for j in range(0, 8):
                            freq[k] = float(line[j])
                            k = k + 1

                    if rncf > 0:
                        line = fp.readline()
                        line = line.strip().split()
                        for i in range(0, rncf):
                            freq[k] = float(line[i])
                            k = k + 1

                    # DF in frequency (dfreq)
                    dfreq = np.zeros(freq.shape[0], 'f')
                    for i in range(0, freq.shape[0]):
                        if i == 0 or i == (freq.shape[0]-1):
                            dfreq[i] = freq[i] * \
                                (1 + (((freq[-1]/freq[-2])-1)/2)) - freq[i]
                        else:
                            dfreq[i] = freq[i] * \
                                (freq[-1]/freq[-2]) - freq[i]

                    # Directions ---------------------
                    ncd = int(np.floor(nd/7))
                    rncd = int(np.round(7*((float(nd)/7)-ncd)))
                    k = 0
                    for i in range(0, ncd):
                        line = fp.readline()
                        line = line.strip().split()
                        for j in range(0, 7):
                            dire[k] = float(line[j])*180/pi
                            k = k+1

                    if rncd > 0:
                        line = fp.readline()
                        line = line.strip().split()
                        for i in range(0, rncd):
                            dire[k] = float(line[i])*180/pi
                            k = k+1

                    nl = int(floor((nf*nd)/7.))
                    rnl = int(np.round(7*((float(nf*nd)/7)-nl)))
                    auxs = np.zeros((nf*nd), 'f')
                    wnds = np.zeros((nt), 'f')
                    wndd = np.zeros((nt), 'f')

                    for t in range(0, nt):

                        cabc = fp.readline()
                        cabc.strip().split()[0]
                        mtime[t] = np.double(timegm(strptime(
                            cabc.strip().split()[0]+cabc.strip().split()[1][0:2], '%Y%m%d%H')))
                        cabc = fp.readline()
                        cabc = cabc.strip().split()
                        print(cabc)
                        if len(cabc) >8:
                            # Format: ["'42085", "'", '17.86', '-66.52', '126.1', '1.12', '146.9']
                            namep = cabc[0][1:]
                            lat = float(cabc[2])
                            lon = float(cabc[3])
                            depth = float(cabc[4])
                            wnds_index = 5
                            wndd_index = 6
                        elif len(cabc) == 8:
                            # Format: ["'46021", "'", '57.70-160.00', '50.8', '5.47', '6.9', '0.00', '270.0']
                            namep = cabc[0][1:]
                            lat_lon_str = cabc[2]
                            lat_lon_str = lat_lon_str.strip("'")
                            lat_lon_parts = lat_lon_str.split('-')
                            print(lat_lon_parts)
                            lat = float(lat_lon_parts[0])
                            lon = -float(lat_lon_parts[1])

                            depth = float(cabc[3])
                            wnds_index = 4
                            wndd_index = 5
                        else:
                            continue  # Skip this file as it cannot be processed
                            # sys.exit('Unrecognized format of cabc')

                        wnds[t] = float(cabc[wnds_index])
                        wndd[t] = float(cabc[wndd_index])

                        k = 0
                        for i in range(0, nl):
                            line = fp.readline()
                            line = line.strip().split()
                            for j in range(0, 7):
                                auxs[k] = float(line[j])
                                k = k+1

                        if rncd > 0:
                            line = fp.readline()
                            line = line.strip().split()
                            for i in range(0, rnl):
                                auxs[k] = float(line[i])
                                k = k+1

                        for ic in range(0, nf):
                            for il in range(0, nd):
                                dspec[t, ic, il] = auxs[il*nf+ic]

                    fp.close()
                    del fp

                    mdate = pd.to_datetime(mtime, unit='s').strftime(
                        '%Y-%m-%dT%H:%M:%S.%f')
                    freq1 = freq*np.nan
                    freq2 = freq*np.nan

                    # ------------------
                    # 1D power spectrum
                    pwst = np.zeros((dspec.shape[0], nf), 'f')
                    for t in range(0, dspec.shape[0]):
                        for il in range(0, nf):
                            pwst[t, il] = sum(
                                dspec[t, il, :]*(2*np.pi)/nd)

                        pwst[t, :] = pwst[t, :]*dfreq[:]

                    # organizing directions  -----
                    adspec = np.copy(dspec)
                    inddire = int(np.where(dire == min(dire))[0][0])
                    for t in range(0, dspec.shape[0]):
                        adspec[t, :, 0:nd-(inddire+1)] = dspec[t,:, (inddire+1):]
                        adspec[t, :, nd-(inddire+1):nd] = dspec[t,:, :(inddire+1)]
                        for i in range(0, nd):
                            dspec[t, :, i] = adspec[t, :, nd-i-1]

                        adspec[t, :, :int(nd/2)] = dspec[t, :, int(nd/2):]
                        adspec[t, :, int(nd/2):] = dspec[t, :, :int(nd/2)]
                        dspec[t, :, :] = adspec[t, :, :]

                    dire = np.sort(dire)

                    # 1D directional spectrum
                    d1sp = np.zeros((dspec.shape[0], nf), 'f')

                    for t in range(0, dspec.shape[0]):
                        for il in range(0, nf):
                            a = np.sum(dspec[t, il, :] * np.array(
                                np.sin((pi*dire)/180.)/np.sum(dspec[t, il, :])))
                            b = np.sum(dspec[t, il, :] * np.array(
                                np.cos((pi*dire)/180.)/np.sum(dspec[t, il, :])))
                            aux = math.atan2(a, b)*(180./pi)
                            if aux < 0:
                                aux = aux+360.

                            d1sp[t, il] = float(aux)
                            del a, b, aux

                    hs = np.sqrt(2 * np.trapz(np.trapz(dspec, x=dire, axis=-1), x=freq, axis=-1))
                    tp = freq[np.argmax(np.max(dspec, axis=-1), axis=-1)]

                    # build dictionary
                    result = {'time': mtime, 'date': mdate, 'latitude': lat, 'longitude': lon, 'depth': depth,
                              'wind_spd': wnds, 'wind_dir': wndd, 'freq': freq, 'freq1': freq1, 'freq2': freq2,
                              'deltafreq': dfreq, 'pspec': pwst, 'theta': dire, 'dmspec': d1sp, 'dirspec': dspec, 'Hs': hs, 'Tp': tp,'station_name': stname}

                    results.append(result)
            except Exception as e:
                print(f"Skipping file {fname} for station {stname}: {str(e)}")
                continue

    return results

# Example usage:
# spec_ww3(['file1.spec', 'file2.spec'], ['station1', 'station2'])

