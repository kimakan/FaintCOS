'''This script subtracts the background using darkframes in
the 'ldark' directory. Additionally, it calculates the confidence limits 
and combines all spectra in the working directory. 

Author: Kirill Makan
'''
version = '1.1'

import os
import sys
import math


from astropy.io import fits
from astropy.table import Table, vstack, hstack, Column
from astropy.stats import poisson_conf_interval
from datetime import datetime, timezone
from scipy.interpolate import interp1d
from scipy.stats import linregress
import numpy as np

from faintcos_config import *




def calc_conf_lim_feldman_cousins(counts, bkg):
    ''' Calculates confidence limits according to Feldman & Cousins 1998.
    
    Parameters:
    ----------
    counts: numpy array, source counts
    
    bkg: numpy array, background counts
    
    
    Returns:
    ---------
    cnt_err_down, cnt_err_up : two numpy arrays with lower and upper confidence
                               limits at 1 sigma
    
    '''
    # in order to install CustomConfLim, run 'python setup.py build' and 
    # 'python setup.py install' in the CustomConfLimits folder
    import CustomConfLim    
    
    cnt_err_up = np.zeros(shape=len(counts), dtype=np.float32)
    cnt_err_down = np.zeros(shape=len(counts), dtype=np.float32)
    
    print("Error calculation (Feldman & Cousins): ")
    for i in range(len(counts)):
        sys.stdout.write("|")
        for p in range(50):
            if i > p*len(counts)/50.:
                sys.stdout.write("=")
            else:
                sys.stdout.write(" ")
        sys.stdout.write("| " + \
                         str(round(100.*float(i)/len(counts))) + \
                         " %\r")
        sys.stdout.flush()
        
        if counts[i] > 150:
            limits = (math.sqrt(counts[i]), \
                      math.sqrt(counts[i]))
        else:
            if bkg[i] > counts[i]:
                limits = CustomConfLim.bootstrap_bkg_conf_lim(counts[i], \
                                                              bkg[i], 10000)         
            elif bkg[i] > 0.0:
                limits = \
                CustomConfLim.feldman_cousins_conf_lim(round(counts[i]), \
                                                              bkg[i])
            else:
                limits = (0.0, 0.0)
                
        cnt_err_down[i] = limits[0]
        cnt_err_up[i] = limits[1]    
        if cnt_err_down[i] < 0:
            cnt_err_down[i] = 0.0
        if cnt_err_up[i] < 0:
            cnt_err_up[i] = 0.0
    
    print("\n")
    return cnt_err_down, cnt_err_up

def calc_lya_scatter_model(gcounts, wave_arr):
    """ Calculates the contamination by Lya according to Worseck et al. 2016
    
    Parameters:
    ----------
    gcounts: numpy array, gcounts containing Lya region
    
    wave_arr: numpy array, wavelength grid of gcounts in angstrom
    
    
    Returns:
    ----------
    Blya: numpy array, contamination by Lya in counts in wave_arr grid
    
    Blya_err_up: numpy array, count upper error of the Blya
    
    Blya_err_down: numpy array, count upper error of the Blya
    """
    
    a = 1.7348e-5
    lam_0 = 1254.6
    b = 100.9
    
    # Relative upper error of the Lya contamination
    rel_err = np.array([[   801.,  853.,  905.,  957., 1009., 1061., 1113.,\
                            1165., 1217., 1269., 1321., 1373., 1425., 1477.,\
                            1529., 1581., 1633., 1685., 1737., 1789., 1841.,\
                            1893., 1945., 1997.],\
                           [2.976524   ,2.029496    ,1.3314583   ,0.8864973 ,\
                            0.5695007  ,0.34901363  ,0.1948071   ,0.09911883,\
                            0.05120203 ,0.04276433  ,0.04487024  ,0.05724196,\
                            0.11888563 ,0.23237737  ,0.40143612  ,0.6552495 ,\
                            1.0261706  ,1.4999093   ,2.2647116   ,3.3137312 ,\
                            4.849577   ,7.3380733   ,11.06501    ,16.91639 ],\
                           [0.7585422  ,0.67552364  ,0.57972205  ,0.4754494 ,\
                            0.37029144 ,0.2608441   ,0.169427    ,0.09658019,\
                            0.05248015 ,0.04616276  ,0.04640767  ,0.06427228,\
                            0.11368275 ,0.19347219  ,0.28874132  ,0.4030617 ,\
                            0.51270425 ,0.613789    ,0.70610535  ,0.78464127,\
                            0.84194565 ,0.8906022   ,0.92548954  ,0.9521525 ]])
    err_up_sys = np.interp(wave_arr, rel_err[0], rel_err[1])
    err_down_sys = np.interp(wave_arr, rel_err[0], rel_err[2])
    
    CountsInLyaRegion = gcounts[(wave_arr <= 1215.97) & (wave_arr >= 1215.37)]
    Clya = np.mean(CountsInLyaRegion)
    
    Blya = a*Clya*np.exp(-np.power(wave_arr - lam_0, 2.)/(2.*b**2))
    
    Clya_err_stat = math.sqrt(Clya)/math.sqrt(len(CountsInLyaRegion))/Clya
    err_up = np.sqrt(err_up_sys**2 + Clya_err_stat**2)
    err_down = np.sqrt(err_down_sys**2 + Clya_err_stat**2)
    
    return Blya, err_up*Blya, err_down*Blya
    



def calc_conf_lim_kraft(counts, bkg):
    ''' Calculates confidence limits according to Kraft et al. 1991.
    
    Parameters:
    ----------
    counts: numpy array, source counts
    
    bkg: numpy array, background counts
    
    
    Returns:
    ---------
    cnt_err_down, cnt_err_up : two numpy arrays with lower and upper confidence
                               limits at 1 sigma
    
    '''
    cnt_err_up = np.zeros(shape=len(counts), dtype=np.float32)
    cnt_err_down = np.zeros(shape=len(counts), dtype=np.float32)

    print("Error calculation (Kraft): ")
    for i in range(len(counts)):
        sys.stdout.write("|")
        for p in range(50):
            if i > p*len(counts)/50.:
                sys.stdout.write("=")
            else:
                sys.stdout.write(" ")
        sys.stdout.write("| " + \
                         str(round(100.*float(i)/len(counts))) + \
                         " %\r")
        sys.stdout.flush()
        
        if counts[i] > 99:
            cnt_err_up[i] = math.sqrt(counts[i] - bkg[i])
            cnt_err_down[i] = cnt_err_up[i]
        else:
            try:
                limits = poisson_conf_interval(counts[i], \
                                               background=bkg[i], \
                                               confidence_level=0.6827, \
                                               interval='kraft-burrows-nousek')
                if counts[i] > bkg[i]:
                    cnt_err_down[i] = counts[i] - bkg[i] - limits[0]
                else:
                    cnt_err_down[i] = 0.0
                cnt_err_up[i] = limits[1] - (counts[i] - bkg[i])
            except ValueError:
                cnt_err_down[i] = 0.0
                cnt_err_up[i] = 0.0
            
        if cnt_err_down[i] < 0:
            cnt_err_down[i] = 0.0
        if cnt_err_up[i] < 0:
            cnt_err_up[i] = 0.0
    
    print("\n")
    return cnt_err_down, cnt_err_up



def createCDRPrimaryHeader(hdul, wave):
    ''' Calculates confidence limits according to Kraft et al. 1991.
    
    Parameters:
    ----------
    hdul: hdu list from the x1d file 
    
    wave: wavelength array of the spectroscopic data
    
    
    Returns:
    ---------
    hdu: custom primary header for the cdr file 
    
    '''
    
    hdu = fits.PrimaryHDU()
    hdu.header['DATE'] = (datetime.utcnow().isoformat(timespec='seconds'),\
                            "File creation date")
    hdu.header['FILETYPE'] = (hdul[0].header['FILETYPE'], 
                              'type of data found in the data file')
    hdu.header['TELESCOP'] = (hdul[0].header['TELESCOP'],
                              'telescope used to acquire data')
    hdu.header['INSTRUME'] = (hdul[0].header['INSTRUME'], 
                              'identifier for instrument used to acquire data')
    hdu.header['EQUINOX'] = (hdul[0].header['EQUINOX'], 
                             'equinox of celestial coord. system')
    
    hdu.header['TARGNAME'] = (hdul[0].header['TARGNAME'], 
                              'proposer\'s target name')
    hdu.header['RA_TARG'] = (hdul[0].header['RA_TARG'], 
                             'right ascention of the target (deg)')
    hdu.header['DEC_TARG'] = (hdul[0].header['DEC_TARG'], 
                              'declination of the target (deg)')
    
    hdu.header['PROPOSID'] = (hdul[0].header['PROPOSID'], 
                              'PEP proposal identifier')
    
    hdu.header['OPUS_VER'] = (hdul[0].header['OPUS_VER'], 
                              'data processing software system version')
    hdu.header['CSYS_VER'] = (hdul[0].header['CSYS_VER'], 
                              'calibration software system version id')
    hdu.header['CAL_VER'] = (hdul[0].header['CAL_VER'], 
                             'CALCOS code version')
    hdu.header['FCOS_VER'] = (version, 
                              'FaintCOS code version')
    
    hdu.header['DETECTOR'] = (hdul[0].header['DETECTOR'],
                              'FUV or NUV')
    hdu.header['CENWAVE'] = (hdul[0].header['CENWAVE'], 
                             'central wavelength of the spectrum')
    
    hdu.header['LIFE_ADJ'] = (hdul[0].header['LIFE_ADJ'], 
                            'Life Time Adjustment Position')
    hdu.header['SEGMENT'] = (hdul[0].header['SEGMENT'], 
                            'FUV detector segment name (FUVA, FUVB or BOTH)')
    hdu.header['OPT_ELEM'] = (hdul[0].header['OPT_ELEM'], 
                              'optical element in use')
    
    
    hdu.header['BANDWID'] = (max(wave)-min(wave), 
                             'bandwidth of the data')
    res = cos_res[(cos_res['LP'] == hdu.header['LIFE_ADJ']) & \
                  (cos_res['OPT_ELEM'] == hdu.header['OPT_ELEM']) &\
                  (cos_res['CENWAVE'] == hdu.header['CENWAVE'])]['R'][0]
    hdu.header['SPECRES'] = (res, 
                             'approx. resolving power at central wavelength')
    hdu.header['CENTRWV'] = ((max(wave)+min(wave))/2.0, 
                             'central wavelength of the data')
    hdu.header['MINWAVE'] = (min(wave), 
                             'minimum wavelength in spectrum')
    hdu.header['MAXWAVE'] = (max(wave), 
                             'maximum wavelength in spectrum')
    
    
    if 'DATE-OBS' in hdul[0].header:
        hdu.header['DATE-OBS'] = (hdul[0].header['DATE-OBS'], 
                                  'UT date of start of observation(yyyy-mm-dd)')
        hdu.header['EXPSTART'] = (hdul[0].header['EXPSTART'], 
                                'observation start time (Modified Julian Date)')
        hdu.header['EXPTIME'] = (hdul[0].header['EXPTIME'], 
                                 'exposure duration (seconds)--calculated')
    else:
        hdu.header['DATE-OBS'] = (hdul[1].header['DATE-OBS'], 
                                  'UT date of start of observation(yyyy-mm-dd)')
        hdu.header['EXPSTART'] = (hdul[1].header['EXPSTART'], 
                                'observation start time (Modified Julian Date)')
        hdu.header['EXPTIME'] = (hdul[1].header['EXPTIME'], 
                                 'exposure duration (seconds)--calculated')
    
    
    hdu.header['ASN_ID'] = (hdul[0].header['ASN_ID'], 
                            'unique identifier assigned to association')
    
    return hdu

def createCoAddPrimaryHeader(headers, wave, exptime):
    ''' Calculates confidence limits according to Kraft et al. 1991.
    
    Parameters:
    ----------
    hdul: hdu list from the x1d file 
    
    wave: wavelength array of the spectroscopic data
    
    exptime: exptime array of the spectroscopic data
    
    Returns:
    ---------
    hdu: custom primary header for the final co-add
    
    '''
    
    unique_pids = []
    unique_cenwaves = []
    unique_optelem = []
    for hdr in headers:
        if hdr['PROPOSID'] not in unique_pids:
            unique_pids.append(hdr['PROPOSID'])
        if hdr['CENWAVE'] not in unique_cenwaves:
            unique_cenwaves.append(hdr['CENWAVE'])
        if hdr['OPT_ELEM'] not in unique_optelem:
            unique_optelem.append(hdr['OPT_ELEM'])
    unique_pids.sort()
    
    hdu = fits.PrimaryHDU()
    hdu.header['DATE'] = (datetime.utcnow().isoformat(timespec='seconds'),\
                            "File creation date")
    
    hdu.header['FILETYPE'] = (headers[0]['FILETYPE'], 
                              'type of data found in the data file')
    hdu.header['TELESCOP'] = (headers[0]['TELESCOP'],
                              'telescope used to acquire data')
    hdu.header['INSTRUME'] = (headers[0]['INSTRUME'], 
                              'identifier for instrument used to acquire data')
    hdu.header['EQUINOX'] = (headers[0]['EQUINOX'], 
                             'equinox of celestial coord. system')
    
    hdu.header['TARGNAME'] = (headers[0]['TARGNAME'], 
                              'proposer\'s target name')
    hdu.header['RA_TARG'] = (headers[0]['RA_TARG'], 
                             'right ascention of the target (deg)')
    hdu.header['DEC_TARG'] = (headers[0]['DEC_TARG'], 
                              'declination of the target (deg)')
    
    
    hdu.header['PROPOSID'] = (unique_pids.pop(0), 
                              'PEP proposal identifier')
    for i, pid in enumerate(unique_pids):
        hdu.header['OTH_PID'+str(i)] = (pid, 
                                        'other PEP proposal identifier')
    
    hdu.header['OPUS_VER'] = (headers[0]['OPUS_VER'], 
                              'data processing software system version')
    hdu.header['CSYS_VER'] = (headers[0]['CSYS_VER'], 
                              'calibration software system version id')
    hdu.header['CAL_VER'] = (headers[0]['CAL_VER'], 
                             'CALCOS code version')
    hdu.header['FCOS_VER'] = (version, 
                              'FaintCOS code version')
    
    hdu.header['DETECTOR'] = (headers[0]['DETECTOR'],
                              'FUV or NUV')
    
    
    hdu.header['CENWAVE'] = (unique_cenwaves.pop(0), 
                             'central wavelength of the spectrum')
    for i, cenwave in enumerate(unique_cenwaves):
        hdu.header['CENWAVE'+str(i)] = (cenwave, 
                                    'other central wavelength of the spectrum')
        
    segments = [hdr['SEGMENT'] for hdr in headers]
    segm = ""
    if 'FUVA' in segments and 'FUVB' in segments:
        segm = 'BOTH'
    else:
        segm = segments[0]
        
    hdu.header['SEGMENT'] = (segm, 
                            'FUV detector segment name (FUVA, FUVB or BOTH)')
    if 'G130M' in unique_optelem and 'G160M' in unique_optelem:
        hdu.header['OPT_ELEM'] = ('FUVM', 
                                'optical element in use')
    else:
        hdu.header['OPT_ELEM'] = (unique_optelem[0], 
                                'optical element in use')
    
    
    hdu.header['BANDWID'] = (max(wave)-min(wave), 
                             'bandwidth of the data')
    hdu.header['SPECRES'] = (min([hdr['SPECRES'] for hdr in headers]), 
                             'smallest SPECRES at CENWAVE from all data sets')
    hdu.header['CENTRWV'] = ((max(wave)+min(wave))/2.0, 
                             'central wavelength of the data')
    hdu.header['MINWAVE'] = (min(wave), 
                             'minimum wavelength in spectrum')
    hdu.header['MAXWAVE'] = (max(wave), 
                             'maximum wavelength in spectrum')
    
    asn_ids = []
    expstart = []
    date_obs = []
    for hdr in headers:
        if hdr['ASN_ID'] not in asn_ids:
            asn_ids.append(hdr['ASN_ID'])
            expstart.append(hdr['EXPSTART'])
            date_obs.append(hdr['DATE-OBS'])
    sorted_dates = [x for _,x in sorted(zip(expstart, date_obs))]
    hdu.header['DATE-OBS'] = (sorted_dates.pop(0), 
                            'UT observation date for the first data set')
    for i, d in enumerate(sorted_dates):
        hdu.header['DATEOB'+str(i)] = (d,
                                    'UT observation date for other data sets')
        
    hdu.header['EXPTIME'] = (int(max(exptime)), 
                            'max exposure duration (seconds)--calculated')
    hdu.header['ETIMEMED'] = (int(np.median(exptime)), 
                            'median exposure duration (seconds)--calculated')
    
    return hdu
    

if __name__ == "__main__": 
    
    path_sci = "."
    
    # path to science was given in the command line?
    if len(sys.argv) > 1:
        path_sci = sys.argv[1]
    
    # find all corrtag files in the directory
    path_corrtag = [f for f in os.listdir(path_sci) if "corrtag" in f]
    if path_sci != ".":
        path_corrtag = [path_sci + s for s in path_corrtag]
    else:
        path_sci = ""
    
    # find reference directory
    try:
        path_ref = os.environ['lref']
    except KeyError:
        print("ERROR: lref is not defined!")
        sys.exit() 
    
    # find darkframes directory
    try:
        path_dark = os.environ['ldark']
    except KeyError:
        print("ERROR: ldark is not defined!")
        sys.exit()   
    
    print("FaintCOSv"+version, flush=True)
    
    # load all darkframes in the 'ldark' directory
    print("Loading darkframes...", end=" ", flush=True)
    path_darkframes = [f for f in os.listdir(path_dark) if "corrtag" in f]
    path_darkframes = [path_dark + s for s in path_darkframes]
    dark_file = []
    dark_expstart = []
    dark_segment = []
    dark_voltage = []
    #darkframes = []
    for d in path_darkframes:
        dark_tmp = fits.open(d)
        dark_file.append(dark_tmp[0].header['FILENAME'])
        dark_expstart.append(dark_tmp[1].header['EXPSTART'])
        dark_segment.append(dark_tmp[0].header['SEGMENT'])
        dark_voltage.append(dark_tmp[1].header["HVLEVEL" + \
                                               dark_segment[-1].split("FUV")[1]])
        dark_tmp.close()
        #darkframes.append(fits.open(d))
    darkframes = Table([dark_file, path_darkframes, \
                        dark_expstart, dark_segment, dark_voltage],\
                        names = ("FILE", "PATH", "EXPSTART", \
                                 "SEGMENT", "VOLTAGE"))
    print("OK")
    print(str(len(darkframes['FILE'])) + " darkframes have been found!")
    # find all datasets in the working directory and sort them
    datasets = Table(names=("ASN_ID", "TARGET", "OBS-DATE", "EXP_TIME", \
                            "OPT_ELEM", "CENTRWV", "FUVA", "FUVB"),\
                 dtype=('S10', 'S25', 'S23', 'i4', 'S5', 'i4', 'i4', 'i4'))
    corrtags = Table(names=("TARGET", "ASN_ID", "CORRTAG_FILE", "START_TIME",\
                                      "OPT_ELEM", "SEGMENT", "CENTRWV",\
                                      "EXP_TIME"), \
                               dtype=('S25', 'S10', 'S30', 'S23', 'S5' , 'S5',\
                                      'i4', 'f4'))
    # load corrtags files in a Table
    for f in path_corrtag:
        h0 = fits.open(f)[0].header
        h1 = fits.open(f)[1].header
        if h1['EXPFLAG'] != 'NORMAL':
            continue
        corrtags.add_row([h0['TARGNAME'], h0['ASN_ID'], f.split('/')[-1], \
                                    h1['DATE-OBS'] + " " + h1['TIME-OBS'],\
                                    h0['OPT_ELEM'], h0['SEGMENT'], \
                                    h0['CENTRWV'], h1['EXPTIME']])
    corrtags.sort(["TARGET", "START_TIME", 'ASN_ID', 'SEGMENT', 'CENTRWV'])
    print("Valid 'corrtag' files in the working directory:")
    corrtags.pprint()
    print("##################################################################")
    print("\n")

    unique_datasets = np.unique(np.array(corrtags['ASN_ID']))

    for asn_id in unique_datasets:
        visit = corrtags[corrtags['ASN_ID'] == asn_id]
        hdul = fits.open(path_sci + visit[0]['CORRTAG_FILE'])
        target = hdul[0].header['TARGNAME']
        obs_date = hdul[1].header['DATE-OBS']
        exp_time_a = np.array(visit[visit['SEGMENT'] == 'FUVA']['EXP_TIME']) 
        exp_time_b = np.array(visit[visit['SEGMENT'] == 'FUVB']['EXP_TIME']) 
        exp_time = max(np.array([np.sum(exp_time_a), np.sum(exp_time_b)]))
        num_fuva = len(np.array(visit[visit['SEGMENT'] == 'FUVA']))
        num_fuvb = len(np.array(visit[visit['SEGMENT'] == 'FUVB']))
        centrwv = visit[0]['CENTRWV']   
        opt_elem = visit[0]['OPT_ELEM']
        datasets.add_row([asn_id, target, obs_date, exp_time, \
                          opt_elem, centrwv, num_fuva, num_fuvb])
    datasets.sort(['TARGET','OBS-DATE'])
    print("Valid datasets in the working directory:")
    datasets.pprint()
    print("###################################################################")
    print("\n")
    
    a = input("Do you wish to proceed?(y/n)")
    if a != 'y' and a != 'Y':
        print("Canceled by user!")
        sys.exit()
    
    # determine background for every corrtag

    corr_files = [path_sci + s for s in list(corrtags['CORRTAG_FILE'])] 
    for c in corr_files:
        corr_hdul = fits.open(c)
        corr_data = Table(corr_hdul[1].data)
        corr_prefix = corr_hdul[0].header['FILENAME'].split('_')[0]
        print("Working on " + str(c.split("/")[-1]))
        segm = corr_hdul[0].header['SEGMENT']
        opt_elem = corr_hdul[0].header['OPT_ELEM']
        cenwave = corr_hdul[0].header['CENTRWV']
        voltage = corr_hdul[1].header["HVLEVEL" + segm.split("FUV")[1]]
        exp_start = corr_hdul[1].header['EXPSTART']
        
        # open xtractab to get apertures
        xtractab = Table(fits.getdata(path_ref + \
                   corr_hdul[0].header['XTRACTAB'].split('$')[-1]))
        xtractab = xtractab[(xtractab['SEGMENT'] == segm) &\
                            (xtractab['OPT_ELEM'] == opt_elem) &\
                            (xtractab['CENWAVE'] == cenwave) &\
                            (xtractab['APERTURE'] == 'PSA')]
        ap_spec = [float(xtractab['B_SPEC']) - float(xtractab['HEIGHT'])/2.,\
                   float(xtractab['B_SPEC']) + float(xtractab['HEIGHT'])/2.]
        if 'B_HGT1' in xtractab.columns:
            ap_bkg = [float(xtractab['B_BKG1']) - float(xtractab['B_HGT1'])/2.,\
                  float(xtractab['B_BKG1']) + float(xtractab['B_HGT1'])/2.,\
                  float(xtractab['B_BKG2']) - float(xtractab['B_HGT2'])/2.,\
                  float(xtractab['B_BKG2']) + float(xtractab['B_HGT2'])/2.]
        else:
            ap_bkg = [float(xtractab['B_BKG1']) - float(xtractab['BHEIGHT'])/2.,\
                  float(xtractab['B_BKG1']) + float(xtractab['BHEIGHT'])/2.,\
                  float(xtractab['B_BKG2']) - float(xtractab['BHEIGHT'])/2.,\
                  float(xtractab['B_BKG2']) + float(xtractab['BHEIGHT'])/2.]
        
        
        # open _pha file to get pha limits
        pha_file = path_ref + corr_hdul[0].header['PHATAB'].split('$')[-1]
        pha_limits = Table(fits.open(pha_file)[1].data)
        pha_max = pha_limits[(pha_limits['OPT_ELEM'] == opt_elem) & \
                             (pha_limits['SEGMENT'] == segm)]['ULT']
        pha_min = pha_limits[(pha_limits['OPT_ELEM'] == opt_elem) & \
                             (pha_limits['SEGMENT'] == segm)]['LLT']
                             
        # open the corresponding 1dx file
        path_x1d = path_sci + corr_prefix + "_x1d.fits"
        hdul_x1d = fits.open(path_x1d)
        data_x1d = Table(hdul_x1d[1].data)
        data_x1d = data_x1d[data_x1d['SEGMENT'] == segm]
                             
        # select only valid darkframes for this corrtag
        val_darks = []
        sel_darks = darkframes[(darkframes['SEGMENT'] == segm) & \
                               (darkframes['VOLTAGE'] == voltage) & \
            (darkframes['EXPSTART'] > exp_start - DARK_EXPSTART_INTERVAL) & \
            (darkframes['EXPSTART'] < exp_start + DARK_EXPSTART_INTERVAL)]
        
        
        for d in sel_darks['PATH']:
            val_darks.append(fits.open(d))
        if len(val_darks) >= MIN_DARKS:
            print(str(len(val_darks)) + " valid darkframes.")
        else:
            print("Not enough darkframes for this _corrtag!") 
            sys.exit()  
        
        # find the dispersion function for xdopp(wavelength)
        wl = np.array(corr_data[(corr_data['WAVELENGTH'] > 1) & \
                                (corr_data['YFULL'] > ap_spec[0]) & \
                                (corr_data['YFULL'] < ap_spec[1]) & \
                                (corr_data['DQ'] == 0)]['WAVELENGTH'])
        xdopp = np.array(corr_data[(corr_data['WAVELENGTH'] > 1) & \
                                   (corr_data['YFULL'] > ap_spec[0]) & \
                                   (corr_data['YFULL'] < ap_spec[1]) & \
                                   (corr_data['DQ'] == 0)]['XDOPP'])
        xdopp_disp = linregress(wl, xdopp)[0:2]

        # find bad regions in xdopp coordinates
        if opt_elem == 'G140L':
            BAD_REGIONS = BAD_REGIONS_G140L
        else:
            BAD_REGIONS = BAD_REGIONS_G130M
            
        bad_reg_darks = []
        for reg in BAD_REGIONS:
            bad_reg_darks.append(np.array(reg)*xdopp_disp[0] + xdopp_disp[1])
        
        # remove bad regions from corrtag
        for reg in BAD_REGIONS:
            corr_data = corr_data[(corr_data['WAVELENGTH'] < reg[0]) | \
                                  (corr_data['WAVELENGTH'] > reg[1])]
            
        # remove counts outside of the detectors active area
        corr_data = corr_data[(corr_data['XDOPP'] < 15000) & \
                                (corr_data['XDOPP'] > 1500)]
        
        # produce a cumulative puls-height distribution for the corrtag
        # select background windows in the corrtag data
        corr_bkg1 = corr_data[(corr_data['YFULL'] < ap_bkg[1]) & \
                              (corr_data['YFULL'] > ap_bkg[0]) & \
                              ((corr_data['DQ'] == 0) | \
                               (corr_data['DQ'] == 512) | \
                               (corr_data['DQ'] == 8192) | \
                               (corr_data['DQ'] == 8704))]
        corr_bkg2 = corr_data[(corr_data['YFULL'] < ap_bkg[3]) & \
                              (corr_data['YFULL'] > ap_bkg[2]) & \
                              ((corr_data['DQ'] == 0) | \
                               (corr_data['DQ'] == 512) | \
                               (corr_data['DQ'] == 8192) | \
                               (corr_data['DQ'] == 8704))]
        corr_bkg = np.append(np.array(corr_bkg1['PHA']), \
                             np.array(corr_bkg2['PHA']))
            
        # cumulative puls-height distribution
        unique, counts = np.unique(corr_bkg, return_counts=True)
        corr_pha_dist = np.zeros(shape=32, dtype=np.int32)
        for i in range(len(unique)):
            corr_pha_dist[unique[i]] = counts[i]
        corr_data_cumsum = np.cumsum(corr_pha_dist)
        corr_max_counts = corr_data_cumsum[-1]
        corr_data_cumsum = corr_data_cumsum / corr_max_counts
        
        # produce cumulative puls_height distribution for every darkframe
        # and compare it to the puls-height distribution of the corrtag
        # with a Kolmogorov–Smirnov test
        KS_values = np.zeros(shape=(len(val_darks)), dtype=np.float32)
        dark_max_counts = np.zeros(shape=(len(val_darks)), dtype=np.float32)
        d = 0
        for i in range(len(val_darks)):
            dark_data = Table(val_darks[i][1].data)
            
            # remove bad regions from the darkframe
            for reg in bad_reg_darks:
               dark_data = dark_data[(dark_data['XDOPP'] < reg[0]) | \
                                     (dark_data['XDOPP'] > reg[1])]
               
            # remove counts outside of the detectors active area
            dark_data = dark_data[(dark_data['XDOPP'] < 15000) & \
                                  (dark_data['XDOPP'] > 1500)]
            
            dark_bkg1 = dark_data[(dark_data['YFULL'] < ap_bkg[1]) & \
                                  (dark_data['YFULL'] > ap_bkg[0]) & \
                                  ((dark_data['DQ'] == 0) | \
                                   (dark_data['DQ'] == 8192))]
            dark_bkg2 = dark_data[(dark_data['YFULL'] < ap_bkg[3]) & \
                                  (dark_data['YFULL'] > ap_bkg[2]) & \
                                  ((dark_data['DQ'] == 0) | \
                                   (dark_data['DQ'] == 8192))]
            dark_bkg = np.append(np.array(dark_bkg1['PHA']), \
                                 np.array(dark_bkg2['PHA']))
        
            # cumulative puls-height distribution
            unique, counts = np.unique(dark_bkg, return_counts=True)
            dark_pha_dist = np.zeros(shape=32, dtype=np.int32)
            for j in range(len(unique)):
                dark_pha_dist[unique[j]] = counts[j]
            dark_data_cumsum = np.cumsum(dark_pha_dist)
            dark_max_counts[i] = dark_data_cumsum[-1]
            dark_data_cumsum = dark_data_cumsum / dark_max_counts[i]
            KS_values[i] = max(abs(dark_data_cumsum - corr_data_cumsum))
        
        # sort K-S statistics
        dark_sorted_KS = KS_values.argsort()
        KS_test = KS_THRESHOLD
        ind_threshold = 0
        while(True):
            for di in range(len(dark_sorted_KS)):
                if KS_values[dark_sorted_KS[di]] < KS_test:
                    ind_threshold = di
            if ind_threshold + 1 < MIN_DARKS:
                KS_test = KS_test + KS_STEP
            else:
                break
        number_of_darks = ind_threshold + 1
        

        dark_total_counts = 0
        
        # save prefixes of the used darks
        dark_prefixes = val_darks[dark_sorted_KS[0]][0].\
                        header['FILENAME'].split("_")[0]
        for i in range(1, number_of_darks):
            dark_prefixes = dark_prefixes + ", " + \
                            val_darks[dark_sorted_KS[i]][0].\
                            header['FILENAME'].split("_")[0]

        # calculate total counts of the used darks
        for i in range(number_of_darks):
            dark_total_counts = dark_total_counts + \
                                dark_max_counts[dark_sorted_KS[i]]
        print("Number of used darks: " + str(number_of_darks))
        print("KS threshold: " + str(KS_test))
        #print("Prefixes: " + dark_prefixes)
        
        # calculate scaling factor between the corrtag and combined darks
        scaling_factor = corr_max_counts / dark_total_counts
        scaling_factor_err = math.sqrt(scaling_factor * \
                            (1. + 1. / dark_total_counts) / dark_total_counts)
        
        print("Average scaling factor: " + str(scaling_factor*number_of_darks))
        
        
        
        # coadd best darkframes
        darks_combined = Table(val_darks[dark_sorted_KS[0]][1].data)
        for i in range(1, number_of_darks):
            darks_combined = vstack([darks_combined, \
                             Table(val_darks[dark_sorted_KS[i]][1].data)])
            
        # extract the spectral window from the combined darks
        dark_psa = darks_combined[(darks_combined['YFULL'] >= ap_spec[0]) &\
                                  (darks_combined['YFULL'] <= ap_spec[1]) & \
                                  (darks_combined['PHA'] >= pha_min) & \
                                  (darks_combined['PHA'] <= pha_max)]
        dark_psa_xfull = np.array(dark_psa['XFULL'])
        
        # shift xfull according to XSHIFT
        xshift = corr_hdul[1].header['SHIFT1' + segm.split("FUV")[1]]
        dark_psa_xshift = dark_psa_xfull - xshift
        
        # collapse the 2d spectrum to 1d
        dark_hist, tmp = np.histogram(dark_psa_xshift, 16384, range=(0, 16384))
        
        dark_hist_mean = np.zeros(shape=16384, dtype=np.float32)
        dark_hist_error = np.zeros(shape=16384, dtype=np.float32)        
        
        science_dq = np.array(data_x1d['DQ'])[0]
        dark_hist_masked = np.ma.masked_where(science_dq != 0, dark_hist)
        unmasked_hist_indices = np.ma.masked_where(science_dq != 0, \
                                np.arange(len(dark_hist_masked))).compressed()
        half_width = int((BKG_AV - 1) / 2)
        first_mean_value = unmasked_hist_indices[0] + half_width
        
        # shift first value, if dq == 0
        last_mean_value = unmasked_hist_indices[-1] - half_width

        # running average between first and last valid points
        for n in range(first_mean_value, last_mean_value + 1):
            win = dark_hist_masked[n - half_width:n + half_width].compressed()
            if len(win) > 0:
                dark_hist_mean[n] = np.sum(win) / len(win)
                dark_hist_error[n] = np.sqrt(np.sum(win)) / len(win)
            else:
                dark_hist_mean[n] = 0.0
                dark_hist_error[n] = 0.0
        
        # extrapolate first and last values to the left and right edges
        for n in range(half_width):
            dark_hist_mean[first_mean_value - half_width + n] = \
            dark_hist_mean[first_mean_value]
            dark_hist_error[first_mean_value - half_width + n] = \
            dark_hist_error[first_mean_value]
            dark_hist_mean[last_mean_value + n] = \
            dark_hist_mean[last_mean_value]
            dark_hist_error[last_mean_value + n] = \
            dark_hist_error[last_mean_value]
        
        # rescale to the corrtag file
        dark_hist_mean = scaling_factor * dark_hist_mean
        
        
        # error propagation due to the scaling factor
        dark_hist_error = np.sqrt(np.power(dark_hist_mean, 2) * \
                                  scaling_factor_err * scaling_factor_err + \
                                  scaling_factor * scaling_factor * \
                                  np.power(dark_hist_error, 2))
        
        calib = np.divide(np.array(data_x1d['NET'])[0], \
                          np.array(data_x1d['FLUX'])[0], \
                          out=np.full(shape=(16384), fill_value=np.nan), \
                          where=np.array(data_x1d['FLUX'])[0] != 0)
        wave = np.array(data_x1d['WAVELENGTH'])[0]
        wave_inter = np.delete(wave, np.where(np.isnan(calib)))
        calib_inter = np.delete(calib, np.where(np.isnan(calib)))
        interp_calib = interp1d(wave_inter, calib_inter, \
                               kind='quadratic', \
                               fill_value=0, \
                               bounds_error=False)
        calib = interp_calib(wave)
        
        
        # include hot spots and low responce regions
        dq = np.array(data_x1d['DQ'])[0]
        dq_wgt = np.zeros(shape=len(dq), dtype=np.int32)
        dq[dq == 1024] = 2
        for q in range(len(dq_wgt)):
            if dq[q] == 0 or dq[q] == 2:
                dq_wgt[q] = 1
            else:
                dq_wgt[q] = 0
        
        gross = np.array(data_x1d['GROSS'][0])
        net = np.array(data_x1d['NET'][0])
        used_px = np.where((gross > 0) & (net > 0))[0]
        valid_gross = gross[used_px]
        valid_net = net[used_px]
        valid_dq = dq[used_px]
        valid_wave = wave[used_px]
        flt_raw = valid_gross/valid_net
        interp_flt = interp1d(valid_wave[np.where(valid_dq == 0)[0]], \
                              flt_raw[np.where(valid_dq == 0)[0]], \
                              kind='linear', \
                              fill_value=1.0, \
                              bounds_error=False)
        flt_curve = interp_flt(wave)

        # calculate contamination by the geocoronal Lya on G140L
        # works only for cenwave 800 and 1105 (see Worseck et al. 2016)
        # otherwise set contamination 0
        
        if opt_elem == 'G140L' and (cenwave == 800 or cenwave == 1105):
            bkg_lya, bkg_lya_err_up, bkg_lya_err_down = \
            calc_lya_scatter_model(data_x1d['GCOUNTS'][0], \
                                   data_x1d['WAVELENGTH'][0])
        
        else:
            bkg_lya = np.zeros(shape = len(gross), dtype = np.float32)
            bkg_lya_err_up = np.zeros(shape = len(gross), dtype = np.float32)
            bkg_lya_err_down = np.zeros(shape = len(gross), dtype = np.float32)
        
        # create a new file to store the results
        
        col1 = fits.Column(name='WAVELENGTH', \
                           format='D', \
                           unit='angstrom',\
                           array=np.array(data_x1d['WAVELENGTH'])[0])
        col2 = fits.Column(name='GCOUNTS', \
                           format='I', \
                           unit='count',\
                           array=np.rint(data_x1d['GCOUNTS'])[0])
        col3 = fits.Column(name='EXPTIME', \
                           format='D', \
                           unit='seconds',\
                           array=np.full(shape = 16384, \
                                         fill_value = data_x1d['EXPTIME']))
        col4 = fits.Column(name='DQ', \
                           format='I', \
                           unit='',\
                           array=dq)
        col5 = fits.Column(name='DQ_WGT', \
                           format='I', \
                           unit='',\
                           array=dq_wgt)
        col6 = fits.Column(name='DARK_CURRENT', \
                           format='D', \
                           unit='counts',\
                           array=dark_hist_mean)
        col7 = fits.Column(name='DARK_CURRENT_ERR', \
                           format='D', \
                           unit='counts',\
                           array=dark_hist_error)
        col8 = fits.Column(name='CALIB', \
                           format='D', \
                           unit='counts cm**2 A / erg',\
                           array=calib)
        col9 = fits.Column(name='FLAT_CORR',\
                           format='D',\
                           unit='',\
                           array = flt_curve)
        col10 = fits.Column(name='LYA_SCATTER', \
                           format='D', \
                           unit='counts',\
                           array=bkg_lya)
        col11 = fits.Column(name='LYA_SCATTER_ERR_UP', \
                           format='D', \
                           unit='counts',\
                           array=bkg_lya_err_up)
        col12 = fits.Column(name='LYA_SCATTER_ERR_DOWN', \
                           format='D', \
                           unit='counts',\
                           array=bkg_lya_err_down)
        
        hdul_x1d[0].header['SEGMENT'] = segm
        hdu = createCDRPrimaryHeader(hdul_x1d, wave)
        
        hdu_binary = fits.BinTableHDU.from_columns([col1, col2, col3, col4,\
                                                    col5, col6, col7, col8,\
                                                    col9, col10, col11, col12])
        hdul = fits.HDUList([hdu, hdu_binary])
        saved_file = path_sci + corr_prefix + "_cdr_" + segm + ".fits"
        hdul.writeto(saved_file, overwrite=True)
        print("Results were written into " + saved_file)
        print("\n")
        # close all darkframes
        for d in val_darks:
            d.close()

    if path_sci != "":
        cdr_files = [f for f in os.listdir(path_sci) if "cdr" in f]
        cdr_files = [path_sci + s for s in cdr_files]
    else:
        cdr_files = [f for f in os.listdir(".") if "cdr" in f]


        
    # coadding exposures for ever dataset
    print("Coadding exposures for every dataset.")
    hdul_ar = []
    visit_data = []
    visit_hdu = []
    for f in cdr_files:
        hdul_ar.append(fits.open(f))
    for asn in unique_datasets:
        targname = ""
        segm = ""
        current_hdul = hdul_ar[0]
        seg_fuva = []
        seg_fuvb = []
        for hdul in hdul_ar:
            if hdul[0].header['ASN_ID'] == asn.decode("utf-8"):
                targname = hdul[0].header['TARGNAME']
                current_hdul = hdul
                if hdul[0].header['SEGMENT'] == 'FUVA':
                    seg_fuva.append(hdul)
                else:
                    seg_fuvb.append(hdul)
        exposures = []
        if len(seg_fuvb) > 0:
            exposures.append(seg_fuvb)
            segm = 'FUVA'
        if len(seg_fuva) > 0:
            exposures.append(seg_fuva)
            segm = 'FUVB'
        if len(seg_fuva) > 0 and len(seg_fuvb) > 0:
            segm = 'BOTH'
        coadded_tab = []
        for expos in exposures:
            if len(expos) > 0:
                wavelength = np.array(Table(expos[0][1].data)['WAVELENGTH'])
                totaltime = np.zeros(shape=len(wavelength), dtype=np.float32)
                total_dq = np.full(shape=len(wavelength), \
                                   fill_value=20000, \
                                   dtype=np.int32)
                total_dq_wgt = np.zeros(shape=len(wavelength), dtype=np.int32)
                totalcounts = np.zeros(shape=len(wavelength), dtype=np.int32)
                totaldark = np.zeros(shape=len(wavelength), dtype=np.float32)
                totaldark_error = np.zeros(shape=len(wavelength), \
                                           dtype=np.float32)
                total_lya = np.zeros(shape=len(wavelength), dtype=np.float32)
                total_lya_err_up = np.zeros(shape=len(wavelength), \
                                           dtype=np.float32)
                total_lya_err_down = np.zeros(shape=len(wavelength), \
                                           dtype=np.float32)
                
                total_flux_calib = np.zeros(shape=len(wavelength), \
                                            dtype=np.float32)
                total_flt = np.array(Table(expos[0][1].data)['FLAT_CORR'])

                for i in range(len(expos)):
                    exp_data = Table(expos[i][1].data)
                    exp_time = np.array(exp_data['EXPTIME'])
                    dq = np.array(exp_data['DQ'])
                    dq_wgt =  np.array(exp_data['DQ_WGT'])

                    # take the maximum of DQ_WGT of the subexposures
                    total_dq_wgt = np.maximum(total_dq_wgt, dq_wgt)
                    totaltime = totaltime + dq_wgt * exp_time
                    # take the minimum of DQ of the subexposures
                    total_dq = np.minimum(dq, total_dq)
                    
                    gcounts = np.array(exp_data['GCOUNTS'])
                    totalcounts = totalcounts + dq_wgt * gcounts
                    
                    darkcurrent = np.array(exp_data['DARK_CURRENT'])
                    totaldark = np.add(totaldark, dq_wgt * darkcurrent)
                           
                    darkcurrent_error = np.array(exp_data['DARK_CURRENT_ERR'])
                    totaldark_error = \
                    np.sqrt(np.add(np.power(totaldark_error, 2), \
                                   np.power(dq_wgt * darkcurrent_error, 2)))
                    
                    lya_scatter = np.array(exp_data['LYA_SCATTER'])
                    total_lya = np.add(total_lya, dq_wgt * lya_scatter)
                    ''''       
                    lya_scatter_err_up = np.array(exp_data['LYA_SCATTER_ERR_UP'])
                    total_lya_err_up = \
                    np.sqrt(np.add(np.power(total_lya_err_up, 2), \
                                   np.power(dq_wgt * lya_scatter_err_up, 2)))
                    
                    lya_scatter_err_down = np.array(exp_data['LYA_SCATTER_ERR_DOWN'])
                    total_lya_err_down = \
                    np.sqrt(np.add(np.power(total_lya_err_down, 2), \
                                   np.power(dq_wgt * lya_scatter_err_down, 2)))
                    '''
                    lya_scatter_err_up = np.array(exp_data['LYA_SCATTER_ERR_UP'])
                    total_lya_err_up = total_lya_err_up + dq_wgt * lya_scatter_err_up
                    
                    lya_scatter_err_down = np.array(exp_data['LYA_SCATTER_ERR_DOWN'])
                    total_lya_err_down = total_lya_err_down + dq_wgt * lya_scatter_err_down
                    
                    np.seterr(divide='ignore')
                    flux_calib = np.array(exp_data['CALIB'])
                    for  n in range(len(total_flux_calib)):
                        if (total_flux_calib[n] == 0) | \
                           (np.isnan(total_flux_calib[n])):
                            total_flux_calib[n] = flux_calib[n]
                calib = total_flux_calib
                flux = np.divide((totalcounts - totaldark - total_lya), \
                                 (calib * totaltime * total_flt), \
                                 out=np.zeros_like(totaldark), \
                                 where=((calib != 0) & (totaltime != 0)))
                tdata = Table([wavelength, flux, totalcounts, totaldark, \
                               totaldark_error, total_dq, total_dq_wgt, \
                               totaltime, calib, total_flt, total_lya, \
                               total_lya_err_up, total_lya_err_down], \
                              names=("WAVELENGTH", "FLUX", "GCOUNTS", \
                                     "DARK_CURRENT", "DARK_CURRENT_ERR", "DQ", \
                                     "DQ_WGT", "EXPTIME", "CALIB",\
                                     "FLAT_CORR", "LYA_SCATTER", \
                                     "LYA_SCATTER_ERR_UP", "LYA_SCATTER_ERR_DOWN"))
                
                coadded_tab.append(tdata)
        
        if len(coadded_tab) > 1:
            max_wl_fuvb = max(np.array(coadded_tab[0]['WAVELENGTH']))
            min_wl_fuva = min(np.array(coadded_tab[1]['WAVELENGTH']))
            seg_overlap = max_wl_fuvb - min_wl_fuva
            seg_cutoff = seg_overlap/2.
            max_wl_fuvb = max_wl_fuvb - seg_cutoff
            min_wl_fuva = min_wl_fuva + seg_cutoff
            coadded_data = \
            vstack([coadded_tab[0][coadded_tab[0]['WAVELENGTH'] \
                                   < max_wl_fuvb],\
                    coadded_tab[1][coadded_tab[1]['WAVELENGTH'] \
                                   > min_wl_fuva]])
        else:
            coadded_data = coadded_tab[0]
        
        current_hdul[0].header['SEGMENT'] = segm
        hdu = createCDRPrimaryHeader(current_hdul, coadded_data['WAVELENGTH'])
        hdu.header['EXPTIME'] = max(coadded_data['EXPTIME'])
        hdu.header['COMMENT'] = "Coadded spectra for a single dataset " + \
                                "with improved calibration"
                            
        visit_data.append(coadded_data)
        visit_hdu.append(hdu)

        binary_hdu = fits.BinTableHDU(coadded_data)
        hdul = fits.HDUList([hdu, binary_hdu])
        saved_file = path_sci + asn.decode("utf-8") + "_dataset_sum.fits"
        hdul.writeto(saved_file, overwrite=True)
        print(asn.decode("utf-8") + " is complete." + \
              " The spectrum is stored under " + saved_file)
    
    # Bin all _cdr_raw.fits files with the binsize = BIN_PX
    if BIN_DATASET:
        print("\n")
        print("Binning every dataset with " + str(BIN_PX) + " pixels.")
        for d in range(len(visit_data)):
            
            asn = visit_hdu[d].header['ASN_ID']
            targname = visit_hdu[d].header['TARGNAME']
            
            cdr_wave = np.array(visit_data[d]['WAVELENGTH'])
            cdr_gcounts = np.array(visit_data[d]['GCOUNTS'])
            cdr_dc = np.array(visit_data[d]['DARK_CURRENT'])
            cdr_dc_err = np.array(visit_data[d]['DARK_CURRENT_ERR'])
            cdr_lya = np.array(visit_data[d]['LYA_SCATTER'])
            cdr_lya_err_up = np.array(visit_data[d]['LYA_SCATTER_ERR_UP'])
            cdr_lya_err_down = np.array(visit_data[d]['LYA_SCATTER_ERR_DOWN'])
            cdr_dq = np.array(visit_data[d]['DQ'])
            cdr_dq_wgt = np.array(visit_data[d]['DQ_WGT'])
            cdr_exptime = np.array(visit_data[d]['EXPTIME'])
            cdr_calib = np.array(visit_data[d]['CALIB'])
            cdr_flt = np.array(visit_data[d]['FLAT_CORR'])
            
            binned_wave = np.ndarray(shape=int(len(cdr_wave)/BIN_PX - 1), \
                                     dtype=np.float32)
            binned_gcounts = np.ndarray(shape=len(binned_wave), dtype=np.float32)
            binned_dc = np.ndarray(shape=(len(binned_wave)), dtype=np.float32)
            binned_dc_err = np.ndarray(shape=(len(binned_wave)), dtype=np.float32)
            binned_lya = np.ndarray(shape=(len(binned_wave)), dtype=np.float32)
            binned_lya_err_up = np.ndarray(shape=(len(binned_wave)), dtype=np.float32)
            binned_lya_err_down = np.ndarray(shape=(len(binned_wave)), dtype=np.float32)
            binned_dq = np.ndarray(shape=len(binned_wave), dtype=np.int32)
            binned_dq_wgt = np.ndarray(shape=len(binned_wave), dtype=np.int32)
            binned_exptime = np.ndarray(shape=len(binned_wave), dtype=np.float32)
            binned_calib = np.ndarray(shape=len(binned_wave), dtype=np.float32)
            binned_flt = np.ndarray(shape=len(binned_wave), dtype=np.float32)
            
            for i in range(len(binned_wave)):
                edge_1 = i*BIN_PX
                edge_2 = (i+1)*BIN_PX
                tmp_wave = cdr_wave[edge_1:edge_2]
                tmp_gcounts = cdr_gcounts[edge_1:edge_2]
                tmp_dc = cdr_dc[edge_1:edge_2]
                tmp_lya = cdr_lya[edge_1:edge_2]
                tmp_lya_err_up = cdr_lya_err_up[edge_1:edge_2]
                tmp_lya_err_down = cdr_lya_err_down[edge_1:edge_2]
                tmp_dq = cdr_dq[edge_1:edge_2]
                tmp_dq_wgt = cdr_dq_wgt[edge_1:edge_2]
                tmp_exptime = cdr_exptime[edge_1:edge_2]
                tmp_calib = cdr_calib[edge_1:edge_2]
                tmp_flt = cdr_flt[edge_1:edge_2]
                tmp_dc_err = cdr_dc_err[edge_1:edge_2]


                binned_wave[i] = np.sum(tmp_wave)/float(BIN_PX)
                binned_gcounts[i] = np.sum(tmp_gcounts*tmp_dq_wgt)
                binned_dc[i] = np.sum(tmp_dc*tmp_dq_wgt)
                
                
                binned_dc_err[i] = np.sum(tmp_dc_err*tmp_dq_wgt)
                binned_lya[i] = np.sum(tmp_lya*tmp_dq_wgt)
                binned_lya_err_up[i] = np.sum(tmp_lya_err_up)
                binned_lya_err_down[i] = np.sum(tmp_lya_err_down)
                binned_dq[i] = min(tmp_dq)
                binned_dq_wgt[i] = max(tmp_dq_wgt)
                binned_exptime[i] = np.sum(tmp_exptime*tmp_dq_wgt)
                binned_calib[i] = np.sum(tmp_calib)/float(BIN_PX)
                binned_flt[i] = np.sum(tmp_flt)/float(BIN_PX)
            binned_flux = np.divide((binned_gcounts - binned_dc - binned_lya), \
                            (binned_calib * binned_exptime * binned_flt), \
                            out=np.zeros_like(binned_dc), \
                            where=((binned_calib != 0) & \
                                   (binned_exptime != 0) & \
                                   (binned_flt != 0)))
            binned_bkg_err_up = np.sqrt(np.power(binned_dc_err, 2.) + \
                                        np.power(binned_lya_err_up, 2.))
            binned_bkg_err_down = np.sqrt(np.power(binned_dc_err, 2.) + \
                                          np.power(binned_lya_err_down, 2.))
        
            if FELDMAN_COUSINS:
                cnt_err = calc_conf_lim_feldman_cousins(binned_gcounts,\
                                                        binned_dc + binned_lya)
            else:
                cnt_err = calc_conf_lim_kraft(binned_gcounts, binned_dc + binned_lya)
            
            
            cnt_err_down = cnt_err[0]
            cnt_err_up = cnt_err[1]
            flux_err_up = np.divide(cnt_err_up, \
                    (binned_calib * binned_exptime * binned_flt), \
                    out=np.zeros_like(binned_dc), \
                    where=((binned_calib != 0) & \
                           (binned_exptime != 0) & \
                           (binned_flt != 0)))
            flux_err_down = np.divide(cnt_err_down, \
                            (binned_calib * binned_exptime * binned_flt), \
                            out=np.zeros_like(binned_dc), \
                            where=((binned_calib != 0) & \
                                   (binned_exptime != 0) & \
                                   (binned_flt != 0)))
          
            col1 = fits.Column(name='WAVELENGTH', \
                               format='D', \
                               unit='angstrom',\
                               array=binned_wave)
            col2 = fits.Column(name='FLUX', \
                               format='D', \
                               unit='erg/s/cm^2/A',\
                               array=binned_flux)
            col3 = fits.Column(name='FLUX_ERR_UP', \
                               format='D', \
                               unit='erg/s/cm^2/A',\
                               array=flux_err_up)
            col4 = fits.Column(name='FLUX_ERR_DOWN', \
                               format='D', \
                               unit='erg/s/cm^2/A',\
                               array=flux_err_down)
            col5 = fits.Column(name='GCOUNTS', \
                               format='I', \
                               unit='count',\
                               array=np.rint(binned_gcounts))
            col6 = fits.Column(name='BACKGROUND', \
                               format='D', \
                               unit='counts',\
                               array=binned_dc + binned_lya)
            col7 = fits.Column(name='BKG_ERR_UP', \
                               format='D', \
                               unit='counts',\
                               array=binned_bkg_err_up)
            col8 = fits.Column(name='BKG_ERR_DOWN', \
                               format='D', \
                               unit='counts',\
                               array=binned_bkg_err_down)
            col9 = fits.Column(name='DARK_CURRENT', \
                               format='D', \
                               unit='counts',\
                               array=binned_dc)
            col10 = fits.Column(name='DARK_CURRENT_ERR', \
                               format='D', \
                               unit='counts',\
                               array=binned_dc_err)      
            col11 = fits.Column(name='EXPTIME', \
                               format='D', \
                               unit='seconds',\
                               array=binned_exptime)
            col12 = fits.Column(name='DQ', \
                               format='I', \
                               unit='',\
                               array=binned_dq)
            col13 = fits.Column(name='CALIB', \
                               format='D', \
                               unit='counts cm**2 A / erg',\
                               array=binned_calib)
            col14 = fits.Column(name='FLAT_CORR',\
                               format='D',\
                               unit='',\
                               array = binned_flt)
            col15 = fits.Column(name='LYA_SCATTER', \
                               format='D', \
                               unit='counts',\
                               array=binned_lya)
            col16 = fits.Column(name='LYA_SCATTER_ERR_UP', \
                               format='D', \
                               unit='counts',\
                               array=binned_lya_err_up)
            col17 = fits.Column(name='LYA_SCATTER_ERR_DOWN', \
                               format='D', \
                               unit='counts',\
                               array=binned_lya_err_down)
            
            #hdu = visit_hdu[d]
            hdu = createCDRPrimaryHeader([visit_hdu[d]], binned_wave)
            
            binned_hdu_binary = fits.BinTableHDU.from_columns([col1, col2, col3, \
                                                        col4, col5, col6, col7, \
                                                        col8, col9, col10, col11,\
                                                        col12, col13, col14, col15,\
                                                        col16, col17])
            hdul = fits.HDUList([hdu, binned_hdu_binary])
            saved_file = path_sci + targname + "_" + asn + "_" + \
                         str(BIN_PX) + "px_bin.fits"
            hdul.writeto(saved_file, overwrite=True)
            print(asn + " is binned and stored in " + saved_file)
            print("\n")

    # Co-adding routine for all datasets in the working folder
    # It works for different setups, but only for the same object and grating!!!
    if COADD_ALL_DATASETS:
        # Coadding all exposures    
        cdr_data = []
        cdr_hdu = []
        for f in cdr_files:
            tmp = fits.open(f)
            t = Table(tmp[1].data)
            h = tmp[0].header
            t = t[t['CALIB'] > 0]
            cdr_data.append(t)
            cdr_hdu.append(h)
        print("Coadding all datasets in the working folder!")
        print("Binning: " + str(BIN_SIZE) + " Angstrom")
        
        # find the minimum and maximum wavelength of all datasets
        wl_min = 10000.
        wl_max = 0.
        if CUSTOM_INTERVAL:
            wl_min = WAVE_MIN
            wl_max = WAVE_MAX
        else:
            for d in cdr_data:
                tmp_min = min(np.array(d['WAVELENGTH']))
                tmp_max = max(np.array(d['WAVELENGTH']))
                if tmp_min < wl_min:
                    wl_min = tmp_min
                if tmp_max > wl_max:
                    wl_max = tmp_max
                
        tot_wavelength = np.arange(wl_min, wl_max, BIN_SIZE)
        print("min wavelength = " + str(wl_min))
        print("max wavelength = " + str(wl_max))
        bins = len(tot_wavelength)
        tot_exptime = np.zeros(shape=len(tot_wavelength), dtype=np.float32)
        norm_exptime = np.zeros(shape=len(tot_wavelength), dtype=np.float32)
        tot_dq = np.full(shape=len(tot_wavelength), \
                         fill_value=20000, dtype=np.int32)
        tot_calib = np.zeros(shape=len(tot_wavelength), dtype=np.float32)
        tot_flt = np.zeros(shape=len(tot_wavelength), dtype=np.float32)
        tot_counts = np.zeros(shape=len(tot_wavelength), dtype=np.float32)
        tot_darks = np.zeros(shape=len(tot_wavelength), dtype=np.float32)
        tot_dc_err = np.zeros(shape=len(tot_wavelength), dtype=np.float32)
        tot_lya = np.zeros(shape=len(tot_wavelength), dtype=np.float32)
        tot_lya_err_up = np.zeros(shape=len(tot_wavelength), dtype=np.float32)
        tot_lya_err_down = np.zeros(shape=len(tot_wavelength), dtype=np.float32)
        
        print("Coadding spectra:")
        for i in range(bins):
            sys.stdout.write("|")
            for p in range(50):
                if i > p*bins/50.:
                    sys.stdout.write("=")
                else:
                    sys.stdout.write(" ")
            sys.stdout.write("| " + str(round(100.*float(i)/bins)) + " %\r")
            sys.stdout.flush()
            # calculate edges of the bin
            edge_min = tot_wavelength[i] - 0.5 * BIN_SIZE
            edge_max = tot_wavelength[i] + 0.5 * BIN_SIZE
            # calculate flux for the bin
    
            data_slices = []
    
            for j in range(len(cdr_data)):
                d_slice = \
                cdr_data[j][(cdr_data[j]['WAVELENGTH'] >= edge_min) & \
                            (cdr_data[j]['WAVELENGTH'] < edge_max)]
                
                # only append bin data from the datasets that contain 
                # the wavelength range in the first place
                if len(np.array(d_slice['WAVELENGTH'])) > 0:
                    data_slices.append(d_slice)
    
            bin_exptime = np.ndarray(shape=len(data_slices), dtype=np.float32)
            bin_normtime = np.ndarray(shape=len(data_slices), dtype=np.float32)
            bin_calib = np.ndarray(shape=len(data_slices), dtype=np.float32)
            bin_flat = np.ndarray(shape=len(data_slices), dtype=np.float32)
            bin_dq = 16384
            bin_dq_wgt = 0
            bin_counts = np.ndarray(shape=len(data_slices), dtype=np.float32)
            bin_darks = np.ndarray(shape=len(data_slices), dtype=np.float32)
            bin_dc_err = np.zeros(shape=len(data_slices), dtype=np.float32)
            bin_lya = np.ndarray(shape=len(data_slices), dtype=np.float32)
            bin_lya_err_up = np.zeros(shape=len(data_slices), dtype=np.float32)
            bin_lya_err_down = np.zeros(shape=len(data_slices), dtype=np.float32)
            bin_calib_wgt = np.ndarray(shape=len(data_slices), \
                                       dtype=np.float32)
            for d in range(len(data_slices)):
                dq_wgt = np.array(data_slices[d]['DQ_WGT'])
                if bin_dq_wgt < max(dq_wgt):
                    bin_dq_wgt = max(dq_wgt)
                dq = np.array(data_slices[d]['DQ'])
                if bin_dq > min(dq):
                    bin_dq = min(dq)
                exptime = np.array(data_slices[d]['EXPTIME'])
                calib = np.array(data_slices[d]['CALIB'])
                flat_corr = np.array(data_slices[d]['FLAT_CORR'])
                gcounts = np.array(data_slices[d]['GCOUNTS'])
                darks = np.array(data_slices[d]['DARK_CURRENT'])
                dc_err = dq_wgt*np.array(data_slices[d]['DARK_CURRENT_ERR'])
                lya = np.array(data_slices[d]['LYA_SCATTER'])
                lya_err_up = dq_wgt*np.array(data_slices[d]['LYA_SCATTER_ERR_UP'])
                lya_err_down = dq_wgt*np.array(data_slices[d]['LYA_SCATTER_ERR_DOWN'])
                calib_wgt = np.array(data_slices[d]['EXPTIME'])
                bin_calib_wgt[d] = np.sum(dq_wgt * calib_wgt)
                bin_exptime[d] = np.sum(dq_wgt * exptime)
                bin_normtime[d] = max(dq_wgt * exptime)
                bin_calib[d] = np.sum(calib) / float(len(calib))
                bin_flat[d] = np.sum(flat_corr) / float(len(flat_corr))
                bin_counts[d] = np.sum(dq_wgt * gcounts)
                bin_darks[d] = np.sum(dq_wgt * darks)
                bin_lya[d] = np.sum(dq_wgt * lya)
                
                # error propagation (the neighboring pixels in the same exposure
                # are highly covariant)
                for n in range(len(dc_err)):
                    bin_dc_err[d] = bin_dc_err[d] + dc_err[n]
                for n in range(len(lya_err_up)):
                    bin_lya_err_up[d] = bin_lya_err_up[d] + lya_err_up[n]
                for n in range(len(lya_err_down)):
                    bin_lya_err_down[d] = bin_lya_err_down[d] + lya_err_down[n]
    
    
            tot_exptime[i] = np.sum(bin_exptime)
            norm_exptime[i] = np.sum(bin_normtime)
            if np.sum(bin_calib_wgt) != 0:
                tot_calib[i] = np.sum(bin_calib_wgt * bin_calib)/\
                               np.sum(bin_calib_wgt)
                tot_flt[i] = np.sum(bin_calib_wgt * bin_flat)/\
                               np.sum(bin_calib_wgt)
            else:
                tot_calib[i] = 0.0
                if i > 0:
                    tot_flt[i] = tot_flt[i-1]
                else:
                    tot_flt[i] = 1.0
            tot_dq[i] = bin_dq
            tot_counts[i] = np.sum(bin_counts)
            tot_darks[i] = np.sum(bin_darks)
            tot_dc_err[i] = np.sqrt(np.sum(np.power(bin_dc_err, 2)))
            tot_lya[i] = np.sum(bin_lya)
            tot_lya_err_up[i] = np.sum(bin_lya_err_up)
            tot_lya_err_down[i] = np.sum(bin_lya_err_down)
        
        print("\n")
        tot_flux = np.divide(tot_counts - tot_darks - tot_lya, \
                          tot_calib*tot_exptime*tot_flt, \
                          out=np.zeros_like(tot_counts - tot_darks), \
                          where=tot_calib*tot_exptime*tot_flt != 0)
        
        if FELDMAN_COUSINS:
            cnt_err = calc_conf_lim_feldman_cousins(tot_counts, tot_darks + tot_lya)
        else:
            cnt_err = calc_conf_lim_kraft(tot_counts, tot_darks + tot_lya)
        cnt_err_down = cnt_err[0]
        cnt_err_up = cnt_err[1]
        
        tot_bkg_err_up = np.sqrt(np.power(tot_dc_err, 2.) + \
                                        np.power(tot_lya_err_up, 2.))
        tot_bkg_err_down = np.sqrt(np.power(tot_dc_err, 2.) + \
                                          np.power(tot_lya_err_down, 2.))
    
        # calculate err_up and err_down for the flux
        flux_err_up = np.divide(cnt_err_up,\
                                tot_calib*tot_exptime*tot_flt,\
                                out = np.zeros_like(tot_bkg_err_up),\
                                where=tot_calib*tot_exptime*tot_flt != 0)
        flux_err_down = np.divide(cnt_err_down,\
                                tot_calib*tot_exptime*tot_flt,\
                                out = np.zeros_like(tot_bkg_err_down),\
                                where=tot_calib*tot_exptime*tot_flt != 0)
    
        col1 = fits.Column(name='WAVELENGTH', \
                           format='D', \
                           unit='angstrom',\
                           array=tot_wavelength)
        col2 = fits.Column(name='FLUX', \
                           format='D', \
                           unit='erg/s/cm^2/A',\
                           array=tot_flux)
        col3 = fits.Column(name='FLUX_ERR_UP', \
                           format='D', \
                           unit='erg/s/cm^2/A',\
                           array=flux_err_up)
        col4 = fits.Column(name='FLUX_ERR_DOWN', \
                           format='D', \
                           unit='erg/s/cm^2/A',\
                           array=flux_err_down)
        col5 = fits.Column(name='GCOUNTS', \
                           format='I', \
                           unit='count',\
                           array=np.rint(tot_counts))
        col6 = fits.Column(name='BACKGROUND', \
                           format='D', \
                           unit='counts',\
                           array=tot_darks + tot_lya)
        col7 = fits.Column(name='BKG_ERR_UP', \
                           format='D', \
                           unit='counts',\
                           array=tot_bkg_err_up)
        col8 = fits.Column(name='BKG_ERR_DOWN', \
                           format='D', \
                           unit='counts',\
                           array=tot_bkg_err_down)
        col9 = fits.Column(name='DARK_CURRENT', \
                           format='D', \
                           unit='counts',\
                           array=tot_darks)
        col10 = fits.Column(name='DARK_CURRENT_ERR', \
                           format='D', \
                           unit='counts',\
                           array=tot_dc_err)      
        col11 = fits.Column(name='EXPTIME', \
                           format='D', \
                           unit='seconds',\
                           array=tot_exptime)
        col12 = fits.Column(name='DQ', \
                           format='I', \
                           unit='',\
                           array=tot_dq)
        col13 = fits.Column(name='CALIB', \
                           format='D', \
                           unit='counts cm**2 A / erg',\
                           array=tot_calib)
        col14 = fits.Column(name='FLAT_CORR',\
                           format='D',\
                           unit='',\
                           array = tot_flt)
        col15 = fits.Column(name='LYA_SCATTER', \
                           format='D', \
                           unit='counts',\
                           array=tot_lya)
        col16 = fits.Column(name='LYA_SCATTER_ERR_UP', \
                           format='D', \
                           unit='counts',\
                           array=tot_lya_err_up)
        col17 = fits.Column(name='LYA_SCATTER_ERR_DOWN', \
                           format='D', \
                           unit='counts',\
                           array=tot_lya_err_down)
        
        #hdu = fits.PrimaryHDU()
        #hdu = cdr_hdu[0]
        #hdu.header['TARGNAME'] = targname
        #hdu.header['OPT_ELEM'] = opt_elem
        #hdu.header["DATE"] = (datetime.utcnow().isoformat(timespec='seconds'),\
        #                  "File creation date")
        hdu = createCoAddPrimaryHeader(cdr_hdu, tot_wavelength, norm_exptime)
        hdu_binary = fits.BinTableHDU.from_columns([col1, col2, col3, col4,\
                                                    col5, col6, col7, col8,\
                                                    col9, col10, col11, col12,\
                                                    col13, col14, col15, col16,\
                                                    col17])
        hdul = fits.HDUList([hdu, hdu_binary])
        saved_file = path_sci + targname + "_spectrum.fits"
        hdul.writeto(saved_file, overwrite=True)
        print("Co-added spectrum is stored in " + saved_file)

    print("DONE!")
    
    
    
    
    
    
    
