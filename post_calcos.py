'''
This is the FaintCOS v1.2 post_calcos script

Purpose:
1. Few percent-level accurate dark current subtraction of a collection of COS FUV exposures
   (Worseck et al. 2016, ApJ, 825, 144; Makan et al. 2021, ApJ, 912, 38)
2. Estimation of scattered light from geocoronal Lyman alpha for G140L (Worseck et al. 2016)
3. Co-addition and rebinning of exposures from multiple COS FUV data sets
   (several visits, both FUV M gratings, several central wavelengths)
4. Calculation of proper Poisson flux error bars (e.g. Worseck et al. 2016)
5. Generation of 1D spectrum with HLSP-compliant metadata (as of Apr 09, 2021,
   https://outerspace.stsci.edu/display/MASTDOCS/HLSP+Contributor+Guide) 

Input:
1. Path to science data (single target, either COS G140L *or* G130M/G160M gratings)
2. Optional path to target-specific configuration file faintcos_config.py,
   otherwise take the default in the FaintCOS installation directory

Output:
1. MEF file with two FITS binary tables, extension 1: co-added and rebinned COS FUV spectrum,
   extension 2: provenance metadata
2. If HLSP output is desired (see config file) the primary output spectrum is copied to an
   hlsp* file following the HLSP naming convention
Note: The choice of the primary output spectrum depends on the number of data sets.
      For a target with a single data set take the combined spectrum rebinned by an integer number of pixels.
      For a target with multiple data sets, the individual data sets are rebinned onto a common wavelength
      grid with constant dispersion, by default assuming that target variability and flux calibration
      errors are negligible. For multiple data sets the rebinned spectra of individual data sets
      are optional (see config file).

Authors: Kirill Makan, Gabor Worseck
'''
version = '1.2'

import os
import sys
import math


from astropy.io import fits
from astropy.table import Table, vstack, hstack, Column
from astropy.stats import poisson_conf_interval
from astropy.time import Time
from datetime import datetime, timezone
from scipy.interpolate import interp1d
from scipy.stats import linregress
import numpy as np
from shutil import copy2


def calc_conf_lim_feldman_cousins(counts, bkg):
    '''
    Purpose:
    Calculation of two-sided equal-tailed 68.26% confidence (1 sigma) statistical uncertainty for Poisson
    counts accounting for known background (frequentist method: Feldman & Cousins 1998, Phys.Rev. D., 57, 3873).
    The confidence interval is for the signal, i.e. counts minus background. If counts < background
    compute a one-sided uncertainty from the one-sided 1 sigma upper limit (84.13% c.l.).
    Feldman & Cousins call this a sensitivity limit (their Section 6).    
    
    Parameters:
    1. counts: numpy array with (integer) Poisson counts (signal + background)
    2. bkg: numpy array with (non-integer) background
    
    Returns:
    cnt_err_down, cnt_err_up: two numpy arrays with lower and upper statistical uncertainty,
                              the sum of both corresponds to a 68.26% confidence level (1 sigma)
    '''
    # in order to install CustomConfLim, run 'python setup.py build' and 
    # 'python setup.py install' in the CustomConfLimits folder
    import CustomConfLim    
    
    cnt_err_up = np.zeros(shape=len(counts), dtype=np.float32)
    cnt_err_down = np.zeros(shape=len(counts), dtype=np.float32)
    
    print("Calculating stat. flux errors (frequentist method: Feldman & Cousins 1998):")
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


def calc_conf_lim_kraft(counts, bkg):
    '''
    Purpose:
    Calculation of two-sided minimal 68.26% integrated posterior density (1 sigma) statistical uncertainty
    for Poisson counts accounting for known background (Bayesian method: Kraft et al. 1991, ApJ, 374, 344;
    implemented in astropy.stats module). The credible interval is for the signal, i.e. counts minus background,
    and corresponds to the posterior maximum. The statistical uncertainty is calculated from the minimal
    68.27% credible interval (as in Kraft et al., hopefully).
    WARNING: Most astronomers prefer/assume equal-tailed credible intervals. The difference is subtle,
    but publications should be specific.
    
    Parameters:
    1. counts: numpy array with (integer) Poisson counts (signal + background)
    2. bkg: numpy array with (non-integer) background
    
    Returns:
    cnt_err_down, cnt_err_up: two numpy arrays with lower and upper statistical uncertainty,
                              the sum of both corresponds to a 68.26% credibility (1 sigma)    
    '''
    cnt_err_up = np.zeros(shape=len(counts), dtype=np.float32)
    cnt_err_down = np.zeros(shape=len(counts), dtype=np.float32)

    print("Calculating stat. flux errors (Bayesian method: Kraft et al. 1991):")
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


def calc_lya_scatter_model(gcounts, wave_arr, dq):
    '''
    Purpose:
    Calculation of scattered geocoronal Lyman alpha emission in COS G140L spectra
    using the simple model (Gaussian profile scaled to maximum Lyman alpha counts)
    from Worseck et al. 2016, ApJ, 825, 144.
    The maximum of the geocoronal Lyman alpha counts is estimated from the mean
    at 1215.37-1215.97A (broad unresolved profile).

    Parameters:
    1. gcounts: numpy array with 1D extracted count spectrum
    2. wave_arr: numpy array with wavelength grid in Angstroem
    
    Returns:
    1. Blya: numpy array with estimated scattered light profile given in counts
    2. Blya_err_up: numpy array with upper error of Blya
    3. Blya_err_down: numpy array with lower error of Blya
    The sum of both errors gives the 1 sigma (68% confidence level) uncertainty
    of the scattered light given the scarce data (estimated by parametric bootstrap
    in Worseck et al. 2016). The error is likely underestimated, because it assumes
    a rigid Gaussian shape for the scattered light profile. When calculating the
    total background uncertainty, this error must be treated as a systematic error.
    '''
    
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

    CountsInLyaRegion = gcounts[(wave_arr <= 1215.97) & (wave_arr >= 1215.37) & (dq == 0)]
    if len(CountsInLyaRegion)<3:
        print('Warning: Grid wire on Lyman alpha. Peak determined from 5A region.')
        CountsInLyaRegion = gcounts[(wave_arr <= 1218.17) & (wave_arr >= 1213.17) & (dq == 0)] 
    Clya = np.mean(CountsInLyaRegion)
    
    Blya = a*Clya*np.exp(-np.power(wave_arr - lam_0, 2.)/(2.*b**2))
    
    Clya_err_stat = math.sqrt(Clya)/math.sqrt(len(CountsInLyaRegion))/Clya
    err_up = np.sqrt(err_up_sys**2 + Clya_err_stat**2)
    err_down = np.sqrt(err_down_sys**2 + Clya_err_stat**2)
    
    return Blya, err_up*Blya, err_down*Blya


def createCDRPrimaryHeader(hdul, wave):
    '''
    Purpose:
    Create custom primary header for FITS file with 1D spectrum of single subexposure,
    i.e. a former CALCOS _x1d file after including the improved dark current estimate
    and scattered light (for G140L). The FITS header is almost HLSP-compliant and
    includes the relevant metadata for a COS FUV 1D spectrum.

    Parameters:
    1. hdul: HDU list from the CALCOS _x1d file 
    2. wave: numpy array with wavelengths of the 1D spectrum in Angstroem
    
    Returns:
    hdu: custom header for primary HDU of the cdr file 
    '''
    
    hdu = fits.PrimaryHDU()
    hdu.header['DATE'] = (datetime.utcnow().isoformat(timespec='seconds'),\
                            "File creation date")
    hdu.header['FILETYPE'] = (hdul[0].header['FILETYPE'], 
                              'type of data found in data file')
    hdu.header['TELESCOP'] = (hdul[0].header['TELESCOP'],
                              'telescope used to acquire data')
    hdu.header['INSTRUME'] = (hdul[0].header['INSTRUME'], 
                              'identifier for instrument used to acquire data')
    hdu.header['EQUINOX'] = (hdul[0].header['EQUINOX'], 
                             'equinox of celestial coord. system')
    hdu.header['TARGNAME'] = (hdul[0].header['TARGNAME'], 
                              'proposer\'s target name')
    hdu.header['RA_TARG'] = (hdul[0].header['RA_TARG'], 
                             '[deg] right ascention of the target')
    hdu.header['DEC_TARG'] = (hdul[0].header['DEC_TARG'], 
                              '[deg] declination of the target')
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
    hdu.header['SEGMENT'] = (hdul[0].header['SEGMENT'], 
                             'FUV detector segment name (FUVA, FUVB or BOTH)')
    hdu.header['LIFE_ADJ'] = (hdul[0].header['LIFE_ADJ'], 
                              'Life Time Adjustment Position')
    hdu.header['APERTURE'] = (hdul[0].header['APERTURE'], 
                              'aperture name')
    hdu.header['OPT_ELEM'] = (hdul[0].header['OPT_ELEM'], 
                              'optical element in use')
    hdu.header['CENWAVE'] = (hdul[0].header['CENWAVE'], 
                             '[Angstrom] grating central wavelength')
    hdu.header['BANDWID'] = (max(wave)-min(wave), 
                             '[Angstrom] bandwidth of the data')
    res = cos_res[(cos_res['LP'] == hdu.header['LIFE_ADJ']) & \
                  (cos_res['OPT_ELEM'] == hdu.header['OPT_ELEM']) &\
                  (cos_res['CENWAVE'] == hdu.header['CENWAVE'])]['R'][0]
    hdu.header['SPECRES'] = (res, 
                             'approx. resolving power at CENWAVE')
    hdu.header['CENTRWV'] = ((max(wave)+min(wave))/2.0, 
                             '[Angstrom] central wavelength of the data')
    hdu.header['MINWAVE'] = (min(wave), 
                             '[Angstrom] minimum wavelength in spectrum')
    hdu.header['MAXWAVE'] = (max(wave), 
                             '[Angstrom] maximum wavelength in spectrum')
    hdu.header['TIMESYS'] = ('UTC     ', 'time scale of time-related keywords')
    hdu.header['DATE-BEG'] = (hdul[1].header['DATE-OBS']+'T'+hdul[1].header['TIME-OBS'], 'exposure start time (ISO-8601 DateTime)')
    t_end = Time(hdul[1].header['EXPEND'], format='mjd', scale='utc', precision=0)
    hdu.header['DATE-END'] = (t_end.isot, 'exposure end time (ISO-8601 DateTime)')
    hdu.header['MJD-BEG'] = (hdul[1].header['EXPSTART'], '[d] exposure start time (Modified Julian Date)')
    hdu.header['MJD-END'] = (hdul[1].header['EXPEND'], '[d] exposure end time (Modified Julian Date)')
    hdu.header['EXPTIME'] = (hdul[1].header['EXPTIME'], '[d] exposure duration (seconds)--calculated')
    hdu.header['ASN_ID'] = (hdul[0].header['ASN_ID'], 
                            'unique identifier assigned to association')
    return hdu


def createDatasetPrimaryHeader(headers, wave):
    '''
    Purpose:
    Create custom primary header for FITS file with coadded calibrated spectrum of COS data set,
    i.e. from a collection of exposures with the same grating setting from one _asn association file.
    The FITS header includes the relevant metadata for a COS FUV 1D spectrum in HLSP-compliant keywords.
    HLSP keywords are added if HLSP output is desired.

    Parameters:
    1. headers: array with primary HDUs of individual exposures from _asn association file
    2. wave: numpy array with wavelengths of the coadded 1D spectrum in Angstroem
    
    Returns:
    hdu: custom header for primary HDU of the coadded spectrum from a single data set.
         This is written into the _dataset_sum file (unbinned spectrum) and the
         _bin file (spectrum binned by given number of original pixels)
    '''

    hdu = fits.PrimaryHDU()
    hdu.header['DATE'] = (datetime.utcnow().isoformat(timespec='seconds'),\
                            "File creation date")
    hdu.header['FILETYPE'] = (headers[0]['FILETYPE'], 
                              'type of data found in data file')
    hdu.header['OBSERVAT'] = (headers[0]['TELESCOP'],
                              'observatory used to obtain observation')
    hdu.header['TELESCOP'] = (headers[0]['TELESCOP'],
                              'telescope used to acquire data')
    hdu.header['INSTRUME'] = (headers[0]['INSTRUME'], 
                              'identifier for instrument used to acquire data')
    hdu.header['EQUINOX'] = (headers[0]['EQUINOX'], 
                             'equinox of celestial coord. system')    
    hdu.header['TARGNAME'] = (headers[0]['TARGNAME'], 
                              'proposer\'s target name')
    hdu.header['RA_TARG'] = (headers[0]['RA_TARG'], 
                             '[deg] right ascention of the target')
    hdu.header['DEC_TARG'] = (headers[0]['DEC_TARG'], 
                              '[deg] declination of the target')
    hdu.header['PROPOSID'] = (headers[0]['PROPOSID'], 
                              'PEP proposal identifier')
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
    hdu.header['SEGMENT'] = (headers[0]['SEGMENT'], 
                             'FUV detector segment name (FUVA, FUVB or BOTH)')
    hdu.header['LIFE_ADJ'] = (headers[0]['LIFE_ADJ'], 
                              'Life Time Adjustment Position')
    hdu.header['APERTURE'] = (headers[0]['APERTURE'], 
                             'aperture name')
    hdu.header['DISPERSR'] = (headers[0]['OPT_ELEM'], 
                             'name of dispersive element used')
    hdu.header['CENWAVE'] = (headers[0]['CENWAVE'], 
                             '[Angstrom] grating central wavelength')
    hdu.header['BANDWID'] = (max(wave)-min(wave), 
                             '[Angstrom] bandwidth of the data')
    res = cos_res[(cos_res['LP'] == headers[0]['LIFE_ADJ']) & \
                  (cos_res['OPT_ELEM'] == headers[0]['OPT_ELEM']) &\
                  (cos_res['CENWAVE'] == headers[0]['CENWAVE'])]['R'][0]
    hdu.header['SPECRES'] = (res, 
                             'approx. resolving power at CENWAVE')
    hdu.header['CENTRWV'] = ((max(wave)+min(wave))/2.0, 
                             '[Angstrom] central wavelength of the data')
    hdu.header['MINWAVE'] = (min(wave), 
                             '[Angstrom] minimum wavelength in spectrum')
    hdu.header['MAXWAVE'] = (max(wave), 
                             '[Angstrom] maximum wavelength in spectrum')
    hdu.header['BINNING'] = (BIN_PX, '[pixel] binning factor w.r.t. orig. pixel size')
    hdu.header['BUNIT'] = ('erg s^-1 cm^-2 A^-1', 'brightness unit')
    
    date_beg = []
    mjd_beg = []
    mjd_end = []
    exptime = 0
    for hdr in headers:
        if hdr['DATE-BEG'] not in date_beg:
            date_beg.append(hdr['DATE-BEG'])
            mjd_beg.append(hdr['MJD-BEG'])
            mjd_end.append(hdr['MJD-END'])
            exptime += hdr['EXPTIME']
    sorted_dates = [x for _,x in sorted(zip(mjd_beg, date_beg))]
    hdu.header['TIMESYS'] = (headers[0]['TIMESYS'], 'time scale of time-related keywords')
    hdu.header['DATE-BEG'] = (sorted_dates.pop(0), 'ISO-8601 DateTime of first exposure start')
    t_end = Time(max(mjd_end), format='mjd', scale='utc', precision=0)
    hdu.header['DATE-END'] = (t_end.isot, 'ISO-8601 DateTime of last exposure end')
    hdu.header['MJD-BEG'] = (min(mjd_beg), '[d] Mod. Julian Date of first exposure start')
    hdu.header['MJD-END'] = (max(mjd_end), '[d] Mod. Julian Date of last exposure end')
    hdu.header['MJD-MID'] = (0.5*(max(mjd_end)+min(mjd_beg)), '[d] MJD mid-time -- likely meaningless')
    hdu.header['XPOSURE'] = (exptime, '[s] observation duration --calculated')    
    hdu.header['FILE_ID'] = (headers[0]['ASN_ID'], 'Dataset identifier')

    if HLSP_write:
      hdu.header['HLSPTARG'] = (headers[0]['TARGNAME'], 'HLSP target designation')
      hdu.header['HLSPID'] = (HLSP_id, 'HLSP identifier (acronym)')
      hdu.header['HLSPNAME'] = (HLSP_name, 'title for HLSP project')
      hdu.header['HLSPLEAD'] = (HLSP_lead, 'full name of HLSP project lead')
      hdu.header['HLSPVER'] = (HLSP_ver, 'version identifier for HLSP product')
      hdu.header['DOI'] = (HLSP_doi, 'HLSP Digital Object Identifier')
      hdu.header['REFERENC'] = (HLSP_referenc, 'bibliogr. identifier (ADS bibcode)')
      hdu.header['LICENSE'] = ('CC BY 4.0', 'license for use of data')
      hdu.header['LICENURL'] = ('https://creativecommons.org/licenses/by/4.0/', 'Data license URL')
    
    return hdu


def createCoAddPrimaryHeader(headers, wave, exptime):
    '''
    Purpose:
    Create custom primary header for FITS file with coadded and rebinned spectrum from several COS data sets,
    i.e. from a collection of exposures from several _asn association files for the same target
    The FITS header includes the relevant metadata for a COS FUV 1D spectrum in HLSP-compliant keywords.
    HLSP keywords are added if HLSP output is desired.

    Parameters:
    1. headers: array with primary HDUs of individual exposures
    2. wave: numpy array with wavelengths of the coadded 1D spectrum in Angstroem
    3. exptime: numpy array with pixel exposure time in seconds
    
    Returns:
    hdu: custom header for primary HDU of the coadded and rebinned spectrum (_spectrum file).
    '''

    unique_pid = []
    unique_lifeadj = []
    unique_aper = []
    unique_disp = []
    unique_cenwave = []
    
    for hdr in headers:
        if hdr['PROPOSID'] not in unique_pid:
            unique_pid.append(hdr['PROPOSID'])
        if hdr['LIFE_ADJ'] not in unique_lifeadj:
            unique_lifeadj.append(hdr['LIFE_ADJ'])
        if hdr['APERTURE'] not in unique_aper:
            unique_aper.append(hdr['APERTURE'])
        if hdr['OPT_ELEM'] not in unique_disp:
            unique_disp.append(hdr['OPT_ELEM'])
        if hdr['CENWAVE'] not in unique_cenwave:
            unique_cenwave.append(hdr['CENWAVE'])

    hdu = fits.PrimaryHDU()
    hdu.header['DATE'] = (datetime.utcnow().isoformat(timespec='seconds'),\
                            "File creation date")
    hdu.header['FILETYPE'] = (headers[0]['FILETYPE'], 
                              'type of data found in data file')
    hdu.header['OBSERVAT'] = (headers[0]['TELESCOP'],
                              'observatory used to obtain observation')
    hdu.header['TELESCOP'] = (headers[0]['TELESCOP'],
                              'telescope used to acquire data')
    hdu.header['INSTRUME'] = (headers[0]['INSTRUME'], 
                              'identifier for instrument used to acquire data')
    hdu.header['EQUINOX'] = (headers[0]['EQUINOX'], 
                             'equinox of celestial coord. system')
    hdu.header['TARGNAME'] = (headers[0]['TARGNAME'], 
                              'proposer\'s target name')
    hdu.header['RA_TARG'] = (headers[0]['RA_TARG'], 
                             '[deg] right ascention of the target')
    hdu.header['DEC_TARG'] = (headers[0]['DEC_TARG'], 
                              '[deg] declination of the target')
    if len(unique_pid) == 1:
        hdu.header['PROPOSID'] = (unique_pid[0], 'PEP proposal identifier')
    else:
        hdu.header['PROPOSID'] = ('MULTI', 'PEP proposal identifier')
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
    segments = [hdr['SEGMENT'] for hdr in headers]
    segm = ""
    if 'FUVA' in segments and 'FUVB' in segments:
        segm = 'BOTH'
    else:
        segm = segments[0]        
    hdu.header['SEGMENT'] = (segm, 'FUV detector segment name (FUVA, FUVB or BOTH)')
    if len(unique_lifeadj) == 1:
        hdu.header['LIFE_ADJ'] = (unique_lifeadj[0], 'Life Time Adjustment Position')
    else:
        hdu.header['LIFE_ADJ'] = ('MULTI', 'Life Time Adjustment Position')
    if len(unique_aper) == 1:
        hdu.header['APERTURE'] = (unique_aper[0], 'aperture name')
    else:
        hdu.header['APERTURE'] = ('MULTI', 'aperture name')
    if len(unique_disp) == 1:
        hdu.header['DISPERSR'] = (unique_disp[0], 'name of dispersive element used')
    else:
        hdu.header['DISPERSR'] = ('MULTI', 'name of dispersive element used')
    if len(unique_cenwave) == 1:
        hdu.header['CENWAVE'] = (unique_cenwave[0], '[Angstrom] grating central wavelength')
    else:
        hdu.header['CENWAVE'] = ('MULTI', 'grating central wavelength')
    hdu.header['BANDWID'] = (max(wave)-min(wave), 
                             '[Angstrom] bandwidth of the data')
    hdu.header['SPECRES'] = (min([hdr['SPECRES'] for hdr in headers]), 
                             'min. res. power at CENWAVE from all data sets')
    hdu.header['CENTRWV'] = ((max(wave)+min(wave))/2.0, 
                             '[Angstrom] central wavelength of the data')
    hdu.header['MINWAVE'] = (min(wave), 
                             '[Angstrom] minimum wavelength in spectrum')
    hdu.header['MAXWAVE'] = (max(wave), 
                             '[Angstrom] maximum wavelength in spectrum')
    hdu.header['BINNING'] = (BIN_SIZE, '[Angstrom] bin size of wavelength axis')
    hdu.header['BUNIT'] = ('erg s^-1 cm^-2 A^-1', 'brightness unit')
    hdu.header['TIMESYS'] = (headers[0]['TIMESYS'], 'time scale of time-related keywords')
    
    date_beg = []
    mjd_beg = []
    mjd_end = []
    tottime = 0
    for hdr in headers:
        if hdr['DATE-BEG'] not in date_beg:
            date_beg.append(hdr['DATE-BEG'])
            mjd_beg.append(hdr['MJD-BEG'])
            mjd_end.append(hdr['MJD-END'])
            tottime += hdr['EXPTIME']
    sorted_dates = [x for _,x in sorted(zip(mjd_beg, date_beg))]
    hdu.header['DATE-BEG'] = (sorted_dates.pop(0), 'ISO-8601 DateTime of first exposure start')
    t_end = Time(max(mjd_end), format='mjd', scale='utc', precision=0)
    hdu.header['DATE-END'] = (t_end.isot, 'ISO-8601 DateTime of last exposure end')
    hdu.header['MJD-BEG'] = (min(mjd_beg), '[d] Mod. Julian Date of first exposure start')
    hdu.header['MJD-END'] = (max(mjd_end), '[d] Mod. Julian Date of last exposure end')
    hdu.header['MJD-MID'] = (0.5*(max(mjd_end)+min(mjd_beg)), '[d] MJD mid-time -- likely meaningless')
    hdu.header['XPOSURE'] = (round(tottime,3), '[s] total exposure time of data sets')
    hdu.header['MXPOSURE'] = (round(float(max(exptime)),3), '[s] max. exposure time in spectral range')

    if HLSP_write:
      hdu.header['HLSPTARG'] = (headers[0]['TARGNAME'], 'HLSP target designation')
      hdu.header['HLSPID'] = (HLSP_id, 'HLSP identifier (acronym)')
      hdu.header['HLSPNAME'] = (HLSP_name, 'title for HLSP project')
      hdu.header['HLSPLEAD'] = (HLSP_lead, 'full name of HLSP project lead')
      hdu.header['HLSPVER'] = (HLSP_ver, 'version identifier for HLSP product')
      hdu.header['DOI'] = (HLSP_doi, 'HSLP Digital Object Identifier')
      hdu.header['REFERENC'] = (HLSP_referenc, 'bibliogr. identifier (ADS bibcode)')
      hdu.header['LICENSE'] = ('CC BY 4.0', 'license for use of data')
      hdu.header['LICENURL'] = ('https://creativecommons.org/licenses/by/4.0/', 'Data license URL')
    
    return hdu
    

if __name__ == "__main__": 

    print("FaintCOS v"+version, flush=True)
    path_sci = "."
    
    # Get paths to science data and FaintCOS config file from command line.
    # If 2nd path to config file is given then frontload it into the search path for imports to force local import.
    # Then import FaintCOS config file
    if len(sys.argv) == 3:
        path_sci = sys.argv[1]
        sys.path.insert(0, sys.argv[2])
    else:
        if len(sys.argv) == 2:
            path_sci = sys.argv[1]
    from faintcos_config import *
    
    # find all corrtag files in the science directory
    path_corrtag = [f for f in os.listdir(path_sci) if "corrtag" in f]
    if path_sci != ".":
        path_corrtag = [path_sci + s for s in path_corrtag]
    else:
        path_sci = ""
        
    # load all relevant metadata from corrtags and datasets into astropy tables and sort them
    # dataset table is built from unique dataset IDs in corrtag files
    datasets = Table(names=("PROPOSID", "FILE_ID", "TARGNAME", "DATE-OBS", "DATE-BEG", "DATE-END", "MJD-BEG", "MJD-END", \
                            "LIFE_ADJ", "APERTURE", "DISPERSR", "CENWAVE", "XPOSURE", "EXPA", "EXPB"),\
                     dtype=('i4', 'S10', 'S25', 'S10', 'S19', 'S19', 'f8', 'f8', 'i4', 'S4', 'S5', 'i4', 'f8', 'i4', 'i4'))
    corrtags = Table(names=("TARGNAME", "FILE_ID", "CORRTAG_FILE", "DATE-BEG",\
                            "DISPERSR", "SEGMENT", "CENWAVE",\
                            "EXPTIME"), \
                     dtype=('S25', 'S10', 'S30', 'S19', 'S5' , 'S5', 'i4', 'f8'))
    datasets['XPOSURE'].format = '4.3f'
    corrtags['EXPTIME'].format = '4.3f'
    for f in path_corrtag:
        h0 = fits.open(f)[0].header
        h1 = fits.open(f)[1].header
        if h1['EXPFLAG'] != 'NORMAL':
            continue
        corrtags.add_row([h0['TARGNAME'], h0['ASN_ID'], f.split('/')[-1], \
                          h1['DATE-OBS'] + "T" + h1['TIME-OBS'],\
                          h0['OPT_ELEM'], h0['SEGMENT'], \
                          h0['CENWAVE'], h1['EXPTIME']])
    corrtags.sort(["TARGNAME", 'FILE_ID', "DATE-BEG", 'SEGMENT'])
    unique_datasets = np.unique(np.array(corrtags['FILE_ID']))
    for asn_id in unique_datasets:
        visit = corrtags[corrtags['FILE_ID'] == asn_id]
        hdul = fits.open(path_sci + visit[len(visit)-1]['CORRTAG_FILE'])
        mjd_end = hdul[1].header['EXPEND']
        t_end = Time(mjd_end, format='mjd', scale='utc', precision=0)
        date_end = t_end.isot
        hdul.close()
        hdul = fits.open(path_sci + visit[0]['CORRTAG_FILE'])
        propid = hdul[0].header['PROPOSID']
        target = hdul[0].header['TARGNAME']
        date_obs = hdul[1].header['DATE-OBS']
        date_beg = hdul[1].header['DATE-OBS']+"T"+hdul[1].header['TIME-OBS']
        mjd_beg = hdul[1].header['EXPSTART']
        life_adj = hdul[0].header['LIFE_ADJ']
        aper = hdul[0].header['APERTURE']
        hdul.close()
        exp_time_a = np.array(visit[visit['SEGMENT'] == 'FUVA']['EXPTIME']) 
        exp_time_b = np.array(visit[visit['SEGMENT'] == 'FUVB']['EXPTIME']) 
        exp_time = max(np.array([np.sum(exp_time_a), np.sum(exp_time_b)]))
        num_fuva = len(np.array(visit[visit['SEGMENT'] == 'FUVA']))
        num_fuvb = len(np.array(visit[visit['SEGMENT'] == 'FUVB']))
        cenwave = visit[0]['CENWAVE']
        opt_elem = visit[0]['DISPERSR']
        datasets.add_row([propid, asn_id, target, date_obs, date_beg, date_end, mjd_beg, mjd_end, \
                          life_adj, aper, opt_elem, cenwave, exp_time, num_fuva, num_fuvb])
    datasets.sort(['PROPOSID','FILE_ID','DATE-OBS'])

    # reduce individual exposures as in Worseck et al. 2016: For every exposure estimate the dark current in
    # science extraction aperture from dark frames with similar pulse height distribution (sensitive to environmental conditions)
    # as in offset calibration windows in science exposure _corrtag file. For G140L exposures also estimate the scattered light
    if REDUCE_EXPOSURES:
        # Get directory with CALCOS calibration files
        try:
            path_ref = os.environ['lref']
        except KeyError:
            print("ERROR: lref is not defined!")
            sys.exit() 
    
        # Get directory with dark frames (_corrtag files, ideally reduced with the same CALCOS version and calibration files)
        try:
            path_dark = os.environ['ldark']
        except KeyError:
            print("ERROR: ldark is not defined!")
            sys.exit()     
    
        # load all dark frames in the 'ldark' directory
        print("Loading dark frames...", end=" ", flush=True)
        path_darkframes = [f for f in os.listdir(path_dark) if "corrtag" in f]
        path_darkframes = [path_dark + s for s in path_darkframes]
        dark_file = []
        dark_expstart = []
        dark_segment = []
        dark_voltage = []
        for d in path_darkframes:
            dark_tmp = fits.open(d)
            dark_file.append(dark_tmp[0].header['FILENAME'])
            dark_expstart.append(dark_tmp[1].header['EXPSTART'])
            dark_segment.append(dark_tmp[0].header['SEGMENT'])
            dark_voltage.append(dark_tmp[1].header["HVLEVEL" + \
                                                   dark_segment[-1].split("FUV")[1]])
            dark_tmp.close()
        darkframes = Table([dark_file, path_darkframes, \
                            dark_expstart, dark_segment, dark_voltage],\
                           names = ("FILE", "PATH", "EXPSTART", \
                                    "SEGMENT", "VOLTAGE"))
        print("OK")
        print(str(len(darkframes['FILE'])) + " dark frames have been found!")

        # print sorted lists of corrtag files and datasets in terminal
        print("Valid 'corrtag' files in the input directory:")
        corrtags['FILE_ID','TARGNAME','CORRTAG_FILE','CENWAVE','EXPTIME'].pprint()
        print("##################################################################")
        print("\n")

        print("Valid datasets in the input directory:")
        datasets['FILE_ID','TARGNAME','DATE-OBS','XPOSURE','DISPERSR','CENWAVE','EXPA','EXPB'].pprint()
        print("###################################################################")
        print("\n")
    
        a = input("Do you wish to proceed? (y/n)")
        if a != 'y' and a != 'Y':
            print("Canceled by user!")
            sys.exit()
    
        # determine background (dark current and scattered light) for every corrtag
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

            # open pulse height calibration file to get pulse height limits
            # these limits should be identical to those given in the FaintCOS config file
            pha_file = path_ref + corr_hdul[0].header['PHATAB'].split('$')[-1]
            pha_limits = Table(fits.open(pha_file)[1].data)
            pha_max = pha_limits[(pha_limits['OPT_ELEM'] == opt_elem) & \
                                 (pha_limits['SEGMENT'] == segm)]['ULT']
            pha_min = pha_limits[(pha_limits['OPT_ELEM'] == opt_elem) & \
                                 (pha_limits['SEGMENT'] == segm)]['LLT']
                             
            # open the uncalibrated 1D spectrum produced by CALCOS (_x1d)
            path_x1d = path_sci + corr_prefix + "_x1d.fits"
            hdul_x1d = fits.open(path_x1d)
            data_x1d = Table(hdul_x1d[1].data)
            data_x1d = data_x1d[data_x1d['SEGMENT'] == segm]
                             
            # select only valid dark frames for this corrtag (taken within defined time frame with the same detector segment at the same detector voltage)
            # if there are not enough matching dark frames it is very likely that the detector voltages do not match, so you have no choice as to accept a
            # systematic error in the dark current (either over- or underestimate depending on the calibration windows and the gain sag state of the detector)
            val_darks = []
            sel_darks = darkframes[(darkframes['SEGMENT'] == segm) & \
                                   (darkframes['VOLTAGE'] == voltage) & \
                                   (darkframes['EXPSTART'] > exp_start - DARK_EXPSTART_INTERVAL) & \
                                   (darkframes['EXPSTART'] < exp_start + DARK_EXPSTART_INTERVAL)]
            for d in sel_darks['PATH']:
                val_darks.append(fits.open(d))
            if len(val_darks) >= MIN_DARKS:
                print(str(len(val_darks)) + " valid dark frames.")
            else:
                print("Not enough dark frames for this exposure taken at voltage level {} !".format(voltage))
                print("This likely means that the science voltage level has not been monitored by STScI.")
                print("A mismatch in the science vs. dark frame voltage will lead to ~10% systematic errors in the estimated dark current (Makan et al. 2021, ApJ, 912, 38).")
                a = input("Do you wish to include dark frames taken at all voltage levels? (y/n)")
                if a == 'y' or a == 'Y':
                    val_darks = []
                    sel_darks = darkframes[(darkframes['SEGMENT'] == segm) & \
                                           (darkframes['EXPSTART'] > exp_start - DARK_EXPSTART_INTERVAL) & \
                                           (darkframes['EXPSTART'] < exp_start + DARK_EXPSTART_INTERVAL)]
                    for d in sel_darks['PATH']:
                        val_darks.append(fits.open(d))
                    if len(val_darks) >= MIN_DARKS:
                        print(str(len(val_darks)) + " valid dark frames.")
                else:
                    print("Canceled by user!")
                    sys.exit()
        
            # find the dispersion function for the Doppler shift xdopp(wavelength), this defines the calibration windows
            wl = np.array(corr_data[(corr_data['WAVELENGTH'] > 1) & \
                                    (corr_data['YFULL'] > ap_spec[0]) & \
                                    (corr_data['YFULL'] < ap_spec[1]) & \
                                    (corr_data['DQ'] == 0)]['WAVELENGTH'])
            xdopp = np.array(corr_data[(corr_data['WAVELENGTH'] > 1) & \
                                       (corr_data['YFULL'] > ap_spec[0]) & \
                                       (corr_data['YFULL'] < ap_spec[1]) & \
                                       (corr_data['DQ'] == 0)]['XDOPP'])
            xdopp_disp = linregress(wl, xdopp)[0:2]

            # find bad regions (with extended geocoronal emission) in xdopp coordinates
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
            
            # remove counts outside of detector's active area
            corr_data = corr_data[(corr_data['XDOPP'] < 15000) & \
                                  (corr_data['XDOPP'] > 1500)]
        
            # produce a cumulative pulse height distribution for the corrtag
            # select background calibration windows in the corrtag data
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
            
            # cumulative pulse height distribution
            unique, counts = np.unique(corr_bkg, return_counts=True)
            corr_pha_dist = np.zeros(shape=32, dtype=np.int32)
            for i in range(len(unique)):
                corr_pha_dist[unique[i]] = counts[i]
            corr_data_cumsum = np.cumsum(corr_pha_dist)
            corr_max_counts = corr_data_cumsum[-1]
            corr_data_cumsum = corr_data_cumsum / corr_max_counts
        
            # produce cumulative pulse height distribution for every dark frame
            # and compare it to the pulse height distribution of the corrtag
            # with a Kolmogorov–Smirnov test
            KS_values = np.zeros(shape=(len(val_darks)), dtype=np.float32)
            dark_max_counts = np.zeros(shape=(len(val_darks)), dtype=np.float32)
            d = 0
            for i in range(len(val_darks)):
                dark_data = Table(val_darks[i][1].data)
                # remove bad regions from the dark frame
                for reg in bad_reg_darks:
                    dark_data = dark_data[(dark_data['XDOPP'] < reg[0]) | \
                                          (dark_data['XDOPP'] > reg[1])]
                # remove counts outside of detector's active area
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
                # cumulative pulse height distribution
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
        
            # calculate scaling factor between the corrtag and combined darks
            scaling_factor = corr_max_counts / dark_total_counts
            scaling_factor_err = math.sqrt(scaling_factor * \
                                           (1. + 1. / dark_total_counts) / dark_total_counts)
        
            print("Average scaling factor: " + str(scaling_factor*number_of_darks))
        
            # coadd best dark frames
            darks_combined = Table(val_darks[dark_sorted_KS[0]][1].data)
            for i in range(1, number_of_darks):
                darks_combined = vstack([darks_combined, \
                                         Table(val_darks[dark_sorted_KS[i]][1].data)])
            
            # in the combined dark extract the 2d dark current spectrum within the science aperture using the science PHA limits
            dark_psa = darks_combined[(darks_combined['YFULL'] >= ap_spec[0]) &\
                                      (darks_combined['YFULL'] <= ap_spec[1]) & \
                                      (darks_combined['PHA'] >= pha_min) & \
                                      (darks_combined['PHA'] <= pha_max)]
            dark_psa_xfull = np.array(dark_psa['XFULL'])
        
            # shift xfull according to XSHIFT (darks are taken at nominal FPPOS whereas science data use multiple offsets)
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

            # Estimate smoothed dark current with running average between first and last valid points
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
            
            # Back out the flux calibration curve used by CALCOS
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
        
            # Reset CALCOS data flags and weights: Good data (DQ=0) and low-response regions (DQ=1024) are included in coadd,
            # all other regions (including grid wires with DQ=4) are excluded from coadd.
            # Experience shows that grid wires cannot be flatfielded out in Poisson-limited data!
            dq = np.array(data_x1d['DQ'])[0]
            dq_wgt = np.zeros(shape=len(dq), dtype=np.int32)
            dq[dq == 1024] = 2
            for q in range(len(dq_wgt)):
                if dq[q] == 0 or dq[q] == 2:
                    dq_wgt[q] = 1
                else:
                    dq_wgt[q] = 0

            # Back out 1D low-order flatfield calibration applied by CALCOS
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

            # calculate scattered light from geocoronal Lyman alpha emission in G140L spectra using model from Worseck et al. 2016,
            # requires Lyman alpha placed on COS detector, so this works for cenwave 800A and 1105A
            # For any other cenwave (1230A, 1280A) set contamination to zero
            # For mixed programs (e.g. 1105A and 1280A) the scattered light in the 1280A spectra may be estimated a posteriori
            # from the 1105A spectra for similar orientations of the orbits (length, solar altitude)       
            if opt_elem == 'G140L' and (cenwave == 800 or cenwave == 1105):
                bkg_lya, bkg_lya_err_up, bkg_lya_err_down = \
                calc_lya_scatter_model(data_x1d['GCOUNTS'][0], \
                                       data_x1d['WAVELENGTH'][0], data_x1d['DQ'][0])
            else:
                bkg_lya = np.zeros(shape = len(gross), dtype = np.float32)
                bkg_lya_err_up = np.zeros(shape = len(gross), dtype = np.float32)
                bkg_lya_err_down = np.zeros(shape = len(gross), dtype = np.float32)
        
            # create FITS file with fully reduced 1D spectrum as a binary table,
            # naming convention is exposure prefix + _cdr_ + detector segment,
            # all coadds and rebins use these individual files, original x1d files are not used any further
            col1 = fits.Column(name='WAVELENGTH', \
                               format='D', \
                               unit='Angstrom',\
                               array=np.array(data_x1d['WAVELENGTH'])[0])
            col2 = fits.Column(name='GCOUNTS', \
                               format='I', \
                               unit='count',\
                               array=np.rint(data_x1d['GCOUNTS'])[0])
            col3 = fits.Column(name='EXPTIME', \
                               format='D', \
                               unit='s',\
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
                               unit='count',\
                               array=dark_hist_mean)
            col7 = fits.Column(name='DARK_CURRENT_ERR', \
                               format='D', \
                               unit='count',\
                               array=dark_hist_error)
            col8 = fits.Column(name='CALIB', \
                               format='D', \
                               unit='count cm^2 A / erg',\
                               array=calib)
            col9 = fits.Column(name='FLAT_CORR',\
                               format='D',\
                               unit='',\
                               array = flt_curve)
            col10 = fits.Column(name='LYA_SCATTER', \
                                format='D', \
                                unit='count',\
                                array=bkg_lya)
            col11 = fits.Column(name='LYA_SCATTER_ERR_UP', \
                                format='D', \
                                unit='count',\
                                array=bkg_lya_err_up)
            col12 = fits.Column(name='LYA_SCATTER_ERR_DOWN', \
                                format='D', \
                                unit='count',\
                                array=bkg_lya_err_down)
            hdul_x1d[0].header['SEGMENT'] = segm
            hdu = createCDRPrimaryHeader(hdul_x1d, wave)
            hdu_binary = fits.BinTableHDU.from_columns([col1, col2, col3, col4,\
                                                        col5, col6, col7, col8,\
                                                        col9, col10, col11, col12])
            hdul = fits.HDUList([hdu, hdu_binary])
            saved_file = path_sci + corr_prefix + "_cdr_" + segm + ".fits"
            hdul.writeto(saved_file, overwrite=True)
            print("Results were written into " + saved_file +"\n")
    
            # close all darkframes
            for d in val_darks:
                d.close()

    # Find all fully calibrated 1D spectra in science directory
    if path_sci != "":
        cdr_files = [f for f in os.listdir(path_sci) if "cdr" in f]
        cdr_files = [path_sci + s for s in cdr_files]
    else:
        cdr_files = [f for f in os.listdir(".") if "cdr" in f]
        
    # coadding exposures for every dataset (unbinned)
    # because all exposures in a dataset share the same wavelength calibration this is simple: stacking of the counts and recalibrating them in flux
    print("Coadding exposures for every dataset ...")
    hdul_ar = []
    visit_data = []
    visit_hdu = []
    for f in cdr_files:
        hdul_ar.append(fits.open(f))
    for asn in unique_datasets:
        targname = ""
        segm = ""
        seg_fuva = []
        seg_fuvb = []
        headers = []
        for hdul in hdul_ar:
            if hdul[0].header['ASN_ID'] == asn.decode("utf-8"):
                headers.append(hdul[0].header)
                targname = hdul[0].header['TARGNAME']
                if hdul[0].header['SEGMENT'] == 'FUVA':
                    seg_fuva.append(hdul)
                else:
                    seg_fuvb.append(hdul)
        exposures = []
        if len(seg_fuvb) > 0:
            exposures.append(seg_fuvb)
            segm = 'FUVB'
        if len(seg_fuva) > 0:
            exposures.append(seg_fuva)
            segm = 'FUVA'
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
                    # sum up dark current, propagate estimated error
                    totaldark = np.add(totaldark, dq_wgt * darkcurrent)
                    darkcurrent_error = np.array(exp_data['DARK_CURRENT_ERR'])
                    totaldark_error = \
                    np.sqrt(np.add(np.power(totaldark_error, 2), \
                                   np.power(dq_wgt * darkcurrent_error, 2)))
                    # sum up scattered light, estimated error accounted as systematic error
                    lya_scatter = np.array(exp_data['LYA_SCATTER'])
                    total_lya = np.add(total_lya, dq_wgt * lya_scatter)   
                    lya_scatter_err_up = np.array(exp_data['LYA_SCATTER_ERR_UP'])
                    total_lya_err_up = total_lya_err_up + dq_wgt * lya_scatter_err_up
                    lya_scatter_err_down = np.array(exp_data['LYA_SCATTER_ERR_DOWN'])
                    total_lya_err_down = total_lya_err_down + dq_wgt * lya_scatter_err_down
                    # calculate total flux calibration curve
                    np.seterr(divide='ignore')
                    flux_calib = np.array(exp_data['CALIB'])
                    for  n in range(len(total_flux_calib)):
                        if (total_flux_calib[n] == 0) | \
                           (np.isnan(total_flux_calib[n])):
                            total_flux_calib[n] = flux_calib[n]
                calib = total_flux_calib
                # calculate total flux and store arrays in astropy table
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
        # store detector segment spectra together in a single table
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

        # Create primary HDU with proper keywords and store the coadded spectrum as binary table
        # visit_hdu and visit_data store primary HDUs and coadds from all datasets of the target
        # the header is for the rebinned spectrum by default, so generate a copy and overwrite the binning
        # _dataset_sum file contains coadded unbinned spectrum for single data set, this is typically for information only
        headers[0]['SEGMENT'] = segm
        hdu = createDatasetPrimaryHeader(headers, coadded_data['WAVELENGTH'])
        hdu.header['COMMENT'] = "Coadded rebinned spectrum for single data set with proper calibration" 
        visit_data.append(coadded_data)
        visit_hdu.append(hdu)
        hdu1 = createDatasetPrimaryHeader(headers, coadded_data['WAVELENGTH'])
        hdu1.header['BINNING'] = (1, '[pixel] binning factor w.r.t. orig. pixel size')
        hdu1.header['COMMENT'] = "Coadded unbinned spectrum for single data set with proper calibration" 
        binary_hdu = fits.BinTableHDU(coadded_data, name='SCI', ver=1)
        hdul = fits.HDUList([hdu1, binary_hdu])
        saved_file = path_sci + asn.decode("utf-8") + "_dataset_sum.fits"
        hdul.writeto(saved_file, overwrite=True)
        print(asn.decode("utf-8") + " is complete." + \
              " The spectrum is stored in " + saved_file)

    # Bin each coadded 1D spectrum (individual data set) with an integer binning factor BIN_PX
    # to reduce oversampling in wavelength and to increase discrete sampling of flux
    # calculate 1 sigma statistical uncertainty of the flux (either frequentist or Bayesian method)
    if BIN_DATASET:
        print("\n")
        print("Binning every dataset by " + str(BIN_PX) + " pixels ...")
        for d in range(len(visit_data)):
            asn = visit_hdu[d].header['FILE_ID']
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
            # Rebinning
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
            # calculate 1 sigma stat. uncertainties, see above descriptions for details
            # convert uncertainties from counts to flux
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

            # Store arrays in astropy table for further processing, sum the background components
            outtab = Table([binned_wave,binned_flux,flux_err_up,flux_err_down,np.rint(binned_gcounts),binned_dc+binned_lya,binned_bkg_err_up,binned_bkg_err_down,binned_dc,binned_dc_err,binned_exptime,binned_dq,binned_calib,binned_flt,binned_lya,binned_lya_err_up,binned_lya_err_down], names=('WAVELENGTH','FLUX','FLUX_ERR_UP','FLUX_ERR_DOWN','GCOUNTS','BACKGROUND','BKG_ERR_UP','BKG_ERR_DOWN','DARK_CURRENT','DARK_CURRENT_ERR','EXPTIME','DQ','CALIB','FLAT_CORR','LYA_SCATTER','LYA_SCATTER_ERR_UP','LYA_SCATTER_ERR_DOWN'))

            # Optionally trim detector edges (shortest and longest wavelengths outside active area) from the rebinned spectrum
            if TRIM_EDGE:
                tmpindex = np.where(outtab['DQ']==0)[0]
                imin = tmpindex[0]
                imax = tmpindex[len(tmpindex)-1]
                outtab = outtab[imin:imax+1]
            # Optionally restrict wavelength range of output spectrum, this is good for blue modes (e.g. G140L/800A) that include poorly calibrated low-sensitivity range <1100A
            if TRIM_WAVE:
                outtab = outtab[(outtab['WAVELENGTH']>=TRIM_MIN) & (outtab['WAVELENGTH']<=TRIM_MAX)]

            # Create output binary FITS table HDU from astropy table columns with specified data types, formatting, and units
            col1 = fits.Column(name='WAVELENGTH', format='D', disp='F10.4', unit='Angstrom', array=outtab['WAVELENGTH'])
            col2 = fits.Column(name='FLUX', format='D', disp='E13.7', unit='erg s^-1 cm^-2 Angstrom^-1', array=outtab['FLUX'])
            col3 = fits.Column(name='FLUX_ERR_UP', format='D', disp='E13.7', unit='erg s^-1 cm^-2 Angstrom^-1', array=outtab['FLUX_ERR_UP'])
            col4 = fits.Column(name='FLUX_ERR_DOWN', format='D', disp='E13.7', unit='erg s^-1 cm^-2 Angstrom^-1', array=outtab['FLUX_ERR_DOWN'])
            col5 = fits.Column(name='GCOUNTS', format='I', disp='I8', unit='count', array=outtab['GCOUNTS'])
            col6 = fits.Column(name='BACKGROUND', format='D', disp='F12.4', unit='count', array=outtab['BACKGROUND'])
            col7 = fits.Column(name='BKG_ERR_UP', format='D', disp='F12.4', unit='count', array=outtab['BKG_ERR_UP'])
            col8 = fits.Column(name='BKG_ERR_DOWN', format='D', disp='F12.4', unit='count', array=outtab['BKG_ERR_DOWN'])
            col9 = fits.Column(name='DARK_CURRENT', format='D', disp='F12.4', unit='count', array=outtab['DARK_CURRENT'])
            col10 = fits.Column(name='DARK_CURRENT_ERR', format='D', disp='F15.4', unit='count', array=outtab['DARK_CURRENT_ERR'])      
            col11 = fits.Column(name='EXPTIME', format='D', disp='F10.3', unit='s', array=outtab['EXPTIME'])
            col12 = fits.Column(name='DQ', format='I', disp='I5', unit='', array=outtab['DQ'])
            col13 = fits.Column(name='CALIB', format='D', disp='E13.7', unit='count cm^2 Angstrom/erg',array=outtab['CALIB'])
            col14 = fits.Column(name='FLAT_CORR', format='D', disp='F9.5', unit='', array=outtab['FLAT_CORR'])
            col15 = fits.Column(name='LYA_SCATTER', format='D', disp='F12.4', unit='count', array=outtab['LYA_SCATTER'])
            col16 = fits.Column(name='LYA_SCATTER_ERR_UP', format='D', disp='F12.4', unit='count', array=outtab['LYA_SCATTER_ERR_UP'])
            col17 = fits.Column(name='LYA_SCATTER_ERR_DOWN', format='D', disp='F12.4', unit='count', array=outtab['LYA_SCATTER_ERR_DOWN'])            
            binned_hdu_binary = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17], name='SCI', ver=1)
            # Create binary FITS table HDU with provenance information (HLSP requirement)
            # Because the output is for a single data set this table just has one row, but this results in consistent output for targets with single or multiple data sets
            # Output is saved to file named target + dataset + binning
            # For COS targets with a single dataset this is the primary output spectrum!
            col1 = fits.Column(name='PROPOSID', format='I4', disp='I8', unit='', array=[datasets[d]['PROPOSID']])
            col2 = fits.Column(name='FILE_ID', format='10A', disp='A10', unit='', array=[datasets[d]['FILE_ID']])
            col3 = fits.Column(name='TARGNAME', format='25A', disp='A25', unit='', array=[datasets[d]['TARGNAME']])
            col4 = fits.Column(name='DATE-BEG', format='19A', disp='A19', unit='', array=[datasets[d]['DATE-BEG']])
            col5 = fits.Column(name='DATE-END', format='19A', disp='A19', unit='', array=[datasets[d]['DATE-END']])
            col6 = fits.Column(name='MJD-BEG', format='D', disp='F12.5', unit='d', array=[datasets[d]['MJD-BEG']])
            col7 = fits.Column(name='MJD-END', format='D', disp='F12.5', unit='d', array=[datasets[d]['MJD-END']])
            col8 = fits.Column(name='LIFE_ADJ', format='I1', disp='I8', unit='', array=[datasets[d]['LIFE_ADJ']])
            col9 = fits.Column(name='APERTURE', format='8A', disp='A8', unit='', array=[datasets[d]['APERTURE']])
            col10 = fits.Column(name='DISPERSR', format='8A', disp='A8', unit='', array=[datasets[d]['DISPERSR']])
            col11 = fits.Column(name='CENWAVE', format='I2', disp='I8', unit='Angstrom', array=[datasets[d]['CENWAVE']])
            col12 = fits.Column(name='XPOSURE', format='E', disp='F9.3', unit='s', array=[datasets[d]['XPOSURE']])
            col13 = fits.Column(name='EXPA', format='I1', disp='I4', unit='', array=[datasets[d]['EXPA']])
            col14 = fits.Column(name='EXPB', format='I1', disp='I4', unit='', array=[datasets[d]['EXPB']])
            provenance = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14], name='PROVENANCE', ver=1)
            hdul = fits.HDUList([visit_hdu[d], binned_hdu_binary, provenance])
            saved_file = path_sci + targname + "_" + asn + "_" + \
                         str(BIN_PX) + "px_bin.fits"
            hdul.writeto(saved_file, overwrite=True)
            print(asn + " is binned and stored in " + saved_file + "\n")

        # If HLSP output is desired and there is just one data set, then copy the output file to a hlsp* file following the HLSP naming convention
        if HLSP_write == True and len(unique_datasets) == 1:
            hdul = fits.open(saved_file,memmap=True)
            h = hdul[0].header
            opt_elem = h['DISPERSR'].lower()
            hlsp_file = path_sci+"hlsp_"+HLSP_id.lower()+"_"+h['OBSERVAT'].lower()+"_"+h['INSTRUME'].lower()+"-"+h['DETECTOR'].lower()+"_"+h['HLSPTARG'].lower()+"_"+opt_elem+"_"+HLSP_ver.lower()+"_sci.fits"
            hdul.close()
            copy2(saved_file,hlsp_file)
            print("Binned spectrum has been copied to HLSP file\n"+hlsp_file)


    # Coadding routine for all data sets in the science directory using a common wavelength grid with constant dispersion
    # It works for different setups, but only for the same object and for similar resolving power (e.g. G130M+G160M is allowed)!
    # The resolving power varies with wavelength and with detector lifetime position, so set the wavelength bin size wisely
    # Recommendations for BIN_SIZE are in the FaintCOS config file
    # Coadd in this way only if corresponding switch is true and if there are multiple data sets, otherwise stick to output with integer binning factor (set BIN_PX wisely)
    if COADD_ALL_DATASETS == True and len(unique_datasets) > 1:
        # Coadding all individual exposures
        cdr_data = []
        cdr_hdu = []
        for f in cdr_files:
            tmp = fits.open(f)
            t = Table(tmp[1].data)
            h = tmp[0].header
            t = t[t['CALIB'] > 0]
            cdr_data.append(t)
            cdr_hdu.append(h)
        print("Coadding all datasets in the input directory ...")
        print("Binning: " + str(BIN_SIZE) + " Angstrom")
        
        # find the minimum and maximum wavelength of all exposures
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
                # error propagation (the neighboring pixels in the same exposure are highly covariant)
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
        # calculate total flux
        tot_flux = np.divide(tot_counts - tot_darks - tot_lya, \
                          tot_calib*tot_exptime*tot_flt, \
                          out=np.zeros_like(tot_counts - tot_darks), \
                          where=tot_calib*tot_exptime*tot_flt != 0)
        # calculate 1 sigma stat. uncertainties, see above descriptions for details
        # convert uncertainties from counts to flux
        # calculate total estimated background error
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
        flux_err_up = np.divide(cnt_err_up,\
                                tot_calib*tot_exptime*tot_flt,\
                                out = np.zeros_like(tot_bkg_err_up),\
                                where=tot_calib*tot_exptime*tot_flt != 0)
        flux_err_down = np.divide(cnt_err_down,\
                                tot_calib*tot_exptime*tot_flt,\
                                out = np.zeros_like(tot_bkg_err_down),\
                                where=tot_calib*tot_exptime*tot_flt != 0)

        # Store arrays in astropy table for further processing, sum the background components
        outtab = Table([tot_wavelength,tot_flux,flux_err_up,flux_err_down,np.rint(tot_counts),tot_darks+tot_lya,tot_bkg_err_up,tot_bkg_err_down,tot_darks,tot_dc_err,tot_exptime,tot_dq,tot_calib,tot_flt,tot_lya,tot_lya_err_up,tot_lya_err_down], names=('WAVELENGTH','FLUX','FLUX_ERR_UP','FLUX_ERR_DOWN','GCOUNTS','BACKGROUND','BKG_ERR_UP','BKG_ERR_DOWN','DARK_CURRENT','DARK_CURRENT_ERR','EXPTIME','DQ','CALIB','FLAT_CORR','LYA_SCATTER','LYA_SCATTER_ERR_UP','LYA_SCATTER_ERR_DOWN'))

        # Optionally trim detector edges (shortest and longest wavelengths outside active area) from the rebinned spectrum
        if TRIM_EDGE:
            tmpindex = np.where(outtab['DQ']==0)[0]
            imin = tmpindex[0]
            imax = tmpindex[len(tmpindex)-1]
            outtab = outtab[imin:imax+1]
        # Optionally restrict wavelength range of output spectrum, this is good for blue modes (e.g. G140L/800A) that include poorly calibrated low-sensitivity range <1100A
        if TRIM_WAVE:
            outtab = outtab[(outtab['WAVELENGTH']>=TRIM_MIN) & (outtab['WAVELENGTH']<=TRIM_MAX)]
            
        # Create output binary FITS table HDU from astropy table columns with specified data types, formatting, and units    
        col1 = fits.Column(name='WAVELENGTH', format='D', disp='F10.4', unit='Angstrom', array=outtab['WAVELENGTH'])
        col2 = fits.Column(name='FLUX', format='D', disp='E13.7', unit='erg s^-1 cm^-2 Angstrom^-1', array=outtab['FLUX'])
        col3 = fits.Column(name='FLUX_ERR_UP', format='D', disp='E13.7', unit='erg s^-1 cm^-2 Angstrom^-1', array=outtab['FLUX_ERR_UP'])
        col4 = fits.Column(name='FLUX_ERR_DOWN', format='D', disp='E13.7', unit='erg s^-1 cm^-2 Angstrom^-1', array=outtab['FLUX_ERR_DOWN'])
        col5 = fits.Column(name='GCOUNTS', format='I', disp='I8', unit='count', array=outtab['GCOUNTS'])
        col6 = fits.Column(name='BACKGROUND', format='D', disp='F12.4', unit='count', array=outtab['BACKGROUND'])
        col7 = fits.Column(name='BKG_ERR_UP', format='D', disp='F12.4', unit='count', array=outtab['BKG_ERR_UP'])
        col8 = fits.Column(name='BKG_ERR_DOWN', format='D', disp='F12.4', unit='count', array=outtab['BKG_ERR_DOWN'])
        col9 = fits.Column(name='DARK_CURRENT', format='D', disp='F12.4', unit='count', array=outtab['DARK_CURRENT'])
        col10 = fits.Column(name='DARK_CURRENT_ERR', format='D', disp='F15.4', unit='count', array=outtab['DARK_CURRENT_ERR'])      
        col11 = fits.Column(name='EXPTIME', format='D', disp='F10.3', unit='s', array=outtab['EXPTIME'])
        col12 = fits.Column(name='DQ', format='I', disp='I5', unit='', array=outtab['DQ'])
        col13 = fits.Column(name='CALIB', format='D', disp='E13.7', unit='count cm^2 Angstrom/erg',array=outtab['CALIB'])
        col14 = fits.Column(name='FLAT_CORR', format='D', disp='F9.5', unit='', array=outtab['FLAT_CORR'])
        col15 = fits.Column(name='LYA_SCATTER', format='D', disp='F12.4', unit='count', array=outtab['LYA_SCATTER'])
        col16 = fits.Column(name='LYA_SCATTER_ERR_UP', format='D', disp='F12.4', unit='count', array=outtab['LYA_SCATTER_ERR_UP'])
        col17 = fits.Column(name='LYA_SCATTER_ERR_DOWN', format='D', disp='F12.4', unit='count', array=outtab['LYA_SCATTER_ERR_DOWN'])
        # Create primary HDU with FITS header using HLSP-compliant keywords, store spectrum in binary FITS table
        hdu = createCoAddPrimaryHeader(cdr_hdu, tot_wavelength, norm_exptime)
        hdu.header['COMMENT'] = "Coadded rebinned spectrum from several data sets with proper calibration"
        hdu_binary = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17], name='SCI', ver=1)

        # Create binary FITS table HDU with provenance information (HLSP requirement)
        # Save primary HDU and both table HDUs to file named target_spectrum.fits
        # For COS targets with multiple datasets this is the primary output spectrum!
        col1 = fits.Column(name='PROPOSID', format='I4', disp='I8', unit='', array=datasets['PROPOSID'])
        col2 = fits.Column(name='FILE_ID', format='10A', disp='A10', unit='', array=datasets['FILE_ID'])
        col3 = fits.Column(name='TARGNAME', format='25A', disp='A25', unit='', array=datasets['TARGNAME'])
        col4 = fits.Column(name='DATE-BEG', format='19A', disp='A19', unit='', array=datasets['DATE-BEG'])
        col5 = fits.Column(name='DATE-END', format='19A', disp='A19', unit='', array=datasets['DATE-END'])
        col6 = fits.Column(name='MJD-BEG', format='D', disp='F12.5', unit='d', array=datasets['MJD-BEG'])
        col7 = fits.Column(name='MJD-END', format='D', disp='F12.5', unit='d', array=datasets['MJD-END'])
        col8 = fits.Column(name='LIFE_ADJ', format='I1', disp='I8', unit='', array=datasets['LIFE_ADJ'])
        col9 = fits.Column(name='APERTURE', format='8A', disp='A8', unit='', array=datasets['APERTURE'])
        col10 = fits.Column(name='DISPERSR', format='8A', disp='A8', unit='', array=datasets['DISPERSR'])
        col11 = fits.Column(name='CENWAVE', format='I2', disp='I8', unit='Angstrom', array=datasets['CENWAVE'])
        col12 = fits.Column(name='XPOSURE', format='E', disp='F9.3', unit='s', array=datasets['XPOSURE'])
        col13 = fits.Column(name='EXPA', format='I1', disp='I4', unit='', array=datasets['EXPA'])
        col14 = fits.Column(name='EXPB', format='I1', disp='I4', unit='', array=datasets['EXPB'])
        provenance = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14], name='PROVENANCE', ver=1)
        hdul = fits.HDUList([hdu, hdu_binary, provenance])
        saved_file = path_sci + targname + "_spectrum.fits"
        hdul.writeto(saved_file, overwrite=True)
        print("Coadded spectrum is stored in " + saved_file)

        # If HLSP output is desired and there are multiple data sets, then copy the output file to a hlsp* file following the HLSP naming convention
        if HLSP_write:
            hdul = fits.open(saved_file,memmap=True)
            h = hdul[0].header
            if h['DISPERSR']=='MULTI':
                opt_elem = "g130m-g160m"
            else:
                opt_elem = h['DISPERSR'].lower()
            hlsp_file = path_sci+"hlsp_"+HLSP_id.lower()+"_"+h['OBSERVAT'].lower()+"_"+h['INSTRUME'].lower()+"-"+h['DETECTOR'].lower()+"_"+h['HLSPTARG'].lower()+"_"+opt_elem+"_"+HLSP_ver.lower()+"_sci.fits"
            hdul.close()
            copy2(saved_file,hlsp_file)
            print("Coadded spectrum has been copied to HLSP file\n"+hlsp_file)
            
#    print("DONE!")
    
    
    
    
    
    
    
