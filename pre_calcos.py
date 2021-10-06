'''This script changes the header of '_rawtag' files to prepare them for the
CALCOS pipeline. Additionally, it adapts the extraction windows in the _1dx 
files and PHA limits in the _pha files. 

Author: Kirill Makan
'''

import os
import sys

from astropy.io import fits
from astropy.table import Table

from faintcos_config import *




def get_files():
    '''Get the directories to rawtag and reference files.
     
    The path to the reference files should be stored in the 
    environment variable "lref". 
    
    Returns: 
    -----------
    path_rawtag: a list with full path to each rawtag file
    
    path_ref: path to the reference files directory
    '''
      
    path_sci = "."
    
    # path to science was given in the command line?
    if len(sys.argv) > 1:
        path_sci = sys.argv[1]
    
    # find all rawtag files in the directory
    path_rawtag = [f for f in os.listdir(path_sci) if "rawtag" in f]
    if path_sci != ".":
        path_rawtag = [path_sci + s for s in path_rawtag]
    print("Found " + str(len(path_rawtag)) + " 'rawtag' files.")
    
    # find reference directory
    try:
        path_ref = os.environ['lref']
    except KeyError:
        print("ERROR: lref is not defined!")
        sys.exit()   
    
    print("Reference files: " + path_ref)
    
    return path_rawtag, path_ref
    

def change_rawtag(path_rawtag, path_ref):
    ''' Adapts the CALIBRATION SWITCHES in the header of the rawtag files. 
    
    Also changes references to some of the calibration files, i.e. 
    gain sag table. Also, checks whether all reference files are available.
    
    Parameters:
    -----------
    path_rawtag: list with the full path to the rawtag files
    
    Returns:
    -----------
    success: 'True' - if successfuly changed rawtag and all reference files
    are available, 'False' otherwise
    '''
    
    success = True
    
    # keywords for the reference filesaccess from your network is temporarily disabled 
    ref_keywords = \
        ['FLATFILE', 'DEADTAB', 'BPIXTAB', 'SPOTTAB', 'GSAGTAB',\
        'HVTAB', 'BRFTAB', 'GEOFILE', 'DGEOFILE', 'TRACETAB', 'PROFTAB',\
        'TWOZXTAB', 'XWLKFILE', 'YWLKFILE', 'PHATAB', 'PHAFILE', 'BADTTAB',\
        'XTRACTAB', 'LAMPTAB', 'DISPTAB', 'IMPHTTAB', 'FLUXTAB', 'WCPTAB',\
        'BRSTTAB', 'TDSTAB', 'SPWCSTAB']
        
    ref_files = os.listdir(path_ref)
    
    # change the calibration switches and reference files
    print("Changing the header of the rawtag files...", end="", flush=True)
    for f in path_rawtag:
        hdr = fits.open(f)[0].header
        # Flatfield correction only for the science data
        if hdr['EXPTYPE'].strip() == 'DARK':
            fits.setval(f, 'FLATCORR', value='OMIT')
        else:
            fits.setval(f, 'FLATCORR', value='PERFORM')
        fits.setval(f, 'RANDSEED', value=0)
        fits.setval(f, 'BACKCORR', value='OMIT')
        fits.setval(f, 'TRCECORR', value='OMIT')
        fits.setval(f, 'ALGNCORR', value='OMIT')
        fits.setval(f, 'XTRCTALG', value='BOXCAR')
        fits.setval(f, 'SPOTTAB', value='N/A')
        #fits.setval(f, 'GSAGTAB', value='lref$'+GSAG_FILE) 
        if fits.open(f)[0].header['EXPTYPE'].strip() == 'DARK':
            # turns off the PHA filtering for darkframes
            fits.setval(f, 'PHACORR', value='OMIT')
    print("OK")  
    
    # check availability of the reference files defined in the header
    print("Checking the availability of the reference files.")
    for f in path_rawtag:
        print(f)
        hdr = fits.open(f)[0].header
        for k in ref_keywords:
            if hdr[k].split('$')[-1] not in ref_files and \
               hdr[k].split('$')[-1] != 'N/A':
                print(str(hdr[k].split('$')[-1]) + " is missing!")
                success = False
    if success:
        print("Reference files: OK")
    
         
    return success

def set_aperture(path_rawtag, path_ref):
    ''' Sets the aperture windows according to 'custom_xtractab'.'''
    
    
    print("Writing new XTRACTAB data into the _1dx file.")
    print("Following data was written:")
    changed_setups = [-1]
    for f in path_rawtag:
        hdr = fits.open(f)[0].header
        
        # check whether it's a darkframe (change of the _1dx file is not
        # possible for the darkframe)
        if hdr['EXPTYPE'].strip() == 'DARK':
            continue
        segment = hdr['SEGMENT']
        life_adj = hdr['LIFE_ADJ']
        opt_elem = hdr['OPT_ELEM']
        cenwave = hdr['CENWAVE']
        aperture = 'PSA'
        path_1dx = path_ref + hdr['XTRACTAB'].split('$')[-1]
        
        # select the HST/COS setup according to the rawtag header
        custom_1dx = custom_xtractab[(custom_xtractab['LP'] == life_adj) & \
                                     (custom_xtractab['SEGMENT'] == segment) &\
                                   (custom_xtractab['OPT_ELEM'] == opt_elem) &\
                                   (custom_xtractab['CENWAVE'] == cenwave)]
        
        # Did we already change the file for this setup?
        if len(custom_1dx['INDEX']) > 1:
            print("SETUP WITH LP=" + str(life_adj) + \
                  ", OPT_ELEM=" + str(opt_elem) + \
                  ", SEGMENT=" + str(segment) + \
                  "and CENWAVE=" + str(cenwave) + \
                  " has multiple defintions for windows!")
            print("Every setup should have just one definition!")
            sys.exit()
        if len(custom_1dx['INDEX']) < 1:
            print("WARNING: " + str(hdr['FILENAME']))
            print("No definition of the windows for the following setup:")
            print("LP=" + str(life_adj) + \
                  ", OPT_ELEM=" + str(opt_elem) + \
                  ", SEGMENT=" + str(segment) + \
                  ", CENWAVE=" + str(cenwave))
            print("The windows for this setup will not be changed!!!")
            print("\n")
            continue
        if int(custom_1dx['INDEX']) in changed_setups:
            continue
        changed_setups.append(int(custom_1dx['INDEX']))
        
        # find the row with the selected setup in the 1dx to change it later
        row = 0
        data_1dx = fits.getdata(path_1dx) 
        for i in range(len(data_1dx)):
            if (data_1dx[i][0] == segment) & \
               (data_1dx[i][1] == opt_elem) & \
               (data_1dx[i][2] == cenwave) & \
               (data_1dx[i][3] == aperture):
                row = i
                break
        
        with fits.open(path_1dx, mode='update') as hdul:
            hdul[1].data[row]['B_SPEC'] = float(custom_1dx['B_SPEC'])
            hdul[1].data[row]['HEIGHT'] = float(custom_1dx['HEIGHT'])  
            hdul[1].data[row]['SLOPE'] = 0.0
            hdul[1].data[row]['B_BKG1'] = float(custom_1dx['B_BKG1']) 
            hdul[1].data[row]['B_BKG2'] = float(custom_1dx['B_BKG2'])
            if 'B_HGT1' in Table(data_1dx).columns: 
                hdul[1].data[row]['B_HGT1'] = float(custom_1dx['B_HGT1']) 
                hdul[1].data[row]['B_HGT2'] = float(custom_1dx['B_HGT2']) 
            else:
                # old 1dx files have same height for both background windows
                hdul[1].data[row]['BHEIGHT'] = float(custom_1dx['B_HGT1'])
            hdul.flush()
        data_1dx = fits.getdata(path_1dx)
        print(data_1dx[row])


def set_pha_limits(path_ref):
    ''' Changes the PHA limits in the _pha file.'''
    
    pha_files = [f for f in os.listdir(path_ref) if "_pha.fits" in f]
    print("Updating PHA limits...", end="", flush=True)
    for f in pha_files:
        pha_data = fits.getdata(path_ref + f)
        for i in range(len(pha_data)):
            segment = pha_data[i]['SEGMENT']
            opt_elem = pha_data[i]['OPT_ELEM']
            with fits.open(path_ref + f, mode='update') as hdul:
                custom_data = custom_pha[(custom_pha['SEGMENT'] == segment) &\
                                         (custom_pha['OPT_ELEM'] == opt_elem)]
                hdul[1].data[i]['LLT'] = custom_data['LLT']
                hdul[1].data[i]['ULT'] = custom_data['ULT']
                hdul.flush()
    print("DONE")        
    print("Updated following files: " + str(pha_files))            
        
if __name__ == "__main__": 
    path_rawtag, path_ref = get_files()
    change_rawtag(path_rawtag, path_ref)
    set_aperture(path_rawtag, path_ref)
    set_pha_limits(path_ref)
    print("Now you can run CALCOS!")

           
