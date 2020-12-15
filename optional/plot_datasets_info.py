''' This code creates the 2d spectra of the stacked corrtags for each dataset.

The code stacks corrtags in the directory (working directory) given by the user 
as a command line argument. The corrtags are stucked for each dataset and
segment. 

Additionaly, it creates PHA histograms for counts in the extraction box.

The output files are put in the "datasets_info/" folder in the 
working directory. The "datasets.txt" file contain a table of all corrtags
in the working directory sorted by their ASN_ID (Dataset ID). 
The "ASN_ID_SEGMENT.pdf" files contain 2d spectra for counts with DQ == 0.
The used extraction box indicated by the horizontal dashed lines. 
It also contains the PHA histogram for all counts with DQ==0 and DQ==512 
in the extraction box. 

Author: Kirill Makan
'''


import os
import sys

from astropy.io import fits
from astropy.table import Table, vstack
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np




path_to_data = sys.argv[1]

DARKS = False

targname = ""

corrtags = [f for f in os.listdir(path_to_data) if "corrtag" in f]
corrtags.sort()


corrtags_t = Table(names=("ASN_ID", "CORRTAG_FILE", "START_TIME",\
                                      "OPT_ELEM", "SEGMENT", "CENTRWV",\
                                      "EXP_TIME", "EXPFLAG", "LP"), \
                               dtype=('S20', 'S30', 'S23', 'S5' , 'S5',\
                                      'i4', 'f4', 'S10', 'i4'))
                        
for f in corrtags:
    hdr = fits.open(path_to_data + f)[0].header
    targname = hdr['TARGNAME']
    if targname == "DARK":
        if "corrtag_b.fits" in f:
            continue
        DARKS = True
    hdr1 = fits.open(path_to_data + f)[1].header
    corrtags_t.add_row([hdr['ASN_ID'], f, \
                    hdr1['DATE-OBS'] + " " + hdr1['TIME-OBS'],\
                    hdr['OPT_ELEM'], hdr['SEGMENT'], \
                    hdr['CENTRWV'], hdr1['EXPTIME'], \
                    hdr1['EXPFLAG'], hdr['LIFE_ADJ']])

if DARKS:
    corrtags_t.sort(['START_TIME'])
else:
    corrtags_t.sort(['ASN_ID'])


# create the directory for the output
os.system("mkdir -p " + path_to_data + "datasets_info")

# create a list for all corrtag files in the working directory
corrtags_t.meta['comments'] = {"OBJECT: " + targname}
corrtags_t.write(path_to_data + "datasets_info/"+ "datasets.txt", \
                format='ascii.fixed_width', overwrite=True)


if not DARKS:
    # MAKE PLOTS FOR SCIENCE DATA SETS ONLY!!!
    asn = np.unique(corrtags_t['ASN_ID'])
    datasets = []
    # create subtables for each dataset and segment
    for a in asn:
        datasets.append(corrtags_t[(corrtags_t['ASN_ID'] == a) & \
            (corrtags_t['SEGMENT'] == 'FUVA')])
        datasets.append(corrtags_t[(corrtags_t['ASN_ID'] == a) & \
            (corrtags_t['SEGMENT'] == 'FUVB')])
    for d in datasets:
        dataset = d[(d['EXPFLAG'] == 'NORMAL')]
        # if there are no good exposures in this dataset then skip
        if len(dataset['ASN_ID']) == 0:
            continue
        
        data = Table(fits.getdata(path_to_data + dataset['CORRTAG_FILE'][0]))
        
        # read the corrtag headers and x1d data
        hdr = fits.open(path_to_data + dataset['CORRTAG_FILE'][0])[0].header
        hdr1 = fits.open(path_to_data + dataset['CORRTAG_FILE'][0])[1].header
        x1d_data = Table(fits.getdata(path_to_data + \
            dataset['CORRTAG_FILE'][0].split("_")[0] + "_x1d.fits"))
        x1d_data = x1d_data[x1d_data['SEGMENT'] == hdr['SEGMENT']]
        
 
        if len(x1d_data['Y_UPPER_OUTER']) == 0:
            bad_data = True
        else:
            bad_data = False
            y_max = x1d_data[x1d_data['SEGMENT'] == \
                hdr['SEGMENT']]['Y_UPPER_OUTER'][0][0]
            y_min = x1d_data[x1d_data['SEGMENT'] == \
                hdr['SEGMENT']]['Y_LOWER_OUTER'][0][0]
            
        # stack all corrtags of the current data set / ASN_ID
        for i in range(1, len(dataset['ASN_ID'])):
            tmp_data = Table(fits.getdata(path_to_data + \
                dataset['CORRTAG_FILE'][i]))
            data = vstack([data, tmp_data])

        
        plt.figure(figsize=(10, 2.5))
        fig = plt.gcf()
        fig.suptitle(targname + " (" + dataset['ASN_ID'][0] + ", " + \
            hdr['SEGMENT'] + ", " + str(int(hdr['CENTRWV'])) + \
            ", EXPTIME = " + str(np.sum(dataset['EXP_TIME'])) + \
            ", EXPOSURES = " + str(len(dataset['ASN_ID'])) + \
            ", LP = " + str(hdr['LIFE_ADJ']) + ")", fontsize=9)
        
        gs = gridspec.GridSpec(1, 4, wspace=0.25)
        ax = plt.subplot(gs[0])
        
        # exclude the geocoronal Lya and OI from the PHA histogram
        pha = data[(data['YFULL'] >= y_min) & (data['YFULL'] <= y_max) & \
                  ((data['DQ'] == 0) | (data['DQ'] == 512))]['PHA']
        if not bad_data:
            plt.hist(pha, bins = 32, range = (0, 32), color = 'black', \
                alpha = 0.3, align = 'left', density = 1)
        plt.xlabel("PHA", fontsize = 6)
        plt.ylabel("PDF", fontsize = 6)
        plt.xticks(fontsize = 7)
        plt.yticks(fontsize = 7)
        ax.xaxis.set_label_coords(0.5, -0.15)
        plt.tick_params(axis='x', which='both', direction = 'in')
        plt.tick_params(axis='y', which='both', direction = 'in')
        plt.minorticks_on()


        ax = plt.subplot(gs[1:])
        data = data[data['DQ'] == 0]
        plt.scatter(data['WAVELENGTH'], data['YFULL'], s = 1.5, \
            color = 'black', alpha = 0.015, rasterized=True)
        plt.minorticks_on()
        if not bad_data:
            plt.axhline(y = y_max, color = 'red', linewidth = 0.8, \
                linestyle = '--', alpha = 0.8)
            plt.axhline(y = y_min, color = 'red', linewidth = 0.8, \
                linestyle = '--', alpha = 0.8)
            plt.ylim(y_min - (60 - (y_max - y_min)/2), \
                y_max + (60 - (y_max - y_min)/2))
        else:
            plt.ylim(250, 800)
        plt.xlabel("WAVELENGTH", fontsize = 6)
        plt.ylabel("YFULL", fontsize = 6)
        plt.xlim(min(x1d_data['WAVELENGTH'][0]), max(x1d_data['WAVELENGTH'][0]))
        plt.tick_params(axis='x', which='both', direction = 'in')
        plt.tick_params(axis='y', which='both', direction = 'in')
        plt.xticks(fontsize = 7)
        plt.yticks(fontsize = 7)
        
        plt.subplots_adjust(left=0.05, right=0.99, top=0.91, bottom=0.15)

        plt.savefig(path_to_data + "datasets_info/"+ dataset['ASN_ID'][0] + \
            "_" + hdr['SEGMENT'] + ".pdf", format='pdf')
        plt.close()
