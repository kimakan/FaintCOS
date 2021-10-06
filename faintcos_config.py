'''Parameters for the pre_calcos.py and post_calcos.py scripts. 

Author: Kirill Makan
'''


from astropy.table import Table



'''
-------------------------------------------------------------------------------
                PIPELINE PARAMETERS
-------------------------------------------------------------------------------
'''

# constants for the selection of darks
DARK_EXPSTART_INTERVAL = 30.# in days; time interval for the slection of the
                            # darks; EXPSTART(science) +/- DARK_EXPSTART_INTERVAL
MIN_DARKS = 5           # min. number of darks needed for the background model
KS_THRESHOLD = 0.03     # initial Kolmogorow-Smirnow statistic threshold
KS_STEP = 0.005         # add this value to KS_THRESHOLD if number of selected
                        # dark is below MIN_DARKS

# window for the running average of the dark current (must be an ODD number!)
BKG_AV = 501
# -----------------------------------------------------------------------------

# binning 
BIN_SIZE = 0.24       # in angstrom for the total co-add
BIN_PX = 7            # in pixels for the co-add per visit


# set a custom wavelength interval for the co-added spectrum
# works only for COADD_ALL_DATASETS = True
CUSTOM_INTERVAL = False
WAVE_MIN = 1100.
WAVE_MAX = 1700.


# Switches for binning and coadding routines
BIN_DATASET = False          # Bins every dataset in the working folder in BIN_PX
                            # pixels. 
                            
COADD_ALL_DATASETS = True   # Co-adds AND bins ALL datasets in the working folder
                            # in BIN_SIZE(in angstroms) bins. MAKE SURE THAT 
                            # EVERY EXPOSURE IN THE FOLDER IS FOR THE SAME 
                            # TARGET AND GRATING BEFORE YOU TURN THIS ON!!!



# bad regions near the geocoronal emission lines (Ly alpha, NI, OI)
BAD_REGIONS_G130M = [[1199.0, 1235.0], \
                     [1195.0, 1200.0], \
                     [1295.0, 1312.0]]    
BAD_REGIONS_G140L = [[1175.0, 1240.0], \
                     [1285.0, 1325.0]]  

''' Calculations of the confidence limits
Standard version of the code calculates confidence limits
according to Kraft et al. 1991
In order to use confidence limits from Feldman & Cousins 1998, you need to 
install additional module CustomConfLim (see README.md)
'''
FELDMAN_COUSINS = False






'''PHA limits for the PHACORR

Please check the PHA distribution of your scientific data to determine the 
limits. WARNING! It will change all _pha files in the 'lref' directory!!!
'''

PHA_G130M_FUVA = [2, 16]
PHA_G130M_FUVB = [2, 16]

PHA_G160M_FUVA = [2, 16]
PHA_G160M_FUVB = [2, 16]

PHA_G140L_FUVA = [2, 16]
PHA_G140L_FUVB = [2, 16]


'''Definitions of the extraction and background windows for different setups.

The default values for the extraction windows refer to the 95% of the total
light. Please change the values ONLY if you know exactly what you are doing!

The background windows were chosen empirically according to the positions of 
the extraction window and calibration spectrum. The default height of the
windows is 40px.


!!!WARNING!!! EVERY SETUP SHOULD HAVE ONLY ONE DEFINITION OF THE WINDOWS!!!!!


The format for the window parameters is following:
opt_elem.append(OPT_ELEM)
cen_wave.append(CENWAVE)
segment.append(SEGMENT)
ext_win.append([YFULL_MIN, YFULL_MAX])
bkg_win.append([YFULL_MIN1, YFULL_MAX1, YFULL_MIN2, YFULL_MAX2])
'''

lp = []
opt_elem = []
cen_wave = []
segment = []
ext_win = []
bkg_win = []

'''
-------------------------------------------------------------------------------
                          LIFETIME POSITION 1
-------------------------------------------------------------------------------
'''

# ---------------- OPTICAL ELEMENT G130M --------------------------------------

lp.append(1)
opt_elem.append("G130M")
cen_wave.append(1291)
segment.append("FUVA")
ext_win.append([475, 500])
bkg_win.append([403, 443, 528, 568])

lp.append(1)
opt_elem.append("G130M")
cen_wave.append(1291)
segment.append("FUVB")
ext_win.append([534, 557])
bkg_win.append([464, 504, 589, 629])

lp.append(1)
opt_elem.append("G130M")
cen_wave.append(1300)
segment.append("FUVA")
ext_win.append([476, 498])
bkg_win.append([403, 443, 528, 568])

lp.append(1)
opt_elem.append("G130M")
cen_wave.append(1300)
segment.append("FUVB")
ext_win.append([534, 555])
bkg_win.append([464, 504, 589, 629])

lp.append(1)
opt_elem.append("G130M")
cen_wave.append(1309)
segment.append("FUVA")
ext_win.append([475, 497])
bkg_win.append([402, 442, 527, 567])

lp.append(1)
opt_elem.append("G130M")
cen_wave.append(1309)
segment.append("FUVB")
ext_win.append([534, 554])
bkg_win.append([463, 503, 588, 628])

lp.append(1)
opt_elem.append("G130M")
cen_wave.append(1318)
segment.append("FUVA")
ext_win.append([476, 496])
bkg_win.append([401, 441, 526, 566])

lp.append(1)
opt_elem.append("G130M")
cen_wave.append(1318)
segment.append("FUVB")
ext_win.append([535, 553])
bkg_win.append([462, 502, 587, 627])

lp.append(1)
opt_elem.append("G130M")
cen_wave.append(1327)
segment.append("FUVA")
ext_win.append([475, 495])
bkg_win.append([401, 441, 526, 566])

lp.append(1)
opt_elem.append("G130M")
cen_wave.append(1327)
segment.append("FUVB")
ext_win.append([534, 553])
bkg_win.append([462, 502, 587, 627])


# ---------------- OPTICAL ELEMENT G140L --------------------------------------

lp.append(1)
opt_elem.append("G140L")
cen_wave.append(1105)
segment.append("FUVA")
ext_win.append([482, 507])
bkg_win.append([409, 449, 546, 586])

lp.append(1)
opt_elem.append("G140L")
cen_wave.append(1230)
segment.append("FUVA")
ext_win.append([483, 507])
bkg_win.append([410, 450, 547, 587])

lp.append(1)
opt_elem.append("G140L")
cen_wave.append(1230)
segment.append("FUVB")
ext_win.append([541, 568])
bkg_win.append([468, 508, 605, 645])


# ---------------- OPTICAL ELEMENT G160M --------------------------------------

lp.append(1)
opt_elem.append("G160M")
cen_wave.append(1577)
segment.append("FUVA")
ext_win.append([471, 491])
bkg_win.append([398, 438, 523, 563])

lp.append(1)
opt_elem.append("G160M")
cen_wave.append(1577)
segment.append("FUVB")
ext_win.append([530, 547])
bkg_win.append([457, 497, 582, 622])

lp.append(1)
opt_elem.append("G160M")
cen_wave.append(1589)
segment.append("FUVA")
ext_win.append([472, 491])
bkg_win.append([402, 442, 517, 557])

lp.append(1)
opt_elem.append("G160M")
cen_wave.append(1589)
segment.append("FUVB")
ext_win.append([530, 546])
bkg_win.append([461, 501, 576, 616])

lp.append(1)
opt_elem.append("G160M")
cen_wave.append(1600)
segment.append("FUVA")
ext_win.append([470, 490])
bkg_win.append([397, 437, 522, 562])

lp.append(1)
opt_elem.append("G160M")
cen_wave.append(1600)
segment.append("FUVB")
ext_win.append([529, 546])
bkg_win.append([456, 496, 581, 621])

lp.append(1)
opt_elem.append("G160M")
cen_wave.append(1611)
segment.append("FUVA")
ext_win.append([468, 489])
bkg_win.append([401, 441, 516, 556])

lp.append(1)
opt_elem.append("G160M")
cen_wave.append(1611)
segment.append("FUVB")
ext_win.append([529, 545])
bkg_win.append([455, 495, 580, 620])

lp.append(1)
opt_elem.append("G160M")
cen_wave.append(1623)
segment.append("FUVA")
ext_win.append([466, 489])
bkg_win.append([395, 435, 520, 560])

lp.append(1)
opt_elem.append("G160M")
cen_wave.append(1623)
segment.append("FUVB")
ext_win.append([528, 545])
bkg_win.append([456, 496, 581, 621])


'''
-------------------------------------------------------------------------------
                          LIFETIME POSITION 2
-------------------------------------------------------------------------------
'''

# ---------------- OPTICAL ELEMENT G130M --------------------------------------

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1055)
segment.append("FUVA")
ext_win.append([501, 551])
bkg_win.append([433, 473, 573, 613])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1055)
segment.append("FUVB")
ext_win.append([559, 612])
bkg_win.append([497, 537, 635, 675])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1096)
segment.append("FUVA")
ext_win.append([499, 553])
bkg_win.append([434, 474, 572, 612])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1096)
segment.append("FUVB")
ext_win.append([560, 612])
bkg_win.append([498, 538, 632, 672])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1222)
segment.append("FUVA")
ext_win.append([510, 541])
bkg_win.append([441, 481, 568, 608])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1222)
segment.append("FUVB")
ext_win.append([570, 600])
bkg_win.append([503, 543, 628, 668])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1291)
segment.append("FUVA")
ext_win.append([511, 534])
bkg_win.append([448, 488, 563, 603])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1291)
segment.append("FUVB")
ext_win.append([574, 596])
bkg_win.append([508, 548, 623, 663])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1300)
segment.append("FUVA")
ext_win.append([512, 534])
bkg_win.append([448, 488, 563, 603])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1300)
segment.append("FUVB")
ext_win.append([574, 594])
bkg_win.append([508, 548, 623, 663])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1309)
segment.append("FUVA")
ext_win.append([513, 532])
bkg_win.append([448, 488, 563, 603])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1309)
segment.append("FUVB")
ext_win.append([575, 594])
bkg_win.append([508, 548, 623, 663])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1318)
segment.append("FUVA")
ext_win.append([513, 532])
bkg_win.append([448, 488, 563, 603])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1318)
segment.append("FUVB")
ext_win.append([575, 593])
bkg_win.append([508, 548, 623, 663])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1327)
segment.append("FUVA")
ext_win.append([513, 531])
bkg_win.append([448, 488, 563, 603])

lp.append(2)
opt_elem.append("G130M")
cen_wave.append(1327)
segment.append("FUVB")
ext_win.append([576, 593])
bkg_win.append([508, 548, 623, 663])


# ---------------- OPTICAL ELEMENT G140L --------------------------------------

lp.append(2)
opt_elem.append("G140L")
cen_wave.append(1105)
segment.append("FUVA")
ext_win.append([522, 552])
bkg_win.append([451, 491, 583, 623])

lp.append(2)
opt_elem.append("G140L")
cen_wave.append(1280)
segment.append("FUVA")
ext_win.append([524, 548])
bkg_win.append([451, 491, 583, 623])


# ---------------- OPTICAL ELEMENT G160M --------------------------------------

lp.append(2)
opt_elem.append("G160M")
cen_wave.append(1577)
segment.append("FUVA")
ext_win.append([509, 526])
bkg_win.append([439, 479, 554, 594])

lp.append(2)
opt_elem.append("G160M")
cen_wave.append(1577)
segment.append("FUVB")
ext_win.append([570, 587])
bkg_win.append([499, 539, 614, 654])

lp.append(2)
opt_elem.append("G160M")
cen_wave.append(1600)
segment.append("FUVA")
ext_win.append([508, 527])
bkg_win.append([439, 479, 554, 594])

lp.append(2)
opt_elem.append("G160M")
cen_wave.append(1600)
segment.append("FUVB")
ext_win.append([570, 586])
bkg_win.append([499, 539, 614, 654])

lp.append(2)
opt_elem.append("G160M")
cen_wave.append(1623)
segment.append("FUVA")
ext_win.append([506, 529])
bkg_win.append([439, 479, 554, 594])

lp.append(2)
opt_elem.append("G160M")
cen_wave.append(1623)
segment.append("FUVB")
ext_win.append([569, 588])
bkg_win.append([499, 539, 614, 654])


'''
-------------------------------------------------------------------------------
                          LIFETIME POSITION 3
-------------------------------------------------------------------------------
'''

# ---------------- OPTICAL ELEMENT G130M --------------------------------------

lp.append(3)
opt_elem.append("G130M")
cen_wave.append(1222)
segment.append("FUVA")
ext_win.append([435, 470])
bkg_win.append([369, 409, 498, 538])

lp.append(3)
opt_elem.append("G130M")
cen_wave.append(1222)
segment.append("FUVB")
ext_win.append([497, 531])
bkg_win.append([430, 470, 559, 599])

lp.append(3)
opt_elem.append("G130M")
cen_wave.append(1291)
segment.append("FUVA")
ext_win.append([442, 466])
bkg_win.append([377, 417, 492, 532])

lp.append(3)
opt_elem.append("G130M")
cen_wave.append(1291)
segment.append("FUVB")
ext_win.append([503, 525])
bkg_win.append([437, 477, 552, 592])

lp.append(3)
opt_elem.append("G130M")
cen_wave.append(1300)
segment.append("FUVA")
ext_win.append([443, 465])
bkg_win.append([376, 416, 491, 531])

lp.append(3)
opt_elem.append("G130M")
cen_wave.append(1300)
segment.append("FUVB")
ext_win.append([503, 524])
bkg_win.append([437, 477, 552, 592])

lp.append(3)
opt_elem.append("G130M")
cen_wave.append(1309)
segment.append("FUVA")
ext_win.append([443, 463])
bkg_win.append([375, 415, 490, 530])

lp.append(3)
opt_elem.append("G130M")
cen_wave.append(1309)
segment.append("FUVB")
ext_win.append([503, 522])
bkg_win.append([436, 476, 551, 591])

lp.append(3)
opt_elem.append("G130M")
cen_wave.append(1318)
segment.append("FUVA")
ext_win.append([442, 462])
bkg_win.append([375, 415, 490, 530])

lp.append(3)
opt_elem.append("G130M")
cen_wave.append(1318)
segment.append("FUVB")
ext_win.append([503, 521])
bkg_win.append([435, 475, 550, 590])

lp.append(3)
opt_elem.append("G130M")
cen_wave.append(1327)
segment.append("FUVA")
ext_win.append([441, 461])
bkg_win.append([374, 414, 489, 529])

lp.append(3)
opt_elem.append("G130M")
cen_wave.append(1327)
segment.append("FUVB")
ext_win.append([503, 520])
bkg_win.append([434, 474, 549, 589])


# ---------------- OPTICAL ELEMENT G140L --------------------------------------

lp.append(3)
opt_elem.append("G140L")
cen_wave.append(1105)
segment.append("FUVA")
ext_win.append([445, 478])
bkg_win.append([374, 414, 505, 545])

lp.append(3)
opt_elem.append("G140L")
cen_wave.append(1280)
segment.append("FUVA")
ext_win.append([448, 477])
bkg_win.append([376, 416, 505, 545])


# ---------------- OPTICAL ELEMENT G160M --------------------------------------

lp.append(3)
opt_elem.append("G160M")
cen_wave.append(1577)
segment.append("FUVA")
ext_win.append([436, 457])
bkg_win.append([375, 415, 480, 520])

lp.append(3)
opt_elem.append("G160M")
cen_wave.append(1577)
segment.append("FUVB")
ext_win.append([497, 515])
bkg_win.append([429, 469, 544, 584])

lp.append(3)
opt_elem.append("G160M")
cen_wave.append(1589)
segment.append("FUVA")
ext_win.append([436, 456])
bkg_win.append([370, 410, 485, 525])

lp.append(3)
opt_elem.append("G160M")
cen_wave.append(1589)
segment.append("FUVB")
ext_win.append([498, 515])
bkg_win.append([429, 469, 544, 584])

lp.append(3)
opt_elem.append("G160M")
cen_wave.append(1600)
segment.append("FUVA")
ext_win.append([435, 458])
bkg_win.append([370, 410, 485, 525])

lp.append(3)
opt_elem.append("G160M")
cen_wave.append(1600)
segment.append("FUVB")
ext_win.append([498, 514])
bkg_win.append([429, 469, 544, 584])

lp.append(3)
opt_elem.append("G160M")
cen_wave.append(1611)
segment.append("FUVA")
ext_win.append([435, 457])
bkg_win.append([470, 410, 485, 525])

lp.append(3)
opt_elem.append("G160M")
cen_wave.append(1611)
segment.append("FUVB")
ext_win.append([499, 514])
bkg_win.append([429, 469, 544, 584])

lp.append(3)
opt_elem.append("G160M")
cen_wave.append(1623)
segment.append("FUVA")
ext_win.append([433, 455])
bkg_win.append([369, 409, 484, 524])

lp.append(3)
opt_elem.append("G160M")
cen_wave.append(1623)
segment.append("FUVB")
ext_win.append([497, 513])
bkg_win.append([427, 467, 542, 582])


'''
-------------------------------------------------------------------------------
                          LIFETIME POSITION 4
-------------------------------------------------------------------------------
'''

# ---------------- OPTICAL ELEMENT G130M --------------------------------------

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1222)
segment.append("FUVA")
ext_win.append([406, 443])
bkg_win.append([360, 400, 455, 495])

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1222)
segment.append("FUVB")
ext_win.append([468, 501])
bkg_win.append([420, 460, 515, 555])

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1223)
segment.append("FUVA")
ext_win.append([406, 441])
bkg_win.append([340, 380, 545, 585])

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1223)
segment.append("FUVB")
ext_win.append([469, 500])
bkg_win.append([403, 443, 605, 645])

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1291)
segment.append("FUVA")
ext_win.append([414, 439])
bkg_win.append([360, 400, 455, 495])

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1291)
segment.append("FUVB")
ext_win.append([474, 497])
bkg_win.append([420, 460, 515, 555])

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1300)
segment.append("FUVA")
ext_win.append([414, 436])
bkg_win.append([350, 390, 459, 499])

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1300)
segment.append("FUVB")
ext_win.append([474, 494])
bkg_win.append([411, 451, 520, 560])

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1309)
segment.append("FUVA")
ext_win.append([414, 435])
bkg_win.append([352, 392, 457, 497])

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1309)
segment.append("FUVB")
ext_win.append([474, 493])
bkg_win.append([412, 452, 517, 557])

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1318)
segment.append("FUVA")
ext_win.append([414, 434])
bkg_win.append([351, 391, 456, 496])

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1318)
segment.append("FUVB")
ext_win.append([475, 492])
bkg_win.append([412, 452, 517, 557])

lp.append(4)
opt_elem.append("G130M")
cen_wave.append(1327)
segment.append("FUVA")
ext_win.append([413, 435])
bkg_win.append([360, 400, 455, 495])


# ---------------- OPTICAL ELEMENT G140L --------------------------------------

lp.append(4)
opt_elem.append("G140L")
cen_wave.append(800)
segment.append("FUVA")
ext_win.append([419, 444])
bkg_win.append([344, 384, 473, 513])

lp.append(4)
opt_elem.append("G140L")
cen_wave.append(1105)
segment.append("FUVA")
ext_win.append([423, 447])
bkg_win.append([346, 386, 567, 607])

lp.append(4)
opt_elem.append("G140L")
cen_wave.append(1280)
segment.append("FUVA")
ext_win.append([425, 448])
bkg_win.append([361, 401, 470, 510])


# ---------------- OPTICAL ELEMENT G160M --------------------------------------

lp.append(4)
opt_elem.append("G160M")
cen_wave.append(1533)
segment.append("FUVA")
ext_win.append([403, 428])
bkg_win.append([340, 380, 460, 500])

lp.append(4)
opt_elem.append("G160M")
cen_wave.append(1533)
segment.append("FUVB")
ext_win.append([465, 486])
bkg_win.append([399, 439, 512, 552])

lp.append(4)
opt_elem.append("G160M")
cen_wave.append(1577)
segment.append("FUVA")
ext_win.append([408, 427])
bkg_win.append([342, 382, 451, 491])

lp.append(4)
opt_elem.append("G160M")
cen_wave.append(1577)
segment.append("FUVB")
ext_win.append([469, 486])
bkg_win.append([405, 445, 510, 550])

lp.append(4)
opt_elem.append("G160M")
cen_wave.append(1589)
segment.append("FUVA")
ext_win.append([407, 427])
bkg_win.append([345, 385, 450, 490])

lp.append(4)
opt_elem.append("G160M")
cen_wave.append(1589)
segment.append("FUVB")
ext_win.append([470, 485])
bkg_win.append([405, 445, 510, 550])

lp.append(4)
opt_elem.append("G160M")
cen_wave.append(1600)
segment.append("FUVA")
ext_win.append([406, 426])
bkg_win.append([347, 387, 451, 491])

lp.append(4)
opt_elem.append("G160M")
cen_wave.append(1600)
segment.append("FUVB")
ext_win.append([469, 484])
bkg_win.append([404, 444, 509, 549])

lp.append(4)
opt_elem.append("G160M")
cen_wave.append(1611)
segment.append("FUVA")
ext_win.append([404, 425])
bkg_win.append([345, 385, 451, 491])

lp.append(4)
opt_elem.append("G160M")
cen_wave.append(1611)
segment.append("FUVB")
ext_win.append([469, 484])
bkg_win.append([404, 444, 509, 549])

lp.append(4)
opt_elem.append("G160M")
cen_wave.append(1623)
segment.append("FUVA")
ext_win.append([403, 425])
bkg_win.append([341, 381, 450, 490])

lp.append(4)
opt_elem.append("G160M")
cen_wave.append(1623)
segment.append("FUVB")
ext_win.append([468, 484])
bkg_win.append([403, 443, 508, 548])



# define a table with the extraction and background windows for the easy access
custom_xtractab = Table(names=("INDEX", "LP", "OPT_ELEM", "SEGMENT", \
                               "CENWAVE","B_SPEC", "HEIGHT", "B_BKG1", \
                               "B_BKG2","B_HGT1", "B_HGT2"), \
                        dtype=('i4', 'i4', 'S7', 'S5', 'i4', 'f4', 'f4',\
                               'f4', 'f4', 'f4', 'f4'))
for i in range(len(opt_elem)):
    custom_xtractab.add_row([i, lp[i], opt_elem[i], \
                            segment[i], int(cen_wave[i]),
                            (ext_win[i][0] + ext_win[i][1])/2.,\
                            ext_win[i][1] - ext_win[i][0],\
                            (bkg_win[i][0] + bkg_win[i][1])/2.,\
                            (bkg_win[i][2] + bkg_win[i][3])/2.,\
                            bkg_win[i][1] - bkg_win[i][0],\
                            bkg_win[i][3] - bkg_win[i][2]])


# define a table with pha limits for the easy access
custom_pha = Table(names=("OPT_ELEM", "SEGMENT", "LLT", "ULT"),\
                   dtype=('S7', 'S5', 'i4', 'i4'))
custom_pha.add_row(["G130M", "FUVA", PHA_G130M_FUVA[0], PHA_G130M_FUVA[1]])
custom_pha.add_row(["G130M", "FUVB", PHA_G130M_FUVB[0], PHA_G130M_FUVB[1]])
custom_pha.add_row(["G160M", "FUVA", PHA_G160M_FUVA[0], PHA_G160M_FUVA[1]])
custom_pha.add_row(["G160M", "FUVB", PHA_G160M_FUVB[0], PHA_G160M_FUVB[1]])
custom_pha.add_row(["G140L", "FUVA", PHA_G140L_FUVA[0], PHA_G140L_FUVA[1]])
custom_pha.add_row(["G140L", "FUVB", PHA_G140L_FUVB[0], PHA_G140L_FUVB[1]])




# correct HST/COS resolution at different LP according to the LSFs
 
cos_res = Table(names=("LP", "OPT_ELEM", "CENWAVE", "R"), \
                dtype=('i4', 'S5', 'i4', 'i4'))

# spectral resolution for LP 1
cos_res.add_row([1, 'G130M', 1222, 13000])
cos_res.add_row([1, 'G130M', 1291, 19000])
cos_res.add_row([1, 'G130M', 1300, 19000])
cos_res.add_row([1, 'G130M', 1309, 19000])
cos_res.add_row([1, 'G130M', 1318, 18500])
cos_res.add_row([1, 'G130M', 1327, 18500])

cos_res.add_row([1, 'G160M', 1577, 19500])
cos_res.add_row([1, 'G160M', 1589, 19500])
cos_res.add_row([1, 'G160M', 1600, 20000])
cos_res.add_row([1, 'G160M', 1611, 20000])
cos_res.add_row([1, 'G160M', 1577, 20000])

cos_res.add_row([1, 'G140L', 1105, 1900])
cos_res.add_row([1, 'G140L', 1230, 2300])
cos_res.add_row([1, 'G140L', 1280, 2400])


# spectral resolution for LP 2
cos_res.add_row([2, 'G130M', 1055, 4500])
cos_res.add_row([2, 'G130M', 1096, 7500])
cos_res.add_row([2, 'G130M', 1222, 13500])
cos_res.add_row([2, 'G130M', 1291, 16500])
cos_res.add_row([2, 'G130M', 1300, 16500])
cos_res.add_row([2, 'G130M', 1309, 16500])
cos_res.add_row([2, 'G130M', 1318, 16500])
cos_res.add_row([2, 'G130M', 1327, 17000])

cos_res.add_row([2, 'G160M', 1577, 18000])
cos_res.add_row([2, 'G160M', 1589, 18000])
cos_res.add_row([2, 'G160M', 1600, 18000])
cos_res.add_row([2, 'G160M', 1611, 18000])
cos_res.add_row([2, 'G160M', 1623, 18000])

cos_res.add_row([2, 'G140L', 1105, 1300])
cos_res.add_row([2, 'G140L', 1280, 1700])


# spectral resolution for LP 3
cos_res.add_row([3, 'G130M', 1222, 12000])
cos_res.add_row([3, 'G130M', 1291, 16500])
cos_res.add_row([3, 'G130M', 1300, 16500])
cos_res.add_row([3, 'G130M', 1309, 16500])
cos_res.add_row([3, 'G130M', 1318, 16500])
cos_res.add_row([3, 'G130M', 1327, 16500])

cos_res.add_row([3, 'G160M', 1577, 18300])
cos_res.add_row([3, 'G160M', 1589, 18000])
cos_res.add_row([3, 'G160M', 1600, 17700])
cos_res.add_row([3, 'G160M', 1611, 17300])
cos_res.add_row([3, 'G160M', 1623, 17000])

cos_res.add_row([3, 'G140L', 1105, 1500])
cos_res.add_row([3, 'G140L', 1280, 1900])


# spectral resolution for LP 4
cos_res.add_row([4, 'G130M', 1222, 13500])
cos_res.add_row([4, 'G130M', 1291, 14500])
cos_res.add_row([4, 'G130M', 1300, 15000])
cos_res.add_row([4, 'G130M', 1309, 15000])
cos_res.add_row([4, 'G130M', 1318, 15000])
cos_res.add_row([4, 'G130M', 1327, 15000])

cos_res.add_row([4, 'G160M', 1533, 16000])
cos_res.add_row([4, 'G160M', 1577, 16500])
cos_res.add_row([4, 'G160M', 1589, 16500])
cos_res.add_row([4, 'G160M', 1600, 16500])
cos_res.add_row([4, 'G160M', 1611, 16000])
cos_res.add_row([4, 'G160M', 1623, 16000])

cos_res.add_row([4, 'G140L', 800, 700])
cos_res.add_row([4, 'G140L', 1105, 900])
cos_res.add_row([4, 'G140L', 1280, 1000])

