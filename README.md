FaintCOS: Improved background subtraction and co-adding code for CALCOS

Authors: Kirill Makan, Gabor Worseck

We used the following version of the software
 - python 3.7.3 with astropy 4.0.1, scipy 1.3.1, numpy 1.17.3
 - CALCOS 3.3.9
 - GCC 7.4.0 (required for the calculation of the Feldman & Cousins uncertainties, see Section 5)

There is no guarantee that FaintCOS will work properly with other versions.

We caution that the default reduction parameters, such as primary science apertures (PSA) and 
pulse height amplitude (PHA) limits, are optimized for faint point sources. For the extended 
sources, you most likely have to adapt these parameters. Please, use "plot_datasets_info.py" to 
check 2D spectra and PHA distributions (see Sections 1 and 6).

CONTENT:

1. INSTRUCTIONS
2. PRODUCED FILES
3. FIT TABLE COLUMNS
4. REDUCTION PARAMETERS
5. CALCULATED UNCERTAINTIES
6. OPTIONAL CODE
-------------------------------------------------------------------------------------
1. INSTRUCTIONS
-------------------------------------------------------------------------------------
 - Download uncalibrated darkframes and associated calibration/reference files
    from https://archive.stsci.edu/hst/search.php
    For this, make the following changes in the standard form:
    
        Target Name: "DARK" 
        Resolver: Don't Resolve
        Select Imagers: COS
        Start Time: "> YYYY-MM-DD" (earliest science data start 
        time - DARK_EXPSTART_INTERVAL, see Section 4).
        Uncheck "Science" and select "Calibration"
        User-specified field 1: select "Start Time" and write "< YYYY-MM-DD" (latest science 
        data start time + DARK_EXPSTART_INTERVAL, see Section 4).
        User-specified field 2: select "Instrument Config" and write "COS/FUV"
        
    Search for the darkframes by clicking "Search".
    All available darks should be now listed in a table. Click "Mark all" and then
    "Submit marked data for retrieval from STDADS". A new window will open and show 
    the retrieval form. Uncheck "Calibrated" and select "Used Reference Files" and 
    "Uncalibrated". Send the retrieval request anonymously or, if you are registered,
    with your STScI ID. 
 - Download uncalibrated science data and associated calibration/reference files
    from https://archive.stsci.edu/hst/search.php
    Make sure to uncheck "Calibrated" and select "Used Reference Files" and 
    "Uncalibrated" in the retrieval form after you have selected the datasets.
 - Put all calibration/reference files from darkframes AND science
    data in one folder (e.g. "calib"). You can easely spot the reference files, 
    because they do not contain a dataset ID in their file name. 
 - Define 'lref' as described in "COS Data Handbook (Version 4.0)" Section 3.6.1 
 
    example:
             export lref="/home/user/Documents/calib/"
             
 - Create a new folder for scientific data for every object and resolution (low and medium 
   resolution data should be ALWAYS in different folders)
 - Check the defined PSA and PHA limits in the faintcos_config.py (add new
    definitions if necessary). Later (after CALCOS), you can run "plot_datasets_info.py" to
    check the PSA position and PHA distribution.
 - Run "pre_calcos.py" with the path to the uncalibrated darkframes (rawtag files) as
    an argument 
    
    example:
              python pre_calcos.py /home/user/Documents/darkframes/
              
   !!!WARNING!!! ALL CODES SHOULD RUN IN THE SAME AstroConda ENVIRONMENT AS CALCOS!!! 
 - Copy "exec_calcos_for_darks.py" into the folder with the uncalibrated darkframes
    and then run it. It will reduce all darkframes and put reduced files in the 
    folder "reduced". (This process can take a while)
 - Put the reduced darkframes files (corrtag files) in a separate folder and define
    'ldark' = path_to_reduced_darkframes, as you defined 'lref'. 
 - Run 'pre_calcos.py" with the path to the folder with uncalibrated data for
    a single target and resolution (M or L).
    
    example:
              python pre_calcos.py /home/user/Documents/data/QSO-231145-141752/
              
 - Reduce the scientific data with CALCOS
 - Set the switches (binning, co-adding, wavelength range etc.) in the 'faintcos_config.py'
    (see Section 4 below)
 - Run 'post_calcos.py' with the path to the REDUCED science data (corrtags and 1dx
    files) as an command-line argument (it is the working folder for the pipeline).
    
    example: 
              python post_calcos.py /home/user/Documents/data/reduced_QSO-231145-141752/
              
 - All visits and exposures will be listed after you run post_calcos.py. If you 
    want to reduce the listed data, then proceed with the reduction.
 - After post_calcos.py is done, the reduced files appear in the working folder
    

    
The default version of the code calculates the statistical count errors according
to Kraft et al. 1991. 
If you wish to use the more correct Feldman & Cousins 1998 confidence limits instead:

 - Install CustomConfLim module by running in the CustomConfLimits folder: 

       "python setup.py build"
       "python setup.py install" 
       
 - Set FELDMAN_COUSINS=True in the faintcos_config.py file

   
   
-------------------------------------------------------------------------------------    
2. PRODUCED FILES
-------------------------------------------------------------------------------------
"DATASET_dataset_sum.fits"
Co-added spectrum for a single data set. Similar to the standard CALCOS x1dsum.fits
file but with the improved data reduction.

"TARGETNAME_DATASET_Npx_binned.fits":
Binned spectrum of a single dataset in N pixels bins. Number of binned pixels 
is defined in post_calcos.py with BIN_PX. This file will only be produced if the 
switch "BIN_VISIT" is set to "True". If there are several data sets in the working
folder then the pipeline produces new file for every data set. 


"TARGETNAME_spectrum.fits":
Co-added and binned spectrum of all data sets in the working folder in BIN_SIZE bins.
The width of the bins in angstroms is defined in post_calcos.py with BIN_SIZE. This 
file will only be produced if the switch "COADD_ALL_VISITS" is set to "True".
WARNING!!! Make sure that every data set in the working folder is for the same
target and resolution (M or L)!!!


------------------------------------------------------------------------------------
3. DESCRIPTION OF THE COLUMNS
------------------------------------------------------------------------------------

|Column Name            |Units              |Desctiption                            |
|:----------------------|:------------------|:--------------------------------------|
|WAVELENGTH             |angstrom           |Wavelength scale                       |
|FLUX                   |ergs/s/cm^2/A      |Calibrated flux                        |
|FLUX_ERR_UP            |ergs/s/cm^2/A      |Upper error estimate                   |
|FLUX_ERR_DOWN          |ergs/s/cm^2/A      |Lower error estimate                   |
|GCOUNTS                |counts             |Gross counts                           |
|BACKGROUND             |counts             |Dark current + Scattered light         |
|BKG_ERR_UP             |counts             |Upper sys. error for BACKGROUND        |
|BKG_ERR_DOWN           |counts             |Lower sys. error for BACKGROUND        |
|DARK_CURRENT           |counts             |Estimated dark current model           |
|DARK_CURRENT_ERR       |counts             |Sys. error of DARK_CURRENT             |
|DQ                     |                   |Data quality flag                      |
|EXPTIME                |seconds            |Pixel exposure time                    |
|CALIB                  |cnts cm^2 A/erg    |Flux calibration factor, lin. interpolated NET/FLUX from 1dx files |
|FLAT_CORR              |                   |Flatfield and deadtime factor, (GROSS - BACKGROUND)/NET from 1dx files |
|LYA_SCATTER            |counts             |Estimated contamination by Lya, according to Worseck et al. 2016  |
|LYA_SCATTER_ERR_UP     |counts             |Upper error for LYA_SCATTER            |
|LYA_SCATTER_ERR_DOWN   |counts             |Lower error for LYA_SCATTER            |





-------------------------------------------------------------------------------------
4. REDUCTION PARAMETERS (faintcos_config.py)
-------------------------------------------------------------------------------------


PHA_OPT_ELEM_SEGMENT: Lower and upper threshold for the valid counts 
                      (see "COS Data Handbook (Version 4.0)" Section 3.4.9)
                      This parameter should be chosen according to the PHA distribution
                      of the detected photons in the PSA, since it changes with time
                      (see Figure 1.9 in "COS Data Handbook (Version 4.0)" 
                      and Worseck et al. 2016).
                      


DARK_EXPSTART_INTERVAL: For every exposure the code searches for darkframes in the 
                        time period +/-DARK_EXPSTART_INTERVAL around the exposure
                        date. The value should be given in days. Standard value is 
                        30 days. Can be increased, if there are not enough darkframes
                        in this time interval.
                        
MIN_DARKS:              Minimum number of darks that will be selected. 

KS_THRESHOLD:           Kolmogorov-Smirnov test threshold for the comparison of the 
                        cumulative PHA distributions (science vs. darkframe). Only
                        darkframes with KS-test lower than KS_THRESHOLD will be selected.
                        If the number of selected darks is lower than MIN_DARKS
                        than the KS_THRESHOLD will be increased by KS_STEP until 
                        sufficient number of darkframes is reached.
                        
KS_STEP:                KS_THRESHOLD will be increased with this value if the number
                        of selected darkframes is lower than MIN_DARKS
                        
BKG_AV:                 The width of the window in pixels for the central moving 
                        average to calculate the dark current model. The central 
                        moving average is applied on the PSA of the stacked darkframes,
                        which were selected with the KS-Test. Since the algorithm uses 
                        central moving average, BKG_AV must be an odd number! 

BIN_SIZE:               Bin size of the co-added spectrum in Angstrom. It is relevant
                        only if COADD_ALL_VISITS = True. 

BIN_PX:                 Bin size of the binned spectrum for every data set in pixels. It
                        is relevant only if BIN_VISIT = True.

CUSTOM_INTERVAL:        If TRUE, the co-add routine will use custom wavelength region
                        for the final co-add of all data sets in the working folder. The 
                        regions is defined by WAVE_MIN and WAVE_MAX in Angstrom. It works
                        only for the total co-add of all data sets (COADD_ALL_VISITS = True).
                        If FALSE, the co-add routine will use max. and min. wavelength of 
                        all available data.
                        
BIN_VISIT:              If TRUE, the spectra for every data set will be binned in BIN_PX pixels
                        and stored in the working folder as TARGETNAME_DATASET_Npx_binned.fits
                        
COADD_ALL_VISITS:       If TRUE, the spectra of all data sets in the working folder will be 
                        co-added to a single spectrum with binning size BIN_SIZE Angstrom.
                        !!!FOR THIS TO WORK PROPERLY, MAKE SURE TO HAVE SPECTRA FOR THE SAME
                        OBJECT AND RESOLUTION (M OR L) IN THE WORKING FOLDER!!!
                      
FELDMAN_COUSINS:        if TRUE, the pipeline will use the algorithm from Feldman & Cousins 1998
                        to calculate the statistical count errors in poisonian regime. 
                        Otherwise it uses the algorithm from Kraft et al. 1991.
                        For this to work, you need to install the custom C library. (see INSTRUCTIONS)





-------------------------------------------------------------------------------------
5. CALCULATED UNCERTAINTIES
-------------------------------------------------------------------------------------

We provide two methods for the calcualtion of the 1\sigma uncertainties.

 - The Bayesian method (Kraft et al. 1991, http://articles.adsabs.harvard.edu/pdf/1991ApJ...374..344K) 
   with the shortest 68.26% confidence interval. Kraft et al. (1991) provide only the confindence 
   limits. We modified it by assuming that the measured signal is N-B (N - measured counts, 
   B - background counts). 
 - The frequentist method (Feldman & Cousins 1998, https://arxiv.org/abs/physics/9711021), 
   which use likelihood ratio ordering for the measured signal (N-B). For the unphysical case N<B, 
   we use Monte Carlo simulations to calculate the limits.
   
Then, the calculated uncertainties are ERR_UP = upper_limit - (N-B) and 
ERR_DOWN = (N-B) - lower_limit. The uncertainties for the unphysical cases where the signal 
is negative (the measured counts N is smaller than the background B) should be treated with caution. 
For this, we suggest to use the "sensitivity" as defined by Feldman & Cousins (1998) which is the upper 
limit on the poisson distributed background with no signal. FaintCOS does not provide sensitivity calcualtions.



-------------------------------------------------------------------------------------
6. OPTIONAL CODE
-------------------------------------------------------------------------------------
"exec_calcos_for_darks.py": executes CALCOS on dark frames and puts the reduced corrtag files
                            into the "/reduced/" folder. Copy this code into the folder with
                            the downloaded dark frame "rawtag" files and run it.
                            
                            
"plot_datasets_info.py": plots the stacked 2D spectra from the "corrtag" files for each dataset
                         and segment (FUVA/FUVB). Additionally, it plots the histogram for the PHA
                         of the counts in the primary science aperture. The primary science aperture
                         is indicated with dashed red lines in the 2D spectrum. 
                         Run the code with the path to the corrtags and 1dx files as command-line argument.
                         The code creates a new subdirectory "/datasets_info" where it stores the resulted 
                         2D spectra. It also creates a .txt file with a table for all corrtags in the diretory.
