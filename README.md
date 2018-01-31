# EchoPlanarReco
Matlab/Ismrmrd script for single-shot EPI reconstruction

## Purpose

This is a simple example of how to reconstruct a single-shot EPI images from data. The code is really similar to theses gadgets (sometimes it purely of c++-> matlab conversion):

* NoiseAdjustGadget.cpp NoiseAdjustGadget.h
* EPIReconXGadget.cpp EPIReconXGadget.h
* EPICorrGadget.cpp and EPICorrGadget.h with the linear fit option
* PCACoilGadget.cpp PCACoilGadget.h
* CoilReductionGadget.cpp CoilReductionGadget.h

The script used some routines or codes from https://github.com/ismrmrd/ismrmrd/tree/master/examples/matlab 

## Experimental data 

* fully sampled dataset: phantom images available from the gadgetron in test/integration using the python script get data.py :
meas_MID517_nih_ep2d_bold_fa60_FID82077.h5 & meas_MID517_nih_ep2d_bold_fa60_NOISE82077.h5
* undersampled dataset with GRAPPA 2 : a gel agar images from my ftp server or send me an email : 00085_epi_RF_1dyn_GRAPPA2_FID10853.h5 & 00085_epi_RF_1dyn_GRAPPA2_NOISEFID10853.h5

The code has only been tested with the two mentionned datasets, minor modification might be necessary to make it works in other cases.

## Running the code

* use main.m
* download the ismrm_sunrise_matlab-master https://github.com/hansenms/ismrm_sunrise_matlab	
* go to line 57 and modify the link addpath

* download the data and change the links line 68 and 69 for fully sampled dataset
* download the data and change the links line 75 and 76 for undersampled dataset

## Capabilities:
   * 2D
   * use noise scans (pre-whitening)
   * remove oversampling in the readout direction
   * regridding of EPI traectory
   * ghost corretion using linear fit of phase diff between +/- readout
   * virtual coil using pca
   * coil reduction
   * magnitude reconstruction     

##  Limitations:
   * only works with a single encoded space
   * no partial fourier
   * could handle repetitions averages, phases, segments and sets 
   * phase reconstruction is creepy
   * only grappa 2 undersampling with a 24 reference/calibration scans 
   

## Restriction
   * working with an existing ISMRMRD data set.
   * required a sequence that fill all EPI flags and header according to the ismsmrd convention
   * 3 navigator lines/scans at the beginning of acquisition with positive and negative readout
   

## Some arbitrary conventions 
Data is a structure that include data from the hdf5 without kspace reordering. The matrix dimensions are [R0 "shots" CHA], "shots" meaning the number of lines or scans that has been acquired successively whatever the encoding. After noise prewhitening and regridding, the data are reordered into an other structure called kspace, the matrix dimensions are [R0 E1 CHA]
Finally reconstruction is done for slice 1=1 , contrast =1, average =1, generalisation is easy, using a loop, you can create a matrix [R0 E1 CHA SLC SET AVE] 
 
* data.kspace.raw is raw data
* data.kspace1D.raw is the raw data after ifft into the readout direction
* data.kspace1D.pre_whitening  is the data after  pre_whitening and with an ifft into the readout direction
* data.kspace1D.regridding  is the data after regridding and with an ifft into the readout direction

* kspace.image is the kspace of the image 
* kspace.calibration is the kspace of the calibration scan (24 lines usign Frappa 2 )   

* kspace.image_ghost_corrected is the kspace of the image after ghost correction 
* kspace.calibrationafter ghost correction  is the kspace of the
* calibration scan (24 lines using % Grappa 2 ) after ghost correction  
* etc

## Figure undersampled dataset

![Figure 2](https://github.com/valeryozenne/EchoPlanarReco/blob/master/figures_grappa_sampling/figure2.png)

![Figure 3](https://github.com/valeryozenne/EchoPlanarReco/blob/master/figures_grappa_sampling/figure3.png)

![Figure 4](https://github.com/valeryozenne/EchoPlanarReco/blob/master/figures_grappa_sampling/figure4.png)
