# eis_slot
Software for processing EIS slot data

This software was developed for a project to measure properties of the 40" slot of the EUV Imaging Spectrometer (EIS) on board the Hinode spacecraft. The results will be published in Young & Ugarte-Urra (2022, in preparation).

The routine eis_fit_slot_exposure performs a 7-parameter fit to the slot emission in the X-direction for each Y-pixel along a single slot exposure. 

An example of how to use the routine is given below

IDL> file=eis_find_file('31-jan-2008 22:35',/level1)
IDL> output=eis_fit_slot_exposure(file,11)

Here the routine processes the 11th exposure of the dataset. The routine automatically looks for the Fe XII 195.12 line and fits this.

The routine can also accept a windata structure. For example:

IDL> file=eis_find_file('31-jan-2008 22:35',/level1)
IDL> wd=eis_getwindata(file,195.12)
IDL> output=eis_fit_slot_exposure(wd,11)

Some options available for the routine:

npix=  

This reduces the size of the wavelength window, with the new window centered on the wavelength specified by wavel=. This is useful for processing full-CCD data. For example:

IDL> output=eis_fit_slot_exposure(wd,11,wavel=195.12,npix=64)


ybin= & ystart=

These enable y-binning to be applied and for the start pixel of the binning to be adjusted. For example, to use a y-binning of 5 and to start from pixel 23, do:

IDL> output=eis_fit_slot_exposure(wd,11,ybin=5,ystart=23)


gwidth= and slot_width=

These specify fixed values for the LSF Gaussian width and/or the slot width. These may be useful if you are struggling to get a good fit (perhaps because the window is too narrow). 

IDL> output=eis_fit_slot_exposure(wd,0,/force,gwidth=0.0285,slot_width=40.949,ybin=10)

The keyword /force is needed if the wavelength windows are only 40 pixels wide, since the routine by default rejects these datasets.
