
FUNCTION eis_fit_slot_exposure, file, iexp, wavel=wavel, npix=npix, quiet=quiet, $
                                outdir=outdir, iexp_prp=iexp_prp, ybin=ybin, $
                                ystart=ystart, force=force, gwidth=gwidth, $
                                slot_width=slot_width, _extra=extra

;+
; NAME:
;     EIS_FIT_SLOT_EXPOSURE
;
; PURPOSE:
;     Takes a single EIS 40" slot exposure and fits the
;     emission in the X-direction with a parametric formula. 
;
; CATEGORY:
;     Hinode; EIS; slot.
;
; CALLING SEQUENCE:
;     Result = EIS_FIT_SLOT_EXPOSURE( File, Iexp )
;
; INPUTS:
;     File:   Either an EIS level-1 filename, or a windata structure. 
;     Iexp:   Index of the exposure to be fit. This is the exposure
;             index as used in the windata structure, which is in
;             reverse time order. The index is counted from zero.
;
; OPTIONAL INPUTS:
;     Wavel: A wavelength used to identify which EIS window to
;            process. If not specified then 195.12 is assumed.
;     Npix:  This allows the number of pixels in the window to be
;            specified. Needed for full CCD data.
;     Outdir: Specifies a directory to which the output data structure
;             will be saved, and the tilt plot file.
;     Iexp_Prp:  Only used if there are multiple exposures per slit
;                position. It is an integer specifying which exposure
;                to use.
;     Ybin:  An integer specifying the binning to be applied in the
;            Y-direction (implemented through eis_bin_windata).
;     Ystart: An integer specifying at which Y-pixel to begin
;             Y-binning. Default is Y-pixel 0. Only active if the YBIN
;             input is set.
;     Gwidth: Specifies a fixed value to use for the Gaussian width of
;             the slot line spread function. Units: Angstroms. Should
;             be around 0.028.
;     Slot_Width:  Specifies a fixed value to use for the slot
;             width. Units: arcsec. Should be around 41.
;
; KEYWORD PARAMETERS:
;     QUIET: If set then information messages and the tilt plot are
;            suppressed.
;     FORCE: The routine will not run if the wavelength window only
;            has 40 pixels as the results will not be reliable. If you
;            want to proceed, then set this keyword.
;	
; OUTPUTS:
;     An IDL structure containing the results. The tags are:
;      AA              FLOAT     Array[7, 512]
;      SIGMAA          FLOAT     Array[7, 512]
;      OUTFIT          FLOAT     Array[48, 512]
;      WAVEL           FLOAT           195.120
;      WVL             DOUBLE    Array[48]
;      INT             FLOAT     Array[48, 512]
;      ERR             FLOAT     Array[48, 512]
;      NMISS           INT       Array[512]
;      FILE            STRING    '/Volumes/PRY_DATA_3/hinode/eis/level1/2008/01/31/eis_l1_20080131_223514.fits.'...
;      IEXP            INT             11
;      NEXP            LONG                16
;      EXP_IND         LONG                 5
;      EXP_IND_STR     STRING    '5 of 16'
;      YIP             LONG               256
;      XCEN            DOUBLE           1061.3451
;      YCEN            DOUBLE           170.41593
;      FID             STRING    '20080131_2235'
;      TILT            FLOAT         -0.113630
;      TILTERR         FLOAT       0.000863131
;
;     AA - 7 fit parameters for each Y-pixel.
;     SIGMAA - 1-sigma errors on the fit parameters.
;     OUTFIT -  fit for each Y-pixel.
;     TILT - slot tilt in milli-Angstroms/pixel.
;
; CALLS:
;     EIS_GETWINDATA, EIS_AIA_OFFSETS, EIS_SLOT_FIT_FUNCTION,
;     EIS_GET_SLOT_PARAMS, EIS_BIN_WINDATA
;
; EXAMPLE:
;     IDL> file=eis_find_file('31-jan-2008 22:35',/level1)
;     IDL> output=eis_fit_slot_exposure(file,11)
;
; MODIFICATION HISTORY:
;     Ver.1, 10-Sep-2021, Peter Young
;     Ver.2, 06-Oct-2021, Peter Young
;       Fixed bug for iexp_prp implementation.
;     Ver.3, 07-Oct-2021, Peter Young
;       Added ybin= optional input; now catches case when fit fails.
;     Ver.4, 09-Dec-2021, Peter Young
;       Added ystart= optional input (passed on to eis_bin_windata);
;       added /force keyword; added GWIDTH= and SLOT_WIDTH= optional
;       inputs; now uses the slot width as a free parameter rather
;       than the right edge of the slot.
;     ver.5, 14-Dec-2021, Peter Young
;       Added int_avg and bad_pix to output structure; check if the
;       bottom slot edge is in the dataset, and set pixels below this
;       to missing.
;     Ver.6, 14-Feb-2022, Peter Young
;       Two bug fixes.
;     Ver.7, 16-Feb-2022, Peter Young
;       Minor bug fix.
;-



IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> output=eis_fit_slot_exposure( file, iexp [, wavel=, npix=, outdir=, iexp_prp= '
  print,'                                         /quiet, ybin=, ystart=, /force, gwidth=, slot_width= ] )'
  print,''
  print,'  iexp - exposure index counted from left-to-right, beginning with zero'
  print,'         (i.e., reverse time order)'
  return,-1
ENDIF 



IF n_tags(file) EQ 0 THEN BEGIN
  IF n_elements(wavel) EQ 0 THEN wavel=195.12
  wd=eis_getwindata(file,wavel)
  filename=file
ENDIF ELSE BEGIN
  wd=file
  filename=wd.filename
ENDELSE

;
; I set the bottom edge of the slit to detector y-pixel 7. All
; y-pixels below this get set to missing.
;
yws=wd.hdr.yws
IF yws LT 7 THEN wd.err[*,*,0:6-yws,*]=wd.missing


IF n_elements(ystart) EQ 0 THEN ystart=0

IF n_elements(ybin) NE 0 THEN BEGIN
  wd=eis_bin_windata(temporary(wd),ybin=ybin,ystart=ystart)
ENDIF ELSE BEGIN
  ybin=1
ENDELSE 
;
nexp=wd.nx
nexp_prp=wd.hdr.nexp_prp
calcheck=strpos(wd.hdr.calstat,'ABS')

s=eis_obs_structure(wd.hdr.date_obs,wd.hdr.date_end,/quiet,count=count)
IF count EQ 0 THEN s=eis_obs_structure(wd.hdr.date_obs,wd.hdr.ccsds_en,/quiet,count=count)
exptimes=str_sep(s[0].exp_times,' ')
IF nexp_prp GT 1 THEN BEGIN
   IF n_elements(iexp_prp) EQ 0 THEN BEGIN
      print,'% EIS_FIT_SLOT_EXPOSURE: there are multiple exposures times per slit position. Please call the routine'
      print,'                         again and set IEXP_PRP to one of the values below.'
      FOR i=0,n_elements(bits)-1 DO BEGIN
         print,format='(10x,i2,". ",a7)',i,bits[i]
      ENDFOR 
      return,-1
   ENDIF ELSE BEGIN
      print,'% EIS_FIT_SLOT_EXPOSURE: selected exposure time is '+exptimes[iexp_prp]+' seconds.'
   ENDELSE 
ENDIF ELSE BEGIN
   iexp_prp=0
ENDELSE


IF n_elements(iexp) EQ 0 THEN BEGIN
   IF nexp EQ 1 THEN BEGIN
      iexp=0
   ENDIF ELSE BEGIN
      print,'% EIS_FIT_SLOT_EXPOSURE: IEXP has not been specified. Please choose an exposure number from 0 to '+trim(nexp-1)+'.'
      return,-1
   ENDELSE
ENDIF 

IF n_elements(wavel) EQ 0 THEN BEGIN
   wavel=195.12
   print,'% EIS_FIT_SLOT_EXPOSURE: WAVEL not specified. Using '+string(format='(f6.2)',wavel)+'.'
ENDIF 




wvl=wd.wvl
nw=n_elements(wvl)
getmin=min(abs(wvl-wavel),imin)

IF nw LT 48 AND NOT keyword_set(force) THEN BEGIN
   print,'% EIS_FIT_SLOT_EXPOSURE: the wavelength window is only '+trim(nw)+' pixels wide. This routine will only'
   print,'                         work accurately if the window is at least 48 pixels wide. If you really want to'
   print,'                         proceed, use the /FORCE keyword.'
   print,'                         Returning...'
   return,-1
ENDIF


IF nw EQ 1024 AND n_elements(npix) EQ 0 THEN BEGIN
   print,'% EIS_FIT_SLOT_EXPOSURE: This is a full-CCD dataset. Please use the input NPIX to reduce the wavelength'
   print,'                         window size. Values of 48 to 64 pixels are recommended. Returning...'
   return,-1
ENDIF



i0=0 & i1=nw-1
IF n_elements(npix) NE 0 THEN BEGIN 
   IF npix LT 48 THEN BEGIN
      print,'% EIS_FIT_SLOT_EXPOSURE: The input NPIX must be at least 48. Returning...'
      return,-1
   ENDIF
  ;
   IF npix LT nw THEN BEGIN
      i0=max([imin-npix/2,0])
      i1=i0+npix-1
   ENDIF
ENDIF 

wvl=wvl[i0:i1]
nw=n_elements(wvl)
int=reform(wd.int[i0:i1,iexp,*,iexp_prp])
err=reform(wd.err[i0:i1,iexp,*,iexp_prp])

;
; If the /noabs keyword was given to eis_prep, then the lines below
; multiply int and err by the calibration factor.
;
IF NOT calcheck THEN BEGIN
   calib=eis_slot_calib_factor(windata=wd)
   int=int*calib
   k=where(err NE wd.missing)
   err[k]=err[k]*calib
ENDIF 


s=size(int,/dim)
nx=s[0]
ny=s[1]

;
; y-pixel indices.
;
ypix=indgen(ny)*round(wd.scale[1])+yws+ystart+round(wd.scale[1])/2


pix_size=wvl[1]-wvl[0]

nparams=7
nfitparams=nparams

output=fltarr(nparams,ny)
error=fltarr(nparams,ny)
outfit=fltarr(nx,ny)
chi2=fltarr(ny)+wd.missing

;
; Set limits on parameters.
; 
parinfo=replicate({fixed: 0, limited: [0,0], limits:[0.d,0.d]},nparams)
parinfo[6].limited=[1,0]  ; Background must be >0
parinfo[6].limits[0]=0.0  ;
parinfo[5].limited=[1,1]       ; Gauss width must be between 0.01 and 0.10.
parinfo[5].limits=[0.01,0.10]  ;
parinfo[4].limited=[1,1]            ; Slot width 
parinfo[4].limits=[39,43]*pix_size  ; must be between 39 and 43"

;
; If gwidth is set, then the Gaussian width is set fixed to this value
; for the fit.
;
IF n_elements(gwidth) NE 0 THEN BEGIN
  parinfo[5].fixed=1
  init_width=gwidth
  nfitparams=nfitparams-1
ENDIF ELSE BEGIN
  init_width=0.03
ENDELSE

;
; If slot_width is specified, then it is assumed a fixed parameter.
; It is specified in arcsec, so has to be converted to angstroms.
;
IF n_elements(slot_width) NE 0 THEN BEGIN 
  parinfo[4].fixed=1
  IF slot_width LT 30. THEN BEGIN
    print,'% EIS_FIT_SLOT_EXPOSURE: The slot_width should be specified in arcsec, and so should have a value '
    print,'                         close to 40 arcsec. Please revise your input.'
    return,-1
  ENDIF 
  init_slot_width=slot_width*pix_size
  nfitparams=nfitparams-1
ENDIF ELSE BEGIN
  init_slot_width=41*pix_size
ENDELSE


;
; This keeps track of number of missing pixels (along the slit), and
; bad pixels (y-pixels where no fit was possible).
;   bad_pix=0  good pixel
;   bad_pix=1  fit failed
;   bad_pix=2  signal too low to perform fit.
;
nmiss=intarr(ny)
bad_pix=bytarr(ny)

FOR i=0,ny-1 DO BEGIN
 ;
  no_fit=0b
 ;
  err_i=err[*,i]
  int_i=int[*,i]
  k=where(err_i NE wd.missing,nk)
  nmiss[i]=nw-nk
 ;
  IF nk LT nw/2 THEN no_fit=1b
 ;
 ; This deals with missing pixels. I require at least half the pixels
 ; to not be missing. 
 ;
  IF nw NE nk AND nk GE nw/2 THEN BEGIN
    x=wvl[k]
      z=err_i[k]
      y=int_i[k]
     ;
      y2=spl_init(x,y)
      int_i=spl_interp(x,y,y2,wvl)
      int[*,i]=int_i
     ;
      z2=spl_init(x,z)
      err_i=spl_interp(x,z,z2,wvl)
      err[*,i]=err_i
   ENDIF  

  ;
  ; Set initial parameters.
  ;  - Gaussian width is set in Angstroms.
  ;  - slot boundaries make use of input WAVEL
  ;
   init=[ max(int_i), 0., 0., $
          wavel-18*pix_size, init_slot_width, $
          init_width, $
          30. ]
  ;
  ; Check signal to noise. Only fit if the intensity in the center of
  ; the slot is x2 above the errror.
  ;
   getmin=min(abs(wvl-(init[3]+init[4]/2.)),imin)
   int_chck=average(int_i[imin-3:imin:3])
   err_chck=average(err_i[imin-3:imin:3])
   IF int_chck LE 2.*err_chck THEN no_fit=1b

   IF NOT no_fit THEN  BEGIN 
   
     aa = MPFITEXPR('eis_slot_fit_FUNCTION(x,p)', wvl, int_i, err_i, init, $
                    perr=sigmaa, /quiet, bestnorm=bestnorm,yfit=yfit, $
                    parinfo=parinfo,status=status)
  ;
  ; If the fit fails, then sigmaa won't be defined. I flag this
  ; by setting nmiss to the maximum value.
  ;
     IF n_elements(sigmaa) NE 0 THEN BEGIN
       output[*,i]=aa
       error[*,i]=sigmaa
       outfit[*,i]=yfit
       chi2[i]=bestnorm/(n_elements(int_i)-nfitparams)
     ENDIF ELSE BEGIN
       nmiss[i]=nw
       bad_pix[i]=1b
     ENDELSE
   ENDIF ELSE BEGIN
     bad_pix[i]=2b
   ENDELSE 
ENDFOR

;
; Get xcen and ycen, correcting for AIA offset.
;
xy=eis_aia_offsets(wd.time_ccsds[iexp])
xcen=wd.solar_x[iexp]+xy[0]
ycen=wd.solar_y[ny/2]+xy[1]


missing=wd.missing
y_arr=make_array(ny,value=missing,/float)

;
; Define output structure.
;
data={aa: output, $
      sigmaa: error, $
      outfit: outfit, $
      chi2: chi2, $
      wavel: wavel, $
      wvl: wvl, $
      int: int, $
      err: err, $
      nmiss: nmiss, $
      file: filename, $
      iexp: iexp, $
      nexp: nexp, $
      exp_ind: nexp-iexp, $
      exp_ind_str: trim(nexp-iexp)+' of '+trim(nexp), $
      yip: yws, $
      xcen: xcen, $
      ycen: ycen, $
      ypix: ypix, $
      fid: time2fid(wd.hdr.date_obs,/full,/time), $
      missing: missing, $
      slot_wid: y_arr, $
      slot_wid_err: y_arr, $
      slot_cen: y_arr, $
      slot_cen_err: y_arr, $
      slot_gwid: y_arr, $
      slot_gwid_err: y_arr, $
      int_avg: y_arr, $
      bad_pix: bad_pix, $
      tilt: 0., $
      tilterr: 0., $
      width: 0., $
      widtherr: 0., $
      gwid: 0., $
      gwiderr: 0., $
      t_obs: wd.time_ccsds[iexp], $
      stud_acr: wd.hdr.stud_acr, $
      pix_size: pix_size, $
      ybin: ybin, $
      ystart: ystart}

IF n_elements(outdir) NE 0 THEN BEGIN
   chck=file_info(outdir)
   IF chck.exists EQ 0 THEN file_mkdir,outdir
  ;
   outfile='data_'+data.fid+'_'+strpad(trim(data.iexp),2,fill='0')+'.save'
   outfile=concat_dir(outdir,outfile)
  ;
   plotfile='plot_'+data.fid+'_'+strpad(trim(data.iexp),2,fill='0')+'.png'
   plotfile=concat_dir(outdir,plotfile)
ENDIF 

;
; This loads the tilt values into data
;
output=eis_get_slot_params(data,quiet=quiet,plotfile=plotfile,_extra=extra)
junk=temporary(output)

IF n_elements(outdir) NE 0 THEN save,file=outfile,data

return,data

END
