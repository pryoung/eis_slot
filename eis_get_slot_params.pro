
function eis_get_slot_params, data, quiet=quiet, plotfile=plotfile, $
                            keep_outliers=keep_outliers, _extra=extra


;+
; NAME:
;     EIS_GET_SLOT_PARAMS()
;
; PURPOSE:
;     Computes parameters for the EIS 40" slot using the output
;     from EIS_FIT_SLOT_EXPOSURE. Note that is called within
;     EIS_FIT_SLOT_EXPOSURE. 
;
; CATEGORY:
;     Hinode; EIS; slot.
;
; CALLING SEQUENCE:
;     Result = EIS_GET_SLOT_PARAMS( Data )
;
; INPUTS:
;     Data:  The structure output by EIS_FIT_SLOT_EXPOSURE. Note that
;            some of the tags are modified by EIS_GET_SLOT_PARAMS.
;
; OPTIONAL INPUTS:
;     Plotfile:  The name of file to which the IDL plot showing the
;                fit to the slot tilt will be sent. If specified, then
;                the plot will not be displayed on the screen.
;	
; KEYWORD PARAMETERS:
;     QUIET:  If set, then no plot is produced and no information
;             messages. 
;     KEEP_OUTLIERS: By default the routine removes outliers from the
;                    slot tilt fit (defined by points that are more
;                    than 5-sigma from the fit). If this keyword is
;                    set, then the outliers are retained.
;
; OUTPUTS:
;     Returns a structure with the following tags:
;      .fid  File ID giving date and time.
;      .yip  Y initial position.
;      .tilt  Tilt in mAng/pixel.
;      .tilterr  Error on tilt.
;
;     However, it is preferable to use DATA, which has several tags
;     updated by this routine including:
;      .slot_wid   Width of slot (angstroms)
;      .slot_cen   Centroid of slot (angstroms)
;      .slot_gwid  Gauss-width of slot (angstroms)
;      .tilt       Slot tilt (mAng/pixel)
;      .width      Average slot width (angstroms)
;      .gwid       Average Gauss-width.
;
;     If the input PLOTFILE is given, then the routine will also
;     produce a plot of the slot tilt fit.
;
; MODIFICATION HISTORY:
;     Ver.1, 13-Sep-2021, Peter Young
;     Ver.2, 07-Oct-2021, Peter Young
;        ypix now taken from input structure.
;     Ver.3, 09-Dec-2021, Peter Young
;        updated so that the fifth fit parameter is now the slot width
;        (not the right edge of the slot).
;     Ver.4, 29-Dec-2021, Peter Young
;        modified how int_avg is calculated for small windows.
;-

IF n_params() LT 1 THEN BEGIN
   print,'Use:  IDL> tilt=eis_get_slot_params(data [,/quiet, plotfile=, /keep_outliers])'
   return,-1
ENDIF

s=size(data.int,/dim)
nx=s[0]
ny=s[1]

ypix=data.ypix

;
; Get positions of left (l) and right (r) edges of slot.
;
lcen=reform(data.aa[3,*])
lcenerr=reform(data.sigmaa[3,*])
rcen=reform(data.aa[3,*]) + reform(data.aa[4,*])
rcenerr=sqrt( reform(data.sigmaa[3,*])^2 + reform(data.sigmaa[4,*])^2 )

;
; Get centroid from mean of lcen and rcen.
;
cen=0.5*(lcen+rcen)
cenerr=sqrt(lcenerr^2 + rcenerr^2)

;
; Get width of slot.
;
wid=reform(data.aa[4,*])
widerr=reform(data.sigmaa[4,*])

;
; Get Gauss width of slot.
;
gwid=reform(data.aa[5,*])
gwiderr=reform(data.sigmaa[5,*])


;
; Get average intensity across slit. Due to fall-off at edge of slot,
; I consider the 35 pixels centered on the centroid.
;
FOR i=0,ny-1 DO BEGIN
  IF data.bad_pix[i] EQ 0 THEN BEGIN 
    getmin=min(abs(data.wvl-cen[i]),imin)
    i0=max([0,imin-17])
    i1=min([imin+17,nx-1])
    data.int_avg[i]=average(data.int[i0:i1,i],missing=data.missing)
  ENDIF 
ENDFOR 

nmiss=data.nmiss
nw=n_elements(data.wvl)

;
; This removes any NaNs from the data, and any points for which there
; are more than 20% missing pixels.
; 
k=where(finite(cen) AND float(nmiss)/float(nw) LT 0.2 AND data.bad_pix EQ 0,n1)
xx=ypix[k]
yy=cen[k]
ee=cenerr[k]

IF NOT keyword_set(quiet) THEN print,'% EIS_GET_SLOT_PARAMS: '+trim(ny-n1)+' pixels were removed because of NaNs or too many missing pixels.'



;
; Initial fit.
;
cc=linfit(xx,yy,measure_errors=ee,sigma=sigma)
fit_all=ypix*cc[1]+cc[0]
fit_vals=xx*cc[1]+cc[0]

;
; Remove outliers and redo fit.
;
IF NOT keyword_set(keep_outliers) THEN BEGIN 
  kk=where(abs(yy-fit_vals) LE 5*ee,n2)
 ;
  IF n2 LT 2 THEN BEGIN
    print,'% EIS_GET_SLOT_PARAMS: not enough points to perform fit. Try setting /keep_outliers. Returning...'
    return,-1
  ENDIF 
 ;
  ind=k[kk]
  IF NOT keyword_set(quiet) THEN print,'% EIS_GET_SLOT_PARAMS: '+trim(n1-n2)+' outliers were removed from fit.'
 ;
  xx=xx[kk]
  yy=yy[kk]
  ee=ee[kk]
  ;
  cc=linfit(xx,yy,measure_errors=ee,sigma=sigma)
ENDIF ELSE BEGIN
  ind=k
ENDELSE 
  
fit_all=ypix*cc[1]+cc[0]
fit_vals=xx*cc[1]+cc[0]



yrange=minmax(fit_all)+[-0.02,0.02]

tilt=cc[1]*1000.
tilterr=sigma[1]*1000.


IF n_elements(plotfile) NE 0 OR keyword_set(quiet) THEN buffer=1 ELSE buffer=0

fs=12
p=errorplot(xx,yy,ee,linestyle='none',dim=[700,450], $
            pos=[0.15,0.12,0.98,0.98], $
            xtitle='Y-pixel', $
            ytitle='Position of slot edge / '+string(197b), $
            font_size=fs,buffer=buffer, $
            xticklen=0.015,yticklen=0.015, $
            xmin=1,/xsty,symbol='+', $
            yrange=yrange,/ysty, $
            xrange=[min(ypix),max(ypix)])
q=plot(/overplot,ypix,fit_all,thick=2)
xr=p.xrange
yr=p.yrange
t=text(0.18,0.18,data.fid,font_size=fs)
IF n_elements(plotfile) NE 0 THEN p.save,plotfile,resolution=192
IF n_elements(plotfile) NE 0 OR keyword_set(quiet) THEN p.close


IF NOT keyword_set(quiet) THEN print,format='("Tilt (mAng/pix): ",f10.4," +- ",f8.4)',tilt,tilterr


output={fid: data.fid, $
        yip: data.yip, $
        iexp: data.iexp, $
        tilt: tilt, $
        tilterr: tilterr }

data.slot_wid[ind]=wid[ind]
data.slot_wid_err[ind]=widerr[ind]
data.slot_cen[ind]=cen[ind]
data.slot_cen_err[ind]=cenerr[ind]
data.slot_gwid[ind]=gwid[ind]
data.slot_gwid_err[ind]=gwiderr[ind]

data.gwid=mean(gwid[ind])
data.gwiderr=stdev(gwiderr[ind])
data.width=mean(wid[ind])
data.widtherr=stdev(widerr[ind])

data.tilt=tilt
data.tilterr=tilterr


return,output

END
