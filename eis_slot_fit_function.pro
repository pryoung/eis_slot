
FUNCTION eis_slot_fit_FUNCTION, x, p


;+
; NAME:
;     EIS_SLOT_FIT_FUNCTION
;
; PURPOSE:
;     A function that describes the emission in the X-direction across
;     an EIS 40" slot exposure.
;
; CATEGORY:
;     Hinode; EIS; slot.
;
; CALLING SEQUENCE:
;     Result = EIS_SLOT_FIT_FUNCTION( X, P )
;
; INPUTS:
;     X:  Array of wavelengths corresponding to the EIS detector
;         pixels.   
;     P:  A seven-element array giving the parameters of the fit
;         function. The parameters are:
;          P[0:2]  Parameters of quadratic function.
;          P[3]    Left edge of boxcar function (angstroms).
;          P[4]    Width of slot (angstroms).
;          P[5]    Gaussian width of Gaussian broadening function.
;          P[6]    Background emission
;
; OUTPUTS:
;     An array of same size as X that gives the fit to the EIS slot
;     emission. 
;
; PROGRAMMING NOTES:
;     The fit function is a quadratic that is convolved with a boxcar
;     (to simulate the fixed width of the slot), and then convolved
;     with a Gaussian to simulate the line spread function. The
;     background is set to a constant.
;
; EXAMPLE:
;     See the routine EIS_FIT_SLOT_EXPOSURE for an example of how to
;     call this routine.
;
; MODIFICATION HISTORY:
;     Ver.1, 10-Sep-2021, Peter Young
;     Ver.2, 09-Dec-2021, Peter Young
;       fixed a problem with edge effects; modified function so that
;       p[4] is the width of the slot (not the right edge).
;     Ver.3, 13-Dec-2021, Peter Young
;       fixed bug for slot width
;-


nx=n_elements(x)
pix_size=x[1]-x[0]

;
; I pad X with 11 pixels on either side to create XX. I then use XX
; for the calculations. This prevents any problems with edge effects
; that may occur if the slit boundary is too close to the edge of X. 
;
x0=(findgen(11)-10.)*pix_size+x[0]-pix_size
x1=findgen(11)*pix_size+x[nx-1]+pix_size
;
xx=[x0,x,x1]
nxx=n_elements(xx)


fn=p[0]+p[1]*xx+p[2]*xx^2   ; quadratic function
bfn=fltarr(nxx)            ; boxcar function

;
; Derive the right edge of the slot (left edge + width)
;
right_edge=p[3]+p[4]

;
; Convolution with boxcar function.
; Note that I account for fractional pixel positions. 
;
getmin=min(abs(xx-p[3]),imin1)
dlambda=(xx[imin1]+pix_size/2.)-p[3]
IF dlambda LT 0. THEN bfn[imin1]=1.0 ELSE bfn[imin1]=dlambda/pix_size
;
getmin=min(abs(xx-right_edge),imin2)
dlambda=right_edge-xx[imin2]+pix_size/2.
IF dlambda GT pix_size THEN bfn[imin2]=1. ELSE bfn[imin2]=dlambda/pix_size
;
bfn[imin1+1:imin2-1]=1.0
;
fn=fn*bfn

;
; Define the Gaussian that is convolved with function.
; It is defined in pixel space, but width (p[5]) is given in
; angstroms. 
;
ng=15
xg=findgen(15)-7
g=exp(-xg^2/2./(p[5]/pix_size)^2)
g=g/total(g)
;
output=convol(fn,g,/edge_constant)

;
; Now trim output back to the original size of X.
;
output=output[11:11+nx-1]

;
; Add constant background.
;
output=output+p[6]

return,output

END
