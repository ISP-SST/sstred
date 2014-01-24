; docformat = 'rst'

;+
; Fit surfaces to xtilts and ytilts, convert results to tilts. 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics, 2012
; 
; 
; :Params:
; 
;    x : in, type=array
; 
;      The X coordinates of the tilts.
; 
;    y : in, type=array
; 
;      The Y coordinates of the tilts.
; 
;    sx : in, type=integer
; 
;      The X dimension of the fitted surface.
; 
;    sy : in, type=integer
; 
;      The Y dimension of the fitted surface.
; 
;    D : in, type=float
; 
;      The telescoep pupil diamter in meters.
; 
;    lambda : in, type=float
; 
;      The wavelength in meters.
; 
;    image_scale : in, type=float
; 
;      The image scale in arcseconds/pixel.  
;
; :Keywords:
; 
;    xtilts : in, optional, type=fltarr
; 
;      The X tilts at the positions x, y
; 
;    ytilts : in, optional, type=fltarr
; 
;      The Y tilts at the positions x, y
; 
;    foc : in, optional
; 
; 
; 
;    xoffs : out, optional
; 
;      The fitted surfaces of X offsets.
; 
;    yoffs : out, optional
; 
;      The fitted surfaces of Y offsets.
; 
;    dfoc : out, optional
; 
; 
; 
;    ddiv : out, optional
; 
; 
; 
; :History:
; 
;   2013-09-10 : MGL. Switch to using mpfit2dfun for fitting and
;                allowing x and y to specify a 2D irregular
;                grid.
;
;   2014-01-24 : MGL. Renamed the dxoffs and dyoffs keywords to xoffs
;                and yoffs. Added documentation. Conversion factor
;                from tilts to shift offsets.
; 
; 
; 
;-
function red_pinh_make_fits_planefunct, x, y, p

  return, p[0] + p[1]*x + p[2]*y

end

pro red_pinh_make_fits, x, y, sx, sy $
                        , D, lambda, image_scale, $
                        , XTILTS = xtilts $ 
                        , YTILTS = ytilts $ 
                        , XOFFS = xoffs $ 
                        , YOFFS = yoffs $ 
                        , FOC = foc $       
                        , DFOC = dfoc $     
                        , DDIV = ddiv       
  
  Nch = (size(xtilts, /dim))[0] 
  Npinh = (size(x, /dim))[0]

  xoffs = fltarr(sx, sy, Nch)
  yoffs = fltarr(sx, sy, Nch)

  ;; Conversion from shift in pixels to tilt parameters in radians,
  ;; multiply by (pi*h)/(2*lambda*F#), where h is the pixel pitch,
  ;; lambda is the wavelength, and F# is the focal ratio. This can be
  ;; rewritten as (pi^2 * D * arcsecperpix)/(lambda*360*3600), where D
  ;; is the pupil diameter.
  conv_fact = (!pi*!pi * D * image_scale)/(lambda*360.*3600.)

  ;; Coordinates for the offset files
  sxx=transpose(reform(rebin(findgen(sx),sx*sy,/samp),sy, sx))
  syy=reform(rebin(findgen(sy),sx*sy,/samp),sx, sy)

  dxtilts = fltarr(Nch, Npinh)
  dytilts = fltarr(Nch, Npinh)
  dfoc = fltarr(Nch, Npinh)
  ddiv = fltarr(Nch)

  for ich = 1, Nch-1 do begin

     ;; X tilts, differential to anchor channel
     dxtilts[ich, *] = xtilts[ich, *]-xtilts[0, *]
;     if use_robust then begin
;        dxtiltsc = robust_planefit(double(xx), double(yy), dxtilts[ich, *, *], dxtiltsfit)
;     endif else begin
;        dxtiltsc = planefit(xx, yy, dxtilts[ich, *, *], 0, dxtiltsfit)
;     endelse
     err = dxtilts[ich, *]*0.+1.
     dxtiltsc = MPFIT2DFUN('red_pinh_make_fits_planefunct', x, y, reform(dxtilts[ich, *]) $
                              , err, [0., 0., 0.]) 
     xoffs[*, *, ich] = red_pinh_make_fits_planefunct(sxx, syy, dxtiltsc)

     ;; Y tilts, differential to anchor channel
     dytilts[ich, *] = ytilts[ich, *]-ytilts[0, *]
;     if use_robust then begin
;        dytiltsc = robust_planefit(double(xx), double(yy), dytilts[ich, *, *], dytiltsfit)
;     endif else begin
;        dytiltsc = planefit(xx, yy, dytilts[ich, *, *], 0, dytiltsfit)
;     endelse 
     err = dytilts[ich, *]*0.+1.
     dytiltsc = MPFIT2DFUN('red_pinh_make_fits_planefunct', x, y, reform(dytilts[ich, *]) $
                              , err, [0., 0., 0.]) 
     yoffs[*, *, ich] = red_pinh_make_fits_planefunct(sxx, syy, dytiltsc)

     if n_elements(foc) ne 0 then begin
        
        ;; This part is not up to date or tested.

        ;; Focus diversity
        ;; The plane fit is just for display, the update is just
        ;; the mean value.
        dfoc[ich, *] = foc[ich, *]-foc[0, *]
;        if use_robust then begin
;           dfocc = robust_planefit(xx, yy, dfoc, dfocfit)
;           ddiv[ich] = biweight_mean(dfoc[ich, *, *]) 
;        endif else begin
;           dfocc = planefit(xx, yy, dfoc, 0, dfocfit)
;           ddiv[ich] = mean(dfoc[ich, *, *]) 
;        endelse
        err = foc[ich, *]*0.+1.
        dfocc = MPFIT2DFUN('red_pinh_make_fits_planefunct', xx, yy, dfoc[ich, *] $
                              , err, [0., 0., 0.]) 
        ddiv[ich] = mean(dfoc[ich, *]) 

;        mnn = min(dfocfit)*3
;        mxx = max(dfocfit)*3
;
;        tv, bytscl([rebin(reform(dfoc), Npinhx*10, Npinhy*10, /samp) $
;                    , rebin(reform(dfocfit, Npinhx, Npinhy), Npinhx*10, Npinhy*10, /samp) $
;                   ], mnn, mxx), 2

     endif                      ; if foc     

  endfor                        ; ich

end
