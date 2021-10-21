; docformat = 'rst'

;+
; Demodulate a single set of four LC images.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
; 
; :Params:
; 
; 
; 
; 
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;   2021-10-15 : MGL. First version.
; 
;-
function red_demodulate_images, ims_lc $
                                , inv_polmat $
                                , isodate $
                                , pref $
                                , tavg

  dims = size(inv_polmat, /dim)
  Nx = dims[2]                  ; Image size
  Ny = dims[3]

  Nlc = 4
  Nstokes = 4
  
  ims_stokes = fltarr(Nx, Ny, 4)  
  
  for ilc = 0L, Nlc-1 do begin
    for istokes = 0L, 3 do begin
      ims_stokes[*,*,istokes] += reform(inv_polmat[ilc,istokes,*,*]) * ims_lc[*,*,ilc]      
    endfor                      ; istokes
  endfor                        ; ilc


;  ;; Combine the T and R Stokes cubes
;  ;;dum = size(ims_stokes_t,/dim)
;  ;;drm = dum / 8.0 
;  xx0 = round(Nx/8. - 1)
;  xx1 = round(Nx - Nx/8. - 1)
;  yy0 = round(Ny/8. - 1)
;  yy1 = round(Ny - Ny/8. - 1)
;  
;  mediant = median(ims_stokes_t[xx0:xx1,yy0:yy1,0])
;  medianr = median(ims_stokes_r[xx0:xx1,yy0:yy1,0])
;  aver = (mediant + medianr) / 2.
;  sct = aver / mediant
;  scr = aver / medianr
;  
;  ims_stokes = (sct * (ims_stokes_t) + scr * (ims_stokes_r)) / 2.
  
  ;; Telescope model 
  ;;line = (strsplit(wbstates[0].fpi_state,'_',/extract))[0]
  ;;print, inam+' : Detected spectral line -> '+line
  
  year = (strsplit(isodate, '-', /extract))[0]
  red_logdata, isodate, tavg, azel = azel,  tilt = tilt

  if min(finite([azel, tilt])) eq 0 then stop

  mtel = red_telmat(pref, {TIME:tavg, AZ:azel[0], ELEV:azel[1], TILT:tilt} $
                    , /no_zero, year=year)
  imtel = invert(mtel) 
  imtel /= imtel[0]             ; Normalized same way as in crisp::demodulate!

  ;; Apply the telescope Muller matrix
  ims_stokes1 = fltarr(Nx, Ny, Nstokes)
  for j=0, Nstokes-1 do begin
    for i=0, Nstokes-1 do begin
      ims_stokes1[*, *, j] += ims_stokes[*, *, i] * imtel[i, j]
    endfor
  endfor
  ;;ims_stokes = temporary(ims_stokes1)
  
  return, ims_stokes1
  
end
