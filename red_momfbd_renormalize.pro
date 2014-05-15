; docformat = 'rst'

;+
; The Zernike and Karhunen-Loeve modes used internally in the MOMFBD
; program are not normalized to unit RMS. This routine
; normalizes the modes and rescales the coefficients to match the
; new modes.
;
; Before this is done, one cannot calculate variances of the
; estimated wavefronts by summing the squares of the coefficients.
;
; :Categories:
;
;   Image restoration, mombfd
;
; :Author:
;
;   Mats Löfdahl, Institute for Solar Physics, mats@astro.su.se.
;
; :Params:
;
;   mr : in, out, type=struct
;
;     A structure as read with the momfbd_read command. On output,
;     same structure but with alpha and modes renormalized.
;
; :Keywords:
;
;   verbose : in, optional, type=boolean
;
;     Set this for more printout.
;
; :History:
;
;   2010-09-21: Written by Mats Löfdahl
;
;   2012-11-21 : Switched to docformat 'rst'.
;
;   2013-08-30 : MGL. Renamed for inclusion in crispred pipeline.
;
;   2014-05-15 : MGL. Don't use non-standard function "dimen". 
;   
;   
;-
pro red_momfbd_renormalize, mr, verbose = verbose

  ;; We want the RMS within the pupil
  pupindx = where(mr.pupil)

  Npass = 3

  ;; Iterate!
  for ip = 0, Npass-1 do begin

     if keyword_set(verbose) then print,'Pass '+strtrim(string(ip), 2)

     ;; Loop over all modes
     for ii=0, (size(mr.mode, /dim))[2]-1 do begin
        
        ;; Calculate the RMS for mode ii:
        modestdev = stdev((mr.mode[*,*,ii])[pupindx])
        
        if keyword_set(verbose) then print, 'RMS of mode '+strtrim(string(ii), 2)+': '+strtrim(string(modestdev))
        
        ;; Normalize modes to unit stdev:
        mr.mode[*,*,ii] = mr.mode[*,*,ii] / modestdev
        
        ;; And compensate the coefficients:
        mr.patch.alpha[ii, *] = mr.patch.alpha[ii, *] * modestdev
        
     endfor

  endfor

end
