pro red_get_fpi_sep,fpi,fpi_wav=fpi_wav
; Adapted from Goran's script
; Calculate the cavity separation so both peaks are
; aligned at fpi.w0
; It replaces the nominal value from the structure fpi
;
; hre alignment
;
  if keyword_set(fpi_wav) then w0=fpi.fpi_wav else w0=fpi.w0
;
  nhr=long(0.5d0+fpi.shr/(w0/2.d0)) 
  fpi.shr=nhr*w0/2.d0
;
; lre alignment
;
  nlr=long(0.5d0+fpi.slr/(w0/2.d0)) 
  fpi.slr=nlr*w0/2.d0
;
  return
end
