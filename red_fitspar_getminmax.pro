; docformat = 'rst'

;+
; Get datamin and datamax values from header or calculate them if
; necessary. 
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
;   datafile : in, type=string
; 
;      The name of the file with the data.
; 
; 
; :Keywords:
; 
;   datamax : out, optional, type=double
;   
;      The largest data value.
; 
; 
;   datamax : out, optional, type=double
;   
;      The largest data value.
; 
;   datamin : out, optional, type=double
;   
;      The smallest data value.
; 
;   force : in, optional, type=boolean
;
;      Calculate the values, regardless of whether they are in the
;      headers. 
; 
;   write : in, optional, type=boolean
;
;      Write the calculated values to the file header.
;
; :History:
; 
;    2020-10-28 : MGL. First version.
; 
;-
pro red_fitspar_getminmax, datafile $
                           , datamax = datamax $
                           , datamin = datamin $
                           , force = force $
                           , write = write

  if n_elements(datafile) eq 0 then return

  if ~arg_present(datamin) and ~arg_present(datamax) then return

  ;; Try the DATA* keywords first.
  hdr = headfits(datafile)
  datamin = fxpar(hdr, 'DATAMIN', count = mincount)
  datamax = fxpar(hdr, 'DATAMAX', count = maxcount)
  
  ;; Are we done?
  if ~keyword_set(force) and mincount eq 1 and maxcount eq 1 then return

  ;; No? Read the file, frame by frame, and calculate the values. 

  dims = fxpar(hdr, 'NAXIS*', Naxis)
  Nwav    = dims[2]
  Nstokes = dims[3]
  Nscans  = dims[4]
  
  Nframes = Nwav * Nstokes * Nscans

  undefine, datamin, datamax
  for iframe = 0, Nframes-1 do begin

    red_progressbar, iframe, Nframes $
                     , /predict $
                     , 'Calculate min/max'
    
    red_fitscube_getframe, datafile, im, iframe = iframe

    red_append, datamin, min(im, max = mx, /nan)
    red_append, datamax, mx
    
  endfor                        ; iframe

  datamin = min(datamin)
  datamax = max(datamax)
  
end

if keyword_set(write) then begin
  fxaddpar, hdr, 'DATAMIN', datamin, 'The minimum data value'
  fxaddpar, hdr, 'DATAMAX', datamax, 'The maximum data value'
  red_fitscube_newheader, datafile, hdr
  endif
  
end
