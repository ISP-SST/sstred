; docformat = 'rst'

;+
; Read data frames from multiple files. Assume all frames have the
; same dimensions.
;
; Keywords will be used when calling red_readdata for the individual
; files. 
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
;    A datacube with the data from the files.
; 
; :Params:
; 
;   fnames : The file names.
; 
; :History:
; 
;   2022-08-22 : MGL. First version.
; 
;-
function red_readdata_multiframe, fnames $
                                  , status = status $
                                  , _ref_extra = extra

  Nfiles = n_elements(fnames)

  status = bytarr(Nfiles)
  
  for ifile =  0, Nfiles-1 do begin
  
    if n_elements(datacube) eq 0 then begin
      tmp = red_readdata(fnames[ifile] $
                         , _strict_extra = extra $  
                         , status = thisstatus)        
      ;; Could check thisstatus here and only create the datacube if
      ;; file read successfully.
      datacube = fltarr([size(tmp, /dim), Nfiles])
      datacube[0, 0, ifile] = tmp
    endif else begin
      datacube[0, 0, ifile] = red_readdata(fnames[ifile] $
                                           , _strict_extra = extra $  
                                           , status = thisstatus)        
    endelse

    status[ifile] = thisstatus
    
  endfor                        ; ifile
  
  return, datacube
  
end
