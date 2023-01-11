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
;   2022-12-22 : MGL. Allow the files to have multiple frames.
; 
;-
function red_readdata_multiframe, fnames $
                                  , status = status $
                                  , _ref_extra = extra

  Nfiles = n_elements(fnames)

  status = bytarr(Nfiles)

  ;; Find out dimensions of datacube
  Nframes_perfile = lonarr(Nfiles)
  for ifile =  0, Nfiles-1 do begin

    red_progressbar, ifile, Nfiles, 'Reading headers from multiple files'

    hdr = red_readhead(fnames[ifile])
    
    if ifile eq 0 then begin
      Nx = red_fitsgetkeyword(hdr,'NAXIS1',count=cnt)
      Ny = red_fitsgetkeyword(hdr,'NAXIS2',count=cnt)
    endif
    
    naxis3 = red_fitsgetkeyword(hdr,'NAXIS3',count=cnt)
    if cnt eq 0 then Nframes_perfile[ifile] = 1 else Nframes_perfile[ifile] = naxis3
    
  endfor                        ; ifile

  Nframes = total(Nframes_perfile)
  datacube = fltarr(Nx, Ny, Nframes)
  
  for ifile =  0, Nfiles-1 do begin

    red_progressbar, ifile, Nfiles, 'Reading data frames from multiple files'

    if ifile eq 0 then iframe = 0 else iframe = total(Nframes_perfile[0:ifile-1])
    
    datacube[0, 0, iframe] = red_readdata(fnames[ifile] $
                                          , _strict_extra = extra $  
                                          , status = thisstatus)        
    
    status[ifile] = thisstatus
    
  endfor                        ; ifile
  
  return, datacube
  
end
