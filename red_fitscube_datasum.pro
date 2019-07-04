; docformat = 'rst'

;+
; Calculate the 32-bit checksum of the data in a fitscube file.
; 
; Do the calculations frame by frame so we don't have to keep
; the whole data cube in memory.
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
;    The 32-bit checksum.
; 
; :Params:
; 
;   fname : in, type=string
; 
;     The file name of the fitscube for which we want the data
;     checksum. 
; 
; 
; :Keywords:
; 
;   checksum_array : out, optional, type=array
;   
;   
; 
; 
; :History:
; 
;    2019-07-03 : MGL. First version.
; 
;-
function red_fitscube_datasum, fname, checksum_array = checksum_array

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
  
  if ~file_test(fname) then begin
    print, inam + ' : File does not exist: '
    print, fname
    return, -1
  endif

  hdr = headfits(fname)
  dims = fxpar(hdr, 'NAXIS*', Naxis)
  Nx      = dims[0]
  Ny      = dims[1]
  Ntuning    = dims[2]
  Nstokes = dims[3]
  Nscans  = dims[4]
  
  undefine, checksum
  Nframes = Nscans*Nstokes*Ntuning
  for iframe = 0, Nframes - 1 do begin
    red_progressbar, iframe, Nframes $
                     , /predict $
                     , 'Calculating checksum of '+file_basename(fname)
    red_fitscube_getframe, fname, frame, iframe = iframe
    checksum32, frame, checksum, /incremental 
  endfor                        ; iframe

  if arg_present(checksum_array) then begin
    iframe = 0
 
    checksum_array = replicate(checksum, Ntuning, Nstokes, Nscans)
    for ituning = 0, Ntuning-1 do begin
      for istokes = 0, Nstokes-1 do begin
        for iscan = 0, Nscans-1 do begin
          red_progressbar, iframe, Nframes $
                           , /predict $
                           , 'Calculating frame by frame checksums of ' $
                           + file_basename(fname)
          red_fitscube_getframe, fname, frame  $
                                 , ituning = ituning $
                                 , istokes = istokes $
                                 , iscan = iscan
          checksum32, frame, checksum_frame
          checksum_array[ituning, istokes, iscan] = checksum_frame
          iframe++
        endfor                  ; iscan
      endfor                    ; istokes
    endfor                      ; ituning
    
  endif

  return, checksum

end

;; Testing (with a small enough cube):

imname = '/scratch/mats/2016.09.19/CHROMIS-dec18/cubes_nb/nb_3950_2016-09-19T09:28:36_scans=67-73_corrected_im.fits'
spname = red_strreplace(imname, '_im','_sp')


check_im = red_fitscube_datasum(imname, checksum_array = charr)
checksum32, red_readdata(spname), check_sp

;; The data checksum should be the same whether calculated frame by
;; frame AND it should be the same in the spectral version of the
;; cube: 
print, check_im, check_sp, check_im-check_sp


end
