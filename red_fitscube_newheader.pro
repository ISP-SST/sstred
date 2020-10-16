; docformat = 'rst'

;+
; Add a new FITS header to an existing fitscube file.
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
;    filename : in, type=string
;
;       The name of the file in which to write the new header.
; 
;    newheader
;
;       The new header.
; 
; :Keywords:
; 
;   Nframes_max : in, optional, type=integer
;   
;     The largest number of frames for which we accept to move the
;     data part on disk, rather than re-writing the file.
; 
; 
; :History:
; 
;   2019-09-12 : MGL. First version.
;
;   2020-10-16 : MGL. Use red_fitscube_copyextensions.
;
;-
pro red_fitscube_newheader, filename, newheader, Nframes_max = Nframes_max

  ;; Name of this subprogram
  inam = red_subprogram(/low, calling = inam1)

  ;; Call it oldfilename
  oldfilename = filename

  print, inam+' : Write a new header to '+oldfilename
  
  if n_elements(Nframes_max) eq 0 then Nframes_max = 10
  
  oldheader = headfits(oldfilename)

  ;; Check that the new header has the right dimensions and data type.
  if fxpar(newheader, 'BITPIX') ne fxpar(oldheader, 'BITPIX') then stop

  NAXIS_old = fxpar(oldheader, 'NAXIS')
  NAXIS_new = fxpar(newheader, 'NAXIS')
  if NAXIS_new ne NAXIS_old then stop

  NAXISi_old = fxpar(oldheader, 'NAXIS*')
  NAXISi_new = fxpar(newheader, 'NAXIS*')
  for i = 0, NAXIS_old-1 do if NAXISi_old[i] ne NAXISi_new[i] then stop

  Nframes = round(product(NAXISi_new[2:*]))
  
  ;; The header is probably OK.

  ;; Check if the header fits in the same number of 2880-byte blocks,
  ;; i.e., 36-line blocks.
  Nlines_old = n_elements(oldheader)
  Nlines_new = n_elements(newheader)

  if Nlines_old/36 eq Nlines_new/36 or Nframes lt Nframes_max then begin
    ;; It's OK to just write the new header in the existing file.
    print, inam+' : Write new header with modfits...(This may take a while)'
    tic
    modfits, oldfilename, 0, newheader
    toc
    print, inam+' : Write new header with modfits...Done!'
    return
  endif 

  ;; Construct a temporary file with the new header and the old data.
  ;; Then copy all extensions. 

  print, inam+' : Write the header to a temporary file and then copy data and extensions.'


  ;; Open and read header from the old file
  red_fitscube_open, oldfilename, oldfileassoc ;, oldfitscube_info

  ;; Initialize the tmp file
  tmpfilename = 'tmp_fitscube_newheader_' + cgTimeStamp(random_digits=6) + '.fits'
  red_fitscube_initialize, tmpfilename, newheader, tmplun, tmpfileassoc, NAXISi_old


  ;; Copy data frames
  for iframe = 0, Nframes-1 do begin
    red_progressbar, iframe, Nframes, 'Copy data frames'
    red_fitscube_getframe, oldfileassoc, frame, iframe = iframe 
    red_fitscube_addframe, tmpfileassoc, frame, iframe = iframe 
  endfor

  ;; Close the files, add the wcs info to the tmp file
  red_fitscube_finish, tmplun   ;, wcs = wcs_coordinates
  red_fitscube_close,  oldfileassoc

  ;; Copy all extensions, including variable keywords, from the old file
  red_fitscube_copyextensions, oldfilename, tmpfilename
  
  ;; Overwrite the old file
  print, inam+' : Move temporary file to the original file name...'
  tic
  file_move, tmpfilename, oldfilename, /overwrite
  toc
  fxhmodify, oldfilename, 'FILENAME', file_basename(oldfilename)
  print, inam+' : Move temporary file to the original file name...Done!'
  
end

;; Testing

origname = '/scratch/mats/2016.09.19/CHROMIS-jan19/cubes_nb/nb_3950_2016-09-19T09:28:36_scans=60-79_corrected_im.fits'

dir = file_dirname(origname)

testname = dir + '/test.fits'
file_copy, origname, testname, /overwrite

;; Make a longer header
newhdr = headfits(testname)
for i = 0, 40 do fxaddpar, newhdr, 'TEST'+string(i, format = '(i02)'), 'test'


red_fitscube_newheader, testname, newhdr, Nframes_max = 3 ;0000

end
