; docformat = 'rst'

;+
; Set pixels with missing data to NaN or median(frame) in a fitscube
; file. 
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
;   filename : in, type=string
; 
;     The path of the file.
; 
; 
; :Keywords:
; 
;   force : in, optional, type=boolean
;
;     Do it again if it's already done.
;
;   missing_type  : in, optional, type=string, default='opposite'
;
;     What value to set missing-data pixels to. One of 'nan' and
;     'median' (of each frame). By default, set it to the opposite of
;     what it appears to be in the input cube.
;
;   noflip  : in, optional, type=boolean
;
;     Don't make a flipped version.
;
;   oname  : in, optional, type=string, default=filename
;
;     The path of the output file.
;
;   overwrite  : in, optional, type=boolean
;   
;     If oname is given and is an existing file, don't
;     overwrite unless this keyword is set.
;   
; 
; 
; :History:
; 
;    2020-06-29 : MGL. First version.
;
;    2020-10-28 : MGL. Remove statistics calculations.
; 
;-
pro red::fitscube_missing, filename $
                           , force = force $
                           , missing_type = missing_type $
                           , noflip = noflip $
                           , oname = oname $
                           , overwrite = overwrite

  inam = red_subprogram(/low, calling = inam1)
  
  ;; Make prpara
  red_make_prpara, prpara, filename
  red_make_prpara, prpara, force
  red_make_prpara, prpara, missing_type
  red_make_prpara, prpara, oname

  hdr = headfits(filename)
  
  ;; Information about processing steps in the formation of the input
  ;; file. 
  prprocs = fxpar(hdr, 'PRPROC*')
  prparas = fxpar(hdr, 'PRPARA*')

  if ~keyword_set(force) then begin
    ;; Check that it is not already done
    pos = where(strmatch(prprocs, inam), Nmatch)
    if Nmatch gt 0 then begin
      print
      print, inam + ' : The input file is already done:'
      print, filename
      print, inam + " : Use /force to do it again."
      return
    endif
  endif

  ;; If we provided oname, then copy the file to there and operate on
  ;; the copy. 
  if n_elements(oname) eq 0 then begin
    fname = filename
  endif else begin
    print, inam+' : Copying original file to output file name.'
    file_copy, filename, oname, overwrite = overwrite
    fname = oname
    fxaddpar, hdr, 'FILENAME', file_basename(fname)
  endelse
  
  red_fitscube_open, fname, fileassoc, fitscube_info, /update

  ;; Cube dimensions
  Nx      = fitscube_info.dimensions[0]
  Ny      = fitscube_info.dimensions[1]
  Ntuning = fitscube_info.dimensions[2]
  Nstokes = fitscube_info.dimensions[3]
  Nscans  = fitscube_info.dimensions[4]

  iprogress = 0L
  Nprogress = long(Ntuning)*Nstokes*long(Nscans)
  old_set_missing_to = ''
  for iscan = 0L, Nscans-1 do begin
    for istokes = 0L, Nstokes-1 do begin
      for ituning = 0L, Ntuning-1 do begin

        if iprogress ne 0 then old_set_missing_to = set_missing_to

        red_progressbar, iprogress, Nprogress, /predict $
                         , 'Missing data in frame '+strjoin(strtrim([iscan, istokes, ituning],2),',')


        ;; Read a frame
        red_fitscube_getframe, fileassoc, frame $
                               , ituning = ituning $
                               , istokes = istokes $
                               , iscan = iscan
        
        ;; Set the padding to the type we want
        red_missing, frame, /inplace $
;                     , image_out = image_out $
                     , missing_type_wanted = missing_type $
                     , missing_type_used = set_missing_to
        
        if iprogress ne 0 && old_set_missing_to ne set_missing_to then begin
          print, inam + ' : set_missing_to changed! (Different padding in different frames.)'
          stop
        endif
        
;        frame = temporary(image_out)
        
        ;; Write the modified frame. We could skip doing this if frame
        ;; hasn't changed.
        red_fitscube_addframe, fileassoc, frame $
                               , ituning = ituning $
                               , istokes = istokes $
                               , iscan = iscan
        
        iprogress++
        
      endfor                    ; ituning
    endfor                      ; istokes
  endfor                        ; iscan
  
  ;; Add info about this step
  case strlowcase(set_missing_to) of
    'nan'    : prref = 'Padding set to NaN'
    'median' : prref = 'Padding set to median'
    else:
  end
  self -> headerinfo_addstep, hdr $
                              , prstep = 'PADDING-CONVERSION' $
                              , prpara = prpara $
                              , prref = prref $
                              , prproc = inam

  ;; Close the file and write the updated header
  red_fitscube_close, fileassoc, fitscube_info, newheader = hdr

  if ~keyword_set(noflip) then begin
    ;; Make a flipped version
    red_fitscube_flip, fname $
                       , flipfile = flipfile $
                       , overwrite = overwrite
  endif

end

dir = '/scratch/mats/2016.09.19/CRISP-aftersummer/cubes_nb/'
filename = dir+'nb_6302_2016-09-19T09:30:20_scans=2-8_stokes_corrected_im.fits'

a = crispred(/dev)
;a -> fitscube_missing, filename, oname = dir+'test.fits', /over, /force
;a -> fitscube_missing, dir+'test.fits', /over, /force
a -> fitscube_missing, dir+'test.fits', /over, /force, missing_type = 'nan'
;a -> fitscube_missing, dir+'test.fits', /over, /force, missing_type = 'median'

end
