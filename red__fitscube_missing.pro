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
;     'median' (of each frame). Byt default, set it to the opposite of
;     what it appears to be in the input cube.
;
;   noflip  : in, optional, type=boolean
;
;     Don't make a flipped version.
;
;   nostatistics  : in, optional, type=boolean
;
;     Don't calculate statistics.
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
;     2020-06-29 : MGL. First version.
; 
;-
pro red::fitscube_missing, filename $
                           , force = force $
                           , missing_type = missing_type $
                           , noflip = noflip $
                           , nostatistics = nostatistics $
                           , oname = oname $
                           , overwrite = overwrite

  inam = red_subprogram(/low, calling = inam1)
  
  ;; Make prpara
  red_make_prpara, prpara, filename
  red_make_prpara, prpara, force
  red_make_prpara, prpara, missing_type
  red_make_prpara, prpara, nostatistics
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
        
        ;; Heuristics for identifying the missing-data pixels. They
        ;; would be connected to the outermost rows and columns and
        ;; they would all have the same value, usually the median of
        ;; the data pixels. We can't simply assume all pixels with
        ;; this value are missing data as there *could* be interior
        ;; pixels with this exact value.

        ;; What kind of padding do we have now?
        
        currently_constant = frame[ 0,  0] eq frame[ 0, -1] and $
                             frame[ 0,  0] eq frame[-1,  0] and $
                             frame[ 0,  0] eq frame[-1, -1]

        currently_nans = ~( finite(frame[ 0,  0]) or $
                            finite(frame[ 0, -1]) or $
                            finite(frame[-1,  0]) or $
                            finite(frame[-1, -1]) $
                          )

        
        if n_elements(missing_type) gt 0 then begin

          ;; missing_type specified
          set_missing_to = missing_type
          
        endif else begin

          ;; missing_type not specified, change to the opposite of
          ;; what we have
          case 1 of

            currently_nans     : set_missing_to = 'median'
          
            currently_constant : set_missing_to = 'nan'
          
            else : begin
              print, inam + ' : Could not identify missing data pixels'
              print, 'Corner pixel intensities: ' $
                     , finite(frame[ 0,  0]) $
                     , finite(frame[ 0, -1]) $
                     , finite(frame[-1,  0]) $
                     , finite(frame[-1, -1])
              stop
            end
            
          endcase
        endelse

        if iprogress ne 0 && old_set_missing_to ne set_missing_to then begin
          print, inam + ' : set_missing_to changed! (Different padding in different frames.)'
          stop
        endif
        
        case strlowcase(set_missing_to) of

          'nan' : begin         ; Assume all corner pixels have the same value.
            
            if currently_nans then begin
              print, inam+' : Padding seems to be NaN already.'
              ;; Close the file and return
              red_fitscube_close, fileassoc, fitscube_info
              return
            endif
            
            ;; Find pixels with the same value as the corners
            mask = frame eq frame[ 0,  0] 

            ;; We want NaNs
            missing_value = !Values.F_NaN
            
          end

          'median' : begin

            if currently_constant then begin
              print, inam+' : Padding seems to be medians (or at least constant values) already.'
              ;; Close the file and return
              red_fitscube_close, fileassoc, fitscube_info
              return
            endif

            ;; Find non-finite pixels
            mask = ~finite(frame)

            ;; We want the median
            indx_data = where(~mask, Ndata)
            if Ndata eq 0 then stop
            missing_value = median(frame(indx_data))
            
          end

          else :  stop
        
        endcase

        ;; Use labal_region to find the pixels that are connected to
        ;; the corners
        mask = red_centerpic(mask, xSize = Nx+8, ySize = ny+8, z = 1) ; Add some rows and columns
        label = label_region(mask)                                    ; Label regions
        label = label[4:-5, 4:-5]                                     ; Remove two rows and columns
        mask = label eq label[0, 0]                                   ; Connected to corner
        
        ;; Change the value of those pixels to NaN
        indx = where(mask, Nwhere)
        if Nwhere gt 0 then frame[indx] = missing_value

        ;; Write the modified frame
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

  if ~keyword_set(nostatistics) then begin
    ;; Calculate statistics
    red_fitscube_statistics, fname, /write
  endif

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
;a -> fitscube_missing, filename, oname = dir+'test.fits', /nostatistics, /over, /force
;a -> fitscube_missing, dir+'test.fits', /nostatistics, /over, /force
a -> fitscube_missing, dir+'test.fits', /nostatistics, /over, /force, missing_type = 'nan'
;a -> fitscube_missing, dir+'test.fits', /nostatistics, /over, /force, missing_type = 'median'

end
