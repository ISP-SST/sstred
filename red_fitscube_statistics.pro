; docformat = 'rst'

;+
; Calculate statistics for the frames in a fitscube file.
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
; :Params:
; 
;   filename : in, type=string
; 
;   frame_statistics : out, optional, type=array
;
;      Frame by frame statistics.
;
;   cube_statistics : out, optional, type=array
;
;      Statistics for the entire cube.
;
; :Keywords:
; 
;   angles : in, optional, type=array
;   
;      Rotation angles for the frames in the cube.
;   
;   cube_comments : out, optional, type=strarr
;   
;      
;   
;   full : in, optional, type=array
;   
;      Parameters that determine the array size of the frames after
;      rotations and shifts.
;   
;   grid : in, optional, type=array
;   
;      Stretch vectors.
;   
;   origNx : in, optional, type=float 
;   
;      Frame size before rotation and shifts.
;   
;   origNy :  in, optional, type=float
;   
;     Frame size before rotation and shifts.
;   
;   percentiles : in, out, optional
;   
;     The percentiles to calculate.
;   
;   write : in, optional, type=boolean
;   
;     Write the statistics to the fitscube. 
;   
;   remove_only : in, optional, type=boolean
;
;     Only remove existing statistics keywords.   
;   
;   shifts : in, optional, type=array
;   
;     Shift vectors.
;
;   update : in, optional, type=boolean
;
;     Calculate (new) statistics only if there are statistics in the
;     header already.
; 
; :History:
; 
;   2019-03-28 : MGL. First version.
; 
;   2019-04-05 : MGL. Use red_fitscube_open and red_fitscube_close. 
; 
;   2019-08-26 : MGL. New keyword remove_only. 
; 
;   2020-01-21 : MGL. New keyword update.
; 
;-
pro red_fitscube_statistics, filename, frame_statistics, cube_statistics $
                             , angles = angles $
                             , full = full $
                             , grid = grid $
                             , origNx = origNx $
                             , origNy = origNy $
                             , percentiles = percentiles $
                             , remove_only = remove_only $
                             , shifts = shifts $
                             , update = update $
                             , write = write $
                             , cube_comments = cube_comments 

  if n_elements(percentiles) eq 0 then percentiles = [.01, .10, .25, .50, .75, .90, .95, .98, .99]
  
  ;; Use angles, shifts, full (if given) to calculate masks that
  ;; select the rotated and shifted area. Can we do something with
  ;; grid as well? Maybe use magnitude of grid shifts to calculate a
  ;; margin around the area?

  ;; Open the file and set up an assoc variable.
  red_fitscube_open, filename, fileassoc, fitscube_info ;$
;                     , lun = lun

  if keyword_set(remove_only) or keyword_set(write) or keyword_set(update) then begin
    ;; Remove any existing statistics keywords from the file header. 
    hdr = fitscube_info.header
    ;; Search for DATA* keywords
    dindx = where(strmid(hdr, 0, 4) eq 'DATA', Ndata)
    ;; Loop through the keywords backwards so we don't move them before
    ;; they are deleted.
    for idata = Ndata-1, 0, -1 do begin
      keyword = strtrim(strmid(hdr[dindx[idata]], 0, 8), 2)
      red_fitsdelkeyword, hdr, keyword
    endfor                      ; idata
    removed_stuff = Ndata gt 0 
  endif else removed_stuff = 0 

  if keyword_set(remove_only) then begin
    ;; We are done. Close if needed and return.
    if removed_stuff then red_fitscube_close, fileassoc, fitscube_info, newheader = hdr
    return
  endif

  ;; If we only want to update existing statistics, then we want to
  ;; return now if there weren't any in the file.
  if keyword_set(update) and ~removed_stuff then return

  ;; If we get this far, we do want to compute statistics!
  
  Nx      = fitscube_info.dimensions[0]
  Ny      = fitscube_info.dimensions[1]
  Ntuning = fitscube_info.dimensions[2]
  Nstokes = fitscube_info.dimensions[3]
  Nscans  = fitscube_info.dimensions[4]

  if n_elements(angles) eq Nscans then begin
    if n_elements(origNx) gt 0 then Nxx = origNx else Nxx = Nx
    if n_elements(origNy) gt 0 then Nyy = origNy else Nyy = Ny
  endif
;  else begin
;    ;; Indices of all image pixels, i.e., no masking.
;    mindx = lindgen(Nx, Ny)
;  endelse

  ;; Calculate statistics for the individual frames
  iprogress = 0
  Nprogress = Nscans * Ntuning * Nstokes
  for iscan = 0L, Nscans - 1 do begin

    if n_elements(angles) eq Nscans then begin
      
      ;; Make a mask by rotating and shifting an image of unit values
      ;; the same way as the images from this scan. 
      
      if n_elements(shifts) eq 0 then begin
        dx = 0
        dy = 0
      endif else begin
        dx = shifts[0,iscan]
        dy = shifts[1,iscan]
      endelse
      
      mask = make_array(Nxx, Nyy, /float, value = 1.) 
      mask = red_rotation(mask, angles[iscan], dx, dy, background = 0, full = full)
      mindx = where(mask gt 0.99)

    endif 

    for ituning = 0L, Ntuning - 1 do begin 
      for istokes = 0L, Nstokes-1 do begin

        red_progressbar, iprogress, Nprogress $
                         , /predict $
                         , 'Calculate frame by frame statistics'

        red_fitscube_getframe, fileassoc, frame $
                               , iscan = iscan, ituning = ituning, istokes = istokes
        
        if (iscan eq 0) and (ituning eq 0) and (istokes eq 0) then begin
          ;; Set up the array if it's the first frame
          if n_elements(mindx) gt 0 then begin
            ;; Calculate statistics based on rotating mask
            frame_statistics = red_image_statistics_calculate(frame[mindx])
          endif else begin
            ;; Calculate statistics based on deselecting missing-data
            ;; pixels in each frame.
            red_missing, frame, indx_data = indx_data
            frame_statistics = red_image_statistics_calculate(frame[indx_data])
          endelse
          frame_statistics = replicate(temporary(frame_statistics), Ntuning, Nstokes, Nscans)
        endif else begin
          if n_elements(mindx) gt 0 then begin
            ;; Calculate statistics based on rotating mask
            frame_statistics[ituning, istokes, iscan] = red_image_statistics_calculate(frame[mindx])
          endif else begin
            ;; Calculate statistics based on deselecting missing-data
            ;; pixels in each frame.
            red_missing, frame, indx_data = indx_data
            frame_statistics[ituning, istokes, iscan] = red_image_statistics_calculate(frame[indx_data])
          endelse
        endelse 

        iprogress++
        
      endfor                    ; istokes
    endfor                      ; ituning
  endfor                        ; iscan

  
  if Nstokes eq 1 then frame_statistics = reform(frame_statistics)
  
  if keyword_set(write) or keyword_set(update) or arg_present(cube_statistics) then begin
    
    ;; Accumulate a histogram for the entire cube, use to calculate
    ;; percentiles.
    cubemin  = min(frame_statistics.datamin)
    cubemax  = max(frame_statistics.datamax)
    Nbins = 2L^16               ; Use many bins!
    binsize = (cubemax - cubemin) / (Nbins - 1.)
    hist = lonarr(Nbins)
    iprogress = 0
    for iscan = 0L, Nscans - 1 do begin


      if n_elements(angles) eq Nscans then begin
        
        ;; Make a mask by rotating and shifting an image of unit values
        ;; the same way as the images from this scan. 
        
        if n_elements(shifts) eq 0 then begin
          dx = 0
          dy = 0
        endif else begin
          dx = shifts[0,iscan]
          dy = shifts[1,iscan]
        endelse
        
        mask = make_array(Nxx, Nyy, /float, value = 1.) 
        mask = red_rotation(mask, angles[iscan], dx, dy, background = 0, full = full)
        mindx = where(mask gt 0.99)
      endif else begin
        mask = 0 * frame+1
        mindx = indgen(n_elements(frame))
      endelse

      for ituning = 0L, Ntuning - 1 do begin 
        for istokes = 0L, Nstokes-1 do begin

          red_progressbar, iprogress, Nprogress $
                           , /predict $
                           , 'Accumulate histogram'

          red_fitscube_getframe, fileassoc, frame $
                                 , iscan = iscan, ituning = ituning, istokes = istokes
          
          if n_elements(mask) eq 0 then begin
            hist += histogram(float(frame), min = cubemin, max = cubemax, Nbins = Nbins, /nan)
          endif else begin
            if ~array_equal(size(mask,/dim),size(frame,/dim)) then stop ; Size mismatch?
            hist += histogram(float(frame[mindx]), min = cubemin, max = cubemax, Nbins = Nbins, /nan)
          endelse
          
          iprogress++
          
        endfor                  ; istokes
      endfor                    ; ituning
    endfor                      ; iscan

    ;; Calculate cube statistics from the histogram and the individual
    ;; frame statistics
    cube_statistics = red_image_statistics_combine(frame_statistics $
                                                   , hist = hist $
                                                   , comments = cube_comments $
                                                   , binsize = binsize)
  endif

  if ~keyword_set(write) and ~keyword_set(update) then begin
    red_fitscube_close, fileassoc, fitscube_info
    return
  endif 

  if removed_stuff then begin
    red_fitscube_close, fileassoc, fitscube_info, newheader = hdr
  endif else begin
    red_fitscube_close, fileassoc, fitscube_info
  endelse
  
  ;; Write the statistics to the fitscube file

  if Ntuning gt 1 then red_append, axis_numbers, 3
  if Nstokes gt 1 then red_append, axis_numbers, 4
  if Nscans gt 1 then red_append, axis_numbers, 5
  
;    if Nstokes gt 1 then begin
;      axis_numbers = [3, 4, 5]  ; (Ntuning, Nstokes, Nscans)
;    endif else begin
;      axis_numbers = [3, 5]     ; (Ntuning, Nscans)
;    endelse

  
  for itag = n_tags(frame_statistics[0])-1, 0, -1 do begin

    itags = where((tag_names(frame_statistics[0]))[itag] eq tag_names(cube_statistics), Nmatch)
    itagc = where((tag_names(frame_statistics[0]))[itag] eq tag_names(cube_comments), Nmatch)

    if Nmatch eq 1 then $
       red_fitscube_addvarkeyword, filename $
                                   , (tag_names(frame_statistics[0]))[itag] $
                                   , frame_statistics.(itag) $
                                   , anchor = anchor $
                                   , keyword_value = cube_statistics.(itags) $
                                   , comment = cube_comments.(itagc) $
                                   , axis_numbers = axis_numbers

  endfor                        ; itag

end
