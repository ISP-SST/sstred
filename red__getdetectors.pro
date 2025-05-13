; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;    dir  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
;
;   2014-11-29 : JdlCR, store cams in a save-file to speed up
;                execution.
;
;   2016-05-22 : THI. Removed CRISP-specific tags.
;
;   2016-05-22 : THI. Bugfix.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.     
;
;   2016-08-31 : MGL. Use red_detectorname instead of red_getcamtag.
;                Rename camtags.idlsave to detectors.idlsave.     
;
;   2017-08-09 : MGL. Better default dir for finding detector tags. 
; 
;   2024-12-18 : MGL. Look for detector info in already summed files.
;
;-
pro red::getdetectors, dir = dir

  ;; Name of this program
  inam = red_subprogram(/low, calling = inam1)

  ptr_free,self.detectors
  tagfil = self.out_dir+'/detectors.idlsave'
  if(file_test(tagfil)) then begin
    restore, tagfil
    if n_elements(detectors) gt 0 then self.detectors = ptr_new(strtrim(detectors, 2), /NO_COPY)
    return
  endif

  if n_elements(dir) eq 0 then begin
    case 1 of
      ;; If any of these directories do not exist you have to do
      ;; something by hand anyway. Like using calibration data from
      ;; another day. Then you can also run getdetectors with the dark
      ;; directory from that other day as the keyword argument.
      ptr_valid(self.dark_dir) : dir = *self.dark_dir
      ptr_valid(self.pinh_dirs) : dir = *self.pinh_dirs
      ptr_valid(self.flat_dir) : dir = *self.flat_dir
      file_test('darks') : dir = 'darks/'
      else : begin
        print, inam + ' : Cannot find detector names.'
        retall
      end
    endcase
  endif

  if n_elements(dir) eq 1 && dir eq 'darks/' then begin

    ;; Summed darks
    files = file_search('darks/*fits', count=nf)

    ;; Both cameras and detectors are in the headers
    file_cams = red_fitsgetkeyword_multifile(files, 'CAMERA', count = cnt)
    file_dets = red_fitsgetkeyword_multifile(files, 'DETECTOR', count = cnt)

    ;; Remove duplicates if any, could be darks with different exposure times.
    file_cams = red_uniquify(file_cams, indx = indx)
    file_dets = file_dets[indx]

    ;; Sort the detectors in the order of self.cameras
    match2, *self.cameras, file_cams, suba, subb
    detectors = file_dets[suba]

    self.detectors = ptr_new(strtrim(detectors, 2), /NO_COPY)
    
  endif else begin
    
    ;; Raw data
    
    for i=0, n_elements(*self.cameras)-1 do begin
      path_spec = dir + '/' + (*self.cameras)[i] + '/*'
      files = file_search(path_spec, count=nf)
      if( nf eq 0 || files[0] eq '' ) then begin
        print, inam + ' : ERROR -> no frames found in [' + dir + '] for ' + (*self.cameras)[i]
        ctag = red_detectorname(files[0])
      endif else begin
        ctag = red_detectorname(files[0])
      endelse
      if ptr_valid(self.detectors) then red_append, *self.detectors, ctag $
      else self.detectors = ptr_new(ctag, /NO_COPY)
    endfor
    
    detectors = *self.detectors

  endelse
  
  save, file=tagfil, detectors

  print, inam + ' : Write to '+tagfil
  print,'   ' + *self.cameras + ' = ' + detectors, format='(a0)'
  
  return
  
end
