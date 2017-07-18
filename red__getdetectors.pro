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
;-
pro red::getdetectors, dir = dir

  inam = 'red::getdetectors : '

  ptr_free,self.detectors
  tagfil = self.out_dir+'/detectors.idlsave'
  if(file_test(tagfil)) then begin
    restore, tagfil
    if n_elements(detectors) gt 0 then self.detectors = ptr_new(detectors, /NO_COPY)
    return
  endif
  
  if(~keyword_set(dir) && ptr_valid(self.dark_dir) ) then dir = *self.dark_dir

  for i=0, n_elements(*self.cameras)-1 do begin
    path_spec = dir + '/' + (*self.cameras)[i] + '/*'
    files = file_search(path_spec, count=nf)
    if( nf eq 0 || files[0] eq '' ) then begin
      print, inam + 'ERROR -> no frames found in [' + dir + '] for ' + (*self.cameras)[i]
      ctag = red_detectorname(files[0])
    endif else begin
      ctag = red_detectorname(files[0])
    endelse
    if ptr_valid(self.detectors) then red_append, *self.detectors, ctag $
    else self.detectors = ptr_new(ctag, /NO_COPY)
  endfor
  
  detectors = *self.detectors

  save, file=tagfil, detectors

  return
  
end
