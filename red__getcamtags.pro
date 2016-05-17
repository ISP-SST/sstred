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
;   2014-11-29 : JdlCR, store cams in a save-file to speed up execution
; 
;-
pro red::getcamtags, dir = dir
  inam = 'red::getcamtags : '

  
  ptr_free,self.cam_tags
  tagfil = self.out_dir+'/camtags.idlsave'
  if(file_test(tagfil)) then begin
     restore, tagfil
     self.camttag = camt
     self.camrtag = camr
     self.camwbtag = camw
     if n_elements(cam_tags) gt 0 then self.cam_tags = ptr_new(cam_tags, /NO_COPY)
     return
  endif
     
  if(~keyword_set(dir)) then dir = self.pinh_dir

  ;; TT cam
  spawn, 'find ' + dir + '/' + self.camt + '/ | grep cam', files
  nf = n_elements(files)
   
  if(files[0] eq '') then begin
     print, inam + 'ERROR -> no frames found in '+dir+' for '+self.camt
     return
  endif
  self.camttag = red_camtag(files[0])
   
  ;; TR cam
  spawn, 'find ' + dir + '/' + self.camr + '/ | grep cam', files
  nf = n_elements(files)
   
  if(files[0] eq '') then begin
     print, inam + 'ERROR -> no frames found in '+dir+' for '+self.camr
     return
  endif
  self.camrtag = red_camtag(files[0])
   
  ;; WB cam
  spawn, 'find ' + dir + '/' + self.camwb + '/ | grep cam', files
  nf = n_elements(files)
   
  if(files[0] eq '') then begin
     print, inam + 'ERROR -> no frames found in '+dir+' for '+self.camwb
     return
  endif
  self.camwbtag = red_camtag(files[0])


  camt = self.camttag
  camr = self.camrtag
  camw = self.camwbtag
  
  FOR i=0, n_elements(*self.cam_channels)-1 DO BEGIN
    spawn, 'find ' + dir + '/' + (*self.cam_channels)[i] + '/ | grep cam', files
    nf = n_elements(files)
    if(files[0] eq '') then begin
        print, inam + 'ERROR -> no frames found in ' + dir + ' for ' + (*self.cam_channels)[i]
        return
    endif
    IF ptr_valid(self.cam_tags) THEN red_append, *self.cam_tags, red_camtag(files[0]) $
    ELSE self.cam_tags = ptr_new(red_camtag(files[0]), /NO_COPY)
  ENDFOR
  cam_tags = *self.cam_tags
  
  save, file=tagfil, camt, camr, camw, cam_tags

  return
end
