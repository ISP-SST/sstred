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
; 
;-
pro red::getcamtags, dir = dir

  if(~keyword_set(dir)) then dir = self.dark_dir

  inam = 'red::getcamtags : '
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
   
  return
end
