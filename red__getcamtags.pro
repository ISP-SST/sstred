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
;-
pro red::getcamtags, dir = dir

    inam = 'red::getcamtags : '

    ptr_free,self.cam_tags
    tagfil = self.out_dir+'/camtags.idlsave'
    if(file_test(tagfil)) then begin
        restore, tagfil
        if n_elements(cam_tags) gt 0 then self.cam_tags = ptr_new(cam_tags, /NO_COPY)
        return
    endif
     
    if(~keyword_set(dir) && ptr_valid(self.dark_dir) ) then dir = *self.dark_dir

    for i=0, n_elements(*self.cam_channels)-1 do begin
        path_spec = dir + '/' + (*self.cam_channels)[i] + '/*'
        files = file_search(path_spec, count=nf)
        if( nf eq 0 || files[0] eq '' ) then begin
            print, inam + 'ERROR -> no frames found in [' + dir + '] for ' + (*self.cam_channels)[i]
            return
        endif
        ctag = red_camtag(files[0])
        if ptr_valid(self.cam_tags) then red_append, *self.cam_tags, ctag $
        else self.cam_tags = ptr_new(ctag, /NO_COPY)
    endfor
    
    cam_tags = *self.cam_tags

    save, file=tagfil, cam_tags

    return
  
end
