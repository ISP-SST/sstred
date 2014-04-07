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
;    overwrite  : 
;   
;   
;   
;    check  : 
;   
;   
;    cams : in, optional, type=strarr
; 
;      A list of cameras (or rather camera subdirs).
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-23 : MGL. Added support for logging.
; 
;   2013-08-27 : MGL. Let the subprogram find out its own name. Added
;                self.done to the selfinfo.
; 
;   2013-09-19 : MGL. Make the sum long integer.
; 
;   2014-03-21 : MGL. Allow for multiple dark_dir. Make it work for
;                blue cameras. New keyword cams.
;
;   2014-04-07 : THI. Bugfix: look for darks in dark_dir.
;
;-
pro red::sumdark, overwrite = overwrite, check = check, cams = cams

  ;; Defaults
  if n_elements(overwrite) eq 0 then overwrite = 0
  if n_elements(check) eq 0 then check = 0

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo1
  help, /struct, self.done, output = selfinfo2 
  red_writelog, selfinfo = [selfinfo1, selfinfo2]

  ;; Default: do CRISP cameras
  if n_elements(cams) eq 0 then begin
     DoCRISP = 1
     Ncams = 3  
  endif else begin
     DoCRISP = 0
     Ncams = n_elements(cams)
  endelse

  if(self.dodark eq 0B) then begin
     print, inam+' : ERROR : undefined dark_dir'
     return
  endif

  ;; create file list
  for ic = 0L, Ncams-1 do begin

     if DoCRISP then begin
        case ic of
           0: begin
              if(self.docamt eq 0B) then begin
                 print, inam+' : Nothing to do for ', self.camt
              endif
              cam = self.camt
              doit = self.docamt
           end
           1: begin
              if(self.docamr eq 0B) then begin
                 print, inam+' : Nothing to do for ', self.camr
              endif
              cam = self.camr
              doit = self.docamr
           end
           2: begin
              if(self.docamwb eq 0B) then begin
                 print, inam+' : Nothing to do for ', self.camwb
              endif
              cam = self.camwb
              doit = self.docamwb
           end
        endcase
     endif else begin
        cam = cams[ic]
        doit = 1
     endelse

     if ~doit then continue

     if strmatch(cam,'Crisp-?') then begin
        ;; Matches Crisp cameras
;        spawn, 'find ' + *self.dark_dir + '/' + cam + '/|grep im.ex', files
        files = file_search(*self.dark_dir + '/' + cam + '/cam*.im.ex.*', count = nf)
    endif else begin
        ;; Matches blue cameras
;        spawn, 'find ' + *self.dark_dir + '/' + cam + '/|grep im.im', files
        files = file_search(*self.dark_dir + '/' + cam + '/cam*.im.im.*', count = nf)
   endelse
;    nf = n_elements(files)

     ;;search_chain = self.dark_dir+'/'+cam+'/*im*'
     ;;print, search_chain
     ;;files = file_search(search_chain, count = nf)
;     if(files[0] eq '') then begin
;        print, inam+' : ERROR, no files found in: '+self.dark_dir
;        return
;     endif

     print, inam+' : Found '+red_stri(nf)+' files in: '
     print, '         '+ *self.dark_dir + '/' + cam + '/'

     if(files[0] eq '') then begin
         print, inam+' : ERROR : '+cam+': no files found in: '
         print, inam+'                      '+*self.dark_dir
         print, inam+' : ERROR : '+cam+': skipping camera!'
        continue
     endif

     
     ;;If everything is ok, sum the darks.
     head = fzhead(files[0])
     dark = red_sumfiles(files, check = check, summed = darksum)

     ;; Normalize dark
     ;;dark = float(tmp / count)
     camtag = (strsplit(file_basename(files[0]),'.',/extract))[0]
     outdir = self.out_dir+'/darks/'
     namout = outdir+camtag+'.dark'
     namout_sum = outdir+camtag+'.summed.0000001'

     outheader = red_dark_h2h(files[nf-1], head, nf)

     file_mkdir, outdir
     print, inam+' : saving ', namout
     fzwrite, dark, namout, 'mean dark'

     print, inam+' : saving ', namout_sum
     fzwrite, long(darksum), namout_sum, outheader

  endfor                        ; (ic loop)

  self.done.sumdark = 1B

  return
end
