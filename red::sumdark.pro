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
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
pro red::sumdark, overwrite = overwrite, check = check
                                ;
  if(self.dodark eq 0B) then begin
     print, 'red::sumdark : ERROR : undefined dark_dir'
     return
  endif
                                ; 
                                ; create file list
  for ic = 0L, 2 do begin
     case ic of
        0: begin
           if(self.docamt eq 0B) then begin
              print, 'red::sumdark : Nothing to do for ', self.camt
           endif
           cam = self.camt
           doit = self.docamt
        end
        1: begin
           if(self.docamr eq 0B) then begin
              print, 'red::sumdark : Nothing to do for ', self.camr
           endif
           cam = self.camr
           doit = self.docamr
        end
        2: begin
           if(self.docamwb eq 0B) then begin
              print, 'red::sumdark : Nothing to do for ', self.camwb
           endif
           cam = self.camwb
           doit = self.docamwb
        end
     endcase
     if(~doit) then continue
                                ;
     spawn, 'find ' + self.dark_dir + '/' + cam + '/|grep im.ex', files
     nf = n_elements(files)
                                ;
                                ;search_chain = self.dark_dir+'/'+cam+'/*im*'
                                ;print, search_chain
                                ;files = file_search(search_chain, count = nf)
     if(files[0] eq '') then begin
        print, 'red::sumdark : ERROR, no files found in: '+self.dark_dir
        return
     endif
     
     ;;If everything is ok, sum the darks.
     head = fzhead(files[0])
     dark = red_sumfiles(files, check = check, summed = tmp)

                                ;
     ;; Normalize dark
                                ;dark = float(tmp / count)
                                ;
     camtag = (strsplit(file_basename(files[0]),'.',/extract))[0]
     outdir = self.out_dir+'/darks/'
     namout = outdir+camtag+'.dark'
     namout_sum = outdir+camtag+'.summed.0000001'
                                ;
     outheader = red_dark_h2h(files[nf-1], head, nf)
                                ;
     file_mkdir, outdir
     print, 'red::sumdark : saving ', namout
     fzwrite, dark, namout, 'mean dark'
                                ;
     print, 'red::sumdark : saving ', namout_sum
     fzwrite, tmp, namout_sum, outheader
                                ;
  endfor                        ; (ic loop)
                                ;
  self.done.sumdark = 1B
                                ;
  return
end
