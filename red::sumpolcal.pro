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
;    remove  : 
;   
;   
;   
;    ucam  : 
;   
;   
;   
;    check : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-09 : MGL. Worked around outdir1, outdir2, and lun when not
;                specifying /old.
; 
; 
; 
;-
pro red::sumpolcal, remove = remove, ucam = ucam, check=check, old = old
  inam = 'red::sumpolcal : '
                                ;
  if(~self.dopolcal) then begin
     print, inam + 'Error -> undefined polcal_dir in '+self.filename
  endif
                                ;
                                ; loop cameras
  for ic = 0, 1 do begin
     firsttime = 1B
     case ic of
        0 : begin
           cam = self.camt
           doit = self.docamt
        end
                                ;
        1 : begin
           cam = self.camr
           doit = self.docamr
        end
     endcase
                                ;
     IF(keyword_set(ucam)) THEN BEGIN
        IF cam NE ucam THEN BEGIN
           print, inam + 'skipping '+cam
           continue
        endif
     endif
     if(~doit) then begin
        print, inam+'nothing to do for '+cam
        continue
     endif
                                ;
                                ; find files
     spawn, 'find ' + self.polcal_dir + '/' + cam + '/ | grep im.ex', files
     
     nt = n_elements(files)
     if(files[0] eq '') then begin
        print, inam+'no files found in '+self.polcal_dir+'/'+cam+', skipping camera!'
        continue
     endif
                                ;
                                ; sort files based on image number
     files = red_sortfiles(files)

                                ;
                                ; get image states
     pstat = red_getstates_polcal(files)
     if(keyword_set(remove)) then pstat.star[*] = red_flagchange(pstat.qw)
                                ;
                                ; unique states
     ustate = pstat.state[uniq(pstat.state, sort(pstat.state))]
     ns = n_elements(ustate)
                                ;
     print, inam+'found '+red_stri(ns)+' states for '+cam
                                ;
                                ; sum images
     outdir = self.out_dir + '/polcal_sums/' + cam+'/'
     file_mkdir, outdir
     if(keyword_set(old)) then begin
        outdir1 = self.out_dir + '/polcal_old/' + cam+'/'
        file_mkdir, outdir1
     endif
     camtag = (strsplit(file_basename(files[0]), '.',/extract))[0]
                                ;
     ddf = self.out_dir+'/darks/'+camtag+'.summed.0000001'
     if(file_test(ddf)) then begin
        ddh = fzhead(ddf) 
        pos = strsplit(ddh,' ')
        ddh = strmid(ddh, pos[1], pos[n_elements(pos)-1])
     endif else ddh = ' '
                                ;
     for ss = 0L, ns - 1 do begin ;ns-1
        pos = where((pstat.state eq ustate[ss]) AND (pstat.star eq 0B), count)
        if(count eq 0) then continue
                                ;
        print, inam+'adding frames for '+cam+' -> '+ustate[ss]
        pcal = red_sumfiles(files[pos], time = time, summed = polcalsum,check=check)
                                ;
        head = fzhead(files[pos[count/2-1]])
        h = head
                                ;
                                ; Save average of the sum (for IDL polcal)
                                ;
        namout = camtag+'.'+ustate[ss]+'.'+pstat.nums[pos[0]]
        print, inam+'saving '+outdir+namout 
        fzwrite, float(pcal), outdir+namout, h
                                ;
                                ; Save un-normalized sums for polcal
                                ;
        if(keyword_set(old)) then begin
           h1 = red_get_polcalheader(files[pos[0]], head, n_elements(pos), pstat, date = date)
           namout1 = camtag + '_im' + date + '.fp0.1234.1234.'+$
                     strcompress(string(float(strmid(pstat.qw[pos[0]],2,3)),format='(F6.2)'),/remove)+'.'+ $
                     pstat.lc[pos[0]]+'.'+pstat.nums[pos[0]]
           print, inam + 'saving ' + outdir1 + namout1
           pref = pstat.pref[pos[0]]
           outdir2 = outdir1 + pref + '/'
           file_mkdir, outdir2 
           fzwrite, polcalsum, outdir2 + namout1, h1
        endif
                                ;
        if(firsttime eq 1B) then begin
           if(keyword_set(old)) then begin
              openw, lun, outdir2+'/'+camtag+'_lg'+date, /get_lun, width = 300
           endif
           firsttime = 0B
        endif
        if(keyword_set(old)) then printf,lun, h1
     endfor                     ; (ss)
     ;; ;
     if(file_test(ddf)) then begin
        if(keyword_set(old)) then begin
           ddlink = camtag+'_dd'+date+'.'+red_stri(long((strsplit(ddh,' ',/extract))[0]), ni='(I07)')
          if(file_test(outdir2+ddlink)) then spawn,'rm '+outdir2+ddlink
          spawn, 'ln -s '+ddf+' '+outdir2+ddlink
       endif
     endif
     if(keyword_set(old)) then begin
        printf, lun, ddh
        free_lun, lun
     endif
  endfor                        ;(ic)
  return
end
