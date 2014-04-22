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
;    ustat  : 
;   
;   
;   
;    old  : 
;   
;   
;   
;    remove  : 
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
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name. Added self.done to the
;                selfinfo.
; 
;   2013-08-29 : MGL. Remove *lcd* files from file list.
; 
;   2013-10-30 : MGL. Get correct prefilter state for WB flats in the
;                presence of the focus file name field. 
; 
;   2013-12-10 : PS  Adapt for multiple flat directories; add lim
;                keyword (passthrough to red_sumfiles)
;
;   2013-12-13 : PS  Only store raw sums if requested.  Save
;                unnormalized sum and add correct number of summed
;                files in header
;
;   2014-01-23 : MGL. Use red_extractstates instead of red_getstates
;                and local extraction of info from file names.
;
;-
PRO red::sumflat, overwrite = overwrite, ustat = ustat, old = old, $
                  remove = remove, cam = ucam, check = check, lim = lim, $
                  store_rawsum = store_rawsum

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo1
  help, /struct, self.done, output = selfinfo2 
  red_writelog, selfinfo = [selfinfo1, selfinfo2]

  if(self.doflat eq 0B) then begin
     print, inam+' : ERROR : undefined flat_dir'
     return
  endif
  outdir = self.out_dir + '/' + 'flats/'
  file_mkdir, outdir
  outdir1 = self.out_dir + '/' + 'flats/summed/'
  IF keyword_set(store_rawsum) THEN file_mkdir, outdir1

  ;; NB cameras
  for ic = 0L, 1 do begin

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
              print, inam+' : Nothing to do for '+ self.camr
           endif
           cam = self.camr
           doit = self.docamr
        end
     endcase
     if(~doit) then continue

     if(n_elements(ucam) ne 0) then begin
        if(ucam ne cam) then begin
           print, inam + 'skipping camera -> '+cam+' (!= '+ucam+')'
           continue
        endif
     endif

 ;    spawn, 'find ' + self.flat_dir + '/' + cam + '/ | grep "camX" | grep -v ".lcd."', files
 ;    nf = n_elements(files)
     
     files = file_search(*self.flat_dir + '/' + cam + '/camX*')
     files = files(where(strpos(files, '.lcd.') LT 0, nf))
     
     print, inam+' : Found '+red_stri(nf)+' files in: '
     print, '         '+ *self.flat_dir + '/' + cam + '/'

     if(files[0] eq '') then begin
         print, inam+' : ERROR : '+cam+': no files found in: '
         print, inam+'                      '+*self.flat_dir
         print, inam+' : ERROR : '+cam+': skipping camera!'
        continue
     endif

     files = red_sortfiles(files)

     ;; We do not actually need all the fields in the stat structure. The
     ;; fullstate and star fields are enough.
     red_extractstates, files, /basename, states = stat, fullstate = fullstate
;     stat = red_getstates(files)
     camtag = (strsplit(file_basename(files[0]), '.', /extract))[0]

     if(~file_test(self.out_dir+'/darks/'+camtag+'.dark')) then begin
        print, inam+' : no darks found for ', camtag
        continue
     endif
     dd = f0(self.out_dir+'/darks/'+camtag+'.dark')

     ;; Unique states
     if(~keyword_set(ustat)) then ustat = fullstate[uniq(fullstate,sort(fullstate))]
;     if(~keyword_set(ustat)) then ustat = stat.state[uniq(stat.state, sort(stat.state))]
     ns = n_elements(ustat)

     ;; Flag first frame after tuning the FPI
;     if(keyword_set(remove)) then red_flagtuning, states
     if(keyword_set(remove)) then red_flagtuning, stat
     
     ;; Loop and sum
     for ss = 0L, ns - 1 do begin
        pos = where((stat.state eq ustat[ss]) AND (stat.star eq 0B), count)
        if(count eq 0) then continue

        outname = camtag + '.' + ustat[ss] + '.flat'
        outname1 = camtag + '.' + ustat[ss] + '.summed.flat'


        ;; If file does not exist, do sum!
        if(keyword_set(overwrite) OR ~file_test(outdir + outname)) then begin
           print, inam+' : adding flats for state -> ' + ustat[ss]
           IF(keyword_set(check)) THEN BEGIN
              openw, lun, self.out_dir + '/flats/'+camtag+'.'+ustat[ss]+'.discarded.txt', width = 500, /get_lun
           endif
           
           ;; Sum files
           IF(keyword_set(check)) THEN BEGIN
                 ;;; if summing from two directories, same frame
                 ;;; numbers from two directories are consecutive ans
                 ;;; produce false drop info.  Re-sort them
               tmplist = files[pos]
               tmplist = tmplist(sort(tmplist))
               flat = red_sumfiles(tmplist, time = time, check = check, $
                                   lun = lun, lim = lim, summed = summed, NSUM = nsum)
           ENDIF ELSE BEGIN 
               flat = red_sumfiles(files[pos], time = time, summed = summed, NSUM = nsum)
           ENDELSE
           
           ;; Remove dark
           flat -= dd
           flat1 = long(temporary(summed))
           
           ;; Output the raw (if requested) and averaged flats

           IF keyword_set(store_rawsum) THEN BEGIN
               headerout = 't='+time+' n_sum='+red_stri(nsum)
               print, inam+' : saving ' + outdir1+outname1
               fzwrite, flat1, outdir1+outname1, headerout
           ENDIF
           
           headerout = 't='+time+' n_aver='+red_stri(nsum)+' darkcorrected'
           print, inam+' : saving ' + outdir+outname
           fzwrite, float(flat), outdir+outname, headerout

           IF(keyword_set(check)) THEN free_lun, lun

        endif else begin
           print, inam+' : file exists: ' + outdir+outname + $
                  ' , skipping! (run sumflat, /overwrite to recreate)'
        endelse
     endfor
  endfor                        ;(ic loop)

  ;; Now WB camera
  if(self.docamwb eq 0B) then begin
     print, inam+' : Nothing to do for '+ self.camwb
     self.done.sumflat = 1B
     return
  endif

  ;; Search files
  cam = self.camwb

  if(n_elements(ucam) ne 0) then begin
     if(ucam ne cam) then begin
        print, inam + 'skipping camera -> '+cam+' (!= '+ucam+')'
        return
     endif
  endif

  files = file_search(*self.flat_dir + '/' + cam + '/camX*')
  
  if(files[0] eq '') then begin
      print, inam+' : ERROR : '+cam+': no files found in:'
      print, inam+'                      '+*self.flat_dir
      print, inam+' : ERROR : '+cam+': skipping camera!'
     self.done.sumflat = 1B
     return
  endif

  files = red_sortfiles(files)
  nt = n_elements(files)

  ;; Get prefilters state
;  wstat = strarr(nt)
;  for jj = 0L, nt - 1 do begin
;     dum = strsplit(file_basename(files[jj]), '.', /extract)
;     Ndum = n_elements(dum)
;     if(Ndum eq 7) then wstat[jj] = dum[Ndum-3] else wstat[jj] = dum[Ndum-5]
;  endfor
  red_extractstates, files, /basename, pref = wstat, cam = camtag
  uwstat = wstat[uniq(wstat, sort(wstat))]
  camtag = camtag[0]

;  camtag = (strsplit(file_basename(files[0]), '.', /extract))[0]
  outdir = self.out_dir + '/' + 'flats/'

  ;; Load dark
  if(~file_test(self.out_dir+'/darks/'+camtag+'.dark')) then begin
     print, inam+' : no darks found for ', camtag
     return
  endif

  dd = f0(self.out_dir+'/darks/'+camtag+'.dark')

  for ss = 0L, n_elements(uwstat) - 1 do begin
     pos = where(wstat eq uwstat[ss], count)

     outname = camtag + '.' + uwstat[ss] + '.flat'
     outname1 = camtag + '.' + uwstat[ss] + '.summed.flat'

     if(file_test(outdir+outname) AND ~keyword_set(overwrite)) then begin
        print, inam+' : file exists: ' + outdir+outname +$
               ' , skipping! (run sumflat, /overwrite to recreate)'
        continue
     endif

     if count eq 0 then begin
        print, inam+' : no files found for WB state -> '+ uwstat[ss]
        continue
     endif

     print, inam+' : adding flats for WB state -> '+uwstat[ss]
     flat = red_sumfiles(files[pos], time = time, check=check, summed = summed, NSUM = nsum)

     ;; Remove dark
     flat1 = long(temporary(summed))
     flat -= dd

     ;; Output the raw (if requested) and averaged flats
     IF keyword_set(store_rawsum) THEN BEGIN
         headerout = 't='+time+' n_sum='+red_stri(nsum)
         print, inam+' : saving ' + outdir1+outname1
         fzwrite, flat1, outdir1+outname1, headerout
     ENDIF
     
     headerout = 't='+time+' n_aver='+red_stri(nsum)+' darkcorrected'
     print, inam+' : saving ' + outdir+outname
     fzwrite, float(flat), outdir+outname, headerout

     
  endfor
  self.done.sumflat = 1B

  return
end  
