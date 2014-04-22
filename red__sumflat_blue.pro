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
;    check : 
;   
;   
;   
;    cams : in, optional, type=strarr
; 
;      A list of cameras (or rather camera subdirs).
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
;   2014-03-21 : MGL. Stripped down version for blue cameras.
;
;
;-
PRO red::sumflat_blue, overwrite = overwrite, ustat = ustat, old = old, $
                  check = check, lim = lim, $
                  store_rawsum = store_rawsum, cams = cams

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

  Ncams = n_elements(cams)
  ;; If cams keyword not given, look in flat_dir for camera directories!

  for ic = 0L, Ncams-1 do begin

     cam = cams[ic]

     files = file_search(*self.flat_dir + '/' + cam + '/cam*', count = Nfiles)
     
     if Nfiles eq 0 then begin
        print, inam+' : ERROR : '+cam+': no files found in: '
        print, inam+'                      '+*self.flat_dir
        print, inam+' : ERROR : '+cam+': skipping camera!'
        continue
     endif else begin
        print, inam+' : Found '+red_stri(Nfiles)+' files in: '
        print, '         '+ *self.flat_dir + '/' + cam + '/'
     endelse
     
     files = red_sortfiles(files)

     ;; We do not actually need all the fields in the stat structure. The
     ;; fullstate and star fields are enough.
     red_extractstates, files, /basename, states = stat, fullstate = fullstate
;     stat = red_getstates(files)
     camtag = (strsplit(file_basename(files[0]), '.', /extract))[0]

     if ~file_test(self.out_dir+'/darks/'+camtag+'.dark') then begin
        print, inam+' : no darks found for ', camtag
        continue
     endif
     dark = f0(self.out_dir+'/darks/'+camtag+'.dark')

     ;; Unique states
     if(~keyword_set(ustat)) then ustat = fullstate[uniq(fullstate,sort(fullstate))]
     Nstates = n_elements(ustat)

     if Nstates gt 1 then begin
        print, inam+' : Nstates > 1'
        help, ustat
        stop
     endif else begin
        print, inam+' : Nstates = 1'
     endelse
     
     
     ;; Loop and sum
     for istate = 0L, Nstates - 1 do begin
        pos = where(stat.state eq ustat[istate], count)
        if count eq 0 then continue

        if Nstates eq 1 then begin
           outname = camtag+'.flat'
           outname1 = camtag+'.summed.flat'
        endif else begin
           outname = camtag + '.' + ustat[istate] + '.flat'
           outname1 = camtag + '.' + ustat[istate] + '.summed.flat'
        endelse

        ;; If file does not exist, do sum!
        if keyword_set(overwrite) OR ~file_test(outdir + outname) then begin
           print, inam+' : adding flats for state -> ' + ustat[istate]
           if keyword_set(check) then begin
              openw, lun, outdir+red_strreplace(outname, 'flat', 'discarded.txt') $
                     , width = 500, /get_lun
           endif
           
           ;; Sum files
           if keyword_set(check) then begin
              ;; if summing from two directories, same frame numbers
              ;; from two directories are consecutive ans produce
              ;; false drop info. Re-sort them
              tmplist = files[pos]
              tmplist = tmplist(sort(tmplist))
              flat = red_sumfiles(tmplist, time = time, check = check, $
                                  lun = lun, lim = lim, summed = summed, nsum = Nsum)
           endif else begin 
               flat = red_sumfiles(files[pos], time = time, summed = summed, nsum = Nsum)
           endelse
           
           ;; Remove dark
           flat -= dark
           
           ;; Output the raw (if requested) and averaged flats
           if keyword_set(store_rawsum) then begin
               headerout = 't='+time+' n_sum='+red_stri(nsum)
               print, inam+' : saving ' + outdir1+outname1
               flat1 = long(temporary(summed))
               fzwrite, flat1, outdir1+outname1, headerout
           endif
           
           headerout = 't='+time+' n_aver='+red_stri(nsum)+' darkcorrected'
           print, inam+' : saving ' + outdir+outname
           fzwrite, float(flat), outdir+outname, headerout

           if keyword_set(check) then free_lun, lun

        endif else begin
           
           print, inam+' : file exists: ' + outdir+outname + $
                  ' , skipping! (run sumflat, /overwrite to recreate)'

        endelse
     endfor                     ; istate
  endfor                        ; ic 

end  
