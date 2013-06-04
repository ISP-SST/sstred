pro red::sumflat, overwrite = overwrite, ustat = ustat, old = old, remove = remove,check=check
                                ;
  if(self.doflat eq 0B) then begin
     print, 'red::sumdark : ERROR : undefined flat_dir'
     return
  endif
  outdir = self.out_dir + '/' + 'flats/'
  file_mkdir, outdir
  outdir1 = self.out_dir + '/' + 'flats/summed/'
  file_mkdir, outdir1

                                ; 
                                ; create file list
  for ic = 0L, 1 do begin
     case ic of

        0: begin
           if(self.docamt eq 0B) then begin
              print, 'red::sumflat : Nothing to do for ', self.camt
           endif
           cam = self.camt
           doit = self.docamt
        end
        1: begin
           if(self.docamr eq 0B) then begin
              print, 'red::sumflat : Nothing to do for '+ self.camr
           endif
           cam = self.camr
           doit = self.docamr
        end
     endcase
     if(~doit) then continue
                                ;
     spawn, 'find ' + self.flat_dir + '/' + cam + '/ | grep "camX"', files
     nf = n_elements(files)
     print, 'red::sumflat : Found '+red_stri(nf)+' files in '+ self.flat_dir + '/' + cam + '/'
                                ;
     if(files[0] eq '') then begin
        print, 'red::sumflat : ERROR : '+cam+': no files found in: '+$
               self.flat_dir+' : skipping camera!'
        continue
     endif
                                ;
     files = red_sortfiles(files)
                                ;
     stat = red_getstates(files)
                                ;
     camtag = (strsplit(file_basename(files[0]), '.', /extract))[0]
                                ;
     if(~file_test(self.out_dir+'/darks/'+camtag+'.dark')) then begin
        print, 'red::sumflat : no darks found for ', camtag
        continue
     endif
     dd = f0(self.out_dir+'/darks/'+camtag+'.dark')
                                ;
                                ; Unique states
     if(~keyword_set(ustat)) then ustat = stat.state[uniq(stat.state, sort(stat.state))]
     ns = n_elements(ustat)

     ;; flag first frame after tunning the FPI
     if(keyword_set(remove)) then red_flagtunning, stat
     
     ;; Loop and sum
     for ss = 0L, ns - 1 do begin
        pos = where((stat.state eq ustat[ss]) AND (stat.star eq 0B), count)
        if(count eq 0) then continue
                                ;
        outname = camtag + '.' + ustat[ss] + '.flat'
        outname1 = camtag + '.' + ustat[ss] + '.summed.flat'
        


                                ;
        ;; If file does not exist, do sum!
        if(keyword_set(overwrite) OR ~file_test(outdir + outname)) then begin
           print, 'red::sumflat : adding flats for state -> ' + ustat[ss]
           IF(keyword_set(check)) THEN BEGIN
              openw, lun, self.out_dir + '/flats/'+camtag+'.'+ustat[ss]+'.discarded.txt', width = 500, /get_lun
           endif
           
           ;; sum files
           IF(keyword_set(check)) THEN flat = red_sumfiles(files[pos], time = time, check = check, lun = lun) ELSE flat = red_sumfiles(files[pos], time = time)
           
           ;; remove dark
           flat1 = flat
           flat-= dd
           
           ;; header for output file
           headerout = 't='+time+' n_aver='+red_stri(count)
           
           ;; output the average flat
           file_mkdir, outdir
           file_mkdir, outdir1

           print, 'red::sumflat : saving ' + outdir1+outname1
           fzwrite, float(flat1), outdir1+outname1, headerout

           headerout+= ' darkcorrected'
           print, 'red::sumflat : saving ' + outdir+outname
           fzwrite, float(flat), outdir+outname, headerout

           IF(keyword_set(check)) THEN free_lun, lun

        endif else begin
           print, 'red::sumflat : file exists: ' + outdir+outname + $
                  ' , skipping! (run sumflat, /overwrite to recreate)'
        endelse
     endfor
  endfor                        ;(ic loop)

  ;; Now WB camera
  if(self.docamwb eq 0B) then begin
     print, 'red::sumflat : Nothing to do for '+ self.camwb
     self.done.sumflat = 1B
     return
  endif

  ;; search files
  cam = self.camwb
  spawn, 'find ' + self.flat_dir + '/' + cam + '/ | grep camX', files

  if(files[0] eq '') then begin
     print, 'red::sumflat : ERROR : '+cam+': no files found in: '+$
            self.flat_dir+' : skipping camera!'
     self.done.sumflat = 1B
     return
  endif

  files = red_sortfiles(files)
  nt = n_elements(files)

  ;; Get prefilters state
  wstat = strarr(nt)
  for jj = 0L, nt - 1 do begin
     dum = strsplit(file_basename(files[jj]), '.', /extract)
     wstat[jj] = dum[4]
  endfor
  uwstat = wstat[uniq(wstat, sort(wstat))]

  camtag = (strsplit(file_basename(files[0]), '.', /extract))[0]
  outdir = self.out_dir + '/' + 'flats/'

  ;; load dark
  if(~file_test(self.out_dir+'/darks/'+camtag+'.dark')) then begin
     print, 'red::sumflat : no darks found for ', camtag
     return
  endif

  dd = f0(self.out_dir+'/darks/'+camtag+'.dark')

  for ss = 0L, n_elements(uwstat) - 1 do begin
     pos = where(wstat eq uwstat[ss], count)

     outname = camtag + '.' + uwstat[ss] + '.flat'
     outname1 = camtag + '.' + uwstat[ss] + '.summed.flat'

     if(file_test(outdir+outname) AND ~keyword_set(overwrite)) then begin
        print, 'red::sumflat : file exists: ' + outdir+outname +$
               ' , skipping! (run sumflat, /overwrite to recreate)'
        continue
     endif

     if count eq 0 then begin
        print, 'red::sumflat : no files found for WB state -> '+ uwstat[ss]
        continue
     endif

     print, 'red::sumflat : adding flats for WB state -> '+uwstat[ss]
     flat = red_sumfiles(files[pos], time = time, check=check)

     ;; remove dark
     flat1 = flat
     flat-= dd

     ;;header for output file
     headerout = 't='+time+' n_aver='+red_stri(count)
     

     
     ;; output the average flat
     file_mkdir, outdir
     
     print, 'red::sumflat : saving ' + outdir+outname
     fzwrite, float(flat), outdir+outname, headerout

     print, 'red::sumflat : saving ' + outdir1+outname1
     fzwrite, float(flat), outdir1+outname1, headerout
  endfor
  self.done.sumflat = 1B

  return
end  
