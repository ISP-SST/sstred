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
;    flatdir  : 
;   
;   
;   
;    no_descatter  : 
;   
;   
;   
;    nthreads  : 
;   
;   
;   
;    cam  : 
;   
;   
;   
;    verbose  : 
;   
;   
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
; 
;   2016-02-15 : MGL. Use loadbackscatter. Remove keyword descatter,
;                new keyword no_descatter.
;
;   2016-03-22 : JLF. Fixed a bug in which red::prepflatcubes_lc4 would
;		 overwrite the results of red::prepflatcubes if there are
;		 both lc4 and lc0-3 datasets at the same prefilter. If
;		 lc0-3 datasets are detected it will output files with
;		 .lc4 added to the state name.
; 
;   2016-05-10 : THI. New keyword pref.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;-
pro red::prepflatcubes_lc4, flatdir = flatdir $
                            , no_descatter = no_descatter $
                            , nthreads = nthreads $
                            , cam = cam $
                            , pref = pref $
                            , verbose = verbose

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  ;; Check keywords
  if(~keyword_set(flatdir)) then flatdir = self.out_dir+'/flats/'
  if(~keyword_set(nthreads)) then nthreads = 2L

  ;; Check flats and get states
  self -> getdetectors, dir = self.data_dir

  if(~keyword_set(cam)) then cam = [self.camttag, self.camrtag]

  nc = n_elements(cam)

  outdir = self.out_dir + '/flats/spectral_flats/'
  file_mkdir, outdir
  for cc = 0, nc -1 do begin
     f = file_search(flatdir+'/'+cam[cc]+'.*lc4*.flat', count = count)
     print, inam + ' : '+red_stri(count)+' wavelenghts found for '+cam[cc]
     if count eq 0 then begin
        print, inam+' : WARNING, skipping camera '+cam[cc]+' -> no files found'
        continue
     endif 

     stat = {state:strarr(count), wav:dblarr(count), f:f, pref:strarr(count), $
             wavs:strarr(count)}

     for ff = 0L, count -1 do begin
        tmp = strsplit(file_basename(f[ff]), '.',/extract)
        stat.state[ff] = tmp[1]+'.'+tmp[2]
        stat.wav[ff] = double((strsplit(tmp[2], '_',/extract))[0]) + $
                       double((strsplit(tmp[2], '_',/extract))[1]) * 0.001d0
        stat.pref[ff] = tmp[1]
        stat.wavs[ff] = tmp[2]
     endfor

     uprefs = stat.pref[uniq(stat.pref, sort(stat.pref))]
     np = n_elements(uprefs)
     uwav = stat.wavs[uniq(stat.wavs,sort(stat.wavs))]

     ;; Loop prefilters
     for pp = 0, np -1 do begin
        if(keyword_set(pref)) then begin
           if uprefs[pp] ne pref then begin
              print, inam + ' : skipping prefilter -> '+uprefs[pp]
              continue
           endif
        endif 
        upref = uprefs[pp]
        pos = where((stat.pref eq upref), nstat)
        
        if ~keyword_set(no_descatter) AND (upref eq '8542' or upref eq '7772') then begin
           self -> loadbackscatter, cam[cc], upref, bg, psf
;           print, inam + ' : loading descatter data for '+cam[cc]
;           bg =  f0(self.descatter_dir + '/' + cam[cc] + '.backgain.f0')
;           psf = f0(self.descatter_dir + '/' + cam[cc] + '.psf.f0')
        endif
        
        for ss = 0L, nstat - 1 do begin
           
           if(ss eq 0) then begin
              fzread, tmp, stat.f[pos[ss]], dump
              dim = size(tmp, /dimension)
              cub = fltarr(nstat, dim[0], dim[1])
              wav = dblarr(nstat)
              tomask = bytarr(dim[0], dim[1])
           endif

           ;; Load flats and (descatter)
           lc4 = file_search(flatdir+'/'+cam[cc]+'.*'+stat.state[pos[ss]]+'.lc4*.flat', count = nlc4)
           if(keyword_set(verbose)) then print, inam + ' : reading -> '+file_basename(lc4)
           tmp = f0(lc4)
           if ~keyword_set(no_descatter) AND (upref eq '8542' or upref eq '7772') then $
              tmp = red_cdescatter(temporary(tmp), bg, psf, /verbose, nthreads = nthreads)

           idx = where(~finite(tmp) OR (tmp LT -0.001), nnan)

           if nnan gt 0 then begin
              tmp[idx] = 0.0
              tomask[idx] = 1B
           endif

           ;; Result
           cub[ss, *, *] = tmp
           wav[ss] = stat.wav[pos[ss]] - double(stat.pref[pos[ss]])

        endfor
        for jj = 0L, dim[1]-1 do for ii = 0L, dim[0]-1 do begin
           if(tomask[ii,jj] eq 1B) then cub[*,ii,jj] = 0.0
        endfor

        ;; Sort states (so far they are sorted as strings -> incorrect order)
        print, inam + ' : sorting wavelengths ... ', FORMAT = '(A,$)'
        ord = sort(wav)
        cub = (temporary(cub))[ord, *,*]
        wav = wav[ord]
        print, 'done'

        ;; are there polarized flat cubes at this prefilter?
        ;; red::prepflatcubes must always be run first!
        ofile = cam[cc]+'.'+ upref +'.flats.sav'
        if file_test(outdir+ofile) then filestem='.lc4' else filestem=''
        
        ;; Save results (separate fits files for old fortran fitgains_ng)
        outname = cam[cc]+ '.'+ upref+filestem+'.flats_data.fits'
        outname1 = cam[cc]+'.'+ upref+filestem+'.flats_wav.fits'
        
        writefits,  outdir + outname, cub
        writefits,  outdir + outname1, wav
        namelist = strarr(nstat)

        ;; Print file names in order into a text file (just in case)
        openw, lun, outdir+cam[cc]+'_filenames.txt', /get_lun
        for ii = 0L, nstat - 1 do begin
           namelist[ii] = cam[cc]+'.'+stat.pref[pos[ord[ii]]]+'.'+ $
                          stat.wavs[pos[ord[ii]]] + '.unpol.flat'
           printf, lun, namelist[ii]
        endfor
        free_lun, lun
        
        ;; Save as structure for new routines (red::fitgain_ng)
        ofile = cam[cc]+'.'+ upref +filestem+'.flats.sav'
        print, inam+' : saving -> '+ofile

        save, file=outdir+ofile, cub, wav, namelist
        
     endfor
  endfor
end
