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
;    nthreads  : 
;   
;   
;   
;    cam  : 
;   
;   
;   
;    pref  : 
;   
;    no_descatter : in, optional, type=boolean 
;   
;      Don't do back-scatter compensation.
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
;                find out its own name.
; 
;   2016-02-15 : MGL. Use loadbackscatter. Remove keyword descatter,
;                new keyword no_descatter.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
; 
; 
;-
pro red::prepflatcubes, flatdir = flatdir, no_descatter = no_descatter, nthreads = nthreads, cam = cam, pref = pref

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

  print, inam + ' : searching for flat images -> '+ flatdir

  for cc = 0, nc -1 do begin
     f = file_search(flatdir+'/'+cam[cc]+'.*lc0*.flat', count = count)
     print, inam + ' : ' + red_stri(count)+' wavelenghts found for '+cam[cc]
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
     for pp = 0, np - 1 do begin
        if(keyword_set(pref)) then begin
           if uprefs[pp] ne pref then begin
              print, inam + ' : skipping prefilter -> '+uprefs[pp]
              continue
           endif
        endif 
                                ; if(np gt 1) then begin
                                ;    print, inam+' : found different spectral lines:'
                                ;    for ss = 0, np - 1 do print, '   '+red_stri(ss, ni='(I3)')+' -> '+ upref[ss]
                                ;    ipref = 0L
                                ;    read, ipref, prompt = inam+' : Select prefilter ID: '
                                ;    upref = upref[ipref]
                                ; endif
                                ;
        upref = uprefs[pp]
                                ;
        mmtf = file_search(self.out_dir+'/polcal/'+cam[cc]+'.'+upref+'.polcal.f0', count = cmmt)
                                ;mmrf = file_search(self.out_dir+'/polcal/'+cam[1]+'.'+upref+'.polcal.f0', count = cmmr)

        if((cmmt ne 1)) then begin
           print, inam + ' : ERROR, '+cam[cc]+'.'+upref+' -> files not found in '+self.out_dir+'/polcal/'
           continue
        endif
        print, inam +' : loading polcal data -> '+mmtf
        immt = f0(mmtf)
                                ;immr = f0(mmrf)
        immt = red_invert_mmatrix(temporary(immt))
                                ;immr = red_invert_mmatrix(temporary(immr))

        pos = where(stat.pref eq upref, nstat)

        ;; Load backscatter data?
        if ~keyword_set(no_descatter) AND (upref eq '8542' or upref eq '7772') then begin
           self -> loadbackscatter, cam[cc], upref, bg, psf
;           print, inam + ' : loading descatter data for '+cam[cc]
;           bg =  f0(self.descatter_dir + '/' + cam[cc] + '.backgain.f0')
;           psf = f0(self.descatter_dir + '/' + cam[cc] + '.psf.f0')
        endif

        ;; Start demodulation
        for ss = 0L, nstat - 1 do begin
           if(ss eq 0) then begin
              fzread, tmp, stat.f[pos[ss]], dump
              dim = size(tmp, /dimension)
              cub = fltarr(nstat, dim[0], dim[1])
              wav = dblarr(nstat)
              tomask = bytarr(dim[0], dim[1])
           endif

           ;; Load flats and demodulate
           lc0 = file_search(flatdir+'/'+cam[cc]+'.*'+stat.state[pos[ss]]+'.lc0*.flat', count = nlc0)
           lc1 = file_search(flatdir+'/'+cam[cc]+'.*'+stat.state[pos[ss]]+'.lc1*.flat', count = nlc1)
           lc2 = file_search(flatdir+'/'+cam[cc]+'.*'+stat.state[pos[ss]]+'.lc2*.flat', count = nlc2)
           lc3 = file_search(flatdir+'/'+cam[cc]+'.*'+stat.state[pos[ss]]+'.lc3*.flat', count = nlc3)

           if((nlc0 NE 1) OR (nlc1 NE 1) OR (nlc2 NE 1) OR (nlc3 NE 1)) then begin
              print, inam+' : ERROR, there are more than 1 flat per state!'
              print, inam+' : lc0 -> '+nlc0
              print, inam+' : lc1 -> '+nlc1
              print, inam+' : lc2 -> '+nlc2
              print, inam+' : lc3 -> '+nlc3
              stop
           endif

           ;; Print info
           print, inam+' : demodulating images: '
           print, '   -> '+lc0
           print, '   -> '+lc1
           print, '   -> '+lc2
           print, '   -> '+lc3

           ;; Load data
           lc0 = f0(lc0)
           lc1 = f0(lc1)
           lc2 = f0(lc2)
           lc3 = f0(lc3)

           ;; Descatter ?
           if(~keyword_set(no_descatter) AND (upref EQ '8542' or upref eq '7772')) then begin

             lc0 = rdx_descatter(temporary(lc0), bg, psf, /verbose, nthreads = nthreads)
             lc1 = rdx_descatter(temporary(lc1), bg, psf, /verbose, nthreads = nthreads)
             lc2 = rdx_descatter(temporary(lc2), bg, psf, /verbose, nthreads = nthreads)
             lc3 = rdx_descatter(temporary(lc3), bg, psf, /verbose, nthreads = nthreads)

           endif

           ;; Demodulate flats
           tmp = reform((red_demodulate_simple(immt, lc0, lc1, lc2, lc3))[*,*,0]) 
                                ;tmp = reform((red_demodulate_simple(immr, lc0, lc1, lc2, lc3))[*,*,0])
                                ;
                                ; nans? They can appear at the borders, especially if 
                                ; there is a mask.
                                ;
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
        ord = sort(wav)
        cub = (temporary(cub))[ord, *,*]
        wav = wav[ord]

        ;; Save results (separate fits files for old fortran fitgains_ng)
        outname = cam[cc]+ '.'+ upref+'.flats_data.fits'
        outname1 = cam[cc]+'.'+ upref+'.flats_wav.fits'
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
        save, file=outdir+cam[cc]+'.'+ upref +'.flats.sav', cub, wav, namelist
     endfor

     cub = 0B
     wav = 0B
     namelist = 0B
     
  endfor

  return
end
