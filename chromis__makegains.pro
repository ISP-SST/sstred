; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CHROMIS pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Keywords:
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
;    pref  : 
;   
;   
;   
;    min  : 
;   
;   
;   
;    max  : 
;   
;   
;   
;    bad : 
;   
;   
;   
;    preserve : 
;   
;   
;   
;    smoothsize  : 
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
;   2017-04-06 : MGL. Make version for CHROMIS without backscatter or
;                zeroing of detector tap borders.
; 
; 
;-
pro chromis::makegains, nthreads = nthreads $
                        , cam = cam $
                        , pref = pref $
                        , min = min $
                        , max = max $
                        , bad=bad $
                        , smoothsize = smoothsize

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  tosearch = self.out_dir+'/flats/*.flat.fits'
  
  files = file_search(tosearch, count = Nfiles)

  if(Nfiles eq 0) then begin
    print, inam+' : No flats found in: ' + tosearch
  endif

  firsttime = 1B
  for ifile = 0L, Nfiles -1 do begin

    tmp = strsplit(file_basename(files[ifile]), '._', /extract)
    if(keyword_set(pref)) then begin
      if(tmp[1] ne pref) then begin
        print, inam+' : skipping prefilter -> '+tmp[1]
        continue
      endif
    endif

    flat = red_readdata(files[ifile])
    
    ;; Only one camera?
    if n_elements(cam) ne 0 then if tmp[0] NE cam then continue

    gain = self->flat2gain(flat, ma=max, mi=min, bad=bad, /preserve, smoothsize=smoothsize)
    
    namout = file_basename(files[ifile], '.flat.fits')+'.gain'
    
    outdir = self.out_dir+'/gaintables/'

    ;; Output gaintable
    file_mkdir, outdir
    print, inam+' : saving '+outdir+namout
    fzwrite, float(gain), outdir+namout, ' '

  endfor                        ;  ifile
  
end