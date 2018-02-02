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
;   2017-12-20 : MGL. Store gains with metadata headers.
; 
; 
;-
pro chromis::makegains, bad=bad $
                        , cam = cam $
                        , max = max $
                        , min = min $
                        , nthreads = nthreads $
                        , pref = pref $
                        , smoothsize = smoothsize 

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
  
  red_make_prpara, prpara, bad
  red_make_prpara, prpara, cam 
  red_make_prpara, prpara, max 
  red_make_prpara, prpara, min 
  red_make_prpara, prpara, pref 
  red_make_prpara, prpara, smoothsize
  
  tosearch = self.out_dir+'/flats/*.flat.fits'
  
  files = file_search(tosearch, count = Nfiles)

  if Nfiles eq 0 then begin
    print, inam+' : No flats found in: ' + tosearch
  endif

  for ifile = 0L, Nfiles -1 do begin
    
    tmp = strsplit(file_basename(files[ifile]), '._', /extract)
    if(keyword_set(pref)) then begin
      if(tmp[1] ne pref) then begin
        print, inam+' : skipping prefilter -> '+tmp[1]
        continue
      endif
    endif

    flat = red_readdata(files[ifile], head = hdr)

    ;; Only one camera?
    if n_elements(cam) ne 0 then if tmp[0] NE cam then continue

    gain = self->flat2gain(flat, ma=max, mi=min, bad=bad, /preserve, smoothsize=smoothsize)
    
    namout = file_basename(files[ifile], '.flat.fits')+'.gain'
    outdir = self.out_dir+'/gaintables/'

    ;; Edit the header
    red_fitsaddkeyword, hdr, 'FILENAME', outdir+namout
    self -> headerinfo_addstep, hdr, prstep = 'Gain making' $
                                , prproc = inam, prpara = prpara
    
    ;; Output gaintable
    file_mkdir, outdir
    print, inam+' : saving '+outdir+namout
    ;;fzwrite, float(gain), outdir+namout, ' '
    overwrite = 1
    red_writedata, outdir+namout, float(gain), header = hdr, filetype='ANA', overwrite = overwrite

  endfor                        ; ifile
                                
end
