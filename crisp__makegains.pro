; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
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
;   2017-12-20 : MGL. Store gains with metadata headers.
; 
;   2018-02-02 : MGL. Rename red::makegains --> crisp::makegains.
;                Adapt to new codebase.
; 
;-
pro crisp::makegains, bad=bad $
                      , cam = cam $
                      ;;, cavityfree=cavityfree
                      , max = max $
                      , min = min $
                      , nthreads = nthreads $
                      , no_descatter = no_descatter $
                      , pref = pref $
                      , preserve=preserve $
                      , smoothsize = smoothsize 
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  red_make_prpara, prpara, bad 
  red_make_prpara, prpara, cam 
  red_make_prpara, prpara, max 
  red_make_prpara, prpara, min 
  red_make_prpara, prpara, pref 
  red_make_prpara, prpara, preserve  
  red_make_prpara, prpara, smoothsize
  red_make_prpara, prpara, no_descatter 
  
  ;;if(keyword_set(cavityfree)) then tosearch = self.out_dir+'/flats/*cavityfree.flat.fits' $
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

    if ~keyword_set(no_descatter) then begin
      if((tmp[1] eq '8542' OR tmp[1] eq '7772') AND self.dodescatter) then begin
        self -> loadbackscatter, tmp[0], tmp[1], bg, psf
;           psff = self.descatter_dir+'/'+tmp[0]+'.psf.f0'
;           bgf = self.descatter_dir+'/'+tmp[0]+'.backgain.f0'
;           if(file_test(psff) AND file_test(bgf)) then begin
;              psf = f0(psff)
;              bg = f0(bgf)
        flat = red_cdescatter(flat, bg, psf, nthreads = nthreads, verbose = 1)
;           endif
      endif
    endif

    gain = self->flat2gain(flat, ma=max, mi=min, bad=bad, preserve=preserve, smoothsize=smoothsize)
    
    namout = file_basename(files[ifile], '.flat.fits')+'.gain'
    outdir = self.out_dir+'/gaintables/'

    ;; Edit the header
    red_fitsaddkeyword, hdr, 'FILENAME', outdir+namout
    self -> headerinfo_addstep, hdr, prstep = 'Gain making' $
                                , prproc = inam, prpara = prpara
    
    ;; Output gaintable
    file_mkdir, outdir
    print, inam+' : saving '+outdir+namout
    ;;fzwrite, float(gain), outdir+namout, h
    overwrite = 1
    red_writedata, outdir+namout, float(gain), header = hdr, filetype='ANA', overwrite = overwrite

  endfor                        ; ifile
                                
end
