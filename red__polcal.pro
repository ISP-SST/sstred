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
; :Returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;    nthreads : 
;   
;   
;   
;    nodual  : 
;   
;   
;   
;    pref  : 
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
;   2018-04-11 : MGL. Adapt to new codebase.
; 
;   2018-04-16 : MGL. Read single FITS files with extensions.
; 
;   2022-09-02 : MGL. Make a spatial mask and use it.
;
;   2024-11-02 : JdlCR. Modifications for new
;                demodulation/flat-fielding scheme.
;
;-
pro red::polcal, offset = offset, nthreads=nthreads, nodual = nodual, pref = pref

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  outdir = self.out_dir + '/polcal/'
  if n_elements(nthreads) eq 0 then nthreads = 7

  ;; Search for cubes
  f = file_search(self.out_dir + '/polcal_cubes/cam*_????_polcalcube.fits', count = Ncubes)
  if Ncubes eq 0 then begin    
    print, inam + ' : ERROR, no valid polcal files found in '+self.out_dir + '/polcal_cubes/'
    return
  endif

  ;; States
  mpref = strarr(Ncubes)
  mcam = strarr(Ncubes)
  for icube = 0, Ncubes-1 do begin
    tmp = strsplit(file_basename(f[icube]),'_',/extract)
    mcam[icube] = tmp[0]
    mpref[icube] = tmp[1]
  endfor
  upref = mpref[uniq(mpref, sort(mpref))]

  select = 1B
  if n_elements(pref) ne 0 then begin
    idx = where(upref eq pref, count)
    if(count eq 0) then begin
      print, inam + ' : ERROR, user-supplied prefilter is not valid -> ' + pref
      select = 1
    endif else select = 0
  endif

  if select then begin
    if n_elements(upref) eq 1 then begin
      pref = upref[0]
    endif else begin
      print, inam + ' : Available prefilters:'
      for ii = 0, n_elements(upref)-1 do print, ii,' -> '+upref[ii],format='(I3,A)'
      idx = 0L
      read, prompt = 'Select prefilter number: ', idx
      pref = upref[idx]
    endelse
  endif
  print, inam + ' : Selected prefilter -> '+pref

  ;; Select cameras
  idx = where(mpref eq pref, count)
  f = f[idx]
  ucam = mcam[idx]

  if(~keyword_set(nodual) AND (count eq 2)) then begin
    print, inam + ' : Found 2 cameras for selected prefilter'
    
    ;; Load data
    dir = self.out_dir + '/polcal_cubes/'
    root = '_'+pref+'_'

    rname = dir+ucam[0]+root+'polcalcube.fits'
    tname = dir+ucam[1]+root+'polcalcube.fits'
    
    r1d = mrdfits(rname,'D1D')
    t1d = mrdfits(tname,'D1D')
    rqw = mrdfits(rname,'QW')
    tqw = mrdfits(tname,'QW')
    rlp = mrdfits(rname,'LP')
    tlp = mrdfits(tname,'LP')
    tpf = mrdfits(tname,'PF')
    rpf = mrdfits(rname,'PF')
    
    ;; Get offset and delta (using QWP angle as reference)
    qlt = (red_polcal_fit(t1d, tqw, tlp, norm=4))[16:17]
    qlr = (red_polcal_fit(r1d, rqw, rlp, norm=4))[16:17]
    ql = (qlt + qlr) * 0.5
    da = ql[1]
    ql[1] = 0.
    print, inam + ' : Detected offset angle -> '+string(da)+' degrees'
    print, inam + ' :          retardance   -> '+string(ql[0])+' degrees'

    ;; Init matrix for each camera
    par_t = red_polcal_fit(t1d, tqw, tlp-da, norm=4, fix=ql)
    par_r = red_polcal_fit(r1d, rqw, rlp-da, norm=4, fix=ql)

    ;; Print efficiency values
    dum = red_eff(invert(reform(par_t[0:15],[4,4])))
    print, inam + ' : Efficiency for transmitted camera -> ' $
           + 'I =' + string(dum[0]*100., format='(F5.1)') + '%, ' $
           + 'Q =' + string(dum[1]*100., format='(F5.1)') + '%, ' $
           + 'U =' + string(dum[2]*100., format='(F5.1)') + '%, ' $
           + 'V =' + string(dum[3]*100., format='(F5.1)') + '%'
    
    dum = red_eff(invert(reform(par_r[0:15],[4,4])))
    print, inam + ' : Efficiency for reflected camera   -> ' $
           + 'I =' + string(dum[0]*100., format='(F5.1)') + '%, ' $
           + 'Q =' + string(dum[1]*100., format='(F5.1)') + '%, ' $
           + 'U =' + string(dum[2]*100., format='(F5.1)') + '%, ' $
           + 'V =' + string(dum[3]*100., format='(F5.1)') + '%'
    
    ;; Do fits and save

    ;; Make a common mask for both R and T cameras.
    rdata = red_readdata(rname)
    tdata = red_readdata(tname)
    totrdata = total(total(total(rdata,1),1),1)
    tottdata = total(total(total(tdata,1),1),1)
    mask = tottdata gt max(tottdata)/10. and totrdata gt max(totrdata)/10.

    file_mkdir, outdir

    mm = red_cpolcal_2d(temporary(rdata), rqw, rlp-da, par_r, nthreads=nthreads, mask = mask)
    oname = outdir + ucam[0]+'_'+pref+'_polcal.fits'
    print, inam + ' : saving '+oname
    dim = size(mm,/dim)

    red_writedata, oname, reform(mm, [16,dim[2],dim[3]], /overwrite), filetype='FITS', /overwrite
    
    mkhdr, ehdr, [rpf], /image
    red_fitsaddkeyword, anchor = anchor, ehdr, 'EXTNAME', 'PF', 'Polcal flatfielding'
    writefits, oname, [rpf], ehdr, /append
      

    
;    data = red_readdata(tname)
;    totdata=total(total(total(data,1),1),1)
;    mask = totdata gt max(totdata)/10. ; Mask that deselects bad pixels and pixels without light
    mm = red_cpolcal_2d(temporary(tdata), tqw, tlp-da, par_t, nthreads=nthreads, mask = mask)
    oname = outdir + ucam[1]+'_'+pref+'_polcal.fits'
    print, inam + ' : saving '+oname
    dim = size(mm,/dim)

    red_writedata, oname, reform(mm, [16,dim[2],dim[3]], /overwrite), filetype='FITS', /overwrite
    
    mkhdr, ehdr, [tpf], /image
    red_fitsaddkeyword, anchor = anchor, ehdr, 'EXTNAME', 'PF', 'Polcal flatfielding'
    writefits, oname, [tpf], ehdr, /append
      
  endif 
  
end
