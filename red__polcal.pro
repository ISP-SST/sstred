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
;    cams  : 
;   
;   
;   
;    offset  : 
;   
;   
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
; 
;-

function red_eff, dmm
  return, 1./sqrt(4*total(dmm^2, 1))
end


pro red::polcal, cams = cams, offset = offset, nthreads=nthreads, nodual = nodual, pref = pref

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])+': '

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  outdir = self.out_dir + '/polcal/'
  if(keyword_set(nthreads)) then threads = nthreads else threads = 7

  ;; Search for cubes
  root = self.out_dir + '/polcal_cubes/camX*.????.3d.f*'
  f = file_search(root, count = ct)
  if(ct eq 0) then begin    
     print, inam + ' : ERROR, no valid polcal files found in '+self.out_dir + '/polcal_cubes/'
     return
  endif

  ;; States
  mpref = strarr(ct)
  mcam = strarr(ct)
  for ii = 0, ct-1 do begin
     tmp = strsplit(file_basename(f[ii]),'.',/extract)
     mcam[ii] = tmp[0]
     mpref[ii] = tmp[1]
  endfor
  upref = mpref[uniq(mpref, sort(mpref))]

  select = 1B
  if(keyword_set(pref)) then begin
     idx = where(upref eq pref, count)
     if(count eq 0) then begin
        print, inam + ' : ERROR, user-supplied prefilter is not-valid -> ' + pref
        select = 1
     endif else select = 0
  endif

  If(select) then begin
     if(n_elements(upref) eq 1) then begin
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
     print, inam + ' : found 2 cameras for selected prefilter'
     
     ;; Load 1D data
     dir = self.out_dir + '/polcal_cubes/'
     root = '.'+pref+'.'
     
     test = file_search(dir+'*.fits', count = ct)
     if(ct gt 0) then begin
        r1d = readfits(dir+ucam[0]+root+'1d.fits')
        t1d = readfits(dir+ucam[1]+root+'1d.fits')
        rqw = readfits(dir+ucam[0]+root+'qw.fits')
        tqw = readfits(dir+ucam[1]+root+'qw.fits')
        rlp = readfits(dir+ucam[0]+root+'lp.fits')
        tlp = readfits(dir+ucam[1]+root+'lp.fits')
        ex = 'fits'
     endif else begin
        r1d = f0(dir+ucam[0]+root+'1d.f0')
        t1d = f0(dir+ucam[1]+root+'1d.f0')
        rqw = f0(dir+ucam[0]+root+'qw.f0')
        tqw = f0(dir+ucam[1]+root+'qw.f0')
        rlp = f0(dir+ucam[0]+root+'lp.f0')
        tlp = f0(dir+ucam[1]+root+'lp.f0')
        ex = 'f0'
     endelse


     ;; Get offset and delta (using QWP angle as reference)
     qlt = (red_polcal_fit(t1d, tqw, tlp, norm=4))[16:17]
     qlr = (red_polcal_fit(r1d, rqw, rlp, norm=4))[16:17]
     ql = (qlt + qlr) * 0.5
     da = ql[1]
     ql[1] = 0.
     print, inam + ' : detected offset angle -> '+string(da)+' degrees'
     print, inam + ' :          retardance   -> '+string(ql[0])+' degrees'


     
     ;; Init matrix for each camera
     par_t = red_polcal_fit(t1d, tqw, tlp-da, norm=4, fix=ql)
     par_r = red_polcal_fit(r1d, rqw, rlp-da, norm=4, fix=ql)

     ;; Print efficiency values
     dum = red_eff(invert(reform(par_t[0:15],[4,4])))
     print, inam+'Efficiency for transmitted camera -> I='+string(dum[0]*100., format='(F6.1)')+'%, Q='+string(dum[1]*100., format='(F6.1)')+'%, U='+string(dum[2]*100., format='(F6.1)')+'%, V='+string(dum[3]*100., format='(F6.1)')+'%'
     
     dum = red_eff(invert(reform(par_r[0:15],[4,4])))
     print, inam+'Efficiency for reflected camera -> I='+string(dum[0]*100., format='(F6.1)')+'%, Q='+string(dum[1]*100., format='(F6.1)')+'%, U='+string(dum[2]*100., format='(F6.1)')+'%, V='+string(dum[3]*100., format='(F6.1)')+'%'
     
     ;; Do fits and save
     if(ex eq 'fits') then begin
        data = readfits(f[0])
     endif else begin
        data = f0(f[0])
     endelse

     mm = red_cpolcal_2d(temporary(data), rqw, rlp-da, par_r, nthreads=threads)
     file_mkdir, outdir
     oname = outdir + ucam[0]+'.'+pref+'.polcal.f0'
     print, inam + ' : saving '+oname
     dim = size(mm,/dim)
     fzwrite, reform(temporary(mm), [16,dim[2],dim[3]]), oname ,' '

     if(ex eq 'fits') then begin
        data = readfits(f[1])
     endif else begin
        data = f0(f[1])
     endelse

     mm = red_cpolcal_2d(temporary(data), tqw, tlp-da, par_t, nthreads=threads)
     oname = outdir + ucam[1]+'.'+pref+'.polcal.f0'
     print, inam + ' : saving '+oname
     dim = size(mm,/dim)
     fzwrite, reform(temporary(mm), [16,dim[2],dim[3]]), oname ,' '
  endif 
  
  return
end
