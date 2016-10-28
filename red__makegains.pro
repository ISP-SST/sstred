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
; 
; 
;-
pro red::makegains, no_descatter = no_descatter, nthreads = nthreads, cam = cam, pref = pref, min = min, max = max, bad=bad, preserve=preserve, smoothsize = smoothsize, cavityfree=cavityfree
                                ;
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo
  
  if(keyword_set(cavityfree)) then tosearch = self.out_dir+'/flats/*cavityfree.flat.fits' $
  else  tosearch = self.out_dir+'/flats/*.flat'
  
  files = file_search(tosearch, count = ct)
                                ;
  if(ct eq 0) then begin
     print, inam+' : No flats found in: ' + tosearch
  endif
                                ;
  firsttime = 1B
  for ii = 0L, ct -1 do begin
     tmp = strsplit(file_basename(files[ii]), '._', /extract)
     if(keyword_set(pref)) then begin
        if(tmp[1] ne pref) then begin
           print, inam+' : skipping prefilter -> '+tmp[1]
           continue
        endif
     endif
     fzread, flat, files[ii], head
                                ;
                                ; Only one camera?
                                ;
     if n_elements(cam) ne 0 then if tmp[0] NE cam then continue
                             ;
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
     
     if(keyword_set(cavityfree)) then namout = file_basename(files[ii], '_cavityfree.flat.fits')+'.gain' $
     else namout = file_basename(files[ii], '.flat')+'.gain'
     
     outdir = self.out_dir+'/gaintables/'
     h = head
                                ;
                                ; Output gaintable
     file_mkdir, outdir
     print, inam+' : saving '+outdir+namout
     fzwrite, float(gain), outdir+namout, h
  endfor
                                ;
  return
end
