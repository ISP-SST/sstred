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
;    step  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
pro red::sumprefilter, step = step
                                ;
  if(~keyword_set(step)) then step = 1L
                                ;
  inam = 'red::sumprefilter : '
                                ;
                                ; Prefilter scan taken usually with T_cam?
                                ;
  cam = self.camt
  
                                ;
                                ; find files
  spawn, 'find ' + self.prefilter_dir + '/' + cam + '/ | grep im', files
  nt = n_elements(files)
  if(files[0] eq '') then begin
     print, inam+'no files found in '+self.polcal_dir+'/'+cam+', skipping camera!'
     return
  endif
                                ;
                                ; sort files based on image number
  files = red_sortfiles(files)
                                ;
                                ; load dark
  camtag =  (strsplit(file_basename(files[0]), '.',/extract))[0]
  df = self.out_dir+'/darks/'+camtag+'.dark'
  if(file_test(df)) then dd = f0(df) else begin
     print, inam+'dark not found in '+ self.out_dir+'/darks/'
     return
  endelse
                                ;
                                ; get states
  pfstat = red_getstates_pref(files)
  red_flagtunning, pfstat
                                ;
                                ; get unique states
  ustate = pfstat.state[uniq(pfstat.state, sort(pfstat.state))]
  ns = n_elements(ustate)
                                ;
  print, inam+' found '+red_stri(ns)+' prefilter states'
                                ;
                                ; Loop states
  outdir = self.out_dir+'/prefilter_scan/'
  file_mkdir, outdir
  prf = fltarr(ns, 1024, 1024)
  wav = lonarr(ns)
                                ;
  for ss = 0L, ns - 1 do begin
                                ;
     pos = where((pfstat.state eq ustate[ss]) AND (pfstat.star eq 0), count)
     if(count eq 0) then continue
     print, inam+'adding images for state-> '+ustate[ss]
                                ;
                                ; sum files and remove dark
     prf[ss,*,*] = red_sumfiles(files[pos], time = time) - dd
                                ;
                                ; Get wavelength from header
     hdr = strsplit(fzhead(files[pos[0]]),/extract)
     hdr = hdr[n_elements(hdr)-1]
     pl = strpos(hdr, 'LRE=u')
     if pl ne -1 then begin
        ph = strpos(hdr, 'HRE=u')
        hex = strmid(hdr, pl+6, 12)+strmid(hdr, ph+6, 12)
        legs = lonarr(6)
        reads, hex, legs, form='(6Z4)'
        lre=legs[0]
        hre=legs[3]
        pfstat.lre[pos] = lre
        pfstat.hre[pos] = hre
        wav[ss] = hre 
     endif
                                ;
  endfor
                                ;
  dw = abs(wav[1] - wav[0])
                                ;
  wav-=wav[0]
  wav*=step
                                ;
  pos = sort(wav)
  wav = wav[pos]
  prf = (temporary(prf))[pos,*,*]
                                ;
                                ; output results
  outname = camtag + '.' + pfstat.pref[0]+'_prefilter.f0'
  outname1 = camtag + '.' + pfstat.pref[0]+'_wavelength.f0'
                                ;
  fzwrite, prf, outdir+outname,' '
  fzwrite, wav, outdir+outname1, 'step='+red_stri(step)+' dHRE='+red_stri(dw)
                                ;
  return
end
