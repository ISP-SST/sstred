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
;    pref  : 
;   
;   
;   
;    files : 
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
pro red::getstats, pref = pref, files=files
  inam = 'red::getstats : '
  if(n_elements(files) eq 0) then begin
     print, inam+'working dir -> '+self.data_dir
     spawn, 'find '+' '+self.data_dir+'/'+self.camwb+'/ | grep camX | grep -v ".lcd."', files 
  endif
  if(n_elements(files) lt 1) then begin
     print, inam + 'ERROR, no files found -> '+self.data_dir
  endif 

  files = red_sortfiles(temporary(files))
  cam = (strsplit(file_basename(files[0]), '.',/extract))[0]

  ddf = self.out_dir + 'darks/'+cam+'.dark'
  if(~file_test(ddf)) then begin
     print, inam + 'error, could not find dark file -> '+ddf
  endif
  dd = f0(ddf)

  stat = red_getstates(files)
  upref = stat.pref[uniq(stat.pref, sort(stat.pref))]

  if(keyword_set(pref)) then begin
     idx = where(stat.pref eq pref, count)
     if(count lt 1) then begin
        print, inam+'error, no images found for prefilter -> '+pref
        return
     endif    
     files = files[idx]
  endif
  nf = n_elements(files)

  stats = fltarr(3,nf)
  
  ntot = 100./(nf-1.0)
  for ii = 0L, nf - 1 do begin
     a = f0(files[ii]) - dd
     me = mean(a)
     st = stdev(a)
     stats[*,ii] = [me, st, st/me]
     print, string(13B), ii*ntot, '%', format='(A,F5.1,A,$)'
  endfor
  print, ' '
  save, file='crispred.imagestats.sav', stat, stats, files

end
