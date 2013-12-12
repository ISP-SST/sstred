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
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
pro red::get_scansquality
  inam = 'red::get_scansquality : '

  ;; files
   
  spawn, 'find ' + (self.data_dir) + '/' + self.camwb + '/ | grep im.ex', wfiles

  ;; sort files by image number
  wfiles = red_sortfiles(temporary(wfiles))
  
  ;; get states
  stat = red_getstates(wfiles)
  
  ;; flag first frame after tuning
  print, 'red::link_data : Flagging first frame after tuning'
  red_flagtuning, stat
   
  camtag = (strsplit(file_basename(wfiles[0]), '.', /extract))[0]
   
  ;; read dark
   
  dd = self.out_dir + 'darks/' + camtag + '.dark'
  if ~file_test(dd) then begin
     print, inam + 'ERROR, file not found -> '+dd
     return
  endif
  dd = f0(dd)
   
  uscan = stat.scan[uniq(stat.scan, sort(stat.scan))]
  ns = n_elements(uscan)
  print, inam + 'Number of scans: '+red_stri(ns)
   
  ;; vars
   
  quality = fltarr(ns)
  scan = strarr(ns)
  nim = lonarr(ns)
   
  ntot = 100. / (ns - 1.0)
  bb = string(13b)
   
  print, bb,inam, 0. * ntot, '%', format = '(A,A,F5.1,A,$)'
  for ii = 0L, ns -1 do begin
     pos = where(stat.scan eq uscan[ii], count)
     if count eq 0 then continue
     aver = 0.0
     for jj = 0L, count - 1 do begin
        tmp = f0(wfiles[pos[jj]]) - dd
        aver+= stdev(tmp) / mean(tmp)
     endfor
     nim[ii] = count
   
     quality[ii] = aver / count
     dum = max(quality[0:ii], p)
   
     print, bb, inam, ii * ntot, '% -> ',red_stri(quality[p])+' -> scan='+red_stri(uscan[p]), format = '(A,A,F5.1,A,A,$)'
  endfor
  print, ' '
   
  pos = reverse(sort(quality))
  openw, lun, self.out_dir+'/scans_quality.txt', /get_lun, width = 200
   
  for ii = 0L, ns -1 do begin
     jj = pos[ii]
     print, 'Scan = '+uscan[pos[ii]] + ' -> qual = '+red_stri(quality[pos[ii]]) + ' -> Nim = '+red_stri(nim[pos[ii]])
     printf, lun, uscan[jj], quality[jj], nim[jj]
  endfor
   
  return
end
