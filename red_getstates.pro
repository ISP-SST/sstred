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
;    files : 
;   
;   
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_getstates, files
  nt = n_elements(files)
  fullstate = strarr(nt)
  nums = strarr(nt)
  wav = strarr(nt)
  lc = strarr(nt)
  pref = strarr(nt)
  scan = strarr(nt)
  dwav = dblarr(nt)
  rscan = strarr(nt)
  hscan =strarr(nt)
                                ;
  for ii = 0L, nt -1 do begin
     tmp = strsplit(file_basename(files[ii]), '._', /extract)
     n = n_elements(tmp)
                                ;
     nums[ii] = tmp[n -1]
     wav[ii] = tmp[5]+'_'+tmp[6]
     pref[ii] = tmp[4]
     lc[ii] = tmp[7]
     scan[ii] = tmp[1]
     rscan[ii] = red_decode_scan(scan[ii], hscan=hs)
     hscan[ii] = hs
     
     dum = strsplit(wav[ii],'_',/extract)
     dwav[ii] = double(dum[0])+double(dum[1])*1.0d-3
                                ;
     fullstate[ii] = pref[ii]+'.'+wav[ii] + '.' + lc[ii]
                                ;
  endfor
                                ;
  stat = {files:files, nums:nums, wav:wav, lc:lc, pref:pref, $
          state:fullstate, star:bytarr(nt), scan:scan, rscan:rscan,$
          hscan:hscan, dwav:dwav}
                                ;
  return, stat
end
