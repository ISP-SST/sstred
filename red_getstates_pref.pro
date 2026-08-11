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
function red_getstates_pref, files
                                ;
  nt = n_elements(files)
                                ;
  prstat = {state:strarr(nt), pref:strarr(nt),$
            wav:strarr(nt), lc:strarr(nt), nums:strarr(nt), files:files, star:bytarr(nt),$
            scan:strarr(nt), lre:lonarr(nt), hre:lonarr(nt)}
                                ;
  for ii = 0L, nt - 1 do begin
     tmp = strsplit(file_basename(files[ii]), '.', /extract)
     prstat.scan[ii] = tmp[1]
     prstat.pref[ii] = tmp[4]
     prstat.wav[ii] = tmp[5]
     prstat.nums[ii] = tmp[7]
     prstat.state[ii] = tmp[4]+'.'+tmp[5]
  endfor
                                ;
  return, temporary(prstat)
end
