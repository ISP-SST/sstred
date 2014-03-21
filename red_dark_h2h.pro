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
;    fname : 
;   
;   
;   
;    head : 
;   
;   
;   
;    count : 
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
;   2014-03-21 : MGL. Use red_extractstates to get file and scan
;                numbers from file name (scan number actually not
;                needed). Substitute real exposure time for hard-coded
;                17.600 ms.
; 
;-
function red_dark_h2h, fname, head, count
                                ;
;  tmp = strsplit(file_basename(fname), '_.', /extract)
;  scan = tmp[1]
;  num = tmp[5]

  red_extractstates,fname,scan=scan,nums=num
  
  tmp = strsplit(head, ' =', /extract)
  nx = tmp[9]                   ; x size
  ny = tmp[11]                  ; y size
  
  time = red_time2double((red_time2double(tmp[18]) + red_time2double(tmp[21]))*0.5d0,$
                     /dir)
  exp = (red_time2double(tmp[21]) - red_time2double(tmp[18])) / 1e-3
  res = '#.' + string(num,format='(I7.7)') $
        + '      ' + num[0] + '  ' + time + '  ' $
        + nx + '  ' + ny + '    ' $
        + string(exp,format='(f6.3)') + '    ' $
        + string(count,format='(I4)') + ' ... SUM ... DD'

  return, res

end
