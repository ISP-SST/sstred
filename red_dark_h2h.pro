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
; 
;-
function red_dark_h2h, fname, head, count
                                ;
  tmp = strsplit(file_basename(fname), '_.', /extract)
                                ;
  scan = tmp[1]
  num = tmp[5]
                                ;
  tmp = strsplit(head, ' =', /extract)
  nx = tmp[9]
  ny = tmp[11]
                                ;
  
  time = red_time2double((red_time2double(tmp[18]) + red_time2double(tmp[21]))*0.5d0,$
                     /dir)
  res = '#.'+string(num,format='(I7.7)') + $
        '      '+num+'  '+time+'  '+nx+'  '+ny+'    17.600    '+$
        string(count,format='(I4)')+' ... SUM ... DD'
                                ;
  return, res
end
