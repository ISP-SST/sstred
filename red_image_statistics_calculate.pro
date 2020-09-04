; docformat = 'rst'

;+
; Calculate statistiscs for an image (or any array).
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
;
;    A struct with statistics.
; 
; :Params:
; 
;    im : in, type=array
; 
;      The data for which the statistics should be calculated.
; 
; 
; :Keywords:
; 
;   percentile_p : in, out, optional, type=fltarr, default="[.01, .02 .05, .10, .25, .50, .75, .90, .95, .98, .99]"
;   
;      The percentile P limits for which the percentile values are to
;      be calculated.
; 
; 
; :History:
; 
;   2018-10-26: MGL. First version.
; 
;   2020-06-26: MGL. Take care of NaNs.
; 
;-
function red_image_statistics_calculate, im, percentile_p = percentile_p 

  if n_elements(percentile_p) eq 0 then $
     percentile_p = [.01, .02, .05, .10, .25, .50, .75, .90, .95, .98, .99]
  
  momnt = moment(im, /double, /nan) ; [mean, variance, skewness, kurtosis]

  indx = where(finite(im), Nwhere)
  if Nwhere gt 0 then perc = double(cgpercentiles(im, percentiles = percentile_p)) $
  else perc = percentile_p + !Values.F_NaN

  output = create_struct('NPIXELS',  long64(n_elements(im)) $
                         , 'DATAMIN',  min(im) $
                         , 'DATAMAX',  max(im) $
                         , 'DATAMEAN', momnt[0] $       
                         , 'DATARMS',  sqrt(momnt[1]) $
                         , 'DATASKEW', momnt[2] $
                         , 'DATAKURT', momnt[3] )

  for i = 0, n_elements(percentile_p)-1 do begin
    case percentile_p[i] of
      0.5  : tag = 'DATAMEDN'
      else : tag = 'DATAP' + string(percentile_p[i]*100, format = '(i02)')
    endcase
    output =  create_struct(tag, perc[i], output)
  endfor                        ; i

  if n_elements(hist) eq 0 then begin

    
  endif
  
  return, output
  
end
