; docformat = 'rst'

;+
; Plot data with vs time with a proper axis.
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
; 
; :Params:
; 
;    t : in, type=array
; 
;       The time coordinate in seconds after midnight.
; 
;    y : in, type=array
; 
;       The dependent coordinate. 
; 
; :Keywords:
; 
;   overplot : in, optional, type=boolean
;
;       Set this keyword if you wish to overplot data on an already
;       exisiting set of axes. It is like calling the IDL OPLOT
;       command.

; 
;   trange : in, optional
;   
;       As xrange.
; 
;   xrange : in, optional
;   
;       See the regular plot command.
; 
; 
; :History:
; 
;   2020-01-14 : MGL. First version.
; 
;-
pro red_timeplot, t, y, overplot = overplot, trange = trange, xrange = xrange, _extra = extra

  if keyword_set(overplot) then begin
    
    cgplot, t, y, _strict_extra = extra, /overplot

  endif else begin
    
    case 1 of
      n_elements(trange) ne 0 : L = red_gen_timeaxis(trange)
      n_elements(xrange) ne 0 : L = red_gen_timeaxis(xrange)
      else : L = red_gen_timeaxis(t)
    endcase
    
    cgplot, t, y, _strict_extra = extra $
            , XTICKV=L.tickv, XTICKS=L.ticks, XMIN=L.minor, XTICKNAM=L.name

  endelse
  
end
