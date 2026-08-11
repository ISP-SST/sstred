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

  if n_elements(trange) gt 0 then begin
    if n_elements(xrange) gt 0 then begin
      print, 'red_timeplot : use only one of the trange and xrange keywords.'
      return
    endif
    xrange = trange
  endif
  
  if keyword_set(overplot) then begin
    
    cgplot, t, y, _strict_extra = extra, /overplot

  endif else begin
    
    if n_elements(xrange) ne 0 then begin
      L = red_gen_timeaxis(xrange)
    endif else begin
      L = red_gen_timeaxis(t)
    endelse
    
    cgplot, t, y, xrange = xrange,_strict_extra = extra $
            , XTICKV=L.tickv, XTICKS=L.ticks, XMIN=L.minor, XTICKNAM=L.name

  endelse
  
end
