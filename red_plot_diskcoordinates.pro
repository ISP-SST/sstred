; docformat = 'rst'

;+
; Plot pointing on the disk.
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
; :Params:
; 
;   diskpos : in, type="fltarr(2,Npoints)"
; 
;      The disk position coordinates X=diskpos[0,*], Y=diskpos[1,*].
;      Assumed to be normalized to unit solar radius, unless Rsun is
;      given as a keyword.
; 
; 
; :Keywords:
; 
;   rplot : in, optional, type=float, default=1.0
; 
;      Sets the plotting range in X and Y to [-rplot,rplot].
; 
;   rsun : in, optional, type=float, default=1.0
; 
;      Normalize diskpos with this number when plotting.
;
;   color : in, optional, type=string, default='black'
; 
;      Plot the positions with this color.
;
;   over :  in, optional, type=boolean
; 
;      Set this to overplot an existing plot.
;
;   indx :  in, optional, type="array(Npoints)"
;   
;      Plot only diskpos[*,indx].
; 
; 
; :History:
; 
;    2017-04-24 : MGL. First version.
; 
;    2017-10-10 : MGL. New keywords title, psym, and symsize. Allow
;                 color to be an array.
; 
;-
pro red_plot_diskcoordinates, diskpos $
                              , rplot = rplot $
                              , rsun = rsun $
                              , color = color $
                              , over = over $
                              , psym = psym $
                              , symsize = symsize $
                              , indx = indx $
                              , title = title $
                              , unit = unit

  if n_elements(unit) eq 0 then unit = 'R$\sun$'
  if n_elements(psym) eq 0 then psym = 3
  if n_elements(symsize) eq 0 then symsize = 1
  if n_elements(rplot) eq 0 then rplot = 1.
  if n_elements(rsun) eq 0 then Rsun = 1.
  if n_elements(color) eq 0 then color = 'black'
  if n_elements(title) eq 0 then title = ''

  
  dims = size(diskpos, /dim)
  if n_elements(indx) eq 0 then indx = indgen(dims[1])

  if ~keyword_set(over) then begin
    
    range = [-1, 1]*rplot
    
    ;; Draw a circle
    x = findgen(201)/100 - 1
    y = sqrt(1-x^2)
    cgplot, x, y, xrange=range, yrange=range, /aspect $
            , xtitle = 'X / '+unit,  ytitle = 'Y / '+unit $
            , title = title
    cgplot, /over, x, -y
    cgplot, /over, !x.crange, [0,0]
    cgplot, /over, [0,0], !y.crange

  endif
  
  ;; Plot the data points
  if n_elements(color) eq 1 then begin
    cgplot, /over, diskpos[0,indx]/Rsun, diskpos[1,indx]/Rsun, psym=psym, symcolor=color, symsize = symsize
  endif else begin
    cgplots, diskpos[0,indx]/Rsun, diskpos[1,indx]/Rsun, psym=psym, symcolor=color, symsize = symsize
  endelse
  
end
