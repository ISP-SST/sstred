;+
; Image displaying routine. 
;
;
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Jaime de la Cruz Rodriguez 2008
; 
; :Params:
;
;     vari : in, type="array or string"
;
;        A 2D image (possibly with color, making the array 3D) or the
;        name of a file with such content.
; 
; :Keywords:
; 
;     wnum : in, out, optional, type=integer
;
;        integer window index
;
;     nowin : in, optional, type=boolean
;
;         Set this to not recreate the window
;
;     offs : in, type=integer, default=57
;
;        Integer removes pixel from the Y-dimension of
;        the screen so compensates for the menubar and so on.
;
;     opt : in, optional, type=boolean 
;
;        Do a histogram optimization on the intensity limits.
;
;
;     reuse : in, optional, type=boolean
;
;        Re-use existing window, if there is one with the appropriate
;        wnum. 
;
;     scroll : in, optional, type=boolean
;
;        Make window with scrollbars for large images instead of
;        shrinking.
;
;     title : in, optional, type=string
;
;        If a window is created, it will have this title.
;
;
; :History:
;
;    2013-07-24 : Renamed red_show for inclusion in crispred pipeline.
;
;    2013-09-20 : MGL. Added title keyword.
;
;    2019-10-17 : MGL. New keyword reuse.
;
;    2020-04-07 : MGL. New keyword scroll.
;
;    20204-07-06 : MGL. Allow file name as input parameter.
;
;-
pro red_show, vari $
              , noscale = noscale $
              , nowin = nowin $
              , offs = offs $
              , opt = opt $
              , reuse = reuse $
              , scroll = scroll $
              , title = title $
              , wnum = wnum 

  ;; Initializes some variables
  if ~keyword_set(wnum) then wnum=0
  if ~keyword_set(offs) then offs=24

  ;; Did we call with a file name?
  if isa(vari, /string) then begin
    var = reform(red_readdata(vari))
    ;; This should be ok for FITS, ANA, MOMFBD format files, as long
    ;; as the contents is 2D.
  endif else begin
    ;; Image
    var=reform(vari)
  endelse
  
  dim=size(var)
  ;;Checks for the right image dimensions
  if dim[0] lt 2 or dim[0] gt 3 then begin
    print,'Wrong dimensions: '+stri(dim[0])+'. Image must be a 2D or 3D array'
    return
  endif
  ;; Color?
  iscolor = dim[0] eq 3
  if iscolor then begin
    xdim = dim[2]
    ydim = dim[3]
  endif else begin
    xdim = dim[1]
    ydim = dim[2]
  endelse
  
  ;; Screen
  sdim=get_screen_size()
  sdim[1]-=offs
  
  if keyword_set(reuse) then begin
    device, window_state=thesewindows
    window_exists = thesewindows[wnum]
    if window_exists then begin
      wset, wnum
      nowin = 1
    endif
  endif
  
  if sdim[0] ge xdim and sdim[1] ge ydim then begin
    ;; Image fits on screen
    if ~keyword_set(nowin) then $
       window, wnum, xsize=xdim, ysize=ydim, title = title
    tvscl, var, true = iscolor
  endif else begin
    ;; Image does not fit on screen
    if keyword_set(scroll) then begin
      scrollwindow, wid = wnum, xs=xdim, ys=ydim $
                    , title = title, sizefraction = 0.75 $
                    , free = ~keyword_set(nowin)
      if iscolor then begin
        tvscl, var, /true
      endif else begin
        if keyword_set(opt) then var = red_histo_opt(var)
        if not keyword_set(noscale) then var = bytscl(var)
        tv,var
      end
    endif else begin
      asp=float(xdim)/float(ydim)
      sasp=sdim[0]/sdim[1]
      if asp ge sasp then begin ;x-dimension bigger than y-dimension
        xsiz=sdim[0]
        ysiz=sdim[0]/asp
      endif else begin          ;y-dimension bigger than x-dimension
        xsiz=sdim[1]*asp
        ysiz=sdim[1]
      endelse
      if not keyword_set(nowin) then $
         window, wnum, xsize=xsiz, ysize=ysiz, title = title
      if iscolor then begin
        tvscl,congrid(var,3,xsiz,ysiz),/true
      endif else begin
        var = congrid(var,xsiz,ysiz)        
        if keyword_set(opt) then var = red_histo_opt(var)
        if not keyword_set(noscale) then var = bytscl(var)
        tv,var
      endelse
    endelse
  endelse
  
end
