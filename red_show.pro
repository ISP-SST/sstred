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
;     vari : in
;
;        2D or 3D image
; 
; 
; 
; :Keywords:
; 
;     wnum : in
;
;        integer window index
;
;     nowin : in, type=boolean
;
;         Set this to not recreate the window
;
;     offs : in, type=integet, default=57
;
;        integer removes pixel from the Y-dimension of
;        the screen so compensates for the menubar and so on.
;
;
;     title : in, optional, type=string
;
;        If a window is created, it will have this title.
;
;
; :history:
;
;    2013-07-24 : Renamed red_show for inclusion in crispred pipeline.
;
;    2013-09-20 : MGL. Added title keyword.
;
;-
pro red_show,vari,wnum=wnum,nowin=nowin,offs=offs,opt=opt,noscale=noscale, title = title

                                ;Initializes some variables
  if ~keyword_set(wnum) then wnum=0
  if ~keyword_set(offs) then offs=24
  var=reform(vari)
  dim=size(var)
  sdim=get_screen_size()
  sdim[1]-=offs
                                ;Checks for the right image dimensions
  if dim[0] lt 2 OR dim[0] gt 3 then begin
    print,'Wrong dimensions: '+stri(dim[0])+'.Image must be a 2D or 3D array'
    return
  endif

                                ;3D Array
  if dim[0] eq 3 then begin
                                ;If both dimensions are smaller than window dimensions
    if sdim[0] ge dim[2] AND sdim[1] ge dim[3] then begin
      if ~keyword_set(nowin) then window,wnum,xsize=dim[2],ysize=dim[3], title = title
      tvscl,var,/true
    endif else begin
      asp=float(dim[2])/float(dim[3])
      sasp=sdim[0]/sdim[1]
      
      if asp ge sasp then begin ;x-dimension bigger than y-dimension
        xsiz=sdim[0]
        ysiz=sdim[0]/asp
      endif else begin          ;y-dimension bigger than x-dimension
        xsiz=sdim[1]*asp
        ysiz=sdim[1]
      endelse
      if not keyword_set(nowin) then window,wnum,xsize=xsiz,ysize=ysiz, title = title
      tvscl,congrid(var,3,xsiz,ysiz),/true
    endelse
  endif
  
                                ;2D Array
  if dim[0] eq 2 then begin
    if sdim[0] ge dim[1] AND sdim[1] ge dim[2] then begin
      if ~keyword_set(nowin) then window,wnum,xsize=dim[1],ysize=dim[2], title = title
      tvscl,var
    endif else begin
      asp=float(dim[1])/float(dim[2])
      sasp=float(sdim[0])/float(sdim[1])
      if asp ge sasp then begin
        xsiz=sdim[0]
        ysiz=sdim[0]/asp
      endif else begin
        xsiz=sdim[1]*asp
        ysiz=sdim[1]
      endelse
      print, xsiz, ysiz, asp, sasp
      if ~keyword_set(nowin) then window,wnum,xsize=xsiz,ysize=ysiz, title = title
      var=congrid(var,xsiz,ysiz)        
      if keyword_set(opt) then var=red_histo_opt(var)
      if not keyword_set(noscale) then var=bytscl(var)
      tv,var
    endelse
  endif

  return
end
