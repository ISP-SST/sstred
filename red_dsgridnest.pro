FUNCTION red_dsgridnest,m1,m2,vg,clips

; 2/10/93	does a series of gridmatch calculations with grid sizes
; specified in vg, returns an offset grid for the last value of vg
; clips to use for each grid are in pixels?
; Grid returned is the cumulative offset of all grids specified in vg.

; Adapted from ANA on April 28, 1999. BDP.
; GRIDMATCH and STRETCH adapted from ANA by L. Strous and coded into
;  an IDL DLM file. 5/2000. 
; Default clip conditions: 12/14/00, TEB.
; Use IDL CONGRID function for faster bilinear interpolation: 12/14/00, TEB.

 bdp_size = SIZE(m1)
 nx = bdp_size(1)
 ny = bdp_size(2)
 nest = N_ELEMENTS(vg)

;Default values for clips:
if N_ELEMENTS(clips) ne nest then begin
  PRINT,'Grid clip array not specified - assigning default value of 20 pixels'
  clips = REPLICATE(20,nest)
end

for k=0,nest-1 do begin 
  n = vg(k)
  stretch_clip = clips(k) 
  ngw = FIX(2.*nx/n)
  nw = FIX(1.25*ngw)
  if(nx gt ny) then begin
     nxg = n
     nyg = LONG(FLOAT(n)*ny/nx +0.5) ;1/3/94 allow for rectangular images
                                    ;note the wx calculation should be FP, 3/18/89
  endif else begin
     nyg = n
     nxg = LONG(FLOAT(n)*nx/ny +0.5) ;1/3/94 allow for rectangular images
                                ;note the wx calculation should be FP, 3/18/89
  endelse

  wx = FLOAT(nx)/nxg
  wy = FLOAT(ny)/nyg
  ;this is supposed to compute the grid the same way as Stu's programs
  gx = FINDGEN(nxg)#(FLTARR(nyg)+1.0) ;integer grid coordinates
  gx = gx*wx+wx/2.-1                  ;in pixels, as required for gridmatch
  gy = (FLTARR(nxg)+1.0)#FINDGEN(nyg)
  gy = gy*wy+wy/2.-1 
  dx = nw
  dy = nw
  gwid = ngw

  if (k eq 0) then begin 
      displ = GRIDMATCH(m1,m2,gx,gy,dx,dy,gwid,stretch_clip)
;      PRINT,'Mean initial = ',STRTRIM(AVG(ABS(displ)),2),$
;            '   Max(abs) = ',STRTRIM(MAX(ABS(displ)),2)
  end else begin   

      if n ne nprev then begin
           ;interpolate the old displ on the new grid and apply to m3
        disx = REFORM(displprev(0,*,*))
        disy = REFORM(displprev(1,*,*))
        fx = CONGRID(disx,nxg,nyg,/INT)
        fy = CONGRID(disy,nxg,nyg,/INT)
        prev = FLTARR(2,N_ELEMENTS(fx(*,0)),N_ELEMENTS(fy(0,*)))
        prev(0,*,*) = fx
        prev(1,*,*) = fy
      end else prev = displprev

      m3=STRETCH(m2,prev)
      displnew = GRIDMATCH(m1, m3, gx, gy, dx, dy, gwid,stretch_clip)
;      PRINT,'Mean incremental = ',STRTRIM(AVG(ABS(displnew)),2),$
;            '   Max(abs) = ',STRTRIM(MAX(ABS(displnew)),2)
      displ = prev + displnew
  end

  if k lt nest-1 then begin  ;save variables for next cycle if needed
    nprev = n
    displprev = displ		
;    gxprev = gx
;    gyprev = gy
  end

end 

;PRINT,'Mean total displacement = ',STRTRIM(AVG(ABS(displ)),2),$
;      '   Max(abs) = ',STRTRIM(MAX(ABS(displ)),2)

;recover some memory
m3=0
disx=0
disy=0
fx=0
fy=0
prev=0
displnew=0
gx = 0
gy = 0
;gxprev=0
;gyprev=0
 
RETURN, displ
END
