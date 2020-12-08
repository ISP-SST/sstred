; docformat = 'rst'

;+
; Does a series of gridmatch calculations for two images, determining
; offset vectors on a grid.
;
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Dick Shine, LMSAL (original ANA version)
; 
; 
; :Returns:
;
;    The cumulative offset of all grids specified in vg. 
; 
; :Params:
; 
;    m1 : in, type="2D array"
;
;       An image.
;
;    m2 : in, type="2D array"
;
;       Another image.
;
;
;    vg : in, type=array
;
;       Subfield grid sizes. The last element specifies the grid of
;       the returned offsets.
;
;    clips : in, type=array
; 
;       Outlier clips in pixels.
;   
;   
;  :Keywords:
;
;    nthreads : in, type=int
;
;        Number of threads to pass to red_stretch_linear
; 
; :History:
; 
;     1993-02-10 : DS. Original ANA version. 
;
;     1994-01-03 : DS. Allow for rectangular images
;
;     1999-04-28 : BDP. Adapted from ANA.
; 
;     2000-05 : L. Strous. GRIDMATCH and STRETCH adapted from ANA and
;               coded into an IDL DLM file.
; 
;     2000-12-14 : TEB. Default clip conditions.
; 
;     2000-12-14 : TEB. Use IDL CONGRID function for faster bilinear
;                  interpolation. 
; 
;     201?-??-?? : JdlCR. Adapted to CRISPRED and Moved to red_
;                  namespace. 
; 
;     2016-11-17 : MGL. Added documentation header. Trust IDL to
;                  destroy local variables at exit.
;
;     2020-10-01 : JdlCR. Use stretch with nearest or bilinear
;                  interpolation. Bugfix: gx and hy must be of type
;                  int32, not float32!
;
;-
function red_dsgridnest, m1, m2, vg, clips, nthreads = nthreads

  bdp_size = size(m1)
  nx = bdp_size(1)
  ny = bdp_size(2)
  nest = n_elements(vg)

  ;; Default values for clips:
  if n_elements(clips) ne nest then begin
    print,'Grid clip array not specified - assigning default value of 20 pixels'
    clips = replicate(20,nest)
  end

  for k=0L, nest-1 do begin 
    n = vg(k)
    stretch_clip = clips(k) 
    ngw = fix(2.*nx/n)
    nw = fix(1.25*ngw)
    if(nx gt ny) then begin
      nxg = n
      nyg = long(float(n)*ny/nx +0.5) ;1/3/94 allow for rectangular images
                                ;note the wx calculation should be FP, 3/18/89
    endif else begin
      nyg = n
      nxg = long(float(n)*nx/ny +0.5) ;1/3/94 allow for rectangular images
      ;; Note the wx calculation should be FP, 3/18/89
    endelse

    wx = float(nx)/nxg
    wy = float(ny)/nyg
    ;; This is supposed to compute the grid the same way as Stu's programs
    gx = findgen(nxg)#(fltarr(nyg)+1.0) ;integer grid coordinates
    gx = long(gx*wx+wx/2.-1)                 ;in pixels, as required for gridmatch
    gy = (fltarr(nxg)+1.0)#findgen(nyg)
    gy = long(gy*wy+wy/2.-1)
    dx = nw
    dy = nw
    gwid = ngw

    if (k eq 0) then begin 
      displ = red_gridmatch(m1,m2,gx,gy,dx,dy,gwid,stretch_clip)
    end else begin   

      if n ne nprev then begin
        ;; Interpolate the old displ on the new grid and apply to m3
        disx = reform(displprev(0,*,*))
        disy = reform(displprev(1,*,*))
        fx = congrid(disx,nxg,nyg,/INT)
        fy = congrid(disy,nxg,nyg,/INT)
        prev = fltarr(2,n_elements(fx(*,0)),n_elements(fy(0,*)))
        prev(0,*,*) = fx
        prev(1,*,*) = fy
      end else prev = displprev

      m3 = red_stretch_linear(m2,prev, nthreads=nthreads)
      displnew = red_gridmatch(m1, m3, gx, gy, dx, dy, gwid, stretch_clip)
      displ = prev + displnew
    end

    if k lt nest-1 then begin   ;save variables for next cycle if needed
      nprev = n
      displprev = displ		
    end

  end 

  return, displ
end
