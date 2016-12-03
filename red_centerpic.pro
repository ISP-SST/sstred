;	   File: 	CenterPic.ANA
;	   Created:	Tue Feb  5 10:00:53 1991
;	   Author: 	LOFDAHL@ROYACS

; Ported from ANA in 2011.

; Return an xSize by ySize version of Pic, that is either cropped or
; framed or both. Background is z.

; If xOffset and/or yOffset are given, the image is shifted by those
; amounts before the image is cropped. If it is to be framed, the
; shift is done afterwards. 

; If the keyword Fourier is set, a Fourier shift (i.e., by half the
; dimensions of the array) is done both before and after the
; cropping/framing. So, the "center" in "centerpic" then refers to the
; origin in a Fourier transform, which is in [0,0] rather than in
; [xsize/2,ysize/2]. 

function red_centerpic, Pic $
                    , Sz = Sz $
                    , xSize = xSize $
                    , ySize = ySize $
                    , z = z $
                    , xOffset = xOffset $
                    , yOffset = yOffset $
                    , Fourier = Fourier
 
  ;; The old dimensions
  xDimen = (size(Pic,/dim))[0]
  yDimen = (size(Pic,/dim))[1]

  if keyword_set(Fourier) then NewPic = shift(Pic, xDimen/2, yDimen/2) else NewPic = Pic

  ;; The new dimensions
  if n_elements(Sz) ne 0 then begin
     xSize = Sz
     ySize = Sz
  end else begin
     if n_elements(xSize) eq 0 then xSize = xDimen
     if n_elements(ySize) eq 0 then ySize = yDimen
  end

  if n_elements(xOffset) eq 0 then xOffset = 0
  if n_elements(yOffset) eq 0 then yOffset = 0

  ;; Background
  if n_elements(z) eq 0 then z = 0.0
  
  ;; Framing
  if (xDimen lt xSize) then begin
     tmp = replicate(Pic[0,0]*0, xSize, yDimen) + z
     tmp[(xSize-xDimen)/2:(xSize-xDimen)/2+xDimen-1, *] = NewPic
     NewPic = tmp
     xDimen = xSize
  endif

  if (yDimen lt ySize) then begin
     tmp = replicate(Pic[0,0]*0, xDimen, ySize) + z
     tmp[*, (ySize-yDimen)/2:(ySize-yDimen)/2+yDimen-1] = NewPic
     NewPic = tmp
     yDimen = ySize
  endif


  ;; Shifting
  NewPic = shift(NewPic, xOffset, yOffset)


  ;; Cropping
  if (xDimen gt xSize) then begin
     xLow = (xDimen-xSize)/2
     xHi  = xLow + xSize - 1
     NewPic = NewPic[xLow:xHi, *]
     xDimen = xSize
  endif

  if (yDimen gt ySize) then begin
     yLow = (yDimen-ySize)/2
     yHi  = yLow + ySize - 1
     NewPic = NewPic[*, yLow:yHi]
     yDimen = ySize
  endif

  if keyword_set(Fourier) then return, shift(NewPic, xDimen/2, yDimen/2)

  return, NewPic
  
END; CenterPic


