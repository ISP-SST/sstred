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
;    Jaime de la Cruz Rodriguez (inspired on Mats Löfdahl's makewindow.pro)
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    arg1 : 
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
function red_taper, arg1

  Sz = long(arg1[0])
  ZM = long(arg1[1])
  SM = long(arg1[2])
                                
  IF 2*(ZM+SM) GT Sz THEN BEGIN
     print,'MakeWindow -- Margins too large for Sz'
     print,Sz,ZM,SM
     RetAll 
  END 
                                
  OneSize = Sz-2*(ZM+SM)
  
  y = FltArr(Sz)
  IF OneSize GT 0 THEN BEGIN
     OnePart = FltArr(OneSize) + 1.0
     y[ZM+SM] = OnePart
  END 
  
  IF SM GT 0 THEN BEGIN
     x = float(IndGen(SM))
     SoftPart = 0.5*(1-cos(2*!pi*x/(2*SM-1)))
     y[ZM] = SoftPart
     y[Sz-(ZM+SM)] = Reverse(SoftPart)
  END 
  
  y2 = FltArr(Sz,Sz)
  
  FOR i = 0,Sz-1 DO y2[0,i] = y
                                
  Return, y2 * transpose(y2)
end
