function red_clipim, pic, cl
                                ;
                                ; Same as momfbd_clip (Pit Sutterlin)
                                ;
  IF cl[0] GT cl[1] THEN BEGIN
     mi_x = cl[1]-1
     ma_x = cl[0]-1
     fx = 1
  ENDIF ELSE BEGIN
     mi_x = cl[0]-1
     ma_x = cl[1]-1
     fx = 0
  ENDELSE
  
  IF cl[2] GT cl[3] THEN BEGIN
     mi_y = cl[3]-1
     ma_y = cl[2]-1
     fy = 1
  ENDIF ELSE BEGIN
     mi_y = cl[2]-1
     ma_y = cl[3]-1
     fy = 0
  ENDELSE
  
  CASE fx+2*fy OF 
     0:  r = 0
     1:  r = 5
     2:  r = 7
     3:  r = 2
  ENDCASE
  
  return, rotate(pic[mi_x:ma_x, mi_y:ma_y], r)
END

;
function red_camtag, file
  return, (strsplit(file_basename(file),'.',/extract))[0]
end
