; Does for data from SectorMean what RoundMatrix does for 
; data from RoundMean


function red_SectorMatrix, arrs, r, angles = angles;, xOffset=xOffset, yOffset=yOffset
  
  na = (size(arrs, /dim))[1]

  ;; Same default angles as red_sectormean
  if n_elements(angles) eq 0 then angles = (indgen(na)/float(na)*!pi)-!pi/2
  
  da = !pi/na

  ang = red_angular_coordinate(r,1) 

  rnd_matrix = fltarr(2*r, 2*r)
  weight_sum = fltarr(2*r, 2*r)

  for j = 0, na-1 do begin

     w1 = ang GT angles(j) - da
     w2 = ang LE angles(j) + da
     sector = w1 AND w2

     weight = (da*1.4-(abs(angles(j)-ang))) > 0
     ;weight = weight^2

     rnd_matrix += (weight*red_roundmatrix(arrs[*, j], r))^2
     weight_sum += weight^2

  endfor

  rnd_matrix = rnd_matrix + shift(reverse(reverse(rnd_matrix,1),2),1,1)
  weight_sum = weight_sum + shift(reverse(reverse(weight_sum,1),2),1,1)
  weight_sum(where(weight_sum eq 0)) =  1.

  rnd_matrix /= weight_sum
  rnd_matrix = sqrt(rnd_matrix)

  return, rnd_matrix

end
