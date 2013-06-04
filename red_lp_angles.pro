function red_lp_angles, time, date
  ave = dblarr(n_elements(time))
  for i = 0L, n_elements(time) - 1 do begin
                                ;head=fzhead(names[i])
                                ;head=strsplit(head,' ,=,-,: ',/extract)
                                ;get_time,head,yy,mm,dd,hh,mi,ss
     tt = double(strsplit(time[i], ':' , /extract))
     dat = double(strsplit(date[i],'-.', /extract))
                                ;
     red_get_sun,dat[0],dat[1],dat[2],red_reform_frac_time(tt[0],tt[1],tt[2]),ha,dec
                                ;
     red_get_azel,ha,dec,az,el
     ave[i]=red_get_rot(az,el,dec)
  endfor
  return,ave
end
