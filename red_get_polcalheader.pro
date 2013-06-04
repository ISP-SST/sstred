function red_get_polcalheader, file, head, nsum, state, date = date
  res = ''
                                ;
                                ; Split header and filename
                                ;
  h1 = strsplit(head, ' =', /extract)
  h2 = strsplit(file_basename(file), '.',/extract)
  n2 = n_elements(h2)
                                ;
                                ; polcal sums have:
                                ; num time nx ny expt qual pref lre hre qwp nsum (LC#) ... POLCALSUM
                                ;
  num = string(long(h2[n2-1]), FORMAT='(I7)')
  time = strsplit(h1[21],'.',/extract)
  time = time[0] + '.' + string(float('0.'+time[1])*1000, format='(I03)')
  nx = h1[9]
  ny = h1[11]
  texp = red_stri((red_time2double(h1[21]) - red_time2double(h1[18])) *1000., ni='(F7.3)')
  qual = '   1000'
  pref = '0'                    ;h2[6]
  hre = '1234'
  lre = '1234'
  qwp = red_stri(float(strmid(h2[5], 2,3)), ni='(F6.2)')
  nim = red_stri(nsum, ni='(I5)')
  lc = '('+strupcase(h2[n2-3])+') ... POLCALSUM'
  date = strsplit(h1[20], '.',/extract)
  month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
  date = date[2]+month[long(date[1])-1]+date[0]
                                ;
  s = ' '
  res = num +s+s+ time +s+s+ nx +s+s+ ny +s+ texp +s+ qual +s+ pref +s+ lre $
        +s+ hre +s+ qwp +s+ nim +s+ lc
                                ;
  
  return, res
end
