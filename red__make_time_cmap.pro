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
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    cmap : 
;   
;   
;   
; 
; :Keywords:
; 
;    rot_dir  : 
;   
;   
;   
;    out_dir  : 
;   
;   
;   
;    result  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
pro red::make_time_cmap, cmap, rot_dir = rot_dir, out_dir = out_dir, result = cm
  inam = 'red::make_time_cmap : '
  if(~keyword_set(rot_dir)) then rot_dir = 0

                                ;
                                ; Load tseries data
                                ;
  file = file_search(self.out_dir+'/calib_tseries/tseries.*.calib.sav', count = ct)
  if(ct eq 0) then begin
     print, inam + 'ERROR -> could not find any calibration file in '+ self.out_dir +'/calib_tseries/'
     return
  endif
  
  print, inam + 'Found states:'
  pref = strarr(ct)
  for ii = 0, ct - 1 do begin
     pref[ii] = (strsplit(file_basename(file[ii]), '.', /extract))[1]
     print, ii, ' -> ' + pref[ii]
  endfor
  
  idx = 0
  if(ct ne 1) then begin
     read, idx, prompt = 'Select prefilter id: ' 
  endif
  print, inam + 'Selected state -> '+pref[idx]

  file = file[idx]
  restore, file
  nt = n_elements(ang)
  dim = size(cmap, /dim)
  nx = dim[0]
  ny = dim[1]




                                ;
                                ; prepare lp_header and open file
                                ;
  typestring = '(integer)'
  datatype = 2
  header = 'datatype='+strtrim(datatype,2)+' '+typestring
  header = header + ', dims='+strtrim(3,2)
  header = header + ', nx='+strtrim(nx,2)
  header = header + ', ny='+strtrim(ny,2)
  header = header + ', nt='+strtrim(nt,2)
  if ((byte(1L, 0, 1))[0] eq 1) then endianstr = 'endian=l'  $ ; little endian
  else endianstr = 'endian=b'                                  ; big endian
  header = header + ', '+endianstr
                                ;
  if(~keyword_set(out_dir)) then out_dir = fold + '/cfg/results/crispex/'
  file_mkdir, out_dir
  ofile = 'crispex.'+pref+'.tseries_cavity_map.icube'
  openw, lun, out_dir+'/'+ofile, /get_lun
  print, inam + 'writing to file: '+out_dir +'/'+ ofile


                                ;
                                ; Write header
                                ;
  hh = bytarr(512)
  hh1 = byte(header)
  hh[0:n_elements(hh1)-1] = hh1
  writeu, lun, hh


                                ;
                                ; create data cube and apply corrections
                                ;
  aa = assoc(lun, intarr(nx,ny), 512)

  for ss = 0L, nt - 1 do begin
     cmap1 = cmap
     aa[ss] = round(rotate(red_stretch(red_rotation(temporary(cmap1), ang[ss], total(shift[0,ss]), total(shift[1,ss])), reform(grid[ss,*,*,*])), rot_dir))
     print, string(13B), inam + 'correcting cmap -> ', 100./(nt-1.0) * ss,'%', FORMAT='(A,A,F5.1,A,$)'
  endfor
  print, ' '


  free_lun, lun

  stop
end
