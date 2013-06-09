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
;    searchroot : 
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
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
pro red::make_time_series, searchroot, rot_dir = rot_dir
  
  inam = 'red::make_time_series : '
                                ;
                                ; files
                                ;
  if(~keyword_set(searchroot)) then begin
     print, inam + 'ERROR, keyword searchroot not set'
     return
  endif
                                ;
  
  files = file_search(searchroot, count = ct)
  if(ct lt 1) then begin
     print, inam + 'Error, no files found using '+searchroot 
     return
  endif
                                ;
  pref = (strsplit(file_basename(files[0]),'.',/extract))[2]
  fold = file_dirname(files[0])
                                ;
                                ; load time-series-calibration data
                                ; vars: tstep, clip, tile, scale, ang, shift, grid, time, date, wfiles, tmean
                                ;
  polish = 0B
  pfile = self.out_dir + '/calib_tseries/tseries.'+pref+'.calib.sav'
  if(file_test(pfile)) then begin
     print, inam + 'Using time-series calibration in ' + pfile
     restore, pfile
     polish = 1B
  endif else print, inam + 'WARNING, no time-series calibration found!'
                                ;
                                ; states
                                ;
  st = red_get_stkstates(files)
                                ;
  nscan = max([st.nscan, n_elements(ang)])
  nwav = st.nwav
                                ;
  dim = size(f0(st.ofiles[0]),/dim)
  nx = dim[0]
  ny = dim[1]
  cub = intarr(nx, ny, st.nwav, nscan)
                                ;
                                ; prepare lp_header and open file
                                ;
  typestring = '(integer)'
  datatype = 2
  header = 'datatype='+strtrim(datatype,2)+' '+typestring
  header = header + ', dims='+strtrim(3,2)
  header = header + ', nx='+strtrim(dim[0],2)
  header = header + ', ny='+strtrim(dim[1],2)
  header = header + ', nt='+strtrim(nscan*nwav,2)
  if ((byte(1L, 0, 1))[0] eq 1) then endianstr = 'endian=l'  $ ; little endian
  else endianstr = 'endian=b'                                  ; big endian
  header = header + ', '+endianstr
                                ;
;  if(~keyword_set(out_dir)) then out_dir = fold + '/cfg/results/crispex/'
  if(~keyword_set(out_dir)) then out_dir = fold + '/crispex'


  file_mkdir, out_dir
  ofile = 'crispex.'+pref+'.tseries_cavity_map.icube'
  openw, lun, out_dir+ofile, /get_lun
                                ;
                                ; Write header
                                ;
  hh = bytarr(512)
  hh1 = byte(header)
  hh[0:n_elements(hh1)-1] = hh1
  writeu, lun, hh
                                ;
                                ; Load and correct images
                                ;
  if(~keyword_set(rot_dir)) then rot_dir=0
  ntot = 100. / (nscan-1.0)
  for ss=0L, nscan-1 do begin
     for ww=0L, st.nwav-1 do begin
        cub[*,*,ww,ss] = fix(round(rotate(round(stretch(red_rotation(f0(st.ofiles[ww,ss])*1.e3, ang[ss], total(shift[0,ss]), total(shift[1,ss])), reform(grid[ss,*,*,*]))), rot_dir)))
        writeu, lun, cub[*,*,ww,ss]
     endfor
     print, string(13b),inam+'creating icube -> ', ntot * ss,'%', format='(A,A,F5.1,A,$)'
  endfor
  print, ' '
  free_lun, lun
  print, inam + 'cube saved to -> '+out_dir+ofile
                                ;
end
