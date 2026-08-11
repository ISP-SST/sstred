; docformat = 'rst'

;+
; Detect raw data images collected while a filter wheel is still moving.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
;    Returns 1 if sync problems detected, 0 otherwise.
; 
; :Params:
; 
;    ims_or_filename : in, type="array[Nx,Ny,Nims] or string"
; 
;      Either an image cube or the path to a FITS file containing a cube.
; 
; 
; :Keywords:
; 
;    c_limit : in, optional, type=float, default=0.1
;   
;      Difference image contrast limit for detection.
;   
;    discard : out, optional, type=integer
; 
;      The number of frames that need to be discarded.
; 
; 
; :History:
; 
;   2023-05-22 : MGL. First version.
; 
;-
function red_detect_wheel_sync, ims_or_filename $
                                , c_limit = c_limit $
                                , discard = discard

  ;; Name of this subprogram
  inam = red_subprogram(/low, calling = inam1)
  
  if n_elements(c_limit) eq 0 then c_limit = 0.1

  if size(ims_or_filename, /tname) eq 'STRING' then begin
    filename = ims_or_filename
    ims = red_readdata(filename)
  endif else begin
    ims = ims_or_filename
  endelse
  
  dims = size(ims, /dim)

  Nx = dims[0]
  Ny = dims[1]
  Nims = dims[2]

  c = fltarr(Nims)
  for iim = 1, (Nims-1) do begin
    diff = ims[*, *, iim] - ims[*, *, iim-1]
    c[iim] = stddev(diff)/median(ims[*, *, iim])
  endfor

  undefine, discard
  if max(c) lt c_limit then return, 0 ; No bad frames detected

  discard = max(where(c ge c_limit)) 

  if n_elements(filename) then title = file_basename(filename)
  
  window, xs = 700, ys = 500, 1
  cgplot, c, psym = 16, color = 'red' $
          , xtitle = 'Frame index', ytitle = 'Diff RMS contrast' $
          , xrange = [-.1, Nims], yrange = [-.01, max(c)+.01], title = title
  cgplot, /over, !x.crange, c_limit*[1, 1]
  cgplot, /over, (discard-.5)*[1, 1], !y.crange

  Nxx = 100
  Nyy = Ny*(Nxx/float(Nx))
  scrollwindow, xs = Nims*Nxx, ys = Nyy, wid = wid
  for i=0,Nims-1 do tvscl, congrid(ims[*, *, i], Nxx, Nyy), i
  tv, bytarr(1, 255)+255b, discard*Nxx, 0

  
  print, inam + ' : Suggested discard '+strtrim(discard, 2)+' frames'

  if discard eq Nims-1 then $
     print, '--> Note that it is difficult to automatically see if the last frame is the (only) ok one.'
  
  ;; Give the user the option to change discard:
  s = ''
  read, 'Number of discarded frames ['+strtrim(discard, 2)+'] : ',s
  if s ne '' then begin
    discard = long(s)
  endif

  wdelete, wid
  wdelete, 1

  return, 1                     ; Bad frames detected
  
 
end

file = "/data/2023/2023-04/2023-04-27/CHROMIS-data/08:01:01/Chromis-N/sst_camXXX_00005_0002595_3950_3999_+0000.fits"
file = "/data/2023/2023-04/2023-04-27/CHROMIS-data/08:01:01/Chromis-N/sst_camXXX_00004_0002160_3950_3999_+0000.fits"

files = file_search('/storage/data/2023/2023-04/2023-04-27/CRISP-data/08:01:20/Crisp-R/sst_camXXXII_*_6173_6173_-0175*.fits', count = Nfiles)
files = file_search('/storage/data/2023/2023-04/2023-04-27/CRISP-data/08:01:20/Crisp-R/sst_camXXXII_*_5876_5876_-1254_*.fits', count = Nfiles)

discard = lonarr(Nfiles)
for ifile = 0, Nfiles-1 do begin
  if red_detect_wheel_sync(files[ifile], discard = nr) then begin
    discard[ifile] = nr
   endif
endfor

indx = where(discard gt 0, cnt)
if cnt gt 0 then hprint, strtrim(discard[indx], 2) + ' ' + files[indx]

end
