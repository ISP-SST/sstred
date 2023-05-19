
pro bug_266, fname

  ims = red_readdata(fname, h = header)

  dims = size(ims, /dim)

  Nx = dims[0]
  Ny = dims[1]
  Nims = dims[2]

  wset, 0
  ;; If you have more than 16 frames per state, you need to change the
  ;; line below and/or the size of window 0 as defined in the main
  ;; program below. The number 15 and/or the number 4 may need
  ;; modification here:
  for iim = 0, (Nims-1) <15 do tvscl, rebin(ims[*, *, iim], Nx/4, Ny/4, /samp), iim
    
  framenum = fxpar(header, 'FRAMENUM', count = cnt)
  if cnt ne 1 then stop
  
  frameinc = fxpar(header, 'FRAMEINC', count = cnt)
  if cnt eq 0 then frameinc = 1

  framenums = indgen(Nims)*frameinc + framenum
  print
  print, fname


  c = fltarr(Nims)
  for iim = 1, (Nims-1) do begin
    diff = ims[*, *, iim] - ims[*, *, iim-1]
    c[iim] = stddev(diff)/median(ims[*, *, iim])
  endfor
  wset, 1
  cgplot, c, psym = 16, xtitle = 'Frame', ytitle = 'Diff RMS contrast'

  ;; You may want to change c_limit
  c_limit = 0.1
  if max(c) gt c_limit then begin
    s = ''
    for i=0,Nims-1,4 do print, framenums[i:(i+3) <(Nims-1)]
    read, s
  endif
  
  wset, 0
  erase
  
end

a = chromisred(/dev, /no)

isodate = file_basename(file_dirname(getenv('PWD')))
date = strsplit(isodate, '-', /extract)

;; Construct the path to your raw data top directory for the day here:
ddir = '/data/' + date[0] + '/' + strjoin(date[0:1], '-') + '/' + strjoin(date, '-') + '/'

;; Change the "*_3999_*" part to examine other tuning points
files = file_search(ddir+'CHROMIS-data/*/Chromis-N/*_3999_*', count = Nfiles)

;; The size of window 0 is based on the number of frames per tuning
;; point being LE 16. If you have more, you need to change this and/or
;; the loop that displays the images in check_frames above.
window, xs = 1920, ys = 1200, 0
window, xs = 700, ys = 500, 1

for ifile = 0, Nfiles-1 do bug_266, files[ifile]


end
