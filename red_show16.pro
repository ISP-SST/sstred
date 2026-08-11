;; Show arrays of dims [16,Nx,Ny] as a grid of 4 by 4 images. Useful
;; for polcal output and similar.
;;
;; 
pro red_show16, ims, fac = fac, mos = mos, title = title, fname = fname $
                , scroll = scroll, pngfile = pngfile
  
  if n_elements(fac) eq 0 then fac = 2.0

  if isa(ims,/string) then begin
    fname = ims
  endif

  if n_elements(fname) gt 0 then begin
    ims = red_readdata(fname)
    if n_elements(title) eq 0 then title = file_basename(fname)
  endif
  
  dims = size(ims, /dim)

  if dims[0] ne 16 then stop

  mos = bytarr(4*dims[1], 4*dims[2])

  i = 0
  for ix = 0, 3 do for iy = 0, 3 do begin

    im = reform(ims[i, *, *])

    ;; Set "missing" pixel values to NaN
    red_missing, im, /inplace

    ;; Set min and max based on biweight statistics so outer, not well
    ;; defined values do not decide the scaling.
    mid = biweight_mean(im(where(finite(im))),sigma)
    mn = mid - sigma*fac
    mx = mid + sigma*fac

    mos[ix*dims[1]:(ix+1)*dims[1]-1, iy*dims[2]:(iy+1)*dims[2]-1] = bytscl(im, mn, mx)

    i++

  endfor                        ; ix,iy

  red_show, mos, title = title $
            , scroll = scroll

  if n_elements(pngfile) ne 0 then write_png, pngfile, bytscl(mos)
  
end

fname = '/scratch_local/mats/2025-08-03/CRISP2-test/polcal/camXXXII_8542_polcal.fits'

red_show16, fname = fname, png = red_strreplace(file_basename(fname), '.fits', '.png')


stop

ims = readfits(fname)


red_show16, ims, title = file_basename(fname)

end
