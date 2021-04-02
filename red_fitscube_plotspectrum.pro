; docformat = 'rst'

;+
; Plot median spectrum of fitscube with atlas.
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
; 
; :Params:
; 
;   filename : in, type=string
; 
;      The name of the file to plot the spectrum of.
; 
; :Keywords:
; 
;   axis_numbers : in, out, optional, type=array
;
;      The fitscube pixel coordinates axis numbers of the
;      frame_statistics array. Used for multiple calls with the same
;      filename.
; 
;   frame_statistics : in, out, optional, type=array
;
;      Intensity statistics for the frames of the fitscube. Used for
;      multiple calls with the same filename.
;
;   nosave : in, optional, type=boolean
;
;      Do not save the plot as a pdf file.
;
;   himargin : in, optional, type=array
;
;      Extend upper end of xrange by this many percent of the
;      wavelength range of the data.
;     
;   lomargin : in, optional, type=array
;
;      Extend lower end of xrange by this many percent of the
;      wavelength range of the data.
;     
;   xrange : in, optional, type=array
;
;      Set xrange of plot explicitly.
; 
; :History:
; 
;   2021-03-02 : MGL. First version.
; 
;-
pro red_fitscube_plotspectrum, filename $
                               , axis_numbers = axis_numbers $
                               , frame_statistics = frame_statistics $
                               , himargin = himargin $
                               , lomargin = lomargin $
                               , nosave = nosave $
                               , xrange = xrange

  if n_elements(lomargin) eq 0 then lomargin = 15. ; Percent of range
  if n_elements(himargin) eq 0 then himargin = 15. ; Percent of range
  
  ;; Get FITS header
  hdr = headfits(filename)

  Ntun    = fxpar(hdr, 'NAXIS3')
  Nstokes = fxpar(hdr, 'NAXIS4')
  Nscans  = fxpar(hdr, 'NAXIS5')

  instrument = strtrim(fxpar(hdr, 'INSTRUME'), 2)
  prefilter  = strtrim(fxpar(hdr, 'FILTER1'), 2)
  units      = strtrim(fxpar(hdr, 'BUNIT'), 2)

  red_fitspar_getdates, hdr $
                        , date_beg = date_beg $
                        , date_end = date_end $
                        , date_avg = date_avg 

  date_split = strsplit(date_avg, 'T', /extract)
  date = date_split[0]
  time = date_split[1]
  
  ;; Get WCS coordinates 
  red_fitscube_getwcs, filename, coordinates = coordinates
  lambda = coordinates[*, 0].wave[0,0] ;Wavelengths in nm
  lambda_min = min(lambda) 
  lambda_max = max(lambda)  
  lambda_delta = lambda_max-lambda_min
  lambda_min -= lambda_delta * lomargin/100.  
  lambda_max += lambda_delta * himargin/100.

  ;; Get statistics
  if n_elements(axis_numbers) gt 0 and n_elements(frame_statistics) gt 0 then begin
    ;; Reuse previously calculated statistics
  end else begin
    red_fitscube_statistics, filename, frame_statistics, axis_numbers = axis_numbers
  endelse
  datamedn = frame_statistics.datamedn
  case 1 of
    array_equal(axis_numbers, [3, 4])    :
    array_equal(axis_numbers, [3, 5])    : 
    array_equal(axis_numbers, [3, 4, 5]) : datamedn = reform(datamedn[*, 0, *])
    else : stop
  endcase
  
  ;; Get mu and zenith angle
  red_logdata, date, time $
               , mu = mu $
               , zenithangle = zenithangle
  
  ;; Get the atlas
  red_satlas, lambda_min, lambda_max, /nm $
              , atlas_lambda, atlas_spectrum $
              , /si, cont = cont 

  ;; We may want to convolve the atlas spectrum here
  ;; Make FPI transmission profile
  dw = atlas_lambda[1] - atlas_lambda[0]
  if instrument eq 'CHROMIS' then begin
    np = round((0.080 * 8) / dw)
    if np/2*2 eq np then np -=1
    tw = (dindgen(np)-np/2)*dw + double(prefilter)
    tr = chromis_profile(tw, erh=-0.09d0)
  endif else begin
    np = long((max(atlas_lambda) - min(atlas_lambda)) / dw) - 2
    if np/2*2 eq np then np -=1
    tw = (dindgen(np)-np/2)*dw                                             
    tr = crisp_fpi_profile(tw, prefilter, erh=-0.01d, /offset_correction)
  endelse
  tr /= total(tr)
  atlas_spectrum_convolved = fftconvol(atlas_spectrum, tr)

  if n_elements(xrange) eq 0 then xrange = [lambda_min, lambda_max]

  ;; Adapt units for cgplot
  plunits = red_strreplace(units, '^-1', '$\exp-1$', n = 3)
  plunits = red_strreplace(plunits, '^-2', '$\exp-2$')
  
  ;; Make the plot
  cgwindow
  cgplot, /add, atlas_lambda/10, atlas_spectrum*1e9 $
          , xtitle = '$\lambda$ / 1 nm', ytitle = 'median(Intensity) / 1 n'+plunits $
          , xrange = xrange
  for iscan = 0, Nscans-1 do cgplot, /add, /over, lambda, datamedn[*, iscan]*1e9, psym = 9, color = 'red'

  plfile = file_dirname(filename) + '/' + file_basename(filename, '.fits') + '.pdf'
  if ~keyword_set(nosave) then cgcontrol, output = plfile
  
end

case 1 of
  0 : begin
    undefine, axis_numbers, frame_statistics
    cd, '/scratch/mats/2016.09.19/CRISP-aftersummer/'
    filename = 'cubes_nb/nb_6302_2016-09-19T09:30:20_scans=0-2_stokes_corrected_im.fits'
    filename = 'cubes_nb_test/nb_6302_2016-09-19T09:30:20_scans=2-8_stokes_corrected_im.fits'
    filename = 'cubes_TEST/nb_6302_2016-09-19T09:30:20_scans=2,3_stokes_corrected_im.fits'
    filename = 'cubes_TEST/nb_6302_2016-09-19T09:30:20_scans=2,3_corrected_im.fits'
    red_fitscube_plotspectrum, filename $
                               , axis_numbers = axis_numbers $
                               , frame_statistics = frame_statistics
  end
  1 : begin
    undefine, axis_numbers, frame_statistics
    cd, '/scratch/olexa/2020-10-16/CHROMIS/'
    filename = 'cubes_nb/nb_3950_2020-10-16T09:11:04_scans=0-3,5,6,8,9_corrected_im.fits'
    red_fitscube_plotspectrum, filename $
                               , axis_numbers = axis_numbers $
                               , frame_statistics = frame_statistics
    red_fitscube_plotspectrum, filename $
                               , axis_numbers = axis_numbers $
                               , frame_statistics = frame_statistics $
                               , xrange=[393,393.7]
  end
  2 : begin
    undefine, axis_numbers, frame_statistics
    cd, '/scratch/mats/2016.09.19/CHROMIS-jan19'
    filename = 'cubes_TEST/nb_3950_2016-09-19T09:28:36_scans=0-3_corrected_im.fits'
    red_fitscube_plotspectrum, filename $
                               , axis_numbers = axis_numbers $
                               , frame_statistics = frame_statistics
    red_fitscube_plotspectrum, filename $
                               , axis_numbers = axis_numbers $
                               , frame_statistics = frame_statistics $
                               , xrange=[393,393.7]
  end
endcase



end

