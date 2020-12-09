; docformat = 'rst'

;+
; Combine statistics for multiple images into statistics for a data
; cube. 
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
;    A struct with statistics.
; 
; 
; :Params:
; 
;    statsarr : in, type=structarr
; 
;       An array of structs as returned by
;       red_image_statistice_calculate.
; 
; :Keywords:
; 
;    hist : in, optional, type=array
; 
;       A histogram for the entire cube. Must be given for the
;       percentiles to be calculated.
;   
;    binsize : in, optional, type=float
;   
;       THe bin size of hist. Must be given for the percentiles to be
;       calculated.
; 
; 
; :History:
; 
;   2018-10-26: MGL. First version.
; 
;   2019-05-29: MGL. Improve percentiles by interpolating in the
;               cumulative histogram.
; 
;   2020-11-30: MGL. Add DATANRMS, remove DATARMS.
; 
;-
function red_image_statistics_combine, statsarr $
                                       , binsize = binsize $
                                       , comments = comments $
                                       , hist = hist 

  tags = tag_names(statsarr[0])
  
  Nframes = n_elements(statsarr)
  Npixels = statsarr[0].Npixels ; Assume same for all frames

  DATARMS = statsarr.DATANRMS * statsarr.DATAMEAN
  
  CUBEMIN  = min(statsarr.datamin)
  CUBEMAX  = max(statsarr.datamax)
  CUBEMEAN = mean(statsarr.datamean)

  CUBERMS  = sqrt(1d/(Nframes*Npixels-1.d) $
                  * ( total( (Npixels-1d) * DATARMS^2 $
                             + Npixels * statsarr.DATAMEAN^2 ) $
                      - Nframes*Npixels * CUBEMEAN^2 ))
  
  CUBESKEW = 1.d/(Nframes*Npixels*CUBERMS^3) $
             * total( Npixels * DATARMS^3 * statsarr.DATASKEW $
                      + Npixels * statsarr.DATAMEAN^3 $
                      + 3*(Npixels-1) * statsarr.DATAMEAN * DATARMS^2) $
             - CUBEMEAN^3 / CUBERMS^3 $
             - 3 * (Npixels*Nframes-1d) * CUBEMEAN / (Npixels*Nframes*CUBERMS)
  
  CUBEKURT = 1.d/(Nframes*Npixels*CUBERMS^4) $
             * total( Npixels * DATARMS^4 * (statsarr.DATAKURT+3) $
                      + 4*Npixels * statsarr.DATAMEAN * DATARMS^3 * statsarr.DATASKEW $
                      + 6*(Npixels-1) * statsarr.DATAMEAN^2 * DATARMS^2 $
                      + Npixels * statsarr.DATAMEAN^4 ) $
             - 4 * CUBEMEAN * CUBESKEW / CUBERMS $
             - 6 * (Npixels*Nframes-1d) * CUBEMEAN^2 / (Npixels*Nframes*CUBERMS^2) $
             - CUBEMEAN^4/CUBERMS^4 - 3 

  
  output = create_struct('NPIXELS' , Nframes*Npixels $
                         , 'DATAMIN' , CUBEMIN $
                         , 'DATAMAX' , CUBEMAX $
                         , 'DATAMEAN', CUBEMEAN $       
                         , 'DATANRMS', CUBERMS/CUBEMEAN $
                         , 'DATASKEW', CUBESKEW $         
                         , 'DATAKURT', CUBEKURT)

  if arg_present(comments) then $
     comments = create_struct('NPIXELS' , 'Number of pixels' $
                              , 'DATAMIN' , 'The minimum data value' $
                              , 'DATAMAX' , 'The maximum data value' $
                              , 'DATAMEAN', 'The average data value' $       
                              , 'DATANRMS', 'The normalized RMS deviation from the mean' $
                              , 'DATASKEW', 'The skewness of the data' $         
                              , 'DATAKURT', 'The excess kurtosis of the data' )
  
  ;; Percentiles from the cumulative histogram but for the P values of
  ;; the percentiles in statsarr.
  
  if n_elements(hist) eq 0 or n_elements(binsize) eq 0 then return, output
  
  indx = where(strmatch(tags,'DATAP*') or strmatch(tags,'DATAMEDN'), Nindx)
  perc = fltarr(Nindx)
  for i = 0, Nindx-1 do begin    
    if tags[indx[i]] eq 'DATAMEDN' then begin
      perc[i] = 0.5
    endif else begin
      perc[i] = strmid(tags[indx[i]], 5)/100.
    endelse
  endfor

  hist_cum = total(hist, /cum, /double) / total(hist, /double)
  binloc = cubemin + ( findgen(n_elements(hist_cum))+1. )*binsize
  p_values = INTERPOL( binloc, hist_cum, perc, /quad )

  ;; p_values will have non-finite values where there are multiple
  ;; data points that are exactly equal to one of the elements of
  ;; perc.
  findx = where(~finite(p_values), Nnonfinite)
  for inonfinite = 0, Nnonfinite-1 do begin
    pindx = where(hist_cum eq perc[findx[inonfinite]], Nmatch)
    ;; Pick the first value
    if Nmatch gt 0 then p_values[findx[inonfinite]] = binloc[pindx[0]]
  endfor                        ; inonfinite

  ;; Apparently some non-finite values are not fixed by the above
  ;; procedure. Linear interpolation seems to fix it.
  findx = where(~finite(p_values), Nnonfinite)
  p_values[findx] = INTERPOL( binloc, hist_cum, perc[findx] )
  
  for i = 0, Nindx-1 do begin
    output = create_struct(tags[indx[i]] $
                           , p_values[i] $
                           , output )
    if arg_present(comments) then $
       comments = create_struct(tags[indx[i]] $
                                , 'The '+string(round(perc[i]*100), format = '(I02)') $
                                + ' percentile of the data' $
                                , comments)
    
  endfor                        ; i

  return, output

end


Nims = 10
ims = randomn(seed, 100, 100, Nims, /double)

cubestats = red_image_statistics_calculate(ims)

imstats = replicate(cubestats, Nims)
for i = 0, Nims-1 do imstats[i] = red_image_statistics_calculate(ims[*, *, i])

combstats = red_image_statistics_combine(imstats)

end

