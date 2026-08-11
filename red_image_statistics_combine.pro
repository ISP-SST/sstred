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
;    binsize : in, optional, type=float
;   
;       THe bin size of hist. Must be given for the percentiles to be
;       calculated.
; 
;    hist : in, optional, type=array
; 
;       A histogram for the entire cube. Must be given for the
;       percentiles to be calculated.
;   
;    old : in, optional, type=boolean
;   
;       Use old, deprecated algorithm.
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
;   2021-04-05: MGL. Implement algorithm by Pébay et al. for combining
;               the moments. 
; 
;-
function red_image_statistics_combine, statsarr $
                                       , binsize = binsize $
                                       , comments = comments $
                                       , hist = hist $
                                       , old = old

  tags = tag_names(statsarr[0])
  
  Nframes = n_elements(statsarr)

  DATARMS = statsarr.DATANRMS * statsarr.DATAMEAN
  CUBEMIN  = min(statsarr.datamin)
  CUBEMAX  = max(statsarr.datamax)
  
  if keyword_set(old) then begin

    ;; Calculation based on my own derivation.

    Npixels = statsarr[0].Ndatapix ; Assume same for all frames

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
    
    output = create_struct('NDATAPIX' , round(total(statsarr.Ndatapix)) $
                           , 'DATAMIN' , CUBEMIN $
                           , 'DATAMAX' , CUBEMAX $
                           , 'DATAMEAN', CUBEMEAN $       
                           , 'DATANRMS', CUBERMS/CUBEMEAN $
                           , 'DATASKEW', CUBESKEW $         
                           , 'DATAKURT', CUBEKURT)
  endif else begin

    ;; Calculation based on Eq. (3.1) of Pébay, P., Terriberry, T. B.,
    ;; Kolla, H., & Bennett, J. 2016, Computational Statistics, 31,
    ;; 1305.

    ;; Note that what we call moments here, [mean, variance, skewness,
    ;; kurtosis] are not the clean 1/n SUM_i (x_i - x_mean)^p of Pébay
    ;; et al., so our values have to be modified into the assumed
    ;; form. And then modified back.

    ;; Initialize with statistics for 0th frame
    Npix = double(statsarr[0].Ndatapix)
    CUBEMEAN = statsarr[0].DATAMEAN
    variance = (statsarr[0].DATANRMS * statsarr[0].DATAMEAN)^2
    Mp = [ variance * (Npix-1)/Npix $                      ; M2 from variance
           , statsarr[0].DATASKEW * variance^(3/2.) $      ; M3 from skewness
           , (statsarr[0].DATAKURT+3) * variance^(4/2.)  $ ; M4 excess kurtosis
         ] * Npix

    ;; Combine pairwise current accumulated frames with next frame
    for iframe = 1, Nframes-1 do begin

      ;; Set A is the so far accumulated frames
      NpixA = Npix
      datameanA = CUBEMEAN
      MpA = Mp

      ;; Set B is the next frame
      NpixB = double(statsarr[iframe].Ndatapix)
      datameanB = statsarr[iframe].DATAMEAN
      variance = (statsarr[iframe].DATANRMS * statsarr[iframe].DATAMEAN)^2
      MpB = [ variance * (NpixB-1)/NpixB $                         ; M2 from variance
              , statsarr[iframe].DATASKEW * variance^(3/2.) $      ; M3 from skewness
              , (statsarr[iframe].DATAKURT+3) * variance^(4/2.)  $ ; M4 from excess kurtosis
            ] * NpixB

      ;; Combined set
      Npix = NpixA + NpixB
      deltaBA = datameanB - datameanA
      CUBEMEAN = datameanA + deltaBA * NpixB / Npix
      for p = 2, 4 do begin     ; Eq. (3.1) 
        Mp[p-2] = MpA[p-2] + MpB[p-2] $
                  + NpixA * (-NpixB/Npix * deltaBA )^p $
                  + NpixB * ( NpixA/Npix * deltaBA )^p 
        for k = 1, p-2 do $
           Mp[p-2] += red_bico(p, k) * deltaBA^k * ( MpA[p-k-2]*(-NpixB/Npix)^k + MpB[p-k-2]*(NpixA/Npix)^k )
      endfor                    ; p

    endfor                      ; iframe

    ;; Construct the IDL-definition moments
    mom =  [CUBEMEAN $
            , Mp[0] / Npix $    ; Variance
            , Mp[1] / Npix $    ; Skewness
            , Mp[2] / Npix $    ; Kurtosis
           ]

    mom[1] *= Npix/(Npix-1.)    ; Adjust Variance
    mom[2] /= mom[1]^(3/2.)     ; Adjust Skewness 
    mom[3] /= mom[1]^(4/2.)     ; Adjust Kurtosis 
    mom[3] -= 3                 ; Excess kurtosis

    ;; And make a struct with the results
    output = create_struct('NDATAPIX' , long64(Npix) $
                           , 'DATAMIN' , CUBEMIN $
                           , 'DATAMAX' , CUBEMAX $
                           , 'DATAMEAN', CUBEMEAN $       
                           , 'DATANRMS', sqrt(mom[1]) / CUBEMEAN $
                           , 'DATASKEW', mom[2] $         
                           , 'DATAKURT', mom[3] $
                          )

    
  endelse


  if arg_present(comments) then $
     comments = create_struct('NDATAPIX' , 'Number of pixels' $
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

