; docformat = 'rst'

;+
; Get calibration curves etc from the prefilter fits.
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
;   states : in, type=structarr
; 
;     File "states" for which we want the curves.
; 
; 
; :Keywords:
; 
;   data_time : in, optional, type=string  
;   
;     Timestamp of the science dataset for which we want the calibration.
; 
;   fitpref_time : in, optional, type=string  
;   
;     Timestamp of the calibration data set.
; 
;   prefilter_curve : out, optional, type=fltarr
;   
;     The prefilter fit curve, possibly the average of T and R cameras.
; 
;   prefilter_wav : out, optional, type=fltarr
;   
;     The prefilter fit wavelengths.
; 
;   r_prefilter_curve : out, optional, type=fltarr
;   
;     The prefilter fit curve of the R camera.
; 
;   t_prefilter_curve : out, optional, type=fltarr
;   
;     The prefilter fit curve of the T camera.
; 
;   units  : out, optional, type=string
;   
;     The units after applying the prefilter curve.
; 
;   wave_shift : out, optional, type=fltarr
;   
;     The wavelength shift from the fit.
; 
; 
; :History:
; 
;   2025-06-08 : MGL. First version.
; 
;-
pro red::get_prefilterfit, states $
                           , data_time = data_time $
                           , fitpref_time = fitpref_time $
                           , prefilter_curve = prefilter_curve $
                           , prefilter_wav = prefilter_wav $
                           , prefilter_wb = prefilter_wb $
                           , prf = prf $
                           , r_prefilter_curve = r_prefilter_curve $
                           , r_prefilter_wb = r_prefilter_wb $
                           , t_prefilter_curve = t_prefilter_curve $
                           , t_prefilter_wb = t_prefilter_wb $
                           , tun_wavelengths = tun_wavelengths $
                           , units = units $
                           , wave_shifts = wave_shifts


  ;; Camera/detector identification
  self -> getdetectors
  wbindx      = where(strmatch(*self.cameras,'*-W'))
  wbcamera    = (*self.cameras)[wbindx[0]]
  wbdetector  = (*self.detectors)[wbindx[0]]
  nbtindx     = where(strmatch(*self.cameras,'*-T')) 
  nbtcamera   = (*self.cameras)[nbtindx[0]]
  nbtdetector = (*self.detectors)[nbtindx[0]]
  nbrindx     = where(strmatch(*self.cameras,'*-R')) 
  nbrcamera   = (*self.cameras)[nbrindx[0]]
  nbrdetector = (*self.detectors)[nbrindx[0]]

  instrument = (strsplit(wbcamera, '-', /extract))[0]

  dims = size(states, /dim)
  
  prefilters = states.prefilter
  
  Ntunings    = n_elements(tun_wavelengths)

  wave_shifts       = fltarr(dims)
  prefilter_curve   = fltarr(dims)
  prefilter_wb      = fltarr(dims)
  r_prefilter_curve = fltarr(dims)
  r_prefilter_wb    = fltarr(dims)
  r_prefilter_wav   = fltarr(dims)
  t_prefilter_curve = fltarr(dims)
  t_prefilter_wb    = fltarr(dims)
  t_prefilter_wav   = fltarr(dims)
  
  uprefs = red_uniquify(states.prefilter, count = Nprefilters)

  ;; What prefilter fits directory to use?
  if ~keyword_set(fitpref_time) then begin
    fitpref_t='_'
    hdr = red_readhead(states[0].filename)
    dt = strtrim(fxpar(hdr, 'DATE-AVG'), 2)
    avg_ts = (strsplit(dt, 'T', /extract))[1]
    avg_time = red_time2double(avg_ts)
    if self -> polarimetric_data() then begin
      pfls = file_search(self.out_dir + '/prefilter_fits/*-T_' + uprefs + $
                         '_[0-9][0-9]:[0-9][0-9]:[0-9][0-9]*save', count=Npfls)
    endif else begin
      pfls = file_search(self.out_dir + '/prefilter_fits/chromis_' + uprefs + $
                         '_[0-9][0-9]:[0-9][0-9]:[0-9][0-9]*save', count=Npfls)
                                ;    if Npfls gt 0 then begin
    endelse
    if Npfls eq Nprefilters then begin
      tt = dblarr(Npfls)
      ts = strarr(Npfls)
      for ii=0,Npfls-1 do begin
        ts[ii] = (strsplit(file_basename(pfls[ii]),'_',/extract))[2]
        tt[ii] = abs(red_time2double(ts[ii]) - avg_time)
      endfor                    ; ii
      mn = min(tt,jj)
      fitpref_t = '_'+ts[jj]+'_'
    endif
  endif else fitpref_t = '_'+fitpref_time+'_' 
  
                                ;: Read the prefilter fit output
  for ipref = 0L, Nprefilters-1 do begin

    idxpref = where(prefilters eq uprefs[ipref], count)

    if count eq 0 then continue
    
    if self -> polarimetric_data() then begin
      
      ;; T camera

      pfile = self.out_dir + '/prefilter_fits/' + instrument + '-T_' $
              + uprefs[ipref] + fitpref_t + 'prefilter.idlsave'
      if ~file_test(pfile) then begin
        red_message, 'Prefilter file not found: '+pfile
        return
      endif
      restore, pfile            ; Restores variable prf which is a struct

      if ipref eq 0 then begin
        units = prf.units
      endif else begin
        if units ne prf.units then begin
          red_message, ['Units in ' + pfile + ' do not match those in earlier read files.' $
                        , 'Please rerun the prefilterfit step for these data.']
          retall
        endif
      endelse

      ;; [m] Shift the wavelengths by this amount
      if n_elements(prf.fitpars) gt 1 then wave_shifts[idxpref] = prf.fitpars[1]/10. else wave_shift = 0.
      if n_elements(prf.wav) eq 1 then begin
        t_prefilter_curve[idxpref] = prf.pref
      endif else begin
        t_prefilter_curve[idxpref] = red_intepf(prf.wav, prf.pref, states[idxpref].tun_wavelength*1.d10)
      endelse
                                ;  nbt_prefilter_curve = prf.pref
      t_prefilter_wav[idxpref] = prf.wav
      t_prefilter_wb[idxpref] = prf.wbint
      
      ;; R camera

      pfile = self.out_dir + '/prefilter_fits/'+instrument + '-R_' $
              + uprefs[ipref] + fitpref_t + 'prefilter.idlsave'
      if ~file_test(pfile) then begin
        red_message, 'Prefilter file not found: '+pfile
        return
      endif
      restore, pfile            ; Restores variable prf which is a struct

      if units ne prf.units then begin
        red_message, ['Units in ' + pfile + ' do not match those in earlier read files.' $
                      , 'Please rerun the prefilterfit step for these data.']
        retall
      endif
      
      
      if n_elements(prf.wav) eq 1 then begin
        r_prefilter_curve[idxpref] = prf.pref
      endif else begin
        r_prefilter_curve[idxpref] = red_intepf(prf.wav, prf.pref, states[idxpref].tun_wavelength*1.d10)
      endelse
                                ;  nbr_prefilter_curve = prf.pref
      r_prefilter_wav[idxpref] = prf.wav
      r_prefilter_wb[idxpref] = prf.wbint

      ;;  if total(this_r_prefilter_wav ne this_t_prefilter_wav) gt 0 then stop
      ;; this_prefilter_wav = this_r_prefilter_wav

      ;;   this_prefilter_curve = (this_t_prefilter_curve + this_r_prefilter_curve) / 2.
      ;;   this_prefilter_wb    = (this_t_prefilter_wb    + this_r_prefilter_wb   ) / 2.

    endif else begin

      ;; N camera (CHROMIS)
      
      ;;pfile = self.out_dir + '/prefilter_fits/chromis_' + unbprefs[inbpref] + fitpref_t + 'prefilter.idlsave'
      pfile = self.out_dir + '/prefilter_fits/chromis_' + uprefs[ipref] + fitpref_t + 'prefilter.idlsave'
      if ~file_test(pfile) then begin
        print, inam + ' : prefilter file not found: '+pfile
        return
      endif
      
      restore, pfile            ; Restores variable prf which is a struct

      if ipref eq 0 then begin
        units = prf.units
      endif else begin
        if units ne prf.units then begin
          red_message, ['Units in ' + pfile + ' do not match those in earlier read files.' $
                        , 'Please rerun the prefilterfit step for these data.']
          retall
        endif
      endelse

      if n_elements(prf.fitpars) gt 1 then begin
        wave_shifts[idxpref] = prf.fitpars[1]/10. ; [nm] Shift the wavelengths by this amount
      endif else begin
        wave_shifts[idxpref] = 0.0
      endelse
    
;      print, 'Wave shifts: ', inbpref, ' ', prefilters[ipref], wave_shift
;      wave_shifts[idxpref] = wave_shift
      
      case !true of

        n_elements(prefilter_curve) eq 1 : prefilter_curve[idxpref] = prf.pref

        array_equal(prf.wav , states[idxpref].tun_wavelength*1.d10) : prefilter_curve[idxpref] = prf.pref

        else :  begin
          prefilter_curve[idxpref] = red_intepf(prf.wav, prf.pref, states[idxpref].tun_wavelength*1.d10)
        endelse

      endcase
      
;      if n_elements(prefilter_curve) eq 1 then begin
;        prefilter_curve[idxpref] = prf.pref
;      endif else begin
;        ;;prefilter_curve[idxpref] = red_intepf(prf.wav, prf.pref, wav[idxpref]*1.d10)
;        prefilter_curve[idxpref] = red_intepf(prf.wav, prf.pref, states[idxpref].tun_wavelength*1.d10)
;      endelse
      
    endelse 

    ;; Both polarimetric and non-polarimetric data:
    
;    red_append, prefilter_curve, this_prefilter_curve
;
;    if n_elements(prefilter_curve) eq 1 then begin
;      red_append, prefilter_wav, prf.wav
;      red_append, prefilter_wb, prf.wbint
;    endif else begin
;      red_append, prefilter_wav, wav[idxpref]*1.d10
;      red_append, prefilter_wb, replicate(prf.wbint, count)
;    endelse
    
  endfor                        ; ipref

  if self -> polarimetric_data() then begin
    prefilter_curve = (t_prefilter_curve + r_prefilter_curve) / 2.
    prefilter_wb    = (t_prefilter_wb    + r_prefilter_wb   ) / 2.
  endif
  
end
