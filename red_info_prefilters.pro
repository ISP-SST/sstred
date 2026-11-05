; docformat = 'rst'

;+
; Return info about (pre)filters of CRISP, CRISP2, and CHROMIS. 
; 
; Source: 
; https://dubshen.astro.su.se/wiki/index.php/CRISP
; https://dubshen.astro.su.se/wiki/index.php/CRISP2
; https://dubshen.astro.su.se/wiki/index.php/CHROMIS
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
;   pref : in, type="string or strarr"
; 
;     The prefilter name, the prefilter CWL in Ångstöm as a 4-digit
;     string. If an array, output keywords will also be arrays of the
;     same length.
; 
; 
; :Keywords:
; 
;   
;   CWL : out, optional, type=float
;   
;     The center wavelength of the tranmission passband in meters.
;     Note: the filter can be tilted so not necessarily the peak of
;     the actual transmission profile.
; 
;   date : in, optional, type=string, default=today
;   
;     The date in ISO yyyy-mm-dd format. Used in case multiple filters
;     have had the same name at different points in time.
; 
;   fwhm : out, optional, type=float
;   
;     The FWHM of the tranmission passband in meters. Note: the filter
;     can be tilted so not necessarily the width of the actual
;     transmission profile.
; 
;   manufacturer : out, optional, type=string
; 
;     The manufacturer of the filter, if known.
; 
;   nb : out, optional, type=boolean
; 
;     True if this filter is used as a NB prefilter.
; 
;   wb : out, optional, type=boolean
; 
;     True if this filter is used as a WB filter.
; 
;   waveband : out, optional, type=string
;   
;     A string describing the main spectral feature(s) in the
;     passband.
; 
;   wavemax : out, optional, type=float
;   
;     CWL + FWHM/2
; 
;   wavemin : out, optional, type=float
;   
;     CWL - FWHM/2
; 
; 
; :History:
; 
;   2026-10-01 : MGL. First version.
; 
;   2026-11-05 : MGL. New keywords nb and wb. Allow array of
;                prefilters.
; 
;-
pro red_info_prefilters, pref $
                         , CWL = CWL $
                         , date = date $
                         , fwhm = fwhm $
                         , manufacturer = manufacturer $
                         , nb = nb $
                         , wb = wb $
                         , waveband = waveband $
                         , wavemax = wavemax $
                         , wavemin = wavemin

  Nprefs = n_elements(pref)
    
  fwhm = fltarr(Nprefs)
  manufacturer = strarr(Nprefs)
  nb = bytarr(Nprefs)
  waveband = strarr(Nprefs)
;  wavemax = fltarr(Nprefs)
;  wavemin = fltarr(Nprefs)
  wb = bytarr(Nprefs)
  cwl = fltarr(Nprefs)

  for ipref = 0, Nprefs-1 do begin
  
    case strtrim(pref[ipref],2) of

      ;; Note: CHROMIS WB and NB filters are not the same but all Ca II
      ;; filters get the same waveband, as do both Hb filters.

      '3925' : begin            ; CHROMIS Ca II K blue wing
        waveband[ipref] = 'Ca II H & K'
        CWL[ipref] = 3925.2e-10
        fwhm[ipref] = 0.41e-9
        manufacturer[ipref] = 'Alluxa'
        nb[ipref] = 1
        wb[ipref] = 0
      end
      '3934' : begin            ; CHROMIS Ca II K core
        waveband[ipref] = 'Ca II H & K'
        CWL[ipref] = 3934.1e-10
        fwhm[ipref] = 0.41e-9
        manufacturer[ipref] = 'Alluxa'
        nb[ipref] = 1
        wb[ipref] = 0
      end
      '3950' : begin            ; CHROMIS Ca II HK wideband
        waveband[ipref] = 'Ca II H & K'
        CWL[ipref] = 3950e-10
        fwhm[ipref] = 1.32e-9
        manufacturer[ipref] = 'Alluxa'
        nb[ipref] = 0
        wb[ipref] = 1
      end
      '3969' : begin            ; CHROMIS Ca II H core
        waveband[ipref] = 'Ca II H & K'
        CWL[ipref] = 3968.8e-10
        fwhm[ipref] = 0.41e-9
        manufacturer[ipref] = 'Alluxa'
        nb[ipref] = 1
        wb[ipref] = 0
      end
      '3978' : begin            ; CHROMIS Ca II H red wing
        waveband[ipref] = 'Ca II H & K'
        CWL[ipref] = 3977.4e-10
        fwhm[ipref] = 0.41e-9
        manufacturer[ipref] = 'Alluxa'
        nb[ipref] = 1
        wb[ipref] = 0
      end
      '3999' : begin            ; CHROMIS Ca II H continuum
        waveband[ipref] = 'Ca II H & K'
        CWL[ipref] = 3999e-10
        fwhm[ipref] = 0.42e-9
        manufacturer[ipref] = 'Alluxa'
        nb[ipref] = 1
        wb[ipref] = 0
      end
      '4862' : begin            ; CHROMIS H-beta core
        waveband[ipref] = 'H-beta'
        CWL[ipref] = 4862.2e-10
        fwhm[ipref] = 0.48e-9
        manufacturer[ipref] = 'Alluxa'
        nb[ipref] = 1
        wb[ipref] = 0
      end
      '4846' : begin            ; CHROMIS H-beta wideband continuum
        waveband[ipref] = 'H-beta'
        CWL[ipref] = 4845.5e-10
        fwhm[ipref] = 0.65e-9
        manufacturer[ipref] = 'Alluxa'
        nb[ipref] = 0
        wb[ipref] = 1
      end

      ;; Below are CRISP and CRISP2 filters. Some passbands have older
      ;; filters manufactured by Barr and newer filters manufactured by
      ;; Alluxa. In the case they have exactly the same name, we check
      ;; the date to pick the version.

      ;; Note that all CRISP/CRISP2 filters are both NB prefilters and
      ;; WB filters.
      
      '5173' : begin            ; CRISP
        waveband[ipref] = 'Mg b ' + strtrim(pref[ipref],2)
        nb[ipref] = 1
        wb[ipref] = 1
        if n_elements(date) gt 0 && date lt red_dates(tag = 'CRISP Alluxa') then begin
          CWL[ipref] = 517.34e-9
          fwhm[ipref] = 0.30e-9
          manufacturer[ipref] = 'Barr'
        endif else begin
          CWL[ipref] = 517.32e-9
          fwhm[ipref] = 0.30e-9
          manufacturer[ipref] = 'Alluxa'
        endelse
      end
      '5250' : begin            ; CRISP
        waveband[ipref] = 'Fe I ' + strtrim(pref[ipref],2)
        CWL[ipref] = 525.0e-9
        fwhm[ipref] = 0.30e-9
        manufacturer[ipref] = 'Barr'
        nb[ipref] = 1
        wb[ipref] = 1
      end
      '5380' : begin            ; CRISP
        waveband[ipref] = 'C I ' + strtrim(pref[ipref],2)
        CWL[ipref] = 538.0e-9
        fwhm[ipref] = 0.33e-9
        manufacturer[ipref] = ''
        nb[ipref] = 1
        wb[ipref] = 1
      end
      '5382' : begin            ; CRISP
        waveband[ipref] = 'C I ' + strtrim(pref[ipref],2)
        CWL[ipref] = 538.2e-9
        fwhm[ipref] = 0.33e-9
        manufacturer[ipref] = 'Barr'
        nb[ipref] = 1
        wb[ipref] = 1
      end
      '5576' : begin            ; CRISP - Is this the same as the 5578 filter?
        waveband[ipref] = 'Fe I ' + strtrim(pref[ipref],2)
        CWL[ipref] = 557.6e-9
        fwhm[ipref] = 0.30e-9
        manufacturer[ipref] = 'Barr'
        nb[ipref] = 1
        wb[ipref] = 1
      end
      '5578' : begin            ; CRISP
        waveband[ipref] = 'Fe I ' + strtrim(pref[ipref],2)
        CWL[ipref] = 557.82e-9
        fwhm[ipref] = 0.30e-9
        manufacturer[ipref] = 'Barr'
        nb[ipref] = 1
        wb[ipref] = 1
      end
      '5876' : begin            ; CRISP
        waveband[ipref] = 'He I D3 ' + strtrim(pref[ipref],2)
        CWL[ipref] = 587.66e-9
        fwhm[ipref] = 0.54e-9
        manufacturer[ipref] = 'Alluxa'
        nb[ipref] = 1
        wb[ipref] = 1
      end
      '5896' : begin            ; CRISP (with new cameras)
        waveband[ipref] = 'Na D ' + strtrim(pref[ipref],2)
        CWL[ipref] = 589.7e-9
        fwhm[ipref] = 0.38e-9
        manufacturer[ipref] = 'Barr'
        nb[ipref] = 1
        wb[ipref] = 1
      end
      '5897' : begin            ; CRISP
        waveband[ipref] = 'Na D ' + strtrim(pref[ipref],2)
        CWL[ipref] = 589.7e-9
        fwhm[ipref] = 0.38e-9
        manufacturer[ipref] = 'Barr'
        nb[ipref] = 1
        wb[ipref] = 1
      end
      '6173' : begin            ; CRISP - Alluxa filter
        waveband[ipref] = 'Fe I ' + strtrim(pref[ipref],2)
        nb[ipref] = 1
        wb[ipref] = 1
                                ; This filter has a story... First the Barr filter with CRISP. A
                                ; new, larger Alluxa filter was installed with CRISP2 in early
                                ; August 2025. But this filter was defective so the old Barr
                                ; filter was reinstalled, probably 29 Aug. And then a new Alluxa
                                ; filter was installed before the 2026 season.
        if n_elements(date) gt 0 then begin
          case !true of
            date gt '2026-01-01': begin
              ;; One of three replacement filters installed before the
              ;; 2026 season.
              CWL[ipref] = 617.38e-9
              fwhm[ipref] = 0.45e-9
              manufacturer[ipref] = 'Alluxa'
            end
            date gt '2025-08-29': begin
              ;; The smaller CRISP Alluxa filter reinstalled because the
              ;; new one was defective..
              CWL[ipref] = 617.38e-9
              fwhm[ipref] = 0.45e-9
              manufacturer[ipref] = 'Alluxa'
            end
            date gt red_dates(tag = 'CRISP2 install') : begin
              ;; CRISP2 original, defective filter. CWL to be confirmed.
              CWL[ipref] = 617.38e-9
              fwhm[ipref] = 0.45e-9
              manufacturer[ipref] = 'Alluxa'
            end
            date gt red_dates(tag = 'CRISP Alluxa') : begin
              ;; The smaller Alluxa filter used with CRISP
              CWL[ipref] = 617.38e-9
              fwhm[ipref] = 0.45e-9
              manufacturer[ipref] = 'Alluxa'
            end
            else : begin
              ;; CRISP - probably identical to the filter labeled '6174'.
              CWL[ipref] = 617.4e-9
              fwhm[ipref] = 0.43e-9
              manufacturer[ipref] = 'Barr'
            end
          endcase
        endif
      end
      '6174' : begin            ; CRISP
        waveband[ipref] = 'Fe I ' + strtrim(pref[ipref],2)
        CWL[ipref] = 617.38e-9
        fwhm[ipref] = 0.45e-9
        manufacturer[ipref] = 'Barr'
        nb[ipref] = 1
        wb[ipref] = 1
      end
      '6302' : begin            ; CRISP
        waveband[ipref] = 'Fe I 6301+6302' 
        nb[ipref] = 1
        wb[ipref] = 1
        if n_elements(date) gt 0 && date lt red_dates(tag = 'CRISP Alluxa') then begin
          CWL[ipref] = 630.26e-9
          fwhm[ipref] = 0.45e-9
          manufacturer[ipref] = 'Barr'
        end else begin
          CWL[ipref] = 630.25e-9
          fwhm[ipref] = 0.45e-9
          manufacturer[ipref] = 'Alluxa'
        endelse
      end
      '6563' : begin            ; CRISP
        waveband[ipref] = 'H-alpha' 
        nb[ipref] = 1
        wb[ipref] = 1
        if n_elements(date) gt 0 && date lt red_dates(tag = 'CRISP Alluxa') then begin
          CWL[ipref] = 656.38e-9
          fwhm[ipref] = 0.49e-9
          manufacturer[ipref] = 'Barr'
        endif else begin
          CWL[ipref] = 656.31e-9
          fwhm[ipref] = 0.49e-9
          manufacturer[ipref] = 'Alluxa'
        endelse
      end
      '7699' : begin            ; CRISP
        waveband[ipref] = 'K I ' + strtrim(pref[ipref],2)
        CWL[ipref] = 769.9e-9
        fwhm[ipref] = 0.6e-9
        manufacturer[ipref] = 'Alluxa'
        nb[ipref] = 1
        wb[ipref] = 1
      end
      '7772' : begin            ; CRISP
        waveband[ipref] = 'O I ' + strtrim(pref[ipref],2)
        CWL[ipref] = 777.25e-9
        fwhm[ipref] = 0.7e-9
        manufacturer[ipref] = 'Barr'
        nb[ipref] = 1
        wb[ipref] = 1
      end
      '8542' : begin            ; CRISP
        waveband[ipref] = 'Ca II ' + strtrim(pref[ipref],2)
        nb[ipref] = 1
        wb[ipref] = 1
        if n_elements(date) gt 0 && date lt red_dates(tag = 'CRISP Alluxa') then begin
          CWL[ipref] = 854.1e-9
          fwhm[ipref] = 0.93e-9
          manufacturer[ipref] = 'Barr'
        endif else begin
          CWL[ipref] = 854.26e-9
          fwhm[ipref] = 0.80e-9
          manufacturer[ipref] = 'Alluxa'
        endelse
      end
      else: begin      
        CWL[ipref] = 0
        fwhm[ipref] = 0
        waveband[ipref] = ''
        manufacturer[ipref] = ''
        nb[ipref] = 1
        wb[ipref] = 1
      end
    endcase

    ;;  CWL = CWL*1e9
    wavemax = (CWL+fwhm/2.)     ;*1e9
    wavemin = (CWL-fwhm/2.)     ;*1e9
  endfor

  if isa(pref, /scalar) then begin
    ;; Output should be scalar
    cwl = cwl[0]
    fwhm = fwhm[0]
    manufacturer = manufacturer[0]
    nb = nb[0]
    waveband = waveband[0]
    wavemax = wavemax[0]
    wavemin = wavemin[0]
    wb = wb[0]
 endif
  
end

pref = '5876'
red_info_prefilters, pref $
                     , fwhm = fwhm $
                     , manufacturer = manufacturer $
                     , waveband = waveband $
                     , CWL = CWL $
                     , wavemax = wavemax $
                     , wavemin = wavemin

waveunit = -9
waveunitfac = 10d^double(waveunit)

print, 'Info about filter '+pref+' [nm] : '
print, 'WAVEBAND = ', waveband
print, 'FWHM = ', fwhm / waveunitfac
print, 'WAVELNTH = ', CWL / waveunitfac
print, 'WAVEMAX = ',  wavemax / waveunitfac
print, 'WAVEMIN = ',  wavemin / waveunitfac
print, 'WAVEUNIT = ', waveunit
print, 'manufacturer = ', manufacturer



end
