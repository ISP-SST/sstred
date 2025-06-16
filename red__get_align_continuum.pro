; docformat = 'rst'

;+
; Get the shifts from the CHROMIS Ca II continuum alignment.
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
;    tun_wavelengths : in, type=fltarr
; 
;      Tuning wavelegths, sorted.
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;   2025-06-06 : MGL. First version.
; 
;-
function red::get_align_continuum, wbgfiles, tun_wavelengths, direction

;;if n_elements(prefilter) eq 0 then prefilter = '3950'

  self -> extractstates, wbgfiles, wbgstates
  wbfilter = wbgstates[0].prefilter
  scannumbers = wbgstates.scannumber

  Ntunings = n_elements(tun_wavelengths)
  Nscans = n_elements(scannumbers)

  ashifts = dblarr(2, Ntunings, Nscans)

  ;; Continuum alignment only done for Ca II scans (so far). H beta is
  ;; not as wide so should be OK.
  if wbfilter ne '3950' then return, ashifts
  ;; Return zero shifts for other wavelength bands.

  ;; Get wavelength-variable shifts based on continuum vs wideband
  ;; alignment.

  ff = strsplit(wbgfiles,'/',/extract)
  dd = (ff.toarray())[*,1]
  timestamps = red_uniquify(dd)

  aligndirs = self.out_dir + '/align/' + timestamps $
              + '/' + WBfilter + '/'

  for jj=0, n_elements(aligndirs)-1 do begin
    
    nname = aligndirs[jj]+'scan_numbers.fz'
    sname = aligndirs[jj]+'continuum_shifts_smoothed.fz'
    
    if ~file_test(nname) or ~file_test(sname) then begin
      print, inam + ' : At least one file missing for aligncont option:'
      print, nname
      print, sname
      retall
    endif
    
    ;; Read the shifts for the continuum images
    fzread, scans, nname
    fzread, shifts, sname, align_header
    red_append, align_scannumbers, scans
    red_append, x_shifts, reform(shifts[0,*])
    red_append, y_shifts, reform(shifts[1,*])
  endfor                        ;aligndirs
  
  ;; Get the wavelengths used for the intra-scan alignment from the
  ;; file header.
  if n_elements(align_header) gt 0 then align_wavelengths = double(strsplit(align_header,/extract))
  if n_elements(align_wavelengths) ne 2 then begin
    ;; Default to WB and NB continuum as used by the align_continuum
    ;; method.
    align_wavelengths = [3.950e-07, 4.000e-07]
  endif
  
  ;; Check that we have alignment for all scan numbers
  match2, scannumbers, align_scannumbers, suba, subb
  missing_indx = where(suba eq -1, Nmissing)
  if Nmissing gt 0 then begin
    print, inam+' : Alignment missing for these scan numbers:'
    print, uscans[missing_indx]
    print, inam+' : Please rerun a -> align_continuum'
    retall
  endif
  
  ;; Select align shifts for the relevant scan numbers.
  nb_shifts = fltarr(2, Nscans)
  ;; The align_shifts are measured with direction=0 so we need to
  ;; take direction into account when interpreting the shifts.
  case direction of
    0 : begin                   ; ( x, y)
      nb_shifts[0, *] =  x_shifts[suba]
      nb_shifts[1, *] =  y_shifts[suba]
    end
    1 : begin                   ; (-y, x)
      nb_shifts[0, *] = -y_shifts[suba]
      nb_shifts[1, *] =  x_shifts[suba]
    end
    2 : begin                   ; (-x,-y)
      nb_shifts[0, *] = -x_shifts[suba]
      nb_shifts[1, *] = -y_shifts[suba]
    end
    3 : begin                   ; ( y,-x)
      nb_shifts[0, *] =  y_shifts[suba]
      nb_shifts[1, *] = -x_shifts[suba]
    end
    4 : begin                   ; ( y, x)
      nb_shifts[0, *] = y_shifts[suba]
      nb_shifts[1, *] = x_shifts[suba]
    end
    5 : begin                   ; (-x, y)
      nb_shifts[0, *] = -x_shifts[suba]
      nb_shifts[1, *] =  y_shifts[suba]
    end
    6 : begin                   ; (-y,-x)
      nb_shifts[0, *] = -y_shifts[suba]
      nb_shifts[1, *] = -x_shifts[suba]
    end
    7 : begin                   ; ( x,-y)
      nb_shifts[0, *] =  x_shifts[suba]
      nb_shifts[1, *] = -y_shifts[suba]
    end
    else : stop
  endcase
  
;    ;; Use interpolation to get the shifts for the selected scans.
;    nb_shifts = fltarr(2, Nscans)
;    for iscan=0L, Nscans-1 do begin
;      pos = where(align_scannumbers eq uscans[iscan], cccc)
;      if cccc eq 1 then nb_shifts[*, iscan] = align_shifts[*, pos] else begin
;        nb_shifts[0, *] = interpol([reform(align_shifts[0, *])] $
;                                   , [float(align_scannumbers)], [float(uscans)])
;        nb_shifts[1, *] = interpol([reform(align_shifts[1, *])] $
;                                   , [float(align_scannumbers)], [float(uscans)])
;      endelse
;    endfor
  pos = where(~finite(nb_shifts), cccc)
  if cccc gt 0 then nb_shifts[pos] = 0
  

  for iscan = 0L, Nscans-1 do begin

;    ts = (strsplit(wbgfiles[iscan],'/',/extract))[1]
    
;    ;; The NB files in this scan, sorted in tuning wavelength order.
;    self -> selectfiles,  files = pertuningfiles, states = pertuningstates $ 
;                                  , cam = nbcamera, scan = uscans[iscan], timestamps = ts $
;                                  , sel = scan_nbindx, count = count
;    scan_nbfiles = pertuningfiles[scan_nbindx]
;    scan_nbstates = pertuningstates[scan_nbindx]
;    sortindx = sort(tun_wavelengths)
;    scan_nbfiles = scan_nbfiles[sortindx]
;    scan_nbstates = scan_nbstates[sortindx]
;
;    icont = where(scan_nbstates.wbfilter eq '3999')

    ;; Interpolate shifts in X and Y
    ashifts[0, *, iscan] = interpol([0., nb_shifts[0, iscan]] $
                                    , align_wavelengths*1e7 $
                                    , tun_wavelengths*1e7)
    ashifts[1, *, iscan] = interpol([0., nb_shifts[1, iscan]] $
                                    , align_wavelengths*1e7 $
                                    , tun_wavelengths*1e7)

    
  endfor                        ; iscan

  return, ashifts

end
