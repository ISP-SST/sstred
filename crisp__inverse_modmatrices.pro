; docformat = 'rst'

;+
; Read inverse modulation matrix maps from disk, make and store them
; if necessary.
; 
; This method does (or at least is meant to do) the same as
; red__polarim.pro in the master (old CRISP) branch, except that it
; stores the inverse modulation matrices to disk rather than pointing
; to them from within the pol class object. (We don't use a pol class
; in the new code base.)
;
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
; 
; :Params:
; 
; 
; 
; 
; 
; 
; :Keywords:
; 
; :History:
; 
; 
;   2018-10-09 : MGL. First version, based on Jaime's red::polarim.
; 
;   2018-12-21 : MGL. Removed keyword smooth.
; 
;   2019-05-29: MGL. Use the cache to avoid reading the inverse
;               modulation matrices over and over.
; 
;-
pro crisp::inverse_modmatrices, prefilter, dir $
                                , camr = camr $
                                , camt = camt $
                                , immr = immr $
                                , immt = immt $
                                , no_ccdtabs = no_ccdtabs $
                                , overwrite = overwrite
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  if n_elements(camr) gt 0 then begin
    red_append, cams, camr
    red_append, prefixes, 'r'
  endif
  
  if n_elements(camt) gt 0 then begin
    red_append, cams, camt
    red_append, prefixes, 't'
  endif

  Ncams = n_elements(cams)

  for icam = 0, Ncams-1 do begin

    detector = self->getdetector( cams[icam] )

    fname = dir + '/imm'+prefixes[icam]+'.fits'

    if file_test(fname) and ~keyword_set(overwrite) then begin

      ;; OK, we should use data from this file. But do we need to read
      ;; it or have we read it already and put it in the cache?

      this_cache = rdx_cacheget(fname, count = cnt)
      if cnt eq 0 then begin
        ;; Get prpara info from hdr and check if they match current
        ;; parameters?
        imm = red_readdata(fname, head = hdr)
        ;; Put imm in the cache so we don't have to read it again.
        rdx_cache, fname, { data:imm }
      endif else begin
        imm = this_cache.data
      endelse
      
    endif else begin

      ;; Polcal output
      pname = self.out_dir+'/polcal/'+detector+'_'+prefilter+'_polcal.fits'
      
      if ~file_test(pname) then begin
        print, inam + ' : ERROR, polcal data not found in ' + self.out_dir + '/polcal/'
        return
      endif
      
      mm = red_readdata(pname)  ; Modulation matrix
      dims = size(mm, /dim)
      Nelements = dims[0]
      Nx = dims[1]
      Ny = dims[2]
      
      ;; Interpolate Sarnoff CCD tabs?
      if strmatch((red_camerainfo(detector)).model,'Sarnoff*') then begin
        if ~keyword_set(no_ccdtabs) then begin
          for ii = 0, Nelements-1 do mm[ii,*,*] = red_mask_ccd_tabs(reform(mm[ii,*,*]))
        endif
      endif

      ;; Mask NaNs
      mask = bytarr(Nx, Ny) + 1B
      for ii = 0, Nelements-1 do mask and= finite(reform(mm[ii,*,*])) 
      if min(mask) lt 1 then for ii = 0, Nelements-1 do mm[ii,*,*] = mask * reform(mm[ii,*,*])

;       ;; Check NaNs
;      for ii = 0, Nelements-1 do begin
;        mask = finite(reform(mm[ii,*,*]))
;        idx = where(mask, count)
;        if count gt 0 and count lt Nx*Ny then $
;           mm[ii,*,*] = red_fillpix(reform(mm[ii,*,*]), mask=mask)
;      endfor                    ; ii

      imm = red_invert_mmatrix(temporary(mm)) ; Inverse modulation matrix

      ;; Mask outlier pixels (but keep any NaN-masking)
      mask and= max(abs(imm),dim=1) le 2

      if min(mask) lt 1 then begin
        ;; Zero before pixe-filling. Should not be needed but
        ;; apparently is.
        for ii = 0, Nelements-1 do imm[ii,*,*] = mask * reform(imm[ii,*,*])
        for ii = 0, Nelements-1 do imm[ii,*,*] = red_fillpix(reform(imm[ii,*,*]), mask=mask)
      endif
      
      ;; Store some prpara info in the header? 
      red_writedata, fname, imm
      
    endelse

    case prefixes[icam] of
      'r'  : immr = imm
      't'  : immt = imm
      else : stop               ; Should not happen!
    endcase
    
  endfor                        ; icam

end
