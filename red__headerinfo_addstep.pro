; docformat = 'rst'

;+
; Add FITS info about a processing step to FITS header.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Params:
; 
;    header : in, type=strarr
;
;      The FITS header to be modified.
; 
; :Keywords:
; 
;    prstep : in, optional, type=string
;
;       The processing step name.
;
;    prproc : in, optional, type=string
;
;       The name of the procedure used.
;
;    prmode : in, optional, type=string
;
;       The processing mode .
;
;    prpara : in, optional, type=string
;
;       List of parameters/options for the procedure.
;
;    addlib : in, optional, type=string
;
;       Library in addition to the default ones.
;
;    level  : in, optional, type=string
;
;       
;
;    version  : in, optional, type=string
;
;       
;
;   
;   
;   
; 
; 
; :History:
;
;   2017-03-16 : MGL. First version.
;
;   2017-03-16 : MGL. Use red_headerinfo_addlib.
; 
; 
;-
pro red::headerinfo_addstep, header $
                             , prstep = prstep $
                             , prproc = prproc $
                             , prmode = prmode $
                             , prpara = prpara $
                             , addlib = addlib $
                             , level = level $
                             , version = version
  
  if n_elements(header) eq 0 then mkhdr, header, 0
  
  prevkey = 'OBS_SHDU';'SOLARNET'

  ;; Look for existing processing steps, set stepnumber to one higher.
  stepnumber = 0
  repeat begin
    stepnumber++
    stp = strtrim(stepnumber, 2)
    tmp = fxpar(header, 'PRSTEP'+stp, count = count)
  endrep until count eq 0

;  if stepnumber gt 1 then prevkey = 

  ;; Add the LONGSTRN keyword, just in case. (Doing it here lets us
  ;; position it a bit better.)
;  fxaddpar_contwarn, header, 'SOLARNET'

  if n_elements(prstep) eq 0 then prstep = 'Unknown'
  fxaddpar, header, 'PRSTEP'+stp, prstep, ' Processing step name', after = prevkey
  prevkey = 'PRSTEP'+stp

  ;; Procedure name
  if n_elements(prproc) ne 0 then begin
    key = 'PRPROC'+stp
    fxaddpar, header, key, prproc, ' Name of procedure used', after = prevkey
    prevkey = key
  endif

  ;; Add headers with library names and versions. (Bug: Should be
  ;; listed in the order they appear in the path!)
  red_headerinfo_addlib, head, 'SSTRED', self.version_pipeline, prevkey = prevkey
  red_headerinfo_addlib, head, 'IDLAstro', self.version_idlastro , prevkey = prevkey
  red_headerinfo_addlib, head, 'Coyote', self.version_coyote, prevkey = prevkey
  red_headerinfo_addlib, head, 'mpfit', self.version_mpfit, prevkey = prevkey
  red_headerinfo_addlib, head, 'reduxdlm', self.version_reduxdlm, prevkey = prevkey
  for ilib = 0, n_elements(addlib)-1 do red_headerinfo_addlib, head, addlib[ilib], prevkey = prevkey
  
  ;; Procedure mode
  if self.developer_mode then red_append, prmode, 'Developer mode' 
  if n_elements(prmode) gt 0 then begin
    prm = strjoin(prmode, ',')
    key = 'PRMODE'+stp
    if prm ne '' then begin
      fxaddpar,header, key, prm, ' Processing mode', after = prevkey
      prekey = key
    end
  endif
  
  ;; Procedure parameters
  if n_elements(prpara) ne 0 then begin
    key = 'PRPARA'+stp
    case typename(prpara) of
      'STRING' : prp = strjoin(strtrim(prpara, 2), ',')
      'DICTIONARY' : begin
        keys = prpara.keys()
        for ipara = 0, n_elements(prpara)-1 do begin
          value = prpara[keys[ipara]]
          if n_elements(value) eq 1 then $
             value = strtrim(value, 2) else value = '['+strjoin(strtrim(value, 2), ',')+']'
          red_append, prp, strtrim(keys[ipara], 2) + '=' + value
        endfor                  ; ipara
        prp = strjoin(strtrim(prp, 2), ',')
      end
      else : stop
    endcase
    fxaddpar, header, key, prp $, after = prevkey $
              , ' List of parameters/options for PRPROC'+stp
    prevkey = key
  endif

  ;; Add headers with other step info.
;  fxaddpar, header, 'VERSION', 0, 'FITS file processing generation/version'
;  fxaddpar, header, 'LEVEL',   0, 'Data level of fits file'

end

