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
  
  ;; Look for existing processing steps, set stepnumber to one higher.
  stepnumber = 0
  repeat begin
    stepnumber++
    stp = strtrim(stepnumber, 2)
    tmp = sxpar(header, 'PRSTEP'+stp, count = count)
  endrep until count eq 0

  ;; Add the LONGSTRN keyword, just in case. (Doing it here lets us
  ;; position it a bit better.)
;  fxaddpar_contwarn, header, 'SOLARNET'

  ;; Add headers with library names and versions. (Bug: Should be
  ;; listed in the order they appear in the path!)
  prlibcomment = ' Software library'
  if n_elements(prproc) ne 0 then prlibcomment += ' containing '+prproc
  fxaddpar, header, 'PRLIB'+stp, after = 'OBS_SHDU' $
            , 'SSTRED', prlibcomment
  fxaddpar, header, 'PRVER'+stp, after = 'PRLIB'+stp $
            , self.version_pipeline, ' Library version/MJD of last update' 
  fxaddpar, header, 'PRLIB'+stp+'A', after = 'PRVER'+stp $
            , 'IDLAstro', ' Additional software library'
  fxaddpar, header, 'PRVER'+stp+'A', after = 'PRLIB'+stp+'A' $
            , self.version_idlastro, ' Library version/MJD of last update'
  fxaddpar, header, 'PRLIB'+stp+'B', after = 'PRVER'+stp+'A' $
            , 'Coyote', ' Additional software library'
  fxaddpar, header, 'PRVER'+stp+'B', after = 'PRLIB'+stp+'B' $
            , self.version_coyote, ' Library version/MJD of last update'
  fxaddpar, header, 'PRLIB'+stp+'C', after = 'PRVER'+stp+'B' $
            , 'mpfit', ' Additional software library'
  fxaddpar, header, 'PRVER'+stp+'C', after = 'PRLIB'+stp+'C' $
            , self.version_mpfit, ' Library version/MJD of last update'
  fxaddpar, header, 'PRLIB'+stp+'D', after = 'PRVER'+stp+'C' $
            , 'reduxdlm', ' Additional software library'
  fxaddpar, header, 'PRVER'+stp+'D', after = 'PRLIB'+stp+'D' $
            , self.version_reduxdlm, ' Library version/MJD of last update'

  for ilib = 0, n_elements(addlib)-1 do begin
    libletter = string(byte(69+ilib)) ; 'D' = 68!
    fxaddpar, header, 'PRLIB'+stp+libletter, after = 'PRVER'+stp+string(byte(68+ilib)) $
              , addlib[ilib], ' Additional software library'
  endfor                        ; ilib
  
  if n_elements(prstep) eq 0 then prstep = 'Unknown'
  fxaddpar, header, 'PRSTEP'+stp, prstep, ' Processing step name', before = 'PRLIB'+stp

;  FXADDPAR_CONTWARN, HEADER, 'PRSTEP'+stp
  
  if self.developer_mode then red_append, prmode, 'Developer mode' 
  if n_elements(prmode) gt 0 then begin
    prm = strjoin(prmode, ',')
    if prm ne '' then $
       fxaddpar,header,'PRMODE'+stp, prm, ' Processing mode', before = 'PRLIB'+stp
  endif
  
  if n_elements(prproc) ne 0 then $
     fxaddpar, header, 'PRPROC'+stp, prproc, ' Name of procedure used', before = 'PRLIB'+stp
  
  if n_elements(prpara) ne 0 then begin
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
  endif else prp = ''
  fxaddpar, header, 'PRPARA'+stp, prp $, before = 'PRLIB'+stp $
            , ' List of parameters/options for PRPROC'+stp


  ;; Add headers with other step info.
;  fxaddpar, header, 'VERSION', 0, 'FITS file processing generation/version'
;  fxaddpar, header, 'LEVEL',   0, 'Data level of fits file'

end

