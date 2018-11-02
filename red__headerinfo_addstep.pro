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
;   2017-12-06 : MGL. Remove trailing blank lines.
;
;   2017-12-20 : MGL. Write prpara dictionaries in json format.
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

  ;; Existing steps
  prsteps_existing = fxpar(header,'PRSTEP*')
  Nexisting = n_elements(prsteps_existing)
  if Nexisting gt 0 then begin
    anchor = 'PRBRA'+strtrim(Nexisting, 2) ; This should be the last existing
  endif else begin
    anchor = 'OBS_HDU'          ;'SOLARNET'
  endelse
  
  ;; Look for existing processing steps, set stepnumber to one higher.
  stepnumber = 0
  repeat begin
    stepnumber++
    stp = strtrim(stepnumber, 2)
    tmp = fxpar(header, 'PRSTEP'+stp, count = count)
  endrep until count eq 0

  if n_elements(prstep) eq 0 then prstep = 'Unknown'
  red_fitsaddkeyword, header, 'PRSTEP'+stp, prstep, 'Processing step name', anchor = anchor

  ;; Procedure name
  if n_elements(prproc) ne 0 then begin
    key = 'PRPROC'+stp
    red_fitsaddkeyword, header, key, prproc, 'Name of procedure used', anchor = anchor
  endif

  ;; Add headers with library names and versions. (Bug: Should be
  ;; listed in the order they appear in the path!)
  red_headerinfo_addlib, header, 'SSTRED', self.version_pipeline, prbranch = self.version_pipeline_branch
  red_headerinfo_addlib, header, 'IDLAstro', self.version_idlastro
  red_headerinfo_addlib, header, 'Coyote', self.version_coyote
  red_headerinfo_addlib, header, 'mpfit', self.version_mpfit
  red_headerinfo_addlib, header, 'reduxdlm', self.version_reduxdlm
  for ilib = 0, n_elements(addlib)-1 do red_headerinfo_addlib, header, addlib[ilib]

  ;; Procedure mode
  if self.developer_mode then red_append, prmode, 'Developer mode' 
  if n_elements(prmode) gt 0 then begin
    prm = strjoin(prmode, ',')
    key = 'PRMODE'+stp
    if prm ne '' then begin
      red_fitsaddkeyword,header, key, prm, 'Processing mode', anchor = anchor
    end
  endif

  ;; Procedure parameters
  if n_elements(prpara) ne 0 then begin
    key = 'PRPARA'+stp
    case typename(prpara) of
      'STRING' : prp = strjoin(strtrim(prpara, 2), ',')
      'DICTIONARY' : prp = json_serialize(prpara)
      else : stop
    endcase
    red_fitsaddkeyword, header, key, prp, anchor = anchor $
                        , 'List of parameters/options for PRPROC'+stp
  endif

  ;; Add headers with other step info.
;  red_fitsaddkeyword, header, 'VERSION', 0, 'FITS file processing generation/version'
;  red_fitsaddkeyword, header, 'LEVEL',   0, 'Data level of fits file'
  
  ;; Remove trailing blank lines
  Nlines = where(strmatch(header, 'END *'), Nmatch)
  if Nmatch eq 0 then stop
  header = header[0:Nlines]

end

;; Mail with list of approved PRSTEP labels from 2018-02.16, Subject:
;; Re: [solarnet-20.3] [solarnet-50.2] Do you sum, add or stack
;; images?

;; Below is our summary of the discussion on what to put in the
;; PRSTEPn keywords. Remember that the *order* of the processing steps
;; is specified by the n. In other words: 

;; PRSTEP1 = ‘BINNING’
;; PRSTEP2 = ‘FIXED-PATTERN-REMOVAL’

;; is different from 

;; PRSTEP1 = ‘FIXED-PATTERN-REMOVAL’
;; PRSTEP2 = ‘BINNING’

;; This makes in unnecessary to have two types of BINNING (on-chip or
;; post-readout).

;; - There was some discussion about the following descriptions: are
;;   they accurate enough? 

;; BINNING               (This is ok, given the explanation above)
;; CALIBRATION           (When necessary, XXX-CALIBRATION or
;;                        CALIBRATION-XXX may be used, with no further
;;                        rules on what XXX might be) 
;; ALIGNMENT             (When necessary, XXX-ALIGNMENT or
;;                        ALIGNMENT-XXX, as for CALIBRATION)  
;; PIXEL-FILLING         (Instead of many other, awkward suggestions) 
;; DISTORTION-CORRECTION (When necessary, XXX-DISTORTION-CORRECTION or
;;                        DISTORTION-CORRECTION-XXX, as for
;;                        CALIBRATION)  


;; - New suggestions for processing steps, any objections?: 

;; SUBTRACTION
;; MULTIPLICATION
;; FILTERING
;; EDGE-DETECTION
;; THRESHOLDING
;; BINARIZATION
;; INTERFEROMETRY-SPECKLING
;; DECONVOLUTION-SPECKLING
;; DEMODULATION
;; DEROTATION  (of rotating FOVs in movies)
;; INVERTING
;; DESTRETCHING
;; LINE-FITTING


;; - For some processing steps we had several suggestions. The
;;   majority of the respondents prefer the following alternatives:  

;; SUMMING               (rather than STACKING or (CO-)ADDING)
;; BIAS-CORRECTION       (rather than BIAS-SUBTRACTION)
;; DARK-SUBTRACTION      (rather than DARK-CORRECTION)
;; FLATFIELDING          (rather than FLAT-CORRECTION or FLAT-DIVISION)

          
;; - There have been no objections to the following descriptions:

;; FIXED-PATTERN-REMOVAL 
;; MOMFBD 
;; DESPIKING   


;; -- 
;; Stein & Terje


