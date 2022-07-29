; docformat = 'rst'

;+
; Class CRISP2
;
; :Author:
; 
;    Tomas Hillberg
; 
; :Keywords:
; 
;    develop : in, optional, type=boolean
; 
;       Run in developer mode.
;
;    no_db : in, optional, type=boolean
;
;       Do not use metadata database.
;
;
; :History:
;
;   2016-04-29 : THI. Split class RED into a base-class (instrument
;                independent parts) and derived classes (CRISP/CHROMIS). 
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding
;                SOLARNET keywords.
;
;   2017-03-09 : MGL. New keyword "develop".
;
;   2017-03-18 : MGL. For now, allow crispred in this branch only in
;                develop mode.
;
;   2017-06-28 : MGL. Removed cam{wb,r,t}-based states.
;
;   2019-04-11 : MGL. Develop mode not needed for crispred anymore.
;
;   2019-05-07 : MGL. Detect directory set up for old pipeline.
;
;   2021-03-03 : MGL. New keyword no_db.
;
;   2022-07-29 : MGL. Define the crisp2red class.
;
;-
pro crisp2::initialize, filename, develop = develop, no_db = no_db

  ;; Call initialize of the base-class first to load common parameters
  self->RED::initialize, filename, develop = develop, no_db = no_db

  spawn, 'grep "camera.*=" '+filename, grepresult1
  if n_elements(grepresult1) eq 1 && grepresult1 eq '' then begin
    spawn, 'grep "cam_wb.*=" '+filename, grepresult2
    print
    if n_elements(grepresult2) eq 1 && grepresult2 eq '' then begin
      print, 'No cameras defined in '+filename
    endif else begin
      print, 'This directory seems to be set up for running CRISP data with the old pipeline.'
      print, 'Please see https://dubshen.astro.su.se/wiki/index.php/Old_CRISPRED for instructions.'
    endelse
    print, 'Returning'
    print
    retall
  endif
  
  return
  
  ;; Then load CRISP2 specific stuff

  self.docamt = 1B
  self.docamr = 1B
  self.docamwb = 1B

  ;; open file and get fields
  openr, lun, filename, /get_lun
  nl = 0L
  while(~eof(lun)) do begin
    line = ''
    readf, lun, line                              ; read line
    field = (strsplit(line,' =',/extract))[0]     ; extract field

    if(strmid(line, 0, 1) eq '#') then continue

    ;; get fields
    case field of
      'cam_t': begin
        tmp = strtrim((strsplit(line, '=', /extract))[1],2)
        IF(strpos(tmp,"'") NE -1) THEN dum = execute('tmp = '+tmp)
        self.camt =  tmp
        IF ptr_valid(self.cameras) THEN red_append, *self.cameras, tmp $
        ELSE self.cameras = ptr_new(tmp, /NO_COPY)
      end
      'cam_r': begin
        tmp = strtrim((strsplit(line, '=', /extract))[1],2)
        IF(strpos(tmp,"'") NE -1) THEN dum = execute('tmp = '+tmp)
        self.camr =  tmp
        IF ptr_valid(self.cameras) THEN red_append, *self.cameras, tmp $
        ELSE self.cameras = ptr_new(tmp, /NO_COPY)
      end
      'cam_wb': begin
        tmp = strtrim((strsplit(line, '=', /extract))[1],2)
        IF(strpos(tmp,"'") NE -1) THEN dum = execute('tmp = '+tmp)
        self.camwb =  tmp
        IF ptr_valid(self.cameras) THEN red_append, *self.cameras, tmp $
        ELSE self.cameras = ptr_new(tmp, /NO_COPY)
      end
      else: begin
      end
    endcase
    nl+= 1
  endwhile
  free_lun, lun

  if(self.camt eq '') then begin
    print, 'red::initialize : WARNING : cam_t is undefined!'
    self.docamt = 0B
  endif
  if(self.camr eq '') then begin
    print, 'red::initialize : WARNING : cam_r is undefined!'
    self.docamr = 0B
  endif
  if(self.camwb eq '') then begin
    print, 'red::initialize : WARNING : cam_wb is undefined!'
    self.docamwb = 0B
  endif
  
  if(self.docamt) then print, 'red::initialize : cam_T = '+ self.camt
  if(self.docamr) then print, 'red::initialize : cam_R = '+ self.camr
  if(self.docamwb) then print, 'red::initialize : cam_wb = '+ self.camwb
end
