; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    filename : 
;   
;   
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-12-10 : PS  adapt for multiple flat_dir
;
;   2014-01-08 : MGL. Initialize new fields: isodate, log_dir, telog, pinhole_spacing.
; 
;   2014-01-10 : PS  New config variable filtype
;
;   2014-01-22 : MGL. Adapt to string functions moved to the str_
;                namespace.
;
;   2014-03-21 : MGL. Allow for multiple dark_dir.
;
;   2014-10-16 : MGL. Set up for Coyote graphics.
;
;   2016-02-15 : MGL. Make descatter_dir local.
;-
pro red::initialize, filename
                                
  ;; Test file
  if(~file_test(filename)) then begin
     print, 'red::initialize : ERROR : file not found: ', filename
     return
  endif
  self.filename = filename
  self.filetype = 'MOMFBD'
  
  ;; Init vars
  ;; self.data_dir = ''
  self.pinh_dir = '' 
  self.polcal_dir = ''
  self.camt = '' 
  self.camr = '' 
  self.camwb = '' 
   
  self.dodata = 1B
  self.doflat = 1B
  self.dopinh = 1B
  self.dodark = 1B
  self.dopolcal = 1B
  self.dodescatter = 1B

  ;; open file and get fields
  openr, lun, filename, /get_lun
  nl = 0L
  while(~eof(lun)) do begin
     line = ''
     readf, lun, line                          ; read line
     field = (strsplit(line,' =',/extract))[0] ; extract field
   
     if(strmid(line, 0, 1) eq '#') then continue
   
     ;; get fields
     case field of
        'root_dir': begin
            self.root_dir = (strsplit(line,' =',/extract))[1] ; extract value
            IF strmid(self.root_dir, strlen(self.root_dir)-1) NE '/' THEN $
              self.root_dir += '/'
        end
;        'descatter_dir': begin
;           self.descatter_dir = (strsplit(line,' =',/extract))[1] ; extract value
;        end
;        'dark_dir': begin
;           self.dark_dir = (strsplit(line,' =',/extract))[1] ; extract value
;        end
        'dark_dir': BEGIN
            tmp = strtrim((strsplit(line, '=', /extract))[1],2)
            IF(strpos(tmp,"'") NE -1) THEN dum = execute('tmp = '+tmp)
            IF ptr_valid(self.dark_dir) THEN red_append, *self.dark_dir, tmp $
            ELSE self.dark_dir = ptr_new(tmp, /NO_COPY)
        end
        'flat_dir': BEGIN
            tmp = strtrim((strsplit(line, '=', /extract))[1],2)
            IF(strpos(tmp,"'") NE -1) THEN dum = execute('tmp = '+tmp)
            IF ptr_valid(self.flat_dir) THEN red_append, *self.flat_dir, tmp $
            ELSE self.flat_dir = ptr_new(tmp, /NO_COPY)
        END
        'out_dir': begin
           self.out_dir = (strsplit(line,' =',/extract))[1] ; extract value
        end
        'pinh_dir': begin
           self.pinh_dir = (strsplit(line,' =',/extract))[1] ; extract value
        end
        'data_dir': begin

           if(strpos(line,"'") ne -1) then begin
              tmp = (strsplit(line,'=',/extract))[1]
              dum = execute('tmp=' + tmp)
           endif else tmp = (strsplit(line,' =',/extract))[1]
           nn = n_elements(tmp)
           self.ndir = nn
           for ii = 0, nn-1 do self.data_list[ii] = tmp[ii]
        end
        'polcal_dir': begin
           self.polcal_dir =  (strsplit(line,' =',/extract))[1] ; extract value
        end
        'filetype': begin
           self.filetype = (strsplit(line,' =',/extract))[1]
        end
        'cam_channel': BEGIN
            tmp = strtrim((strsplit(line, '=', /extract))[1],2)
            IF(strpos(tmp,"'") NE -1) THEN dum = execute('tmp = '+tmp)
            IF ptr_valid(self.cam_channels) THEN red_append, *self.cam_channels, tmp $
            ELSE self.cam_channels = ptr_new(tmp, /NO_COPY)
        END
        'prefilter_dir': begin
           self.prefilter_dir =  (strsplit(line,' =',/extract))[1] ; extract value
        end
        'telescope_d':begin
           self.telescope_d =  (strsplit(line,' =',/extract))[1] ; extract value
        end
        'image_scale':begin
           self.image_scale =  (strsplit(line,' =',/extract))[1] ; extract value
        end
        'pixel_size':begin
           self.pixel_size =  (strsplit(line,' =',/extract))[1] ; extract value
        end
        'camsz':begin
           self.camsz =  (strsplit(line,' =',/extract))[1] ; extract value
        end
        'isodate':begin
           self.isodate =  (strsplit(line,' =',/extract))[1] ; extract value
        end
        'pinhole_spacing': begin
            self.pinhole_spacing = float((strsplit(line,' =',/extract))[1]) ; extract value
        end
        else: begin
           print, 'red::initialize : Skipping line ',$
                  strcompress(string(nl),/remove_all), ': ',line
        end
     endcase
     nl+= 1
  endwhile
  free_lun, lun

  if(strlen(self.root_dir) gt 0) then begin
     print, 'red::initialize : root_dir = '+ self.root_dir
     IF ptr_valid(self.flat_dir) THEN FOR i = 0, n_elements(*self.flat_dir)-1 DO $
       IF strmid((*self.flat_dir)[i], 0, 1) NE '/' AND strlen((*self.flat_dir)[i]) GT 0 $
       THEN (*self.flat_dir)[i] = self.root_dir + (*self.flat_dir)[i]
     IF ptr_valid(self.dark_dir) THEN FOR i = 0, n_elements(*self.dark_dir)-1 DO $
       IF strmid((*self.dark_dir)[i], 0, 1) NE '/' AND strlen((*self.dark_dir)[i]) GT 0 $
       THEN (*self.dark_dir)[i] = self.root_dir + (*self.dark_dir)[i]
;     if(strmid(self.dark_dir, 0, 1) NE '/' AND strlen(self.dark_dir) gt 0) then self.dark_dir = self.root_dir + self.dark_dir
     if(strmid(self.pinh_dir, 0, 1) NE '/' AND strlen(self.pinh_dir) gt 0) then self.pinh_dir = self.root_dir + self.pinh_dir
     if(strmid(self.prefilter_dir, 0, 1) NE '/' AND strlen(self.prefilter_dir) gt 0) then self.prefilter_dir = self.root_dir + self.prefilter_dir
     if(strmid(self.data_list[self.ndir-1], 0, 1) NE '/' AND strlen(self.data_list[self.ndir-1]) gt 0) then self.data_list[0:self.ndir-1] = self.root_dir + self.data_list[0:self.ndir-1]
     if(strmid(self.polcal_dir, 0, 1) NE '/' AND strlen(self.polcal_dir) gt 0) then self.polcal_dir = self.root_dir + self.polcal_dir
;     if(strmid(self.descatter_dir, 0, 1) NE '/' AND strlen(self.descatter_dir) gt 0) then self.descatter_dir = self.root_dir + self.descatter_dir
  endif
  if(strlen(self.telescope_d) eq 0) then self.telescope_d = '0.970'
  print, 'red::initialize : telescope_d = '+self.telescope_d 

  if(strlen(self.image_scale) eq 0) then self.image_scale = '0.0592'
  print, 'red::initialize : image_scale = '+self.image_scale 

  if(strlen(self.pixel_size) eq 0) then self.pixel_size = '16.0E-6'
  print, 'red::initialize : pixel_size = '+self.pixel_size 

  if(strlen(self.camsz) eq 0) then self.camsz = '1024'
  print, 'red::initialize : camsz = '+self.camsz

  if self.pinhole_spacing eq 0.0 then begin
     ;; Default value (5.12 arcseconds between pinholes) was
     ;; measured in 2013 by comparison with SDO/HMI images.
     self.pinhole_spacing = 5.12
  endif 
  print, 'red::initialize : pinhole_spacing = '+strtrim(self.pinhole_spacing, 2)

  ;; check available fields
;  if(self.descatter_dir eq '') then begin
;     print, 'red::initialize : WARNING : descatter_dir is undefined!'
;     self.dodescatter = 0B
;  endif
;  if(self.dark_dir eq '') then begin
  IF ~ptr_valid(self.dark_dir) then begin
     print, 'red::initialize : WARNING : dark_dir is undefined!'
     self.dodark = 0B
  endif
  IF ~ptr_valid(self.flat_dir) then begin
     print, 'red::initialize : WARNING : flat_dir is undefined!'
     self.doflat = 0B
  endif
  if(self.out_dir eq '') then begin
     print, 'red::initialize : ERROR : out_dir is undefined!'
     return
  endif
  ;;if(self.data_dir eq '') then begin
  ;;   print, 'red::initialize : WARNING : data_dir is undefined!'
  ;;   self.dodata = 0B
  ;;endif
  if(self.pinh_dir eq '') then begin
     print, 'red::initialize : WARNING : pinh_dir is undefined!'
     self.dopinh = 0B
  endif
  if(self.polcal_dir eq '') then begin
     print, 'red::initialize : WARNING : polcal_dir is undefined!'
     self.dopolcal = 0B
  endif
  if(self.isodate eq '') then begin
     print, 'red::initialize : WARNING : isodate is undefined!'
     print, '                  Try to get it from PWD!'
     date = stregex(getenv('PWD'),'[12][0-9][0-9][0-9][-.][0-1][0-9][-.][0-3][0-9]',/extr)
     if date eq '' then begin
        print, 'red::initialize : WARNING : No recognizable date in PWD. Giving up.'
     endif else begin
        ;; Do a red_strreplace in case the found date uses dots rather
        ;; than dashes.
        self.isodate = red_strreplace(date, '.', '-', n = 2)
     endelse
  endif
 
  ;; Fields that depend on fields defined above:
  self.log_dir = self.out_dir+'/downloads/sstlogs/'
  self.telog = self.log_dir+'positionLog_'+red_strreplace(self.isodate, '-', '.', n = 2)+'_final'
  self.descatter_dir = self.out_dir+'/downloads/backscatter/'



  ;; print fields
  IF ptr_valid(self.dark_dir) THEN BEGIN
  if(self.dopinh) then print, 'red::initialize : pinh_dir = '+ self.pinh_dir
  if(self.dopolcal) then print, 'red::initialize : polcal_dir = '+ self.polcal_dir
     nn = n_elements(*self.dark_dir)
     IF(nn EQ 1) THEN BEGIN
        print, 'red::initialize : dark_dir = '+ (*self.dark_dir)[0]
     ENDIF ELSE BEGIN
        print, 'red::initialize : dark_dirs :'
        FOR k = 0, nn-1 DO print, string(k, format='(I5)') + ' -> ' + (*self.dark_dir)[k]
     ENDELSE
  ENDIF
  IF ptr_valid(self.flat_dir) THEN BEGIN
     nn = n_elements(*self.flat_dir)
     IF(nn EQ 1) THEN BEGIN
        print, 'red::initialize : flat_dir = '+ (*self.flat_dir)[0]
     ENDIF ELSE BEGIN
        print, 'red::initialize : flat_dirs :'
        FOR k = 0, nn-1 DO print, string(k, format='(I5)') + ' -> ' + (*self.flat_dir)[k]
     ENDELSE
  ENDIF
  IF ptr_valid(self.cam_channels) THEN BEGIN
     nn = n_elements(*self.cam_channels)
     IF(nn EQ 1) THEN BEGIN
        print, 'red::initialize : cam_dir = '+ (*self.cam_channels)[0]
     ENDIF ELSE BEGIN
        print, 'red::initialize : cam_channels :'
        FOR k = 0, nn-1 DO print, string(k, format='(I5)') + ' -> ' + (*self.cam_channels)[k]
     ENDELSE
  ENDIF
  ;if(self.doflat) then print, 'red::initialize : flat_dir = '+ strjoin(*self.flat_dir, '  ')
  if(self.dodata) then begin
     nn = self.ndir
     if(nn eq 1) then begin
        print, 'red::initialize : data_dir = '+ self.data_list[0]
        self.data_dir = self.data_list[0]
     endif else begin
        print, 'red::initialize : data_dirs :'
        for k = 0, nn-1 do print, string(k, format='(I5)') + ' -> ' + self.data_list[k]
        id = 0
        ;; read, id, prompt = "red::initialize : select default folder's id :"
        self.data_dir = self.data_list[id]
     endelse
  endif
  if(self.dodescatter) then print, 'red::initialize : descatter_dir = '+ self.descatter_dir
  if(self.filetype) then print, 'red::initialize : filetype = '+ self.filetype
  print, 'red::initialize : out_dir = '+ self.out_dir

  cgWindow_SetDefs, PS_Decomposed=1
                   
  return
end
