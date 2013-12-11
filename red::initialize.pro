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
;-
pro red::initialize, filename
                                
  ;; Test file
  if(~file_test(filename)) then begin
     print, 'red::initialize : ERROR : file not found: ', filename
     return
  endif
  self.filename = filename
  
  ;; Init vars
  self.dark_dir = '' 
  self.flat_dir = ptr_new('') 
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
  self.docamt = 1B
  self.docamr = 1B
  self.docamwb = 1B
  self.dopolcal = 1B
  self.dodescatter = 1B
   
  ;; open file and get fields
  openr, lun, filename, /get_lun
  nl = 0L
  while(~eof(lun)) do begin
     line = ''
     readf, lun, line                          ; read line
     field = (strsplit(line,' =',/extract))[0] ; extract field
   
     if(strmid(line, 0, 1) eq '#') then begin
        print, 'red::initialize : Skipping commented-out line'
        continue
     endif
   
     ;; get fields
     case field of
        'root_dir': begin
            self.root_dir = (strsplit(line,' =',/extract))[1] ; extract value
            IF strmid(self.root_dir, strlen(self.root_dir)-1) NE '/' THEN $
              self.root_dir += '/'
        end
        'descatter_dir': begin
           self.descatter_dir = (strsplit(line,' =',/extract))[1] ; extract value
        end
        'dark_dir': begin
           self.dark_dir = (strsplit(line,' =',/extract))[1] ; extract value
        end
        'flat_dir': BEGIN
            if(strpos(line,"'") ne -1) THEN $
              dum = execute('tmp = '+(strsplit(line, ' =', /extract))[1]) $
            ELSE $
              tmp = (strsplit(line,' =',/extract))[1]
            ptr_free, self.flat_dir
            self.flat_dir = ptr_new(tmp, /NO_COPY)
        end
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
        'cam_t': begin
           self.camt =  (strsplit(line,' =',/extract))[1] ; extract value
        end
        'cam_r': begin
           self.camr =  (strsplit(line,' =',/extract))[1] ; extract value
        end
        'cam_wb': begin
           self.camwb =  (strsplit(line,' =',/extract))[1] ; extract value
        end
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
        else: begin
           print, 'red::initialize : Skipping line ',$
                  strcompress(string(nl),/remove_all), ': ',line
        end
     endcase
     nl+= 1
  endwhile
  free_lun, lun


  if(strlen(self.root_dir) gt 0) then begin
     if(self.docamt) then print, 'red::initialize : root_dir = '+ self.root_dir
     FOR i = 0, n_elements(*self.flat_dir)-1 DO $
       IF strmid((*self.flat_dir)[i], 0, 1) NE '/' AND strlen((*self.flat_dir)[i]) GT 0 $
       THEN (*self.flat_dir)[i] = self.root_dir + (*self.flat_dir)[i]
     if(strmid(self.dark_dir, 0, 1) NE '/' AND strlen(self.dark_dir) gt 0) then self.dark_dir = self.root_dir + self.dark_dir
     if(strmid(self.pinh_dir, 0, 1) NE '/' AND strlen(self.pinh_dir) gt 0) then self.pinh_dir = self.root_dir + self.pinh_dir
     if(strmid(self.prefilter_dir, 0, 1) NE '/' AND strlen(self.prefilter_dir) gt 0) then self.prefilter_dir = self.root_dir + self.prefilter_dir
     if(strmid(self.data_list[self.ndir-1], 0, 1) NE '/' AND strlen(self.data_list[self.ndir-1]) gt 0) then self.data_list[0:self.ndir-1] = self.root_dir + self.data_list[0:self.ndir-1]
     if(strmid(self.polcal_dir, 0, 1) NE '/' AND strlen(self.polcal_dir) gt 0) then self.polcal_dir = self.root_dir + self.polcal_dir
     if(strmid(self.descatter_dir, 0, 1) NE '/' AND strlen(self.descatter_dir) gt 0) then self.descatter_dir = self.root_dir + self.descatter_dir
  endif
  if(strlen(self.telescope_d) eq 0) then self.telescope_d = '0.970'
  print, 'red::initialize : telescope_d = '+self.telescope_d 

  if(strlen(self.image_scale) eq 0) then self.image_scale = '0.0592'
  print, 'red::initialize : image_scale = '+self.image_scale 

  if(strlen(self.pixel_size) eq 0) then self.pixel_size = '16.0E-6'
  print, 'red::initialize : pixel_size = '+self.pixel_size 

  if(strlen(self.camsz) eq 0) then self.camsz = '1024'
  print, 'red::initialize : camsz = '+self.camsz

  ;; check available fields
  if(self.descatter_dir eq '') then begin
     print, 'red::initialize : WARNING : descatter_dir is undefined!'
     self.dodescatter = 0B
  endif
  if(self.dark_dir eq '') then begin
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

  ;; print fields
  if(self.dodark) then print, 'red::initialize : dark_dir = '+ self.dark_dir
  if(self.doflat) then print, 'red::initialize : flat_dir = '+ strjoin(*self.flat_dir, '  ')
  if(self.dopinh) then print, 'red::initialize : pinh_dir = '+ self.pinh_dir
  if(self.dopolcal) then print, 'red::initialize : polcal_dir = '+ self.polcal_dir
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
  if(self.docamt) then print, 'red::initialize : cam_T = '+ self.camt
  if(self.docamr) then print, 'red::initialize : cam_R = '+ self.camr
  if(self.docamwb) then print, 'red::initialize : cam_wb = '+ self.camwb
  if(self.dodescatter) then print, 'red::initialize : descatter_dir = '+ self.descatter_dir

  print, 'red::initialize : out_dir = '+ self.out_dir
                                
  return
end
