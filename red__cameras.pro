; docformat = 'rst'

;+
; Identify cameras and detectors for WB, NB, PD.
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
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;    2025-09-19 : MGL. First version.
; 
;-
pro red::cameras $
   , nb_cameras = nb_cameras,  nb_detectors = nb_detectors, Nnb = Nnb $
   , pd_camera  = pd_camera,   pd_detector  = pd_detector,  Npd = Npd $
   , wb_camera  = wb_camera,   wb_detector  = wb_detector,  Nwb = Nwb

  self -> getdetectors

  wb_indx = where(strmatch(*self.cameras, '*-W'), Nwb)
  nb_indx = where(strmatch(*self.cameras, '*-[NTR]'), Nnb)
  pd_indx = where(strmatch(*self.cameras, '*-D'), Npd)
  
  wb_camera   = (*self.cameras)[wb_indx[0]]
  wb_detector = (*self.detectors)[wb_indx[0]]

  nb_cameras   = (*self.cameras)[nb_indx]
  nb_detectors = (*self.detectors)[nb_indx]

  if Npd gt 0 then begin
    pd_camera   = (*self.cameras)[pd_indx[0]]
    pd_detector = (*self.detectors)[pd_indx[0]]
  endif

end
