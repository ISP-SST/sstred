; docformat = 'rst'

;+
; Identify cameras and detectors for WB, NB, PD, as well as the instrument.
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
;    instrument : out, optional, type=string
; 
;       The name of the instrument, with first letter in caps, like
;       Crisp.
; 
;    nb_cameras : out, optional, type=string
;   
;       NB camera names.
; 
;    nb_detectors  : out, optional, type=string
;   
;       NB detector names.
;   
;    Nnbcams : out, optional, type=integer
;   
;       Number of NB cameras.
;   
;    nbr_cameras : out, optional, type=string
;   
;       NB reflected camera names.
;   
;    nbr_detectors : out, optional, type=string
;   
;       NB reflected detector names.
;   
;    Nnbrcams : out, optional, type=integer
;   
;       Number of NB reflected cameras.
;   
;    nbt_cameras : out, optional, type=string
;   
;       NB transmitted camera names.
;   
;    nbt_detectors : out, optional, type=string
;   
;      NB transmitted detector names.
;   
;    Nnbtcams  : out, optional, type=integer
;   
;      Number of NB transmitted cameras.
;   
;    pd_camera   : out, optional, type=string
;   
;      PD camera names.
;   
;    pd_detector : out, optional, type=string
;   
;      PD detector names.
;   
;    Npdcams : out, optional, type=string
;   
;      Number of PD cameras.
;   
;    wb_camera  : out, optional, type=string
;   
;      WB camera names.
;   
;    wb_detector : out, optional, type=string
;   
;      WB detector names.
;   
;    Nwbcams : out, optional, type=string
;   
;      Number of WB cameras.
;   
; 
; :History:
; 
;    2025-09-19 : MGL. First version.
;
;    2025-11-05 : MGL. New keywords nb[rt]_cameras, nb[rt]_detectors,
;                 Nnb[rt], instrument.
; 
;-
pro red::cameras $
   , instrument = instrument $
   , nb_cameras  = nb_cameras, nb_detectors = nb_detectors, Nnbcams  = Nnb $
   , nbr_camera = nbr_camera,  nbr_detector = nbr_detector, Nnbrcams = Nnbr $
   , nbt_camera = nbt_camera,  nbt_detector = nbt_detector, Nnbtcams = Nnbt $
   , pd_camera   = pd_camera,  pd_detector  = pd_detector,  Npdcams  = Npd $
   , wb_camera   = wb_camera,  wb_detector  = wb_detector,  Nwbcams  = Nwb

  self -> getdetectors

  wb_indx  = where(strmatch(*self.cameras, '*-W'), Nwb)
  pd_indx  = where(strmatch(*self.cameras, '*-D'), Npd)
  nbt_indx = where(strmatch(*self.cameras, '*-T'), Nnbt)
  nbr_indx = where(strmatch(*self.cameras, '*-R'), Nnbr)
  nb_indx  = where(strmatch(*self.cameras, '*-[NTR]'), Nnb)
  
  wb_camera   = (*self.cameras)[wb_indx[0]]
  wb_detector = (*self.detectors)[wb_indx[0]]
  instrument = (strsplit(wb_camera, '-', /extract))[0]
  
  nb_cameras   = (*self.cameras)[nb_indx]
  nb_detectors = (*self.detectors)[nb_indx]

  if Npd gt 0 then begin
    pd_camera   = (*self.cameras)[pd_indx[0]]
    pd_detector = (*self.detectors)[pd_indx[0]]
  endif

  if Nnbt gt 0 then begin
    nbt_camera   = (*self.cameras)[nbt_indx[0]]
    nbt_detector = (*self.detectors)[nbt_indx[0]]
  endif

  if Nnbr gt 0 then begin
    nbr_camera   = (*self.cameras)[nbr_indx[0]]
    nbr_detector = (*self.detectors)[nbr_indx[0]]
  endif

;  stop
  
end
