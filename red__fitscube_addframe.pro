; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
; 
; 
; :Params:
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
;    2016-03-24 : MGL. First version.
; 
; 
;-
pro red::fitscube_addframe, fileassoc, frame $
                              , ituning = ituning $
                              , istokes = istokes $
                              , iscan = iscan $
                              , Ntuning = Ntuning $
                              , Nstokes = Nstokes $
                              , Nscan = Nscan 

  if n_elements(Ntuning) eq 0 then Ntuning = 1L
  if n_elements(Nstokes) eq 0 then Nstokes = 1L
  if n_elements(Nscan)   eq 0 then Nscans  = 1L

  if n_elements(ituning) eq 0 then ituning = 0L
  if n_elements(istokes) eq 0 then istokes = 0L
  if n_elements(iscan)   eq 0 then iscan   = 0L

  iframe = long(ituning) + long(istokes)*long(Ntuning) $
           + long(iscan)*long(Ntuning)*long(Nstokes)
;  print, 'Adding frame ', iscan, iframe

  fileassoc[iframe] = frame
  
end
