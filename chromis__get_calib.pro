; docformat = 'rst'

;+
; From state data, return file names for calibration data and the data
; itself. 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Params:
; 
;    states : in, type=structarr
;
;       The states, the calibrations of which the keywords refer to. 
; 
; :Keywords:
; 
;    darkname : out, optional, type=strarr 
; 
;       The name(s) of the dark file(s) appropriate to the state(s).
; 
;    darkdata : out, optional, type=array 
; 
;        The data in the dark file(s) appropriate to the state(s).    
; 
;    flatname : out, optional, type=strarr 
; 
;       The name(s) of the flat file(s) appropriate to the state(s).
; 
;    flatdata : out, optional, type=array 
; 
;        The data in the flat file(s) appropriate to the state(s).
; 
;    sflatname : out, optional, type=strarr 
; 
;       The name(s) of the summed flat file(s) appropriate to the state(s).
; 
;    sflatdata : out, optional, type=array 
; 
;        The data in the summed flat file(s) appropriate to the state(s).
; 
;    gainname : out, optional, type=strarr 
; 
;       The name(s) of the gain file(s) appropriate to the state(s).
; 
;    gaindata : out, optional, type=array 
; 
;        The data in the gain file(s) appropriate to the state(s).
; 
;    pinhname : out, optional, type=strarr 
; 
;        The name(s) of the pinhole file(s) appropriate to the state(s).
; 
;    pinhdata : out, optional, type=array
;   
;        The data in the pinhole file(s) appropriate to the state(s).  
; 
;    status : out, optional, type=integer
; 
;        The status of the operation, 0 for success.
; 
;    timestamp : in, optional, type=string
; 
;        Look for darks and flats in timestamped directories below the
;        regular darks/, flats/ subdirectories.
; 
; :History:
; 
;    2016-05-26 : MGL. Initial version including darks, flats, and
;                 pinholes. Added support for summed flats.
; 
;    2016-05-27 : MGL. Bugfixes.
; 
;    2016-06-02 : MGL. Added gains.
;
;    2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                 so the names match those of the corresponding SolarNet
;                 keywords.
;
;    2016-09-04 : THI. Wideband pinhole names without narrowband
;                 tuning info. 
; 
;    2018-04-18 : MGL. Rewrote using filenames() method. Removed
;                 keywords for returning file headers.
; 
;    2018-07-11 : MGL. New keywords, individual status for the
;                 different kinds of data.
; 
;    2018-07-25 : Sum darks and flats, and make gains as needed.
; 
;    2025-03-28 : MGL. New keyword timestamp.
; 
;-
pro chromis::get_calib, states $
                        , no_fits = no_fits $
                        , status = status $
                        , darkstatus  = darkstatus,  darkname  = darkname,  darkdata  = darkdata   $
                        , flatstatus  = flatstatus,  flatname  = flatname,  flatdata  = flatdata   $
                        , gainstatus  = gainstatus,  gainname  = gainname,  gaindata  = gaindata   $
                        , pinhstatus  = pinhstatus,  pinhname  = pinhname,  pinhdata  = pinhdata   $
                        , sflatstatus = sflatstatus, sflatname = sflatname, sflatdata = sflatdata  $
                        , timestamp = timestamp 
  

  Nstates = n_elements(states)

  if Nstates eq 0 then begin
    status = -1
    return
  endif

  if arg_present(darkname)  or arg_present(darkdata) then $
     darkname = self -> filenames('dark'   , states, no_fits = no_fits, timestamp = timestamp)
  if arg_present(flatname)  or arg_present(flatdata) or $
     arg_present(gainname)  or arg_present(gaindata) then $
        flatname = self -> filenames('flat'   , states, no_fits = no_fits, timestamp = timestamp)
  if arg_present(gainname)  or arg_present(gaindata) then $
     gainname = self -> filenames('gain'   , states, no_fits = no_fits)
  if arg_present(pinhname)  or arg_present(pinhdata) then $
     pinhname = self -> filenames('pinh'   , states, no_fits = no_fits)
  if arg_present(sflatname) or arg_present(sflatdata) then $
     sflatname = self -> filenames('sumflat', states, no_fits = no_fits)

  ;; Assume this is all for the same camera type, at least for the
  ;; actual data. Otherwise we cannot return the actual data in a
  ;; single array.
;  detector = states[0].detector
;  caminfo = red_camerainfo(detector)
  
;  if arg_present(darkdata)  then darkdata  = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
;  if arg_present(flatdata)  then flatdata  = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
;  if arg_present(gaindata)  then gaindata  = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
;  if arg_present(pinhdata)  then pinhdata  = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
;  if arg_present(sflatdata) then sflatdata = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
  

  status = 0

  ;; Darks
  if arg_present(darkdata) then begin
    if n_elements(darkname) ne 0  then begin

      for istate = 0, Nstates-1 do begin        
        if ~file_test(darkname[istate]) then begin
          ;; Try summing darks for this camera
          self -> sumdark, /check, /sum_in_rdx $
             , cams = states[istate].camera
        endif
      endfor                    ; istate

      darkdata = red_readdata_multiframe(darkname, status = darkstatus, /silent)
      status = min([status, darkstatus])

    endif else status = -1
  endif                         ; Darks
  

  ;; Flats
  if arg_present(flatdata) then begin
    if n_elements(flatname) ne 0 then begin

      for istate = 0, Nstates-1 do begin
        if ~file_test(flatname[istate]) then begin
          ;; Try summing flats for this state
          self -> sumflat, /check, /sum_in_rdx $
             , cams = states[istate].camera $
             , ustat = states[istate].fullstate
        endif 
      endfor                    ; istate

      flatdata = red_readdata_multiframe(flatname, status = flatstatus, /silent)
      status = min([status, flatstatus])

    endif else status = -1
  endif                         ; Flats

  ;; Summed flats
  if arg_present(sflatdata) then begin
    if n_elements(sflatname) ne 0 then begin

      sflatdata = red_readdata_multiframe(sflatname, status = sflatstatus, /silent)
      status = min([status, sflatstatus])
      
    endif else status = -1        
  endif                         ; Summed flats

  ;; Gains
  if arg_present(gaindata) then begin
    if  n_elements(gainname) ne 0 then begin

      for istate = 0, Nstates-1 do begin
        if ~file_test(gainname[istate]) then begin
          ;; Try summing flats for this state and then making gains
          if ~file_test(flatname[istate]) then begin
            self -> sumflat, /check, /sum_in_rdx $
               , cams = states[istate].camera $
               , ustat = states[istate].fullstate
          endif 
          self -> makegains, smooth=3.0, files = flatname[istate]
        endif 
      endfor                    ; istate

      gaindata = red_readdata_multiframe(gainname, status = gainstatus, /silent)
      status = min([status, gainstatus])

    endif else status = -1
  endif                         ; Gains

  
  ;; Pinholes
  if arg_present(pinhdata) then begin

    if n_elements(pinhname) ne 0 then begin
      
      pinhdata = red_readdata_multiframe(pinhname, status = pinhstatus, /silent)
      status = min([status, pinhstatus])

    endif else status = -1
  endif                         ; Pinholes

  ;; Reduce dimensions if possible
  if Nstates eq 1 then begin

    if arg_present(darkname)  then darkname = darkname[0]
    if arg_present(flatname)  then flatname = flatname[0]
    if arg_present(gainname)  then gainname = gainname[0]
    if arg_present(pinhname)  then pinhname = pinhname[0]
    if arg_present(sflatname) then sflatname = sflatname[0]

    if arg_present(darkdata)  then darkdata = darkdata[*, *, 0]
    if arg_present(flatdata)  then flatdata = flatdata[*, *, 0]
    if arg_present(gaindata)  then gaindata = gaindata[*, *, 0]
    if arg_present(pinhdata)  then pinhdata = pinhdata[*, *, 0] 
    if arg_present(sflatdata) then sflatdata = sflatdata[*, *, 0]
    
  endif

end
