; docformat = 'rst'

;+
; From state data, return file names for calibration data and the data
; itself. 
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
;    cflatname : out, optional, type=strarr 
; 
;       The name(c) of the cavityfree flat file(s) appropriate to the state(s).
; 
;    cflatdata : out, optional, type=array 
; 
;        The data in the cavityfree flat file(s) appropriate to the state(s).
; 
;    gainname : out, optional, type=strarr 
; 
;       The name(s) of the gain file(s) appropriate to the state(s).
; 
;    gaindata : out, optional, type=array 
; 
;        The data in the gain file(s) appropriate to the state(s).
; 
;    sgainname : out, optional, type=strarr 
; 
;       The name(s) of the scan gain file(s) appropriate to the state(s).
; 
;    sgaindata : out, optional, type=array 
; 
;        The data in the scan gain file(s) appropriate to the state(s).
; 
;    pinhname : out, optional, type=strarr 
; 
;        The name(s) of the pinhole file(s) appropriate to the state(s).
; 
;    pinhdata : out, optional, type=array
;   
;        The data in the pinhole file(s) appropriate to the state(s).  
;   
;    polcname : out, optional, type=strarr 
; 
;        The name(s) of the polcal file(s) appropriate to the prefilter(s).
; 
;    polcdata : out, optional, type=array
;   
;        The data in the polcal file(s) appropriate to the prefilter(s).  
; 
;    polsname : out, optional, type=strarr 
; 
;        The name(s) of the polcal sum file(s) appropriate to the state(s).
; 
;    polsdata : out, optional, type=array
;   
;        The data in the polcal sum file(s) appropriate to the state(s).  
; 
;    status : out, optional, type=integer
; 
;        The status of the operation, 0 for success.
; 
; :History:
; 
;    2018-04-20 : MGL. New version based on chromis::get_calib. 
; 
;    2018-07-20 : MGL. New keywords, individual status for the
;                 different kinds of data.
; 
;    2018-07-25 : MGL. Sum darks and flats, and make gains as needed.
; 
;    2018-11-12 : MGL. New keywords cflatname, cflatstatus, cflatdata.
; 
;    2019-10-10 : MGL. New keywords cgainname, cgainstatus, cgaindata.
; 
;    2021-10-06 : MGL. New keywords sgainname, sgainstatus, sgaindata.
; 
;    2022-08-02 : MGL. Change from a CRISP:: method to a RED:: method
;                 to prepare for CRISP camera upgrade.
; 
;-
pro red::get_calib, states $
                    , no_fits = no_fits $
                    , status = status $
                    , darkstatus  = darkstatus,  darkname  = darkname,  darkdata  = darkdata  $
                    , flatstatus  = flatstatus,  flatname  = flatname,  flatdata  = flatdata  $
                    , gainstatus  = gainstatus,  gainname  = gainname,  gaindata  = gaindata  $
                    , pinhstatus  = pinhstatus,  pinhname  = pinhname,  pinhdata  = pinhdata  $
                    ,                            polcname  = polcname,  polcdata  = polcdata  $
                    ,                            polsname  = polsname,  polsdata  = polsdata  $
                    , sflatstatus = sflatstatus, sflatname = sflatname, sflatdata = sflatdata $
                    , cflatstatus = cflatstatus, cflatname = cflatname, cflatdata = cflatdata $
                    , cgainstatus = cgainstatus, cgainname = cgainname, cgaindata = cgaindata $
                    , sgainstatus = sgainstatus, sgainname = sgainname, sgaindata = sgaindata   

  Nstates = n_elements(states)

  if Nstates eq 0 then begin
    status = -1
    return
  endif
  
  if arg_present(darkname)  or arg_present(darkdata) then $
     darkname = self -> filenames('dark'   , states, no_fits = no_fits)
  if arg_present(flatname)  or arg_present(flatdata) or $
     arg_present(gainname)  or arg_present(gaindata) then $
        flatname = self -> filenames('flat'   , states, no_fits = no_fits)
  if arg_present(gainname)  or arg_present(gaindata) then $
     gainname = self -> filenames('gain'   , states, no_fits = no_fits)
  if arg_present(pinhname)  or arg_present(pinhdata) then $
     pinhname = self -> filenames('pinh'   , states, no_fits = no_fits)
  if arg_present(polsname)  or arg_present(polsdata) then $
     polsname = self -> filenames('pols'   , states, no_fits = no_fits)
  if arg_present(polcname)  or arg_present(polcdata) then $
     polcname = self -> filenames('polc'   , states, no_fits = no_fits)
  if arg_present(sflatname) or arg_present(sflatdata) then $
     sflatname = self -> filenames('sumflat', states, no_fits = no_fits)
  if arg_present(cflatname) or arg_present(cflatdata) then $
     cflatname = self -> filenames('cavityflat', states, no_fits = no_fits)
  if arg_present(cgainname) or arg_present(cgaindata) then $
     cgainname = self -> filenames('cavityfree_gain', states, no_fits = no_fits)
  if arg_present(sgainname) or arg_present(sgaindata) then $
     sgainname = self -> filenames('scangain', states, no_fits = no_fits)

  ;; Assume this is all for the same camera type, at least for the
  ;; actual data. Otherwise we cannot return the actual data in a
  ;; single array.
;  detector = states[0].detector
;  caminfo = red_camerainfo(detector)
  
;  if arg_present(darkdata)  then darkdata  = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
;  if arg_present(flatdata)  then flatdata  = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
;  if arg_present(gaindata)  then gaindata  = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
;  if arg_present(pinhdata)  then pinhdata  = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
;  if arg_present(polcdata)  then polcdata  = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
;  if arg_present(polsdata)  then polsdata  = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
;  if arg_present(sflatdata) then sflatdata = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
;  if arg_present(cflatdata) then cflatdata = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
;  if arg_present(cgaindata) then cgaindata = fltarr(caminfo.xsize, caminfo.ysize, Nstates)   
;  if arg_present(sgaindata) then sgaindata = fltarr(caminfo.xsize, caminfo.ysize, Nstates)   


  status = 0

  ;; Darks    
  if arg_present(darkdata) then begin    
    if n_elements(darkname) ne 0 then begin      

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
  
  ;; Cavityfree flats    
  if arg_present(cflatdata) then begin    
    if n_elements(cflatname) ne 0 then begin
      
      cflatdata = red_readdata_multiframe(cflatname, status = cflatstatus, /silent)
      status = min([status, cflatstatus])
      
    endif else status = -1        
  endif                         ; Cavityfree flats
  
  ;; Gains    
  if arg_present(gaindata) then begin    
    if n_elements(gainname) ne 0 then begin      

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

  ;; Cavityfree gains    
  if arg_present(cgaindata) then begin
    if n_elements(cgainname) ne 0 then begin
      
      cgaindata = red_readdata_multiframe(cgainname, status = cgainstatus, /silent)
      status = min([status, cgainstatus])

    endif else status = -1
  endif                         ; Cavityfree gains             
  
  ;; Scan gains    
  if arg_present(sgaindata) then begin    

    if n_elements(sgainname) ne 0 then begin
      
      sgaindata = red_readdata_multiframe(sgainname, status = sgainstatus, /silent)
      status = min([status, sgainstatus])

    endif else status = -1
  endif                    

    
  
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

a = crispred(/dev, /no)

files = file_search('/data/2016/2016.09/2016.09.19/Flats/11:21:22/Crisp-T/*', count = Nfiles)

a -> extractstates, files[0:10], states

a -> get_calib, states, darkdata = darkdata, flatdata = flatdata

end
