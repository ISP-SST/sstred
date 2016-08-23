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
;    darkhead : out, optional, type=array 
; 
;        The header of the dark file(s) appropriate to the state(s).      
; 
;    flatname : out, optional, type=strarr 
; 
;       The name(s) of the flat file(s) appropriate to the state(s).
; 
;    flatdata : out, optional, type=array 
; 
;        The data in the flat file(s) appropriate to the state(s).
; 
;    flathead : out, optional, type=array 
; 
;        The header of the flat file(s) appropriate to the state(s).
; 
;    sflatname : out, optional, type=strarr 
; 
;       The name(s) of the summed flat file(s) appropriate to the state(s).
; 
;    sflatdata : out, optional, type=array 
; 
;        The data in the summed flat file(s) appropriate to the state(s).
; 
;    sflathead : out, optional, type=array 
; 
;        The header of the summed flat file(s) appropriate to the state(s).
; 
;    gainname : out, optional, type=strarr 
; 
;       The name(s) of the gain file(s) appropriate to the state(s).
; 
;    gaindata : out, optional, type=array 
; 
;        The data in the gain file(s) appropriate to the state(s).
; 
;    gainhead : out, optional, type=array 
; 
;        The header of the gain file(s) appropriate to the state(s).
; 
;    pinhname : out, optional, type=strarr 
; 
;        The name(s) of the pinhole file(s) appropriate to the state(s).
; 
;    pinhdata : out, optional, type=array
;   
;        The data in the pinhole file(s) appropriate to the state(s).  
; 
;    pinhhead : out, optional, type=array
;   
;        The header of the pinhole file(s) appropriate to the state(s).  
;   
;    status : out, optional, type=integer
; 
;        The status of the operation, 0 for success.
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
; 
;-
pro chromis::get_calib, states $
                        , status = status $
                        , darkname = darkname, darkdata = darkdata, darkhead = darkhead $
                        , flatname = flatname, flatdata = flatdata, flathead = flathead $
                        , gainname = gainname, gaindata = gaindata, gainhead = gainhead $
                        , pinhname = pinhname, pinhdata = pinhdata, pinhhead = pinhhead $
                        , sflatname = sflatname, sflatdata = sflatdata, sflathead = sflathead 

  Nstates = n_elements(states)

  if Nstates eq 0 then begin
     status = -1
     return
  endif

  if arg_present(darkname)  or arg_present(darkdata)  then darkname = strarr(Nstates)   
  if arg_present(flatname)  or arg_present(flatdata)  then flatname = strarr(Nstates) 
  if arg_present(gainname)  or arg_present(gaindata)  then gainname = strarr(Nstates) 
  if arg_present(pinhname)  or arg_present(pinhdata)  then pinhname = strarr(Nstates) 
  if arg_present(sflatname) or arg_present(sflatdata) then sflatname = strarr(Nstates) 

  if arg_present(darkdata) $
     or arg_present(flatdata) $
     or arg_present(gaindata) $
     or arg_present(sflatdata) $
     or arg_present(pinhdata) then begin

     ;; Assume this is all for the same camera type, at least for the
     ;; actual data. Otherwise we cannot return the actual data in a
     ;; single array.
     detector = states[0].detector
     caminfo = red_camerainfo(detector)

     if arg_present(darkdata) then darkdata = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
     if arg_present(flatdata) then flatdata = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
     if arg_present(gaindata) then gaindata = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
     if arg_present(pinhdata) then pinhdata = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
     if arg_present(sflatdata) then sflatdata = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
  endif

  Nheadlines = 100              ; Assume max numer of header lines
  if arg_present(darkhead) then darkhead = arr(Nheadlines,Nstates) 
  if arg_present(flathead) then flathead = arr(Nheadlines,Nstates) 
  if arg_present(gainhead) then gainhead = arr(Nheadlines,Nstates) 
  if arg_present(pinhhead) then pinhhead = arr(Nheadlines,Nstates) 
  if arg_present(sflathead) then sflathead = arr(Nheadlines,Nstates) 

  status = 0

  for istate = 0, Nstates-1 do begin

     detector = states[istate].detector

     ;; Darks
     if arg_present(darkname) or arg_present(darkdata) or arg_present(darkhead) then begin

        darktag = detector
        if( states[istate].cam_settings ne '' ) then begin
            darktag += '_' + states[istate].cam_settings
        endif

        dname = self.out_dir+'/darks/' + darktag + '.dark'
        if arg_present(darkname) then begin
            darkname[istate] = dname
        endif

        if file_test(dname) then begin
            if arg_present(darkdata) then begin
                darkdata[0, 0, istate] = red_readdata(dname, header = darkhead $
                                                      , status = darkstatus, /silent)
                if status eq 0 then status = darkstatus
            endif else if arg_present(darkhead) then begin
                darkhead[0, istate] = red_readhead(dname, status = darkstatus, /silent)
                if status eq 0 then status = darkstatus
            endif
        endif else begin
            if( arg_present(darkdata) || arg_present(darkhead) ) then status = -1
        endelse

     endif                      ; Darks

     ;; Flats
     if arg_present(flatname) or arg_present(flatdata) or arg_present(flathead) $
        or arg_present(sflatname) or arg_present(sflatdata) or arg_present(sflathead) then begin

        flattag = detector
        if( states[istate].fullstate ne '' ) then begin
            flattag += '_' + states[istate].fullstate
        endif

        fname = self.out_dir+'/flats/' + flattag + '.flat'
        if arg_present(flatname) then begin
            flatname[istate] = fname
        endif

        if file_test(fname) then begin
            if arg_present(flatdata) then begin
                flatdata[0, 0, istate] = red_readdata(fname, header = flathead $
                                                      , status = flatstatus, /silent)
                if status eq 0 then status = flatstatus
            endif else if arg_present(flathead) then begin
                flathead[0, istate] = red_readhead(fname, status = flatstatus, /silent)
                if status eq 0 then status = flatstatus
            endif
        endif else begin
            if( arg_present(flatdata) || arg_present(flathead) ) then status = -1
        endelse

        sfname = self.out_dir+'/flats/' + flattag + '_summed.flat'
        if arg_present(sflatname) then begin
            sflatname[istate] = sfname
        endif
           
        if file_test(sfname) then begin
            if arg_present(sflatdata) then begin
                sflatdata[0, 0, istate] = red_readdata(sfname, header = sflathead $
                                                      , status = sflatstatus, /silent)
                if status eq 0 then status = sflatstatus
            endif else if arg_present(flathead) then begin
                sflathead[0, istate] = red_readhead(sfname, status = sflatstatus, /silent)
                if status eq 0 then status = sflatstatus
            endif
        endif else begin
            if( arg_present(sflatdata) || arg_present(sflathead) ) then status = -1
        endelse

     endif                      ; Flats

     ;; Gains
     if arg_present(gainname) or arg_present(gaindata) or arg_present(gainhead) then begin

        gaintag = detector
        if( states[istate].fullstate ne '' ) then begin
           gaintag += '_' + states[istate].fullstate
        endif

        fname = self.out_dir+'/gaintables/' + gaintag + '.gain'
        if arg_present(gainname) then begin
           gainname[istate] = fname
        endif

        if file_test(fname) then begin
           if arg_present(gaindata) then begin
              gaindata[0, 0, istate] = red_readdata(fname, header = gainhead $
                                                    , status = gainstatus, /silent)
              if status eq 0 then status = gainstatus
           endif else if arg_present(gainhead) then begin
              gainhead[0, istate] = red_readhead(fname, status = gainstatus, /silent)
              if status eq 0 then status = gainstatus
           endif
        endif else begin
           if( arg_present(gaindata) || arg_present(gainhead) ) then status = -1
        endelse

     endif                      ; Gains

     
     ;; Pinholes
     if arg_present(pinhname) or arg_present(pinhdata) or arg_present(pinhhead) then begin

        pinhtag = detector
        if( states[istate].prefilter ne '' ) then begin
            pinhtag += '_' + states[istate].prefilter
        endif
        if( states[istate].tuning ne '' ) then begin
            pinhtag += '_' + states[istate].tuning
        endif

        pname = self.out_dir+'/pinhs/' + pinhtag + '.pinh'
        if arg_present(pinhname) then begin
            pinhname[istate] = pname
        endif

        if file_test(pname) then begin
            if arg_present(pinhdata) then begin
                pinhdata[0, 0, istate] = red_readdata(pname, header = pinhhead $
                                                     , status = pinhstatus, /silent)
                if status eq 0 then status = pinhstatus
            endif else if arg_present(pinhhead) then begin
                pinhhead[0, istate] = red_readhead(pname, status = pinhstatus, /silent)
                if status eq 0 then status = pinhstatus
            endif
        endif else begin
            if( arg_present(pinhdata) || arg_present(pinhhead) ) then status = -1
        endelse
        
     endif                      ; Pinholes

  endfor                        ; istate

  ;; Reduce dimensions if possible
  if Nstates eq 1 then begin

     if arg_present(darkname)  then darkname = darkname[0]
     if arg_present(flatname)  then flatname = flatname[0]
     if arg_present(gainname)  then gainname = gainname[0]
     if arg_present(pinhname)  then pinhname = pinhname[0]
     if arg_present(sflatname) then sflatname = sflatname[0]

     if arg_present(darkdata)  then darkdata = darkdata[*, *, 0]
     if arg_present(flatdata)  then flatdata = flatdata[*, *, 0]
     if arg_present(GAINdata)  then gaindata = gaindata[*, *, 0]
     if arg_present(pinhdata)  then pinhdata = pinhdata[*, *, 0] 
     if arg_present(sflatdata) then sflatdata = sflatdata[*, *, 0]
     
     if arg_present(darkhead)  then darkhead = darkhead[*, 0]
     if arg_present(flathead)  then flathead = flathead[*, 0]
     if arg_present(gainhead)  then gainhead = gainhead[*, 0]
     if arg_present(pinhhead)  then pinhhead = pinhhead[*, 0]
     if arg_present(sflathead) then sflathead = sflathead[*, 0]
     
  endif

end
