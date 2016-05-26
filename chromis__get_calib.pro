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
; 
;-
pro chromis::get_calib, states $
                        , status = status $
                        , darkname = darkname $
                        , darkdata = darkdata $
                        , darkhead = darkhead $
                        , flatname = flatname $
                        , flatdata = flatdata $
                        , flathead = flathead $
                        , sflatname = sflatname $
                        , sflatdata = sflatdata $
                        , sflathead = sflathead $
                        , pinhname = pinhname $
                        , pinhdata = pinhdata $
                        , pinhhead = pinhhead 

  Nstates = n_elements(states)

  if Nstates eq 0 then begin
     status = -1
     return
  endif

  if arg_present(darkname) then darkname = strarr(Nstates)   
  if arg_present(flatname) then flatname = strarr(Nstates) 
  if arg_present(pinhname) then pinhname = strarr(Nstates) 
  if arg_present(sflatname) then sflatname = strarr(Nstates) 

  if arg_present(darkdata) $
     or arg_present(flatdata) $
     or arg_present(sflatdata) $
     or arg_present(pinhdata) then begin

     ;; Assume this is all for the same camera type, at least for the
     ;; actual data. Otherwise we cannot return the actual data in a
     ;; single array.
     camtag = states[0].camtag
     caminfo = red_camerainfo(camtag)

     if arg_present(darkdata) then darkdata = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
     if arg_present(flatdata) then flatdata = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
     if arg_present(pinhdata) then pinhdata = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
     if arg_present(sflatdata) then sflatdata = fltarr(caminfo.xsize, caminfo.ysize, Nstates) 
  endif

  Nheadlines = 100              ; Assume max numer of header lines
  if arg_present(darkhead) then darkhead = arr(Nheadlines,Nstates) 
  if arg_present(flathead) then flathead = arr(Nheadlines,Nstates) 
  if arg_present(pinhhead) then pinhhead = arr(Nheadlines,Nstates) 
  if arg_present(sflathead) then sflathead = arr(Nheadlines,Nstates) 

  status = 0

  for istate = 0, Nstates-1 do begin

     camtag = states[istate].camtag

     ;; Darks

     if arg_present(darkname) or arg_present(darkdata) or arg_present(darkhead) then begin

        darktag = camtag $
                  + '_' + string(states[istate].exposure*1000, format = '(f4.2)') + 'ms' $
                  + '_' + 'G' + string(states[istate].gain, format = '(f05.2)')
        dname = self.out_dir+'/darks/' + darktag + '.dark'
        darkname[istate] = dname

        if n_elements(dname) ne 0 then begin
           
           if arg_present(darkdata) then begin
              darkdata[0, 0, istate] = red_readdata(dname, header = darkhead, status = darkstatus)
              if status eq 0 then status = darkstatus
           endif else if arg_present(darkhead) then begin
              darkhead[0, istate] = red_readhead(dname, status = darkstatus)
              if status eq 0 then status = darkstatus
           endif

        endif else status = -1

     endif

     ;; Flats

     if arg_present(flatname) or arg_present(flatdata) or arg_present(flathead) $
        or arg_present(sflatname) or arg_present(sflatdata) or arg_present(sflathead) then begin

        flattag = camtag $
                  + '_' + string(states[istate].exposure*1000, format = '(f4.2)') + 'ms' $
                  + '_' + 'G' + string(gain, format = '(f05.2)') $
                  + '_' + states[istate].fullstate
        fname = self.out_dir+'/flats/' + flattag + '.flat'
        flatname[istate] = fname

        if n_elements(fname) ne 0 then begin
           
           if arg_present(flatdata) then begin
              flatdata[0, 0, istate] = red_readdata(fname, header = flathead, status = flatstatus)
              if status eq 0 then status = flatstatus
           endif else if arg_present(flathead) then begin
              flathead[0, istate] = red_readhead(fname, status = flatstatus)
              if status eq 0 then status = flatstatus
           endif

        endif else status = -1

        sfname = self.out_dir+'/flats/' + flattag + '_summed.flat'
        sflatname[istate] = sfname

        if n_elements(sfname) ne 0 then begin
           
           if arg_present(sflatdata) then begin
              sflatdata[0, 0, istate] = red_readdata(sfname, header = sflathead, status = sflatstatus)
              if status eq 0 then status = sflatstatus
           endif else if arg_present(flathead) then begin
              sflathead[0, istate] = red_readhead(sfname, status = sflatstatus)
              if status eq 0 then status = sflatstatus
           endif

        endif else status = -1

     endif

     ;; Pinholes

     if arg_present(pinhname) or arg_present(pinhdata) or arg_present(pinhhead) then begin

        pinhtag = camtag + '_' + states[istate].fullstate
        pname = self.out_dir+'/pinhs/' + pinhtag + '.pinh'
        pinhname[istate] = pname

        if n_elements(pname) ne 0 then begin
           
           if arg_present(pinhdata) then begin
              pinhdata[0, 0, istate] = red_readdata(pname, header = pinhhead, status = pinhstatus)
              if status eq 0 then status = pinhstatus
           endif else if arg_present(pinhhead) then begin
              pinhhead[0, istate] = red_readhead(pname, status = pinhstatus)
              if status eq 0 then status = pinhstatus
           endif

        endif else status = -1
        
     endif

  endfor                        ; istate

  ;; Reduce dimensions if possible
  if Nstates eq 1 then begin

     if arg_present(darkname) then darkname = darkname[0]
     if arg_present(flatname) then flatname = flatname[0]
     if arg_present(pinhname) then pinhname = pinhname[0]

     if arg_present(darkdata) then darkdata = darkhead[*, *, 0]
     if arg_present(flatdata) then flatdata = flathead[*, *, 0]
     if arg_present(pinhdata) then pinhdata = pinhhead[*, *, 0] 
     
     if arg_present(darkhead) then darkhead = darkhead[*, 0]
     if arg_present(flathead) then flathead = flathead[*, 0]
     if arg_present(pinhhead) then pinhhead = pinhhead[*, 0]
     
  endif

end