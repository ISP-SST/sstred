; docformat = 'rst'

;+
; Return various info about a momfbd manager and its slaves.
;
; :Categories:
;
;   Image restoration, mombfd
;
; :Author:
;
;   Mats Löfdahl, Institute for Solar Physics, mats@astro.su.se.
;
; :Keywords:
;
;   port : in, optional, type=integer
;
;     A port number. If given, all output is for a momfbd manager
;     running with -p PORT.
;
;   nslaves : out, optional, type=integer
;
;     The number of momfbd_slave processes.
;
;   njobs : out, optional, type=integer
;
;     The number of momfbd jobs in the queue.
;
;   nmanagers : out, optional, type=integer
;
;     The number of managers running.
;
;   managerpid : out, optional, type=string
;
;     The PID of the manager process as a string.
;
;   slaveids : out, optional, type=array
;
;     Array with momfbd slave ids.
;
;   jobids : out, optional, type=array
;
;     Array with momfbd job ids.
;
; :History:
;
;   2010-09-21: Written by Mats Löfdahl
;
;   2012-11-21 : Switched to docformat 'rst'.
;
;   2013-08-30 : MGL. Renamed for inclusion in crispred pipeline.
;
;   2014-01-24 : TH. Added HOST keyword to allow for using (already running)
;                remote managers.
;
;-
pro red_momfbd_check, HOST = host $
                    , PORT = port $
                    , NSLAVES = nslaves $
                    , NJOBS = njobs $
                    , NMANAGERS = nmanagers $
                    , MANAGERPID = managerpid $
                    , SLAVEIDS = slaveids $
                    , JOBIDS = jobids

  hoststring=''
  portstring=''
  
  if n_elements(host) ne 0 then begin   ; use a remote manager
      hoststring=' -m '+host
      if n_elements(port) ne 0 then begin
          portstring=' -p '+strtrim(string(port),2)
      end
  end else begin        ; no host given -> check for local manager
      info = get_login_info()
      spawn, 'ps aux | grep "[[:space:]]manager[[:space:]]" | grep '+info.user_name, psout
      if n_elements(port) ne 0 then begin
          indices = where(strmatch(psout,'*-p '+strtrim(string(port),2)+'*') eq 1, count)
          if count gt 0 then begin
              psout = psout[indices]
              portstring=' -p '+strtrim(string(port),2)
          end else psout=''     ; no match to specified port, zero out the list
      endif
      ;help, psout
      ;print, psout
      
      nmanagers = (size(psout, /dim))[0] gt 0
      ;print, nmanagers

      if nmanagers eq 0 then begin
         print,'No local manager found'
         nslaves=0
         njobs=0
         slaveids=0
         return
      end 
      
      ;; Manager
      managerpid = (token(psout[0]))[1]
      if n_elements(port) eq 0 then begin        ; no port specified, get it from the pid
          spawn, 'netstat -atnpl | grep LISTEN | grep '+managerpid+'/', netstat
          if strmatch(netstat,'') eq 0 then begin
              port = long((strsplit((strsplit(netstat,/EXTRACT))[3],':',/EXTRACT))[1])      ; returns the port number if an empty port variable was provided
              portstring = ' -p '+strtrim(string(port),2)
          end
      end


  end

     ;; Do we really want to return output?

     no_args = ~(arg_present(nslaves) $ 
                 or arg_present(slaveids) $
                 or arg_present(njobs) $
                 or arg_present(jobids) $
                 or arg_present(managerpid) $
                 or arg_present(nmanagers) $
                )

     ;; Slaves
     if no_args or arg_present(nslaves) or arg_present(slaveids) then begin
        spawn, 'jstat '+hoststring+' '+portstring+' -s', jstats
     end

     ;if arg_present(nslaves) then
     nslaves = (size(jstats,/dim))[0] - 2
     nmanagers = 1    ; should be set if the reply is ok, popssibly use a timeout-check with spawn ??
     
     if arg_present(slaveids) and  arg_present(nslaves) then begin
        if nslaves gt 0 then begin
           slaveids = strarr(nslaves)
           for i=0,nslaves-1 do slaveids[i] = (token(jstats[i+2]))[1]
        endif
     endif

     ;; Jobs

     if no_args or arg_present(njobs) or arg_present(jobids) then begin
        spawn, 'jstat '+hoststring+' '+portstring+' -j', jstatj
     end

     if arg_present(njobs) then njobs = (size(jstatj,/dim))[0] - 1

     if arg_present(jobids) then begin     
        jobids = strarr(njobs)
        for i=0,njobs-1 do jobids[i] = (token(jstatj[i+2]))[1]
     end

     if no_args then begin
        ;; We apparently didn't want any data returned. Let's print
        ;; something out instead.
        print,'--------------------------MOMFBD SLAVES-----------------------------------------------------'
        print,jstats
        print,'--------------------------MOMFBD JOB QUEUE--------------------------------------------------'
        print,jstatj
        print,'--------------------------------------------------------------------------------------------'
     end


end
