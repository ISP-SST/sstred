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
;-
pro red_momfbd_check, port = port $
                      , nslaves = nslaves $
                      , njobs=njobs $
                      , nmanagers=nmanagers $
                      , managerpid=managerpid $
                      , slaveids=slaveids $
                      , jobids=jobids
  
  if n_elements(port) eq 0 then begin
     spawn, 'ps aux | grep "[[:space:]]manager[[:space:]]" | grep -v "\-p"', psout
     portstring=' '
  endif else begin
     spawn, 'ps aux |grep "[[:space:]]manager[[:space:]].*-p[[:space:]]'+strtrim(string(PORT), 2)+'\$"', psout
     portstring=' -p '+strtrim(string(port),2)+' '
  endelse

;  help, psout
;  print, psout

  nmanagers = (size(psout, /dim))[0] ne 0

  if nmanagers eq 0 then begin

     print,'No manager is running'

     nslaves=0
     njobs=0
     slaveids=0

  end else begin

     ;; Do we really want to return output?

     no_args = ~(arg_present(nslaves) $ 
                 or arg_present(slaveids) $
                 or arg_present(njobs) $
                 or arg_present(jobids) $
                 or arg_present(managerpid) $
                 or arg_present(nmanagers) $
                )

     ;; Manager

     if arg_present(managerpid) then managerpid = (token(psout))[1]
     
     ;; Slaves

     if no_args or arg_present(nslaves) or arg_present(slaveids) then begin
        spawn, 'jstat  '+portstring+' -s', jstats
     end

     if arg_present(nslaves) then nslaves = (size(jstats,/dim))[0] - 2

     if arg_present(slaveids) and  arg_present(nslaves) then begin
        if nslaves gt 0 then begin
           slaveids = strarr(nslaves)
           for i=0,nslaves-1 do slaveids[i] = (token(jstats[i+2]))[1]
        endif
     endif

     ;; Jobs

     if no_args or arg_present(njobs) or arg_present(jobids) then begin
        spawn, 'jstat  '+portstring+' -j', jstatj
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

end
