; docformat = 'rst'

;+
; Setup a momfbd manager and slaves.
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
;   port : in, out, type=integer
;
;     A port number. If given, all actions are for a momfbd manager
;     running with -p PORT.
;
;   nslaves : in, optional, type=integer, default=0
;
;     Start or stop slaves so this number of slaves are running.
;
;   nmanagers : in, optional, type=boolean, default=0
;
;     If set, make sure a manager is running. If not set, stop any
;     running manager.
;
;   free : in, optional, type=boolean, default=0
;
;     If set, find a free port, use it and return in the PORT keyword.
;
;   nthreads : in, optional, type=integer, default=1
;
;     If starting new slaves, start them with this number of threads.
;
; :History:
;
;   2010-09-21 : Written by Mats Löfdahl
;
;   2012-11-21 : Switched to docformat 'rst'.
;
;   2013-08-30 : MGL. Renamed for inclusion in crispred pipeline.
;
;   2013-09-04 : MGL. Use red_momfbd_check, not momfbd_check.
;
;   2014-01-24 : TH. Added HOST keyword to allow for using (already running)
;                remote managers. 
;
;-
pro red_momfbd_setup, HOST = host, $
                      PORT = port, $
                      NSLAVES = nslaves, $
                      NMANAGERS = nmanagers, $
                      FREE = free, $
                      NTHREADS = nthreads
  
  hoststring=''
  portstring=''
  if keyword_set(host) then begin
      if n_elements(port) eq 0 then begin
          print, 'You have to specify a port when using a remote manager.'
          return
      end
      hoststring='-m '+host
  end

  ; find a free port to use
  if n_elements(free) ne 0 then begin
     port = red_freeport(port=32765)
  end

  if n_elements(nthreads) eq 0 then nthreads=1
     
  red_momfbd_check, HOST=host, PORT=port, NSLAVES=nslavesrunning, NMANAGERS=nmanagersrunning, MANAGERPID=managerpid, SLAVEIDS=slaveids

  ; if the check didn't return a port-number, find a free one
  if n_elements(port) eq 0 then begin
     port = red_freeport(port=32765)
  end
  
  portstring='-p '+strtrim(string(port),2)
  
  if n_elements(nslavesrunning) ne 0 then begin
      if n_elements(nslaves) eq 0 then nslaves = nslavesrunning
  end else begin
      if n_elements(nslaves) eq 0 then nslaves = 0
  end

  if keyword_set(nmanagers) then begin ; We want a manager running
     
     if nmanagersrunning eq 0 then begin
        print,'Start a manager'
        if n_elements(port) ne 0 then print, '  on port ', port
        spawn,'nohup manager -v '+portstring+' &'
        hoststring=''
     end

     if nslaves gt nslavesrunning then begin
        print,'Start ',nslaves - nslavesrunning,' slaves with ',nthreads,' threads each.'
        for i=0,nslaves - nslavesrunning - 1 do begin
           spawn, 'momfbd_slave '+hoststring+' '+portstring+' -n '+strtrim(string(nthreads), 2)+' &'
        end
     end else if nslaves lt nslavesrunning then begin
        print,'Stop ',nslavesrunning - nslaves,' slaves.'
        print,'Slaveids:',slaveids
        for i=0,nslavesrunning - nslaves - 1 do begin
           ;; Kill from end of list so slave IDs don't grow unnecessarily
           print,'Kill slave # ',slaveids[nslavesrunning-1-i]
           spawn, 'sdel '+hoststring+' '+portstring+' '+slaveids[nslavesrunning-1-i]
        end
     end

  end else begin                ; We don't want a manager running

     if n_elements(nslavesrunning) ne 0 and nslavesrunning gt 0 then begin
        print,'Stop ',nslavesrunning,' slaves'
        print,'Slaveids:',slaveids
        for i=0,nslavesrunning - 1 do begin
           print,'Kill slave # ',slaveids[i]
           spawn, 'sdel '+portstring+' '+slaveids[i]
        end
     end

     if n_elements(nmanagersrunning) ne 0 and nmanagersrunning ne 0 and n_elements(managerpid) ne 0 then begin
        print,'Kill the manager with PID=',managerpid
        spawn,'kill '+managerpid
     end

  end

end
