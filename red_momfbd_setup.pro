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
;-
pro red_momfbd_setup, port = port, nslaves = nslaves, nmanagers=nmanagers, free=free, nthreads=nthreads
  
  if keyword_set(free) then begin
     ;; Find a free port
     spawn,'netstat -atn',netstat
     Nn = dimen(netstat, 0)
     ports = lonarr(Nn)
     for i = 2, Nn-1 do ports[i] = long(last(strsplit(token(netstat[i],4),':',/extract)))
     ports = ports[uniq(ports,sort(ports))] ; Uniquify
     freeport = 32765
     while inset(freeport, ports) do freeport += -1 
     port = strtrim(string(freeport), 2)
  end

  if n_elements(nslaves) eq 0 then nslaves=0
  if n_elements(nthreads) eq 0 then nthreads=1
     
  if n_elements(port) eq 0 then begin
     momfbd_check, NSLAVES = nslavesrunning, NMANAGERS=nmanagersrunning,MANAGERPID=managerpid, SLAVEIDS=slaveids
     portstring=' '
  end else begin
     momfbd_check, PORT = port, NSLAVES = nslavesrunning, NMANAGERS=nmanagersrunning,MANAGERPID=managerpid, SLAVEIDS=slaveids
     portstring=' -p '+strtrim(string(port),2)+' '
  end

  if keyword_set(nmanagers) then begin ; We want a manager running
     
     if nmanagersrunning eq 0 then begin
        print,'Start a manager'
        if n_elements(port) ne 0 then print, '  on port ', port
        spawn,'nohup manager -v '+portstring+' &'
     end

     if nslaves gt nslavesrunning then begin
        print,'Start ',nslaves - nslavesrunning,' slaves with ',nthreads,' threads each.'
        for i=0,nslaves - nslavesrunning - 1 do begin
           spawn, 'momfbd_slave '+portstring+' -n '+strtrim(string(nthreads), 2)+' &'
        end
     end else if nslaves lt nslavesrunning then begin
        print,'Stop ',nslavesrunning - nslaves,' slaves.'
        print,'Slaveids:',slaveids
        for i=0,nslavesrunning - nslaves - 1 do begin
           ;; Kill from end of list so slave IDs don't grow unnecessarily
           print,'Kill slave # ',slaveids[nslavesrunning-1-i]
           spawn, 'sdel '+portstring+' '+slaveids[nslavesrunning-1-i]
        end
     end

  end else begin                ; We don't want a manager running

     if nslavesrunning gt 0 then begin
        print,'Stop ',nslavesrunning,' slaves'
        print,'Slaveids:',slaveids
        for i=0,nslavesrunning - 1 do begin
           print,'Kill slave # ',slaveids[i]
           spawn, 'sdel '+portstring+' '+slaveids[i]
        end
     end

     if nmanagers eq 0 then begin
        print,'Kill the manager with PID=',managerpid
        spawn,'kill '+managerpid
     end

  end

  print

end
