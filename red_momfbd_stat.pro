; docformat = 'rst'

;+
; MOMFBD info from rdx_stat.
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
; :Returns:
; 
; 
; :Params:
; 
;   host : in, optional, type=string, default=localhost
; 
; 
; 
;   port : in, optional, type=string
; 
; 
; 
; :Keywords:
; 
;   status : out, optional
;   
;     0 for success, 1 for failure (no manager running).
; 
;   makestruct : in, optional, type=boolean
; 
;     Return the info as a struct.
;
; :History:
; 
;   2025-02-07 : MGL. First version.
; 
;-
function red_momfbd_stat, host, port, status = status, struct = makestruct

  if n_elements(host) eq 0 then host = 'localhost'
  hoststring=' -m '+host

  if n_elements(port) then portstring=' -p '+strtrim(port, 2) else portstring = ''

  spawn, 'rdx_stat -j ' + hoststring + portstring, jstat
  spawn, 'rdx_stat -s ' + hoststring + portstring, sstat
  
  if n_elements(sstat) eq 1 && strmatch(output, 'Connection failed*') then begin
    status = 1                  ; Failed
    return, ''
  endif

  status = 0

  if ~keyword_set(makestruct) then return, [sstat, jstat]

  ;; Manager info
  splt = strsplit(sstat[1], /extract)
  ;; Note that the Runtime column (7) is empty when the status (8)
  ;; is idle. Strsplit messes this up so we need to parse this
  ;; carefully.
  if splt[7] eq 'idle' then begin
    runtime = ''
    slave_status = 'idle'
  endif else begin
    runtime = splt[7]
    slave_status = splt[8]
  endelse

  manager_info = {host:splt[1] $
                  , pid:splt[2] $
                  , threads:splt[3] $
                  , version:splt[4] $
                  , uptime:splt[6] $
                  , runtme:runtime $
                  , status:slave_status $
                 }
  
  return_struct = { manager:manager_info }
  
  ;; Slave info
  slave_strings = sstat[2:*]
  for islave = 0, n_elements(slave_strings)-1 do begin

    splt = strsplit(slave_strings[islave], /extract)
    
    ;; Note that the Runtime column (8) is empty when the status (9)
    ;; is idle. Strsplit messes this up so we need to parse this
    ;; carefully.
    if splt[8] eq 'idle' then begin
      runtime = ''
      slave_status = 'idle'
    endif else begin
      runtime = splt[8]
      slave_status = splt[9]
    endelse

    red_append, slave_info, {id:splt[0] $
                             , host:splt[2] $
                             , pid:splt[3] $
                             , threads:splt[4] $
                             , version:splt[5] $
                             , uptime:splt[7] $
                             , runtme:runtime $
                             , status:slave_status $
                            }
    
  endfor                        ; islave

  return_struct = create_struct(return_struct, 'slaves', slave_info)
  


  ;; Job info
  if n_elements(jstat) lt 2 then begin
    id    = ''
    time  = ''
    name  = ''
    user  = ''
    prio  = ''
    state = ''
    job_info = {id:id $
                , time:time $
                , name:name $
                , user:user $
                , prio:prio $
                , state:state $
               }              
  endif else begin
    job_strings = jstat[1:*]
    for ijob = 0, n_elements(job_strings)-1 do begin
      splt = strsplit(job_strings[ijob], /extract)
      id    = splt[1] 
      time  = splt[3] 
      name  = splt[4] 
      user  = splt[5] 
      prio  = splt[6] 
      if n_elements(splt) ge 8 then state = strjoin(splt[7:*], ' ') else state = ''
      red_append, job_info, {id:id $
                             , time:time $
                             , name:name $
                             , user:user $
                             , prio:prio $
                             , state:state $
                            }
    endfor                      ; ijob
  endelse

  return_struct = create_struct(return_struct, 'jobs', job_info)
  

  ;; Return the struct output
  return, return_struct
  
end

port = '31000'
host = 'freija'

out = red_momfbd_stat( host, port, status = status, /struct )

end
