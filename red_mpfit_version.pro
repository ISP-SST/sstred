; docformat = 'rst'

;+
; Find version information for mpfit.
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
; 
; :Params:
; 
;   mpfitdir : in, optional, type=string, default='./'
; 
;     The directory where the mpfit is installed locally.
; 
; 
; :Keywords:
; 
;   local_version : out, optional, type=string
;
;     The version number (date) of the locally installed mpfit library.
;
;   latest_version : out, optional, type=string
;   
;     The version number (date) of the latest version available from
;     the download site.
; 
; 
; :History:
; 
;   2017-04-26 : MGL. First version. 
; 
; 
; 
; 
;-
pro red_mpfit_version, mpfitdir, local_version = local_version, latest_version = latest_version

  if arg_present(local_version) then begin
    
    if n_elements(mpfitdir) eq 0 then mpfitdir = './'
  
    ;; As the local version, we use the latest date within any of the
    ;; files.
    
    spawn, 'grep "\$Id" '+mpfitdir+'/*.pro', mpfit_spawnoutput
    timestamps = STREGEX(mpfit_spawnoutput,'[0-9][0-9][0-9][0-9]/[0-9][0-9]/[0-9][0-9]' $
                         ,/EXTRACT)
    for i = 0, n_elements(timestamps)-1 do timestamps[i] = red_strreplace(timestamps[i], '/', '-', n = 2)
    for i = 0, n_elements(timestamps)-1 do timestamps[i] = red_strreplace(timestamps[i], ' ', 'T')
    local_version = (timestamps(sort(timestamps)))[n_elements(timestamps)-1]
  endif

  if arg_present(latest_version) then begin

    ;; The latest version is given on the web page, next to the
    ;; archive file name.
    
    cmds = ['cd /tmp'  $
            , 'wget wget http://cow.physics.wisc.edu/~craigm/idl/fitting.html' $
            , 'grep mpfit.tar.gz fitting.html' $
           ]
    cmd = strjoin(cmds, ' ; ')
    
    spawn, cmd, spawnoutput
    timestamp = (stregex(spawnoutput, '[A-Z][a-z][a-z] [0-9][0-9] [0-9][0-9][0-9][0-9]', /extract))[1]
    latest_version = date_conv(strjoin((strsplit(strtrim(timestamp,2),' ',/extract))[[1,0,2]],'-'),'FITS')

  endif
  
end
