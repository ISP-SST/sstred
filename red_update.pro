; docformat = 'rst'

;+
; Update crispred and the libraries that it depends on.
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
; 
; 
; 
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;    2017-04-26 : MGL. First version. Does the IDL libraries but not
;                 the rdx DLM.
; 
;    2017-04-27 : MGL. Check if the rdx DLM is the latest version.
; 
;    2017-05-04 : MGL. Do make in the creduc/ subdirectory.
; 
;    2017-05-17 : MGL. Exit IDL after updating.
; 
;    2017-05-18 : MGL. Don't require specific directory names for the
;                 extra libraries. Exit only if something was actually
;                 updated. 
; 
;    
; 
; 
;-
pro red_update

  paths = strsplit(!path,":",/extract) 
  
  print, 'red_update : Update the current crispred branch.'
  srcdir = file_dirname( routine_filepath("red_update"), /mark )
  spawn, 'cd '+srcdir+'; git pull', spawnoutput, /stderr
  print, spawnoutput, format = '(a0)'
  if max(strmatch(spawnoutput, 'Already up-to-date.')) eq 0 then doexit = 1
  spawn, 'cd '+srcdir+'/creduc ; make'
  
  print, 'red_update : Update the coyote library.'
  coyotepaths = file_dirname( File_Search(paths+'/cgcolor.pro', COUNT=Nwhere), /mark )
  case Nwhere of
    0: print, 'The Coyote library does not seem to be installed.'
    1: begin
      spawn, 'cd '+coyotepaths+'; git pull', spawnoutput, /stderr
      print, spawnoutput, format = '(a0)'
      if max(strmatch(spawnoutput, 'Already up-to-date.')) eq 0 then doexit = 1
      if strmatch(spawnoutput[0],'fatal*') then $
         print, 'red_update : The Coyote library not updated, does not seem to be under git control.'
    end
    else: begin
      print, 'red_update : Multiple Coyote directories.'
      print, '             Please make sure there is only one and that it is under git control.'
      doexit = 1
    end
  endcase

  print, 'red_update : Update the IDLAstro library.'
  idlastropaths = file_dirname( File_Search(paths+'/astrolib.pro', COUNT=Nwhere), /mark )
  case Nwhere of
    0: print, 'The IDLAstro library does not seem to be installed.'
    1: begin
      spawn, 'cd '+idlastropaths+'; git pull', spawnoutput, /stderr
      print, spawnoutput, format = '(a0)'
      if max(strmatch(spawnoutput, 'Already up-to-date.')) eq 0 then doexit = 1
      if strmatch(spawnoutput[0],'fatal*') then $
         print, 'red_update : The IDLAstro library not updated, does not seem to be under git control.'
    end
    else: begin
      print, 'red_update : Multiple IDLAstro directories.'
      print, '             Please make sure there is only one and that it is under git control.'
      doexit = 1
    end
  endcase

  print, 'red_update : Update the mpfit library.'
  mpfitpaths = file_dirname( File_Search(paths+'/mpfit.pro', COUNT=Nwhere), /mark )
  case Nwhere of
    0: print, 'The mpfit library does not seem to be installed.'
    1: begin
      red_mpfit_version, mpfitpaths, local_version = local_version, latest_version = latest_version
      print, 'MPFIT local version: ',  local_version
      print, 'MPFIT latest version: ',  latest_version
      if latest_version gt local_version then begin
        cmd = strjoin(['cd '+mpfitpaths $
                       , 'rm -f mpfit.tar.gz' $
                       , 'wget --quiet http://cow.physics.wisc.edu/~craigm/idl/down/mpfit.tar.gz' $
                       , 'tar xvzf mpfit.tar.gz' $
                      ], ' ; ')
        spawn, cmd
        doexit = 1
      endif
    end
    else: begin
      print, 'red_update : Multiple mpfit directories.'
      print, '             Please make sure there is only one.'
      doexit = 1
    end
  endcase

  print, 'red_update : Check version of rdx DLMs.'
  ;; The rdx DLM knows its own git version. However, we can't go to
  ;; the source git repository because the !dlm_path variable only
  ;; tells us where the compiled DLMs are, so we can't use git locally
  ;; to find out if it's the latest version. However, the repo on
  ;; dubshen will always be the latest version, so we can use
  ;; git-ls-remote to find out its version.
  spawn, 'git ls-remote -h git://dubshen.astro.su.se/hillberg/redux master', rdx_latest_hash
  rdx_latest_hash = strmid(rdx_latest_hash, 0, 12)
  help,/dlm,'rdx', output = rdx_dlm_version
  dlmpos = strpos(rdx_dlm_version[1], '-g')  
  rdx_local_hash = strmid(rdx_dlm_version[1],dlmpos+2,12)

  print
  print, 'red_update : IDL libraries updated (if needed and possible).'
  print
  if rdx_latest_hash eq rdx_local_hash then begin
    print, 'red_update : The rdx DLMs are up to date.' 
  endif else begin
    print, 'red_update : Local and latest git hashes for the rdx DLMs:'
    print, '             Local  '+rdx_local_hash
    print, '             Latest '+rdx_latest_hash
    print, '             You may want to update the rdx DLMs (git pull, compile, and install).' 
    print, '             This might mean that also your momfbd-redux code needs to be updated.'
  endelse

  if keyword_set(doexit) then begin
    ;; Exit if there was something that was actually updated, so user
    ;; will start with a fresh session.
    print
    print, 'Will exit IDL now so you can start with the fresh versions of the software.'
    print
    exit
  endif
  
end
