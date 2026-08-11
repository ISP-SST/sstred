; docformat = 'rst'

;+
; Find version info about sstred and other libraries.
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
; :Keywords:
; 
;    branch_pipeline, out, optional, type=string
;   
;      Git repo branch of the SSTRED pipeline library.
; 
;    version_coyote, out, optional, type=string
;   
;      Version of the coyote library.
; 
;    version_idlastro, out, optional, type=string
;   
;      Version of the coyote library.
; 
;    version_mpfit, out, optional, type=string
;   
;      Version of the coyote library.
; 
;    version_pipeline, out, optional, type=string
;   
;      Version of the SSTRED pipeline itself.
; 
;    version_problems, out, optional, type=string
;   
;      Problems with the installation.
; 
;    version_reduxdlm, out, optional, type=string
;   
;      Version of the redux DLMs.
; 
; 
; :History:
; 
;    2026-07-01 : MGL. First version.
; 
;-
pro red_version_info, branch_pipeline    = version_pipeline_branch $
                      , version_coyote   = version_coyote $
                      , version_idlastro = version_idlastro $
                      , version_mpfit    = version_mpfit $
                      , version_pipeline = version_pipeline $
                      , version_problems = version_problems $
                      , version_reduxdlm = version_reduxdlm

  paths = strsplit(!path,":",/extract) 
  git_describe_command = 'git describe --always --abbrev=12 --long --dirty=\ \(Modified\)'
  git_count_command = 'git rev-list HEAD --count'
  git_diff_command = 'git diff HEAD'
  git_status_command = 'git status'

  
  version_problems = ''         ; Initialize

  ;; Pipeline version
  srcdir = file_dirname( routine_filepath("red::initialize"), /mark )
  if srcdir eq './' then begin
    indx = where(file_basename(paths) eq 'sstred', Nwhere)
    if Nwhere eq 0 then stop
    srcdir = paths[indx[0]]
  endif
  spawn, 'cd '+srcdir+'; ' + git_describe_command, pipeline_gitoutput
  pipeline_gitoutput = red_strreplace(pipeline_gitoutput, 'release/', '')
  spawn, 'cd '+srcdir+'; ' + git_status_command, pipeline_status 
  if strmatch(pipeline_gitoutput, '*(Modified)') then $
     version_problems += 'The pipeline is modified. '
  version_pipeline = strjoin((strsplit(pipeline_gitoutput, '-', /extract))[0:1], '-')

  ;; Pipeline branch
  spawn, 'cd '+srcdir+'; ' + git_status_command, pipeline_status  
  version_pipeline_branch = (strsplit(pipeline_status[0],/extract))[-1]
  if strmatch(version_pipeline_branch,'*_dev') then  version_problems += 'Running the development branch ' $
     + version_pipeline_branch
  
  ;; Redux dlm version. We require that the ANA and MOMFBD dlms are
  ;; part of the rdx dlm and the same version.
  help,/dlm,'rdx', output = rdx_dlm_version
  dlmpos = strpos(rdx_dlm_version[1], 'release/')
  version_reduxdlm = strjoin((strsplit(strmid(rdx_dlm_version[1],dlmpos+8), '-', /extract))[0:1], '-')
  help,/dlm,'ana' , output = ana_dlm_version
  dlmpos = strpos(ana_dlm_version[1], 'release/')
  ana_dlm_version = strjoin((strsplit(strmid(ana_dlm_version[1],dlmpos+8), '-', /extract))[0:1], '-')
  help,/dlm,'momfbd', output = momfbd_dlm_version
  dlmpos = strpos(momfbd_dlm_version[1], 'release/')
  momfbd_dlm_version = strjoin((strsplit(strmid(momfbd_dlm_version[1],dlmpos+8), '-', /extract))[0:1], '-')

  if ana_dlm_version ne version_reduxdlm then version_problems += 'ANA DLM not identical to redux DLM'
  if momfbd_dlm_version ne version_reduxdlm then version_problems += 'MOMFBD DLM not identical to redux DLM'
  
  rdx_dlm_required_versionstr = '1.0.0-1'
  rdx_dlm_required_version = fix(strsplit(rdx_dlm_required_versionstr, '-.', /EXTRACT))
  rdx, version=rdx_dlm_version
  rdx_dlm_version = fix(strsplit(rdx_dlm_version, '-.', /EXTRACT))
  n_elem = min([n_elements(rdx_dlm_required_version),n_elements(rdx_dlm_version)])
  for i=0,n_elem-1 do begin
    if rdx_dlm_version[i] gt rdx_dlm_required_version[i] then break
    if rdx_dlm_version[i] lt rdx_dlm_required_version[i] then begin
      version_problems += 'The RDX DLM should be updated to version >= ' + rdx_dlm_required_versionstr
      break
    endif
  endfor

  ;; Coyote library version
  coyotepaths = paths(where(strmatch(paths,'*coyote'), Ncoyote))
  case Ncoyote of
    0: begin
      print, 'The Coyote library does not seem to be in your !path.'
      stop
    end
    1: begin
      ;coyotedir = file_dirname( filepath(root_dir = coyotepaths[0], "cgcolor"), /mark )
      spawn, 'cd '+coyotepaths[0]+'; ' + git_count_command, coyote_gitoutput
      coyote_gitoutput = red_strreplace(coyote_gitoutput, 'release/', '')
      ;; if strmatch(coyote_gitoutput, '(Modified)') then $
      ;;   version_problems += 'The Coyote library is modified. '
      ;; version_coyote = strjoin((strsplit(coyote_gitoutput, '-', /extract))[0:1], '-')
      version_coyote = coyote_gitoutput
      spawn, 'cd '+coyotepaths[0]+'; ' + git_diff_command, coyote_gitoutput
      if size(coyote_gitoutput, /n_dim) gt 0 then version_problems $
         += 'The Coyote library is not the latest version. '
    end
    else: begin
      version_coyote = 'Undefined'
      version_problems += 'Multiple Coyote directories. '
    end
  endcase


  ;; IDLastro library version
  idlastropaths = paths(where(strmatch(paths, '*IDLAstro/pro'), Nwhere))
  case Nwhere of
    0: begin
      print, 'The IDLAstro library does not seem to be in your !path.'
      stop
    end
    1: begin
      spawn, 'cd '+idlastropaths[0]+'; ' + git_count_command, idlastro_gitoutput
      idlastro_gitoutput = red_strreplace(idlastro_gitoutput, 'release/', '')
      ;; if strmatch(idlastro_gitoutput, '(Modified)') then $
      ;;  version_problems += 'The IDLAstro library is modified. '
      ;; version_idlastro = strjoin((strsplit(idlastro_gitoutput, '-', /extract))[0:1], '-')
      version_idlastro = idlastro_gitoutput
      spawn, 'cd '+idlastropaths[0]+'; ' + git_diff_command, idlastro_gitoutput
      if size(idlastro_gitoutput, /n_dim) gt 0 then version_problems $
         += 'The IDLAstro library is not the latest version. '
    end
    else: begin
      version_idlastro = 'Undefined.'
      version_problems += 'Multiple IDLAstro directories. '
    end
  endcase


  ;; mpfit version  
  mpfitpaths = paths(where(strmatch(paths, '*mpfit'), Nwhere))
  case Nwhere of
    0: begin
      print, 'The mpfit library does not seem to be in your !path.'
      stop
    end
    1: begin
      red_mpfit_version, mpfitpaths, local_version = local_version
      version_mpfit = local_version 
    end
    else: begin
      version_mpfit = 'Undefined'
      version_problems += 'Multiple mpfit directories. '
    end
  endcase

end

red_version_info, branch_pipeline = version_pipeline_branch $
                  , version_coyote = version_coyote $
                  , version_idlastro = version_idlastro $
                  , version_mpfit = version_mpfit $
                  , version_pipeline = version_pipeline $
                  , version_problems = version_problems $
                  , version_reduxdlm = version_reduxdlm

print, 'Coyote: ', version_coyote
print, 'IDLAstro: ', version_idlastro
print, 'Mpfit: ', version_mpfit
print, 'SSTRED version: ', version_pipeline
print, 'SSTRED branch: ', version_pipeline_branch
print, 'redux DLMs: ', version_reduxdlm

print, 'Problems: ', version_problems

end
