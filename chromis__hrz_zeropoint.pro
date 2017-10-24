; docformat = 'rst'

;+
; Find the zero point of CHROMIS scans with wavelength given in HRE
; digital units.
;
; The zero points are stored in files to be used by the extractstates
; method. 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl
; 
; 
; 
; :History:
;
;    2016-09-18 : MGL. First version.
;
;    2016-09-21 : MGL. Work in meters, not Å or mÅ. Bugfix in
;                 extracting the hrz tunings. Typos hzr --> hrz.
;                 Remove previously generated hrz_zeropoint files.
;                 Change filter tags to four-digit tags representing
;                 approximate filter wavelength.
;    
;    2017-08-18 : THI. Use the FPI calibration values to get hrz instead
;                 of parsing files.
; 
;-
function range, first, last, step
  ; this is just a dummy function to prevent errors since some linedef.py uses range()
  return, 0
end

pro chromis::hrz_zeropoint

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  infodir = self.out_dir + 'info/'
  file_mkdir, infodir

  ;; Remove previously generated hrz_zeropoint files if any. 
  spawn, 'rm -f '+ infodir + 'hrz_zeropoint_*'

  ;; Get du_ref from the unique list of hrz tunings in raw flats taken
  ;; with the nb camera, list them in all flats directories.

  linedefdir = self.out_dir + 'downloads/'
  ldfiles = file_search(linedefdir + '/linedef.*', count = nLD)
  
  if nLD eq 1 then begin
    openr, lun, ldfiles[0], /GET_LUN
    while not EOF(lun) do begin
      line = ''
      READF, lun, line
      if strmid(line,0,1) ne '#' then begin
        line = red_strreplace( line, '[:]', '', n=100 )
        line = STRSPLIT( line, ';', /EXTRACT, count=nlines )
        for i=0,nlines-1 do begin
          assign_pos = strpos(line[i],'=')
          if assign_pos gt 0 then begin
            ok = execute( line[i] )
            varname = strtrim(strmid(line[i], 0, assign_pos), 2)
            if strpos(varname,'[') eq -1 then begin
              red_append, varnames, varname
            endif
          endif
        endfor
      endif
    endwhile
    free_lun, lun
  endif
  
  idx = where(strmatch(varnames, '*calib*', /FOLD_CASE) EQ 1, compl=cidx)
  varnames = varnames[cidx]               ; skip calib-parameters

  lambda_list = [ 0.0d, $                           ; 0-position in the wheel is open
                  3925.0d-10, $                     ; Ca II K blue wing
                  3933.7d-10, $                     ; Ca II K core
                  3968.5d-10, $                     ; Ca II H core
                  3978.0d-10, $                     ; Ca II H red wing
                  (3998.640d-10 + 1.258d-10), $     ; Clean reference line within the passband
                  4861.0d-10 $                      ; H-beta core
                  ]
  prefilter_list = [ '0000', '3925', '3934', '3969', '3978', '3999', '4862' ]
  disp_list = -1.53d-13*lambda_list/3998.6d-10

  linedef = { name:'', prefilter:'', fpi:lonarr(6), wheel:0b, wheel_tilt:0L, lambda_ref:0.0d, disp:0.0d }
  for i=0,n_elements(varnames)-1 do begin
    vardata = SCOPE_VARFETCH(varnames[i])
    if n_elements(vardata) ge 8 then begin      ;  this is a linedef
      linedef.name = varnames[i]
      linedef.fpi = vardata[0:5]
      linedef.wheel = vardata[6]
      linedef.wheel_tilt = vardata[7]
      if n_elements(vardata) eq 9 then linedef.disp = vardata[8]*1.0d-13 $
      else linedef.disp = 0.0
      if (linedef.wheel gt 0) && (linedef.wheel lt 7) then begin
        linedef.prefilter = prefilter_list[linedef.wheel]
        linedef.lambda_ref = lambda_list[linedef.wheel]
      endif else begin
        ;print, inam+' : Unrecognized linedef for wheel-position '+strtrim(string(linedef.wheel,/print),2)+' will be ignored.'
        continue
      endelse
      red_append, linedefs, linedef
    endif else begin        ; other variables we recognize
      case varnames[i] of
        ;'disp_cak': disp_list[1:2] = vardata           ; TBD: should we always override the auto-generated
        ;'disp_cah': disp_list[3:4] = vardata           ; disp_list with the value from the file?
        ;'disp_hb' : disp_list[6] = vardata
        ;'cfak' : do we need this for anything ?
        else: 
      endcase 
    endelse
  endfor
  
  cont_idx = where(linedefs.wheel eq 5)
  for i=0,n_elements(linedefs)-1 do begin
    if (linedefs[i].wheel eq 5) && (n_elements(cont_idx) gt 1) then begin
      ; if both cac_line & cac_cont are defined, pick the line (highest hrz)
      if linedefs[i].fpi[2] ne max(linedefs[cont_idx].fpi[2]) then continue  ; skip
    endif
    ;print, linedefs[i].name, linedefs[i].fpi, linedefs[i].wheel, linedefs[i].wheel_tilt, linedefs[i].lambda_ref, linedefs[i].disp
     zfile = infodir + 'hrz_zeropoint_' +  linedefs[i].prefilter + '.fz'
     if linedefs[i].disp eq 0.0d then begin ; disp not in linedef, get it from the list
       linedefs[i].disp = disp_list[linedefs[i].wheel]
     endif
     hrz_zero = linedefs[i].fpi[2]
     ;print,'new ', zfile, [linedefs[i].lambda_ref, hrz_zero, linedefs[i].disp]
     fzwrite, [linedefs[i].lambda_ref, hrz_zero, linedefs[i].disp], zfile, ' '
  endfor

end