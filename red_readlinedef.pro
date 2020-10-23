function red_readlinedef, filename

  spawn, 'cat '+filename, lines

  ;; Remove commented lines
  indx = where(strmid(lines, 0, 1) ne '#', Nlines)

  if Nlines eq 0 then return, []
  
  lines = lines[indx]

  ;; Execute python lines as IDL
  for iline = 0, Nlines-1 do begin

    ;; Prevent future conflicts with additions to the linedef.py file
    if strmatch( lines[iline], "end_of_lines", /fold_case ) then break  
    
    ;; python --> IDL
    lines[iline] = red_strreplace(lines[iline], ':', '*')
    lines[iline] = red_strreplace(lines[iline], ';', '&')
    status = execute(lines[iline])

    varname = (strsplit(lines[iline], '[ [=]', /extract))[0]
    if varname ne '' then red_append, varnames, varname
  endfor

  ;; What variables did the script set?
  varnames = varnames[uniq(varnames, sort(varnames))]
  Nvars = n_elements(varnames)
  
  ;; Make a struct with the variables
  for ivar = 0, Nvars-1 do begin
    if ivar eq 0 then begin
      status = execute('struc = create_struct("'+varnames[ivar]+'",'+varnames[ivar]+')')
    endif else begin
      status = execute('struc = create_struct(struc, "'+varnames[ivar]+'",'+varnames[ivar]+')')
    endelse
  endfor                        ; ivar

  ;; Old linedef files may not set all variables. So add some defaults.
  if max(varnames eq 'nb_wl') eq 0 then $
     struc = create_struct(struc, 'nb_wl', [ 0, 3925, 3934, 3969, 3978, 3999, 4862, 0, 0, 0 ])
  if max(varnames eq 'wb_wl') eq 0 then $
     struc = create_struct(struc, 'wb_wl', [ 0, 0, 0, 0, 0, 3950, 4846, 0, 0, 0 ])
  if max(varnames eq 'wb_map') eq 0 then $
     struc = create_struct(struc, 'wb_map', [ 0, 5, 5, 5, 5, 5, 6, 7, 8, 9 ])
  
  return, struc
  
end

filename = '/scratch/mats/2019-11-11/CHROMIS/downloads/linedef.py-2019.11.11-09:00:17'       ; With WB H&K filter
filename = '/scratch/mats/2016.09.19/CHROMIS-jan19/downloads/linedef.py-2016.09.19-07:39:19' ; Without wb_map etc.
filename = '/scratch/mats/2020-08-18/CHROMIS/downloads/linedef.py-2020.08.18-08:04:57'       ; With WB limb filter
linedef = red_readlinedef(filename)

print, linedef

if max(tag_names(linedef) eq 'WB_MAP') gt 0 then begin
  nb = 3934
  nb = 3969
  iwb = linedef.wb_map[where(linedef.nb_wl eq nb, Nmatch)]
  if Nmatch gt 0 then print,  'WB filter matching NB'+strtrim(nb, 2)+': ' + strtrim(linedef.wb_wl[iwb], 2)
endif


end
