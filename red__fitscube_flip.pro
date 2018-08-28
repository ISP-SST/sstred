; docformat = 'rst'

;+
; Make a flipped version of a fitscube, as needed by CRISPEX.
;
; The original file should have the dimensions in order [Nx, Ny,
; Ntuning, Nstokes, Nscans], the flipped version has them in order
; [Ntuning, Nstokes, Nscans, Nx, Ny].
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
;    filename : in, type=string
; 
;      The file containing the fitscube to be flipped.
; 
; 
; :Keywords:
; 
;    flipfile : out, optional, type=string
;   
;      The name of the flipped file is returned in the keyword.
; 
;    overwrite : in, optional, type=boolean
;
;      Don't care if cube is already on disk, overwrite it with a new
;      version.
;   
;   
; 
; 
; :History:
; 
;   2018-08-28 : MGL. First version based on code from
;                red::fitscube_finish. 
; 
; 
; 
; 
;-
pro red::fitscube_flip, filename, flipfile = flipfile, overwrite = overwrite

  inam = red_subprogram(/low, calling = inam1)
  
  ;; Original file header
  him = headfits(filename)
  Naxis = fxpar(him,'NAXIS')
  dimensions = fxpar(him,'NAXIS*')

  bitpix = fxpar(him,'BITPIX')
  
  Nx      = long(dimensions[0])
  Ny      = long(dimensions[1])
  Ntuning = long(dimensions[2])
  Nstokes = long(dimensions[3])
  Nscans  = long(dimensions[4])

  ;; Some condition here for when the file does not need to be flipped!
  
  ;; Construct the name of the flipped file.
  if strmatch(filename, '*_im.fits') then begin
    flipfile = red_strreplace(filename, '_im.fits', '_sp.fits')
  endif else begin
    flipfile = red_strreplace(filename, '.fits', '_sp.fits')
  endelse

  if file_test(flipfile) then begin
    if keyword_set(overwrite) then begin
      print, inam + ' : Overwriting flipped file ' + flipfile
    endif else begin
      print, inam + ' : File exists ' + flipfile
      print, inam + ' : Use /overwrite if you want to replace it.'
      return
    endelse 
  endif else print, inam + ' : Making flipped file ' + flipfile

  reorder = [2, 4, 3, 0, 1]
  if n_elements(reorder) ne Naxis then stop

  print, 'Number of X-Y frames: '+strtrim(Ntuning*Nstokes*Nscans)
  print, 'Number of lambda-scan frames: '+strtrim(Nx*Ny*Nstokes)

  ;; Header for flipped file
  hsp = him
  red_fitsaddkeyword, hsp, 'FILENAME', file_basename(flipfile)
  ;; No variable-keywords
  red_fitsdelkeyword, hsp, 'VAR_KEYS'


  ;; Change keywords that are affected by the axis reordering.
  ;; Keywords for one axis should be renamed so they are for the
  ;; reordered axis. Because axis numbers start at 1, while IDL
  ;; indices start at 0, the operation will be written as rename from
  ;; i+1 --> reorder[i]+1, where i is the IDL index.
  keywords = strmid(hsp, 0, 8)  ; keywords ordered before flipping

  for iax = 0, Naxis-1 do begin
    ;; Reorder NAXISi
    red_fitsaddkeyword, hsp, 'NAXIS'+strtrim(iax+1, 2), dimensions[reorder[iax]]
    ;; Reorder WCS keywords
;    ckeywords = keywords[where(strmatch(keywords.trim(),'C*'+strtrim(iax+1, 2)), Nc)]
;    pkeywords = keywords[where(strmatch(keywords.trim(),'P[SV]*'+strtrim(iax+1, 2)+'_*'), Np)]
    ckeywords = keywords[where(strmatch(keywords.trim(),'C*'+strtrim(reorder[iax]+1, 2)), Nc)]
    pkeywords = keywords[where(strmatch(keywords.trim(),'P[SV]*'+strtrim(reorder[iax]+1, 2)+'_*'), Np)]
    for ikey = 0, Nc-1 do begin
      iline = where(ckeywords[ikey] eq keywords) ; pos in im header
      theline = hsp[iline]
;     cvalue = red_fitsgetkeyword(him, ckeywords[ikey], comment = ccomment)
;      ckeyword = red_strreplace(ckeywords[ikey], strtrim(iax+1, 2), strtrim(reorder[iax]+1, 2))
      ckeyword = red_strreplace(ckeywords[ikey], strtrim(reorder[iax]+1, 2), strtrim(iax+1, 2))
;      red_fitsaddkeyword, hsp, ckeyword, cvalue, ccomment
      strput, theline, ckeyword
      hsp[iline] = theline
      print, 'Changed keyword '+ckeywords[ikey]+' to '+ckeyword
    endfor                      ; ikey
    for ikey = 0, Np-1 do begin
      iline = where(pkeywords[ikey] eq keywords)
      theline = hsp[iline]
;      pvalue = red_fitsgetkeyword(him, pkeywords[ikey], comment = pcomment)
      pkeyword = red_strreplace(pkeywords[ikey], strtrim(reorder[iax]+1, 2)+'_', strtrim(iax+1, 2)+'_')
;      red_fitsaddkeyword, hsp, pkeyword, pvalue, pcomment
      strput, theline, pkeyword
      hsp[iline] = theline
      print, 'Changed keyword '+pkeywords[ikey]+' to '+pkeyword
    endfor                      ; ikey
  endfor                        ; iax

  ;; Transfer data to flipped file
  
  ;; Open files (Might be OK to skip /swap_if_little_endian on both?)
  openr, ilun, filename, /get_lun, /swap_if_little_endian ; Input file
  openw, flun, flipfile, /get_lun, /swap_if_little_endian ; Output file

  ;; Write the header of the flipped cube
  hsp = hsp[0:where(strmid(hsp,0,3) eq 'END')] ; Strip trailing empty lines
  Nlines = n_elements(hsp)     
  bhdr=reform(byte(hsp),80L*Nlines) ; Byte version of header
;  bhdr = replicate(32B, 80L*Nlines)
;  for n = 0L, Nlines-1 do bhdr[80*n] = byte(hsp[n])
  writeu, flun, bhdr
  
  ;; Make assoc variables
  Nlines = n_elements(him)     
  Nblock = (Nlines-1)*80/2880+1 ; Number of 2880-byte blocks
  offset_in = Nblock*2880       ; Offset to start of data
  print, 'offset_in=', offset_in

  Nlines = n_elements(hsp)     
  Nblock = (Nlines-1)*80/2880+1 ; Number of 2880-byte blocks
  offset_out = Nblock*2880      ; Offset to start of data
  print, 'offset_out=', offset_out
  
  case bitpix of
    16 : begin
      cube_in  = assoc(ilun, intarr(Nx, Ny, /nozero), offset_in)
      cube_out = assoc(flun, intarr(Ntuning, /nozero), offset_out)
      subcube  = intarr(Nx, Ny, Ntuning)
    end
    -32 : begin
      cube_in  = assoc(ilun, fltarr(Nx, Ny, /nozero), offset_in)
      cube_out = assoc(flun, fltarr(Ntuning, /nozero), offset_out)
      subcube  = fltarr(Nx, Ny, Ntuning)
    end
    else : stop
  endcase
  
  ;; Loop over subcubes
  iprogress = 0
  for istokes = 0L, Nstokes-1 do begin
    for iscan = 0L, Nscans-1 do begin

      ;; Read a subcube frame by frame
      red_progressbar, iprogress, Nstokes*Nscans, /predict $
                       , 'Flip the cube, istokes,iscan='+strtrim(istokes, 2)+','+strtrim(iscan, 2)+' - read'
      for ituning = 0, Ntuning-1 do begin
        iframe = ituning + Ntuning*(istokes + Nstokes*iscan)
        subcube[0, 0, ituning] = cube_in[iframe]
      endfor                    ; ituning

      ;; Write the subcube spectrum by spectrum
      red_progressbar, iprogress, Nstokes*Nscans, /predict $
                       , 'Flip the cube, istokes,iscan='+strtrim(istokes, 2)+','+strtrim(iscan, 2)+' - write'
      for ix = 0L, Nx-1 do begin
        for iy = 0L, Ny-1 do begin
          ispectrum = iscan + Nscans*(istokes + Nstokes*(ix + Nx*iy))
          cube_out[ispectrum] = reform(subcube[ix, iy, *])
        endfor                  ; iy
      endfor                    ; ix

      iprogress++
      
    endfor                      ; iscan
  endfor                        ; istokes

  ;; Close the files
  free_lun, ilun, flun              
  
  if n_elements(wcs) gt 0 then begin
    ;; Copy WCS extension to flipped file
    red_fits_copybinext, filename, flipfile, 'WCS-TAB'
  endif
  stop
  ;; Copy some variable-keywords from the ordinary nb cube to the
  ;; flipped version.
  self -> fitscube_addvarkeyword, flipfile, 'SCANNUM',  old_filename = filename, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'ATMOS_R0', old_filename = filename, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'DATE-BEG', old_filename = filename, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'DATE-AVG', old_filename = filename, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'DATE-END', old_filename = filename, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'XPOSURE',  old_filename = filename, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'TEXPOSUR', old_filename = filename, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'NSUMEXP',  old_filename = filename, /flipped

end
