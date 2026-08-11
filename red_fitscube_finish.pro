; docformat = 'rst'

;+
; Close the fitscube file, optionally make a "flipped" version.
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
;    lun : in, type=integer
; 
;      The logical unit number of the fitscube file to finish.
; 
; 
; :Keywords:
; 
;    flipfile : out, optional, type=string
;   
;      The precense of this keyword causes a "flipped" version of the
;      [Nx, Ny, Ntuning, Nstokes, Nscans] fitscube to be produced, in
;      which the dimensions are [Ntuning, Nstokes, Nscans, Nx, Ny].
;      The name of the flipped file is returned in the keyword.
; 
;    wcs : in, optional, type=struct
;
;      WCS info to be written as a binary extension.
; 
; :History:
; 
;   2017-09-11 : MGL. First version.
; 
;   2017-11-16 : MGL. Now works with integer cubes.
;
;   2019-09-13 : MGL. A version that is not a class method.
; 
;-
pro red_fitscube_finish, lun, flipfile = flipfile, wcs = wcs

  inam = red_subprogram(/low, calling = inam1)
  
  filename = (fstat(lun)).name  ; Get the file name while the file is still open.
  print, inam + ' : Closing fitscube file ' + filename

  ;; Close the file
  free_lun, lun              

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

  if n_elements(wcs) gt 0 then begin
    red_fitscube_addwcs, filename, wcs, dimensions = dimensions
  endif
  
  if ~arg_present(flipfile) then begin
    ;; If we don't want to flip, we are done now.
    return
  endif

  ;; Update header after WCS was added
  him = headfits(filename)

  ;; Construct the name of the flipped file.
  if strmatch(filename, '*_im.fits') then begin
    flipfile = red_strreplace(filename, '_im.fits', '_sp.fits')
  endif else begin
    flipfile = red_strreplace(filename, '.fits', '_sp.fits')
  endelse

  print, inam + ' : Making flipped file ' + flipfile

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
  offset_in = Nblock*2880          ; Offset to start of data
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

end
