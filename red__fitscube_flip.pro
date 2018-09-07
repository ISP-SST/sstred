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
;    openclose : in, optional, type=boolean
;
;      Open and close the files repeatedly to avoid having two open
;      units at the same time. 
;
;    overwrite : in, optional, type=boolean
;
;      Don't care if cube is already on disk, overwrite it with a new
;      version.
; 
; 
; :History:
; 
;   2018-08-28 : MGL. First version based on code from
;                red::fitscube_finish. 
; 
;   2018-08-30 : MGL. New keyword openclose, prevents input and output
;                files from being open at the same time.
; 
; 
;-
pro red::fitscube_flip, filename $
                        , flipfile = flipfile $
                        , openclose = openclose $
                        , overwrite = overwrite

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
    ckeywords = keywords[where(strmatch(keywords.trim(),'C*'+strtrim(reorder[iax]+1, 2)), Nc)]
    pkeywords = keywords[where(strmatch(keywords.trim(),'P[SV]*'+strtrim(reorder[iax]+1, 2)+'_*'), Np)]
    for ikey = 0, Nc-1 do begin
      iline = where(ckeywords[ikey] eq keywords) ; pos in im header
      theline = hsp[iline]
      ckeyword = red_strreplace(ckeywords[ikey], strtrim(reorder[iax]+1, 2), strtrim(iax+1, 2))
      strput, theline, ckeyword
      hsp[iline] = theline
      print, 'Changed keyword '+ckeywords[ikey]+' to '+ckeyword
    endfor                      ; ikey
    for ikey = 0, Np-1 do begin
      iline = where(pkeywords[ikey] eq keywords)
      theline = hsp[iline]
      pkeyword = red_strreplace(pkeywords[ikey], strtrim(reorder[iax]+1, 2)+'_', strtrim(iax+1, 2)+'_')
      strput, theline, pkeyword
      hsp[iline] = theline
      print, 'Changed keyword '+pkeywords[ikey]+' to '+pkeyword
    endfor                      ; ikey
  endfor                        ; iax

  ;; Transfer data to flipped file
  
  ;; Create flipped file and write header.
  openw, flun, flipfile, /get_lun, /swap_if_little_endian ; Output file

  ;; Write the header of the flipped cube
  hsp = hsp[0:where(strmid(hsp,0,3) eq 'END')] ; Strip trailing empty lines
  Nlines = n_elements(hsp)     
  bhdr = reform(byte(hsp),80L*Nlines) ; Byte version of header
  writeu, flun, bhdr
  free_lun, flun
  
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
    16  : subcube = intarr(Nx, Ny, Ntuning)
    -32 : subcube = fltarr(Nx, Ny, Ntuning)
    else : stop
  endcase
  
  ;; Loop over subcubes
  iprogress = 0
  for istokes = 0L, Nstokes-1 do begin
    for iscan = 0L, Nscans-1 do begin

      ;; Input
      
      if keyword_set(openclose) || (iprogress eq 0) then begin
        openr, ilun, filename, /get_lun, /swap_if_little_endian 
        case bitpix of
          16  : cube_in  = assoc(ilun, intarr(Nx, Ny, /nozero), offset_in) 
          -32 : cube_in  = assoc(ilun, fltarr(Nx, Ny, /nozero), offset_in)
          else : stop
        endcase
      endif
      
;      red_progressbar, iprogress, Nstokes*Nscans, /predict $
;      , 'Flip the cube, istokes,iscan=' + strtrim(istokes, 2) + $
;      ',' + strtrim(iscan, 2) + ' - read'
      for ituning = 0, Ntuning-1 do begin
        ;; Read a subcube frame by frame
        iframe = ituning + Ntuning*(istokes + Nstokes*iscan)
        subcube[0, 0, ituning] = cube_in[iframe]
      endfor                    ; ituning

      if keyword_set(openclose) then free_lun, ilun

      ;; Output
      
      if keyword_set(openclose) || (iprogress eq 0) then begin
        openu, flun, flipfile, /get_lun, /swap_if_little_endian 
        case bitpix of
          16  : cube_out = assoc(flun, intarr(Ntuning, /nozero), offset_out)
          -32 : cube_out = assoc(flun, fltarr(Ntuning, /nozero), offset_out)
          else : stop
        endcase
      endif

;      red_progressbar, iprogress, Nstokes*Nscans, /predict $
;      , 'Flip the cube, istokes,iscan=' + strtrim(istokes, 2) $
;      + ',' + strtrim(iscan, 2) + ' - write'
      for ix = 0L, Nx-1 do begin

        red_progressbar, iprogress, Nstokes*Nscans*Nx, /predict $
                         , 'Write the flipped cube, istokes,iscan,ix=' + strtrim(istokes, 2) $
                         + ',' + strtrim(iscan, 2) + ',' + strtrim(ix, 2) 

        for iy = 0L, Ny-1 do begin
          ;; Write the subcube spectrum by spectrum
          ispectrum = iscan + Nscans*(istokes + Nstokes*(ix + Nx*iy))
          cube_out[ispectrum] = reform(subcube[ix, iy, *])
        endfor                  ; iy

        iprogress++
        
      endfor                    ; ix

      if keyword_set(openclose) then free_lun, flun
      
    endfor                      ; iscan
  endfor                        ; istokes

  ;; Close the files
  if ~keyword_set(openclose) then free_lun, ilun, flun      

  ;; Copy WCS extension
  red_fits_copybinext, filename, flipfile, 'WCS-TAB'

  ;; Copy the variable-keywords from the regular nb cube to the
  ;; flipped version.
  var_keys = red_fits_var_keys(him)
  for ikey = 0, n_elements(var_keys)-1 do begin
    self -> fitscube_addvarkeyword, flipfile, var_keys[ikey] $
                                    ,  old_filename = filename, /flipped
  endfor                        ; ikey ;
  
end
