; docformat = 'rst'

;+
; Make a flipped version of a fitscube, as needed by CRISPEX.
;
; The original file should have the dimensions in order [Nx, Ny,
; Ntuning, Nstokes, Nscans], the flipped version has them in order
; [Ntuning, Nscans, Nstokes, Nx, Ny].
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
;    maxmemory : in, optional, type="scalar number", default="free memory/4"
; 
;      The amount of memory in bytes that the cube/chunk is compared
;      to in order to select the transposing method.
; 
;    method : in, optional, type=string, default='whole'
; 
;      The wanted method for transposing. One of 'whole' (fastest) or
;      'chunks' (also fairly fast). If the selected method is not
;      possible due to the available memory slower method is selected.
;      The fallback is the much slower method of transposing each
;      [Nx,Ny,Ntuning] subcube separately.
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
;   2018-11-05 : MGL. Implement two faster flipping methods, new
;                keywords method and maxmemory.
; 
;-
pro red::fitscube_flip, filename $
                        , flipfile = flipfile $
                        , maxmemory = maxmemory $
                        , method = method $
                        , openclose = openclose $
                        , overwrite = overwrite

  inam = red_subprogram(/low, calling = inam1)

  ;; Use the fastest, whole-cube method by default. 
  if n_elements(method) eq 0 then method = 'whole' 
  
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

  permutation = [2, 4, 3, 0, 1]
  if n_elements(permutation) ne Naxis then stop

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
  ;; i+1 --> permutation[i]+1, where i is the IDL index.
  keywords = strmid(hsp, 0, 8)  ; keywords ordered before flipping

  for iax = 0, Naxis-1 do begin
    ;; Reorder NAXISi
    red_fitsaddkeyword, hsp, 'NAXIS'+strtrim(iax+1, 2), dimensions[permutation[iax]]
    ;; Reorder WCS keywords
    ckeywords = keywords[where(strmatch(keywords.trim() $
                                        ,'C*'+strtrim(permutation[iax]+1, 2)) $
                               , Nc)]
    pkeywords = keywords[where(strmatch(keywords.trim() $
                                        ,'P[SV]*'+strtrim(permutation[iax]+1, 2)+'_*') $
                               , Np)]
    for ikey = 0, Nc-1 do begin
      iline = where(ckeywords[ikey] eq keywords) ; pos in im header
      theline = hsp[iline]
      ckeyword = red_strreplace(ckeywords[ikey] $
                                , strtrim(permutation[iax]+1, 2) $
                                , strtrim(iax+1, 2))
      strput, theline, ckeyword
      hsp[iline] = theline
      print, 'Changed keyword '+ckeywords[ikey]+' to '+ckeyword
    endfor                      ; ikey
    for ikey = 0, Np-1 do begin
      iline = where(pkeywords[ikey] eq keywords)
      theline = hsp[iline]
      pkeyword = red_strreplace(pkeywords[ikey] $
                                , strtrim(permutation[iax]+1, 2)+'_' $
                                , strtrim(iax+1, 2)+'_')
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


  ;; Do some checking of the amount of free memory here
  meminfo = red_meminfo()
  memsize = meminfo.memfree
  cubsize = product(dimensions) * abs(bitpix)/8.d
  if n_elements(maxmemory) eq 0 then maxmemory = memsize/4 
  
  ;; Can we do it in chunks? Compare chunk sizes to maxmemory
  factor, Ny, p, Np, /quiet     ; Need integer subsize of Ny
  Nfactors = n_elements(p)      ; Repeated factors are indicated by elements(Np) gt 1.
  for ifactor = 0, Nfactors-1 do red_append, factors, replicate(p[ifactor], Np[ifactor])
  Nfactors = n_elements(factors) ; Total number of factors, including those repeated.
  Nsubsets = 2^Nfactors          ; Calculate the possible subsizes
  for isubset = 1, Nsubsets-1 do $
     red_append, pp, round(product(factors[where(long((binary(isubset))[-Nfactors:-1]))]))
  chunksizes = pp * abs(bitpix)/8.d * Nx*Nstokes*Nscans*Ntuning ; And the corresponding chunk sizes
  Ny1 = max(pp[where(chunksizes le maxmemory)]) ; Select the largest subsize allowed
  chunksize = Ny1*Nx*Nstokes*Nscans*Ntuning     ; And the corresponding chunk size

  ;; Select method
  if method eq 'whole' and cubsize le maxmemory then begin
    usemethod = 'whole' 
  endif else if method eq 'chunks' and chunksize le maxmemory then begin
    usemethod = 'chunks'
  endif
  
  ;; Should we do it the slow and secure way?
  if usemethod eq 'whole' then begin

    ;; We can fit the whole cube in memory (and then some) so read it
    ;; all in, transpose it, and write it in a single operation. This
    ;; should be the fastest method.

    red_progressbar, 0, 2, 'Transpose the whole cube'

    cube = red_readdata(filename)
    cube = transpose(cube, permutation)

    openu, flun, flipfile, /get_lun, /swap_if_little_endian 
    cube_out = assoc(flun, cube, offset_out)
  
    cube_out[0] = cube

    red_progressbar, 0, 1, 'Transpose the whole cube'

  endif else if usemethod eq 'chunks' then begin
    
    ;; Adapted from Jaime's flipthecube program.
    
    ;; A chunk is a subcube with all the dimensions, except that the Y
    ;; dimension is split into integer subsizes. Each chunk is
    ;; transposed and written in a single operation. This is a fairly
    ;; fast method.

    ;; Possible problem: If Ny is a prime number, no integer subsizes
    ;; are possible as long as we require them to be same size. But we
    ;; could go for as large chunks as possible and then a smaller
    ;; final chunk. ToDo!
    
    Ntimes = Ny / Ny1
    k = 0LL

    openw, llun, /get_lun, 'flip.log'
    
    for itime = 0L, Ntimes - 1 do begin

      case bitpix of
        16  : subcube = intarr(Nx,Ny1,Ntuning,Nstokes,Nscans)
        -32 : subcube = fltarr(Nx,Ny1,Ntuning,Nstokes,Nscans)
        else : stop
      endcase
      
      ;; Open for input
      if keyword_set(openclose) || (itime eq 0) then begin
        openr, ilun, filename, /get_lun, /swap_if_little_endian 
        case bitpix of
          16  : cube_in = assoc(ilun, intarr(Nx, Ny1, /nozero), offset_in) 
          -32 : cube_in = assoc(ilun, fltarr(Nx, Ny1, /nozero), offset_in)
          else : stop
        endcase
      endif

      ;; Read
      for iscan = 0L, Nscans-1 do begin
        for istokes = 0L, Nstokes-1 do begin
          for ituning = 0L, Ntuning-1 do begin
            ele = iscan * Nstokes * Ntuning * Ntimes $
                  + istokes * Ntuning * Ntimes $
                  + ituning * Ntimes $
                  + itime
            printf, llun, ele
            subcube[*,*,ituning,istokes,iscan] = cube_in[ele]             
          endfor                ; ituning
        endfor                  ; istokes
      endfor                    ; iscan

      ;; Close input
      if keyword_set(openclose) then free_lun, ilun

      ;; Transpose
      subcube = transpose(subcube, permutation)

      ;; Open for output
      if keyword_set(openclose) || (itime eq 0) then begin
        openu, flun, flipfile, /get_lun, /swap_if_little_endian 
        case bitpix of
          16  : cube_out = assoc(flun, intarr(Ntuning, Nscans, Nstokes, /nozero), offset_out)
          -32 : cube_out = assoc(flun, fltarr(Ntuning, Nscans, Nstokes, /nozero), offset_out)
          else : stop
        endcase
      endif
      
      ;; Write
      for iy = 0L, Ny1 - 1 do for ix=0L, Nx-1 do begin
        cube_out[k] = subcube[*,*,*,ix,iy]
        k++
      endfor

      ;; Close output
      if keyword_set(openclose) then free_lun, flun
      
      red_progressbar, itime, Ntimes, /predict $
                       , 'Transpose the cube in '+strtrim(Ntimes, 2)+' chunks'

    endfor                      ; itime

    ;; Close the files
    if ~keyword_set(openclose) then free_lun, ilun, flun      

    free_lun, llun
    
  endif else begin

    ;; If we get to this point, we either want to or have to do it the
    ;; slow way.

    print, inam+' : Use the slow but memory-conservative "X/Y/tuning" method.'
    
    case bitpix of
      16  : subcube = intarr(Nx, Ny, Ntuning)
      -32 : subcube = fltarr(Nx, Ny, Ntuning)
      else : stop
    endcase
    
    ;; Loop over subcubes
    iprogress = 0
    for iscan = 0L, Nscans-1 do begin
      for istokes = 0L, Nstokes-1 do begin

        ;; Input
        
        if keyword_set(openclose) || (iprogress eq 0) then begin
          openr, ilun, filename, /get_lun, /swap_if_little_endian 
          case bitpix of
            16  : cube_in = assoc(ilun, intarr(Nx, Ny, /nozero), offset_in) 
            -32 : cube_in = assoc(ilun, fltarr(Nx, Ny, /nozero), offset_in)
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
        endfor                  ; ituning

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
          endfor                ; iy

          iprogress++
          
        endfor                  ; ix

        if keyword_set(openclose) then free_lun, flun
        
      endfor                    ; istokes
    endfor                      ; iscan

    ;; Close the files
    if ~keyword_set(openclose) then free_lun, ilun, flun      

  endelse

  
  ;; Copy WCS extension
  red_fits_copybinext, filename, flipfile, 'WCS-TAB'

  ;; Copy the variable-keywords from the regular nb cube to the
  ;; flipped version.
  var_keys = red_fits_var_keys(him, count = Nkeys)
  for ikey = 0, Nkeys-1 do begin
    self -> fitscube_addvarkeyword, flipfile, var_keys[ikey] $
                                    ,  old_filename = filename, /flipped
  endfor                        ; ikey ;
  
end
