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
;      The wanted method for transposing. One of 'whole' (fastest),
;      'chunks' (also fairly fast if the number of FOV rows is not a
;      prime), or "slow". If the selected method is not possible due
;      to lack of available memory, a slower method is selected. The
;      fallback is the "slow" method of transposing each
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
;   2018-11-05 : MGL. Two faster flipping methods, new keywords method
;                and maxmemory.
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

  ;; The permutation or the axis order
  permutation = [2, 4, 3, 0, 1] ; ['Tuning','Scan','Stokes','X','Y']
  if n_elements(permutation) ne Naxis then stop

  ;; Calculate the inverse permutation
  diag = diag_matrix(replicate(1, Naxis))
  perm_matrix = diag[permutation,*]
  inv_perm_matrix = round(invert(perm_matrix))
  inv_permutation = inv_perm_matrix # indgen(Naxis)

  orig_axes_types =  ['X','Y', 'Tuning','Stokes','Scan']
  print, 'Original axis order: ', orig_axes_types
  print, 'Flipped order: ', reform(orig_axes_types[permutation])
  
  
  
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
    orig_axis = strtrim(iax+1, 2)
    red_fitsaddkeyword, hsp, 'NAXIS'+orig_axis, dimensions[permutation[iax]]    
  endfor                        ; iax
  
  anchor = 'PC5_5'              ; The C* and P[SV]* keywords should come immediately after the PC keywords.
  for iax = 0, Naxis-1 do begin

    ;; Loop is over the permuted axis numbers
    reord_axis = strtrim(iax+1, 2)
    orig_axis  = strtrim(permutation[iax]+1, 2)

    ;; Reorder WCS keywords
    cpkeywords = keywords[where(strmatch(keywords.trim() $
                                         ,'C*'+orig_axis) $
                                or strmatch(keywords.trim() $
                                            ,'P[SV]*'+orig_axis+'_*') $
                                , Nkeys)]
    
    for ikey = 0, Nkeys-1 do begin
      ;;cpkeyword = red_strreplace(cpkeywords[ikey], reord_axis,
      ;;orig_axis)
      if strmatch(cpkeywords[ikey], '*_*') then begin
        cpkeyword = strtrim(red_strreplace(cpkeywords[ikey], orig_axis+'_', reord_axis+'_'), 2)
      endif else begin
        cpkeyword = strtrim(red_strreplace(cpkeywords[ikey], orig_axis, reord_axis), 2)
      endelse
      ;; Read from original header, write to the new header
      value = red_fitsgetkeyword(him, cpkeywords[ikey], comment = comment)
      red_fitsaddkeyword, anchor = anchor, hsp, cpkeyword, value, comment
      print, 'Changed keyword '+cpkeywords[ikey]+' to '+cpkeyword
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
  Ntimes = Ny / Ny1

;  Nchunks = ceil(cubsize/maxmemory)
;  Ny1values = replicate(Ny/Nchunks, Nchunks)
;  remainder = round(Ny-total(Ny1values))
;  if remainder gt 0 then begin
;    Ny1extra = replicate(0, Nchunks)
;    Ny1extra[0] = replicate(1, remainder)
;    Ny1values += Ny1extra
;  endif
;
;  Ny1values += [1, -2, 1]
;
;  print, inam+' : Could do the '+strtrim(Ny, 2)+' rows in chunks '+strjoin(strtrim(Ny1values,2),' + ')+'.'

  altchunk = 0
  
  ;; Select method
  if method eq 'whole' and cubsize le maxmemory then begin
    usemethod = 'whole' 
  endif else if method eq 'chunks' and chunksize le maxmemory then begin
    usemethod = 'chunks'
  endif else if method eq 'slow' then begin
    usemethod = 'slow'
  endif else begin
    ;; OK, we hadn't specified the method keyword, so just try
    ;; them in order.
    if cubsize le maxmemory then begin
      usemethod = 'whole' 
    endif else if chunksize le maxmemory then begin
      usemethod = 'chunks'
    endif else begin
      usemethod = 'slow'
    endelse
  endelse
  
  ;; Should we do it the slow and secure way?
  if usemethod eq 'whole' then begin

    print, inam + ' : Transpose the whole cube in one go.'

    ;; We can fit the whole cube in memory (and then some) so read it
    ;; all in, transpose it, and write it in a single operation. This
    ;; should be the fastest method.

    cube = red_readdata(filename)
    cube = transpose(cube, permutation)

    openu, flun, flipfile, /get_lun, /swap_if_little_endian 
    cube_out = assoc(flun, cube, offset_out)
  
    cube_out[0] = cube

    red_progressbar, 0, 1, 'Transpose the whole cube'

  endif else if usemethod eq 'chunks' then begin

    if altchunk then begin

      print, inam + ' : Transpose the cube in '+strtrim(Nchunks, 2) + ' chunks, ' $
             + strjoin(strtrim(Ny1values,2),' + ')+'.'

      ;; Adapted from Jaime's flipthecube program.
      
      ;; A chunk is a subcube with all the dimensions, except that the Y
      ;; dimension is split into integer subsizes. Each chunk is
      ;; transposed and written in a single operation. This is a fairly
      ;; fast method.

      
      k = 0LL

      for ichunk = 0L, Nchunks - 1 do begin

        red_progressbar, ichunk, Nchunks, /predict $
                         , 'Transpose the cube in chunks'

        Ny1 = Ny1values[ichunk]
        if ichunk eq 0 then begin
          i0 = 0
        endif else begin
          i0 = Ny1values[ichunk-1]
        endelse
        i1 = i0 + Ny1-1
        
        case bitpix of
          16  : subcube = intarr(Nx,Ny1,Ntuning,Nstokes,Nscans)
          -32 : subcube = fltarr(Nx,Ny1,Ntuning,Nstokes,Nscans)
          else : stop
        endcase
        
        ;; Open for input
        if keyword_set(openclose) || (ichunk eq 0) then begin
          openr, ilun, filename, /get_lun, /swap_if_little_endian 
          case bitpix of
            16  : cube_in = assoc(ilun, intarr(Nx, Ny, /nozero), offset_in) 
            -32 : cube_in = assoc(ilun, fltarr(Nx, Ny, /nozero), offset_in)
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
                    + ichunk
              tmp = cube_in[ele]             
              subcube[*,*,ituning,istokes,iscan] = tmp[*, i0:i1]
            endfor              ; ituning
          endfor                ; istokes
        endfor                  ; iscan

        ;; Close input
        if keyword_set(openclose) then free_lun, ilun

        ;; Transpose
        subcube = transpose(subcube, permutation)

        ;; Open for output
        if keyword_set(openclose) || (ichunk eq 0) then begin
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
        
      endfor                    ; ichunk
      
    endif else begin    
      print, inam + ' : Transpose the cube in '+strtrim(Ntimes, 2)+' chunks.'

      ;; Adapted from Jaime's flipthecube program.
      
      ;; A chunk is a subcube with all the dimensions, except that the Y
      ;; dimension is split into integer subsizes. Each chunk is
      ;; transposed and written in a single operation. This is a fairly
      ;; fast method.

      ;; Possible problem: If Ny is a prime number, no integer subsizes
      ;; are possible as long as we require them to be same size. But we
      ;; could go for as large chunks as possible and then a smaller
      ;; final chunk. ToDo!
      
      k = 0LL

;      openw, llun, /get_lun, 'flip.log'
      
      for itime = 0L, Ntimes - 1 do begin

        red_progressbar, itime, Ntimes, /predict $
                         , 'Transpose the cube in chunks'

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
;              printf, llun, ele
              subcube[*,*,ituning,istokes,iscan] = cube_in[ele]             
            endfor              ; ituning
          endfor                ; istokes
        endfor                  ; iscan

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
        
      endfor                    ; itime 
    endelse
    
    ;; Close the files
    if ~keyword_set(openclose) then free_lun, ilun, flun      

;    free_lun, llun
    
  endif else begin

    print, inam+' : Use the slow but memory-conservative "X/Y/tuning" method.'

    ;; If we get to this point, we either want to or have to do it the
    ;; slow way.
    
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


  ;; Copy the variable-keywords from the regular nb cube to the
  ;; flipped version.
  var_keys = red_fits_var_keys(him, count = Nkeys)
  for ikey = 0, Nkeys-1 do begin
    self -> fitscube_addvarkeyword, flipfile, var_keys[ikey] $
                                    ,  old_filename = filename, /flipped
  endfor                        ; ikey ;
  
  ;; Copy WCS extension
  red_fits_copybinext, filename, flipfile, 'WCS-TAB'

  ;; Copy cavity maps
  cmaps = mrdfits(filename, 'WCSDVARR', chdr, status = status, /silent)
  if status eq 0 then begin
    writefits, flipfile, cmaps, chdr, /append
    ;; The CWERRj, CWDISj, and DWj keywords should already be in the
    ;; header. We just need to copy the WCSDVARR (image) extension.
  endif else begin
    print, inam+' : This file does not seem to have cavity maps.'
  endelse

end
