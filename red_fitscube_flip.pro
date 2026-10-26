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
;    method : in, optional, type=string, default='chunks'
; 
;      The wanted method for transposing. One of 'chunks' (fastest) or
;      "slow". If the 'chunks' is not possible due to lack of
;      available memory, the "slow" method of transposing each
;      [Nx,Ny,Ntuning] subcube separately is used.
; 
;    no_metadata : in, optional, type=boolean
; 
;      Don't copy lots of unnecessary metadata to the flipped cube.
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
;   2019-04-04 : MGL. Removed the "whole" method because it caused
;                problems with integer cubes. And "chunks" with a
;                single chunk should be equivalent.
;
;   2019-10-17 : MGL. Make it a regular subroutine rather than a
;                method. 
;
;   2026-09-21 : MGL. New keyword no_metadata. 
;
;   2026-09-22 : MGL. Rewrote the chunk processing.
; 
;-
pro red_fitscube_flip, filename $
                       , flipfile = flipfile $
                       , maxmemory = maxmemory $
                       , method = method $
                       , no_metadata = no_metadata $
                       , openclose = openclose $
                       , overwrite = overwrite

  inam = red_subprogram(/low, calling = inam1)

  ;; Use the fastest, method by default. 
  if n_elements(method) eq 0 then method = 'chunks' 
  
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

  ;; =========================================================================
  ;; Transfer data to flipped file 
  ;; =========================================================================
  
  ;; 1. Initialize output file and pre-allocate full size on disk to
  ;; prevent fragmentation
  openw, flun, flipfile, /get_lun, /swap_if_little_endian 
  hsp = hsp[0:where(strmid(hsp,0,3) eq 'END')] ; Strip trailing empty lines
  Nlines = n_elements(hsp)     
  bhdr = reform(byte(hsp), 80L*Nlines) ; Byte version of header
  writeu, flun, bhdr
  
  ; Calculate exact data size and offsets using 64-bit integers for large files
  typesize = abs(bitpix) / 8L
  total_data_bytes = Long64(Nx) * Ny * Ntuning * Nstokes * Nscans * typesize
  Nlines_hsp = n_elements(hsp)
  Nblock_out = (Nlines_hsp - 1) * 80 / 2880 + 1
  offset_out = Nblock_out * 2880
  
  ; Force the OS to allocate the full file space continuously
  point_lun, flun, offset_out + total_data_bytes - 1
  writeu, flun, 0b 
  free_lun, flun

  ;; Calculate input file data offset
  Nlines_him = n_elements(him)     
  Nblock_in = (Nlines_him - 1) * 80 / 2880 + 1
  offset_in = Nblock_in * 2880

  ;; 2. Determine optimal chunk size based on available memory (Safe
  ;; for prime numbers)
  meminfo = red_meminfo()
  memsize = meminfo.memfree
  if n_elements(maxmemory) eq 0 then maxmemory = memsize / 4 
  
  ;; Calculate memory footprint of a single Y-row across all other dimensions
  bytes_per_Y_row = Double(Nx) * Ntuning * Nstokes * Nscans * typesize
  
  ;; Determine how many Y-rows (Ny1) can fit into maxmemory
  Ny1 = Long(floor(maxmemory / bytes_per_Y_row))
  if Ny1 lt 1 then Ny1 = 1
  if Ny1 gt Ny then Ny1 = Ny
  
  ;; Fallback to 1 row if 'slow' method is explicitly requested
  if method eq 'slow' then Ny1 = 1
  
  ;; Calculate the total number of chunks needed
  Nchunks = Long(ceil(Float(Ny) / Ny1))
  
  print, inam + ' : Transposing the cube in ' + strtrim(Nchunks, 2) + ' chunks.'

  ;; 3. Main I/O Loop over Y-chunks
  for ichunk = 0L, Nchunks - 1 do begin
    
;    red_progressbar, ichunk, Nchunks, /predict, 'Transpose the cube in chunks'
    
    ; Define the Y-index boundaries for the current chunk (handles the final remainder chunk)
    i0 = ichunk * Ny1
    i1 = (i0 + Ny1 - 1) < (Ny - 1)
    current_Ny1 = i1 - i0 + 1
    
    ; Allocate subcube for this specific chunk
    case bitpix of
      16  : subcube = intarr(Nx, current_Ny1, Ntuning, Nstokes, Nscans, /nozero)
      -32 : subcube = fltarr(Nx, current_Ny1, Ntuning, Nstokes, Nscans, /nozero)
      else : stop
    endcase
    
    ; Open input file and map it using ASSOC for fast 2D frame access
    openr, ilun, filename, /get_lun, /swap_if_little_endian 
    case bitpix of
      16  : cube_in = assoc(ilun, intarr(Nx, Ny, /nozero), offset_in) 
      -32 : cube_in = assoc(ilun, fltarr(Nx, Ny, /nozero), offset_in)
      else : stop
    endcase
    
    ;; Read data slice by slice. Original structure is [Nx, Ny,
    ;; Ntuning, Nstokes, Nscans]. Extract the [Nx, i0:i1] region for
    ;; each combination of the outer dimensions
    for iscan = 0L, Nscans - 1 do begin

      red_progressbar, /add_loopindex, /predict $
                       , ichunk, Nchunks, 'Chunk' $
                       , iscan, Nscans, 'Scan' 
                       
      for istokes = 0L, Nstokes - 1 do begin
        for ituning = 0L, Ntuning - 1 do begin
          ;; Retrieve the full 2D frame from disk
          tmp_frame = cube_in[ituning + Ntuning * (istokes + Nstokes * iscan)]
          ;; Extract only the Y-rows belonging to the current chunk
          subcube[*, *, ituning, istokes, iscan] = tmp_frame[*, i0:i1]
        endfor                  ; ituning
      endfor                    ; istokes
    endfor                      ; iscan
    free_lun, ilun              ; Close input file immediately to free system resources
    
    ;; Transpose the entire chunk in memory to order: [Ntuning,
    ;; Nscans, Nstokes, Nx, Ny1]
    subcube = transpose(subcube, permutation)
    
    ;; Open output file and write the entire chunk sequentially 
    openu, flun, flipfile, /get_lun, /swap_if_little_endian
    
    ;; Compute the exact byte offset for this chunk in the output
    ;; file. In the flipped file, each step in the Y-dimension (i0)
    ;; shifts the position by Ntuning * Nscans * Nstokes * Nx
    ;; elements.
    chunk_byte_offset = offset_out + (Long64(i0) * Nx * Ntuning * Nstokes * Nscans * typesize)
    
    point_lun, flun, chunk_byte_offset
    writeu, flun, subcube
    free_lun, flun

  endfor ; ichunk
;  red_progressbar, Nchunks, Nchunks, 'Chunks';, Nscans, Nscans, 'Scans'
  red_progressbar, /done, message_done = 'Chunks'


  print, inam + ' : Data transposition completed successfully.'
  ;; =========================================================================

  if ~keyword_set(no_metadata) then begin
  
    ;; Copy the variable-keywords from the regular nb cube to the
    ;; flipped version.
    var_keys = red_fits_var_keys(him, count = Nkeys)
    for ikey = 0, Nkeys-1 do begin
      red_progressbar, ikey, Nkeys, 'Add variable keywords'
      red_fitscube_addvarkeyword, flipfile, var_keys[ikey] $
                                  ,  old_filename = filename, /flipped
    endfor                      ; ikey
;    red_progressbar, Nkeys, Nkeys, 'Add variable keywords'
    red_progressbar, /done
 
    ;; Copy WCS extension
    print, inam+' : Copy the WCS extension...'
    tic
    red_fits_copybinext, filename, flipfile, 'WCS-TAB'
    toc
    print, inam+' : Copy the WCS extension... Done!'

    ;; Copy cavity maps
    fits_open, filename, fcb
    free_lun, fcb.unit
    if total(fcb.extname eq 'WCSDVARR') eq 1 then begin
      cmaps = mrdfits(filename, 'WCSDVARR', chdr, status = status, /silent)
      if status ne 0 then stop
      writefits, flipfile, cmaps, chdr, /append
      ;; The CWERRj, CWDISj, and DWj keywords should already be in the
      ;; header. We just need to copy the WCSDVARR (image) extension.
    endif else begin
      print, inam + ' : No cavity maps to copy.'
    endelse

  endif else begin

    red_message, 'We should really delete some more keywords from the header before initializing the file.'
    
  endelse
  
end

;a = chromisred(/no, /dev)

filename = 'cubes_nb/nb_3950_2026-08-25T09:23:49_09:23:49=1-4_corrected_im.fits'
filename = 'cubes_nb_mats/mats_nb_cube.fits'
filename = 'cubes_nb_mats/nb_3950_2026-08-25T08:45:10_08:45:10=0-46_08:59:22=0-12_09:06:15=0-15_09:16:38=0-16_09:23:49=0-160_09:57:59=0-105_corrected_im.fits'
red_fitscube_flip, filename, /no_meta, /over

end
