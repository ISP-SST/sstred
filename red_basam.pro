; docformat = 'rst'

;+
; Calculate BaSAM (Denker & Verma 2019, Denker et al. 2023).
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
;    The BaSAMs as a (Nx, Ny, Nwave, Nstokes) array.
; 
; :Params:
; 
;   filename : in, type=string 
; 
;     The path to a fitscube file.
; 
; 
; :Keywords:
; 
;   dt : in, optional, type=float
;   
;     The length of a sliding temporal window in seconds. Default is
;     to use the entire datacube.
;
;   s : out, optional, type=fltarr
;
;     The S array.
;
;   spectral : in, optional, type=boolean
; 
;     Calculate spectral BaSAMs.
;
;   use_median : in, optional, type=boolean
; 
;     Use median rather than mean.
; 
; :History:
; 
;   2024-01-13 : MGL. First version.
; 
;-
function red_basam, filename, dt = dt, s = s, use_median = use_median, spectral = spectral

  h = headfits(filename)

  stokesname = ['I', 'Q', 'U', 'V']
  
  naxis = fxpar(h, 'NAXIS*')
  Nx      = naxis[0]
  Ny      = naxis[1]
  Nwave   = naxis[2]
  Nstokes = naxis[3]
  Nrepeat = naxis[4]

  if keyword_set(spectral) then begin
    
    BaSAM = fltarr(Nx, Ny, Nstokes,  Nrepeat)
    S = fltarr(Nx, Ny, Nstokes, Nrepeat)
    
    for istokes = 0, Nstokes-1 do begin
      for irepeat = 0, Nrepeat-1 do begin
 
        red_progressbar, irepeat, Nrepeat, 'Processing Stokes '+stokesname[istokes]
        S = fltarr(Nx, Ny, Nwave)

        for iwave = 0, Nwave-1 do begin
          red_fitscube_getframe, filename, frame $
                                 , iscan = irepeat $
                                 , istokes = istokes $
                                 , ituning = iwave
          S[*, *, iwave] = frame
        endfor                  ; iwave

        
        if keyword_set(use_median) then begin
          Savg = median(S, dim = 3)
          BaSAM[*, *, istokes, irepeat] = median(abs(S-rebin(Savg, Nx, Ny, Nrepeat, /samp)), dim = 3)
        endif else begin
          Savg = mean(S, dim = 3)
          BaSAM[*, *, istokes, irepeat] = mean(abs(S-rebin(Savg, Nx, Ny, Nrepeat, /samp)), dim = 3)
        endelse
        
      endfor                    ; irepeat
    endfor                      ; istokes

    basam = reform(basam, Nx, Ny, 1, Nstokes,  Nrepeat)
    
  endif else begin
    
    BaSAM = fltarr(Nx, Ny, Nwave, Nstokes)

    for istokes = 0, Nstokes-1 do begin
      for iwave = 0, Nwave-1 do begin

        red_progressbar, iwave, Nwave, 'Processing Stokes '+stokesname[istokes]
        S = fltarr(Nx, Ny, Nrepeat)

        for irepeat = 0, Nrepeat-1 do begin
          red_fitscube_getframe, filename, frame $
                                 , iscan = irepeat $
                                 , istokes = istokes $
                                 , ituning = iwave
          S[*, *, irepeat] = frame
        endfor                  ; irepeat

        if keyword_set(use_median) then begin
          Savg = median(S, dim = 3)
          BaSAM[*, *, iwave, istokes] = median(abs(S-rebin(Savg, Nx, Ny, Nrepeat, /samp)), dim = 3)
        endif else begin
          Savg = mean(S, dim = 3)
          BaSAM[*, *, iwave, istokes] = mean(abs(S-rebin(Savg, Nx, Ny, Nrepeat, /samp)), dim = 3)
        endelse
        
      endfor                    ; iwave
    endfor                      ; istokes

    basam = reform(basam, Nx, Ny, Nwave, Nstokes, 1)

  endelse
  
  return, BaSAM
  
end

dir = '/scratch/mats/2016.09.19/CRISP-aftersummer/cubes_nb/'
filename = 'nb_6302_2016-09-19T09:30:20_scans=1-43_stokes_corrected_im.fits'
basam = red_basam(dir+filename)
sbasam = red_basam(dir+filename, /spectral)

end
