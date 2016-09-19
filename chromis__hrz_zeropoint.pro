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
;-

pro chromis::hrz_zeropoint

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  infodir = self.out_dir + 'info/'
  file_mkdir, infodir

  ;; Get du_ref from the unique list of hrz tunings in raw flats taken
  ;; with the nb camera, list them in all flats directories.

  ;; Get all raw flat file names
  fnames = file_search(*self.flat_dir + '/Chromis-N/*' + '*.fits', count = Nff)

  self -> extractstates, fnames, states

  upref = states[uniq(states.prefilter, sort(states.prefilter))].prefilter
  Npref = n_elements(upref)

  for ipref = 0, Npref-1 do begin
  
     zfile = infodir + 'hzr_zeropoint_' + upref[ipref] + '.fz'

     ;; Just a single continuum point, set the reference
     ;; wavelength to that wavelength and the reference du to
     ;; the actual du, so the tuning is zero.
     case upref[ipref] of
        'Hb-core'  : lambda_ref = 4861d-10                         ; Hbeta line center
        'CaK-core' : lambda_ref = 3933.7d-10                       ; Ca II K line center
        'CaH-core' : lambda_ref = 3968.5d-10                       ; Ca II H line center
        'CaH-cont' : lambda_ref = (3998.640d-10 + 1.258d-10)       ; Clean reference line within the passband
        else: begin
           print, inam+' : Please add reference wavelength for prefilter '+states[ifile].prefilter
           stop
        end
     endcase 

     convfac = -1.53d                 ; [mÅ/du] for CaH-cont 3998.6 Å
     convfac *= lambda_ref/3998.6e-10 ; Scales with lambda.
     
     indx = where(states.prefilter eq upref[ipref])
     prefstates = states[indx]

     hrz = long((stregex(prefstates.fullstate, 'hrz([0-9]*)', /extract, /subexpr))[1])
     uhrz = hrz[uniq(hrz, sort(hrz))]
     Nhrz = n_elements(uhrz)
     
     ;; The scans are supposed to be symmetric around the reference
     ;; wavelength, so the number of tunings should be odd.
     if Nhrz mod 2 ne 0 then begin

        hzr_zero = uhrz[Nhrz/2]
        
        fzwrite, [lambda_ref, hzr_zero, convfac], zfile, ' '

     endif else begin

        print, inam + ' : Warning, even number of hrz states for prefilter ' + upref[ipref] + '.'
        print, inam + ' : No hrz zero point found.'
        
     endelse

     
  endfor                        ; ipref

  
end
