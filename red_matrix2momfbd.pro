; docformat = 'rst'

;+
; Smooth an inverse modulation matrix map with PSFs from the image
; restoration process.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    lc : 
;   
;       The momfbd data for the different LC states.
;   
;    imodmat : 
;   
;       The inverse modulation matrix.
;   
; 
; :Keywords:
; 
;    amap : in, optional, type=fltarr(3,3)
; 
;       A projective transform matrix.
;
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2014-01-10 : PS. Cleanup (remove unused code).
;
;   2018-10-18 : MGL. Cleanup and change input to array of LC state
;                structs. 
;
;   2018-10-19 : MGL. Adapt to using the projective transforms rather
;                than the old align_clips and xoffs,yoffs mechanism.
;
;-
function red_matrix2momfbd, lc, imodmat $
                            , amap = amap

  ;; Name of this function
  inam = red_subprogram(/low, calling = inam1)      

  Nstokes = 4
  Nlc = n_elements(lc)

  imm = replicate({patch:lc[0].patch, clip:lc[0].clip},[4,4])

  dim=size(imodmat, /dimension)

  mny = dim[2] - 1
  mnx = dim[1] - 1

  dim = size(lc[0].patch, /dimension)

  ;; Figure out if we need to average over PSFs
  psfdim = size(lc[0].patch[0, 0].psf,/dim)
  if n_elements(psfdim) eq 2 then Npsfs = 1 else Npsfs = psfdim(2)
  
  for ix = 0, dim[0]-1 do begin   ; Loop over momfbd subfields
    for iy = 0, dim[1]-1 do begin ; with ix and iy
      
      pnx = lc[0].patch[ix,iy].xh - lc[0].patch[ix,iy].xl + 1
      pny = lc[0].patch[ix,iy].yh - lc[0].patch[ix,iy].yl + 1
      cpmm = dindgen(4,4,pnx,pny)      
      pmm = dindgen(pnx,pny)

      for ilc = 0, Nlc-1 do begin
      
        clip = red_getclips(lc[ilc],ix,iy)

        if n_elements(amap) gt 0 then begin
          ;; Apply the projective transform to clip, the patch corner
          ;; coordinates. 
          c = [mean(clip[0:1]), mean(clip[2:3]), 1]     ; Projective coordinate vector for patch center
          c1 = amap ## c                                ; Projected center
          c1 = c1/c1[2]                                 ; Normalize
          dc = round(c1-c)                              ; Center displacement
          clip = [clip[0:1] + dc[0], clip[2:3] + dc[1]] ; Projected corner coordinates
        endif

        clp = [(clip[0]>0)<mnx, (clip[1]<mnx), (clip[2]>0)<mny, (clip[3]<mny)]
        if clp[0] ne clp[1] and clp[2] ne clp[3] then begin
          tpmm=imodmat[*,clp[0]:clp[1], clp[2]:clp[3]]
;        tpmm=imodmat[*,(clip[0]>0)<mnx:(clip[1]<mnx),(clip[2]>0)<mny:(clip[3]<mny)]

          for istokes=0,Nstokes-1 do begin
            pmm[*,*]=median(tpmm[ilc*4 + istokes,*,*]) ; Fill extra bits with the median
            pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[istokes + ilc*4,*,*]) ; Copy what we do have
            if n_elements(amap) gt 0 then begin
              ;; The projective transform might require a mirroring of
              ;; the patch, in addition to a transpose that the momfbd
              ;; program does (for no good reason).
              case 1 of
                amap[0, 0] lt 0 && amap[1, 1] lt 0 : pmm = rotate(pmm, 6) ; Transpose after mirroring in X and Y
                amap[0, 0] lt 0                    : pmm = rotate(pmm, 3) ; Transpose after mirroring in X
                amap[1, 1] lt 0                    : pmm = rotate(pmm, 1) ; Transpose after mirroring in Y
                else                               : pmm = rotate(pmm, 4) ; Transpose after no mirroring
              endcase 
            endif
            ;; This psf should be the single psf you get when running
            ;; with GET_PSF_AVG (or whatever the keyword is called), or
            ;; the average PSF over the ones corresponding to the
            ;; multiple raw images.
            ;;
            ;; Or, in the case we want SIMPLE, the spatial average of
            ;; PSFs (potentially as well as over raw images). Unless we
            ;; want to do that operation without making.
            case Npsfs of
              1    : psf =       lc[ilc].patch[ix,iy].psf             ; The average PSF delivered by momfbd
              else : psf = total(lc[ilc].patch[ix,iy].psf, 3) / Npsfs ; The average PSF calculated here
            endcase
            imm[istokes,ilc].patch[ix,iy].img = red_convolve(red_fillnan(pmm), psf)
          endfor                ; istokes

        endif else stop
        
      endfor                    ; ilc

    endfor                      ; iy
  endfor                        ; ix

  return, imm
  
end
