; docformat = 'rst'

;+
; Neckel's limb darkening coefficients for 300 nm <= lambda <
; 1100 nm. 
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
;    lambda : in, type=float
; 
;       The wavelength in units of meters.
; 
; :Keywords:
; 
;     analytical : in, optional, type=boolean
; 
;        Get coefficients from Neckel (2005) instead of default Neckel
;        & Labs (1994). 
;   
; :History:
; 
;    2021-09-06 : MGL. First version, based on private IDL library
;                 function neckel_coefficients().
; 
;-
function red_neckel_coefficients, lambda, analytical = analytical

  A = fltarr(6)
  
  if keyword_set(analytical) then begin
    ;; Neckel (2005)
    ;; Lambda input in [m], Neckel uses lambda in microns.
    lambda1 = (double(lambda)*1d6)^(-1)
    lambda5 = (double(lambda)*1d6)^(-5)
    
    if lambda ge 300.00e-9 and lambda lt 372.98e-9 then begin
      print, lambda*1E9, lambda1, ' 1st range'
      a00 =  0.35601d
      a01 = -0.085217d
      a20 = -0.67237d
      a30 =  0.18696d
      a40 =  0.00038d
      a50 =  0.01373d
      a25 =  0.003589d
      a35 = -0.002415d
      a45 =  0.000897d
      a55 = -0.000200d
    endif else if lambda ge 385e-9 and lambda lt 422.57e-9 then begin
      print, lambda*1E9, lambda1, ' 2nd range'
      a00 =  0.09900d
      a01 =  0.010643d
      a20 = -2.80511d
      a30 =  3.32780d
      a40 = -2.17878d
      a50 =  0.58825d
      a25 =  0.024873d
      a35 = -0.029317d
      a45 =  0.018273d
      a55 = -0.004663d
    endif else if lambda ge 422.57e-9 and lambda lt 1100e-9 then begin
      print, lambda*1E9, lambda1, ' 3rd range'
      a00 =  0.75267d
      a01 = -0.265577d
      a20 = -1.89287d
      a30 =  2.42234d
      a40 = -1.71150d
      a50 =  0.49062d
      a25 =  0.012582d
      a35 = -0.017117d
      a45 =  0.011977d
      a55 =  0.003347d
    endif else begin
      print, 'neckel_coefficients: Out of range lambda = ', strtrim(lambda*1e9, 2), ' nm'
      return, 0
    endelse

    ;; if lambda1 gt 3.33 or lambda1 lt 0.91 then begin
    ;;    print, 'Out of range'
    ;;    return, 0
    ;; endif
    ;; if lambda1 gt 2.68 then begin
    ;;    print, lambda*1E9, lambda1, ' 1st range'
    ;;    a00 =  0.35601d
    ;;    a01 = -0.085217d
    ;;    a20 = -0.67237d
    ;;    a30 =  0.18696d
    ;;    a40 =  0.00038d
    ;;    a50 =  0.01373d
    ;;    a25 =  0.003589d
    ;;    a35 = -0.002415d
    ;;    a45 =  0.000897d
    ;;    a55 = -0.000200d
    ;; endif else if lambda1 lt 2.37 then begin
    ;;    print, lambda*1E9, lambda1, ' 3rd range'
    ;;    a00 =  0.75267d
    ;;    a01 = -0.265577d
    ;;    a20 = -1.89287d
    ;;    a30 =  2.42234d
    ;;    a40 = -1.71150d
    ;;    a50 =  0.49062d
    ;;    a25 =  0.012582d
    ;;    a35 = -0.017117d
    ;;    a45 =  0.011977d
    ;;    a55 =  0.003347d
    ;; endif else begin
    ;;    print, lambda*1E9, lambda1, ' 2nd range'
    ;;    a00 =  0.09900d
    ;;    a01 =  0.010643d
    ;;    a20 = -2.80511d
    ;;    a30 =  3.32780d
    ;;    a40 = -2.17878d
    ;;    a50 =  0.58825d
    ;;    a25 =  0.024873d
    ;;    a35 = -0.029317d
    ;;    a45 =  0.018273d
    ;;    a55 = -0.004663d
    ;; endelse

    a10 = 1.d - (a00 + a20 + a30 + a40 + a50)
    a11 = -a01
    a15 = -(a25 + a35 + a45 + a55)
    
    A[0] = a00 + a01*lambda1
    A[1] = a10 + a11*lambda1 + a15*lambda5
    A[2] = a20 + a25*lambda5
    A[3] = a30 + a35*lambda5
    A[4] = a40 + a45*lambda5
    A[5] = a50 + a55*lambda5

  endif else begin
    
    ;; Neckel & Labs 1994, Solar Physics, 153(1-2):91-114. (skip
    ;; columns 2 and 3), lambda given in [nm].
    table1 = fltarr(30, 8)
    table1[ 0,*] = [303.327,0.08011,0.70695,0.49910,-0.31080,-0.02177,0.04642,0.6826]
    table1[ 1,*] = [310.843,0.08160,0.71609,0.69685,-0.87703,0.47008,-0.08760,0.6883]
    table1[ 2,*] = [320.468,0.08833,0.77285,0.65382,-1.04647,0.72921,-0.19775,0.6985]
    table1[ 3,*] = [329.897,0.09188,0.92459,0.19604,-0.39546,0.23599,-0.05303,0.7116]
    table1[ 4,*] = [349.947,0.11012,1.02168,-0.10924,-0.00055,-0.08688,0.06487,0.7260]
    table1[ 5,*] = [365.875,0.12828,1.04969,0.17482,-1.13371,1.23882,-0.45990,0.7445]
    table1[ 6,*] = [374.086,0.12579,0.85402,0.54601,-1.15048,0.88928,-0.26462,0.7288]
    table1[ 7,*] = [390.915,0.12995,0.91836,-0.07566,0.19149,-0.28712,0.12298,0.7204]
    table1[ 8,*] = [401.970,0.12323,1.08648,-0.43974,0.45912,-0.32759,0.09850,0.7303]
    table1[ 9,*] = [416.319,0.12814,1.19947,-0.84407,1.07201,-0.79537,0.23982,0.7380]
    table1[10,*] = [427.930,0.14249,1.28796,-1.19564,1.68603,-1.36658,0.44572,0.7495]
    table1[11,*] = [443.885,0.16220,1.24893,-0.92165,0.89978,-0.50148,0.11220,0.7588]
    table1[12,*] = [445.125,0.15248,1.38517,-1.49615,1.99886,-1.48155,0.44119,0.7596]
    table1[13,*] = [457.345,0.16604,1.38544,-1.52275,2.00232,-1.45969,0.42864,0.7651]
    table1[14,*] = [477.427,0.19571,1.30551,-1.25845,1.50626,-1.05472,0.30570,0.7751]
    table1[15,*] = [492.905,0.20924,1.30798,-1.20411,1.21505,-0.67196,0.14381,0.7823]
    table1[16,*] = [519.930,0.23695,1.29927,-1.28034,1.37760,-0.85054,0.21706,0.7925]
    table1[17,*] = [541.760,0.26073,1.27428,-1.30352,1.47085,-0.96618,0.26384,0.8002]
    table1[18,*] = [559.950,0.26892,1.34319,-1.58427,1.91271,-1.31350,0.37295,0.8061]
    table1[19,*] = [579.880,0.28392,1.36896,-1.75998,2.22154,-1.56074,0.44630,0.8125]
    table1[20,*] = [610.975,0.30854,1.36620,-1.83572,2.33221,-1.63082,0.45959,0.8221]
    table1[21,*] = [640.970,0.33644,1.30590,-1.79238,2.45040,-1.89979,0.59943,0.8290]
    table1[22,*] = [669.400,0.34685,1.37539,-2.04425,2.70493,-1.94290,0.55999,0.8360]
    table1[23,*] = [700.875,0.37885,1.25553,-1.70908,2.19647,-1.59554,0.47378,0.8434]
    table1[24,*] = [748.710,0.40627,1.22842,-1.67877,2.05535,-1.39972,0.38845,0.8524]
    table1[25,*] = [811.760,0.42977,1.25182,-1.85164,2.31949,-1.59101,0.44155,0.8621]
    table1[26,*] = [869.600,0.45396,1.25101,-2.02958,2.75410,-2.02287,0.59338,0.8701]
    table1[27,*] = [948.850,0.47855,1.19813,-1.86296,2.36939,-1.64367,0.46056,0.8773]
    table1[28,*] = [1046.60,0.49870,1.21429,-2.06976,2.80703,-2.05247,0.60221,0.8841]
    table1[29,*] = [1098.95,0.51149,1.19354,-2.00174,2.66936,-1.94981,0.57715,0.8890]
    
    A[0] = interpol(table1[*, 1], table1[*, 0], lambda*1e9)
    A[1] = interpol(table1[*, 2], table1[*, 0], lambda*1e9)
    A[2] = interpol(table1[*, 3], table1[*, 0], lambda*1e9)
    A[3] = interpol(table1[*, 4], table1[*, 0], lambda*1e9)
    A[4] = interpol(table1[*, 5], table1[*, 0], lambda*1e9)
    A[5] = interpol(table1[*, 6], table1[*, 0], lambda*1e9)

  endelse

  
  A = A/total(A)                ; Unity at mu = 1

  return, A

end
