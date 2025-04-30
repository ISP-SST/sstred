; docformat = 'rst'

;+
; Provide information about SST science cameras.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author: 
;
;    Mats Löfdahl, ISP, 2016-05-01
; 
; 
; :Returns:
;
;    A struct with info about the camera.
; 
; 
; :Params:
; 
;    x : in, type="string or integer"
;   
;      The number of the camera. If x is a string, assume it is the
;      camera number as a roman numeral, possibly prepended with the
;      string "cam".
; 
; :History:
; 
;   2016-05-01 : MGL. First version.
; 
;   2016-05-04 : MGL. Added PointGrey cameras. Allow arrays of camera
;                numbers as input.
; 
;-
function red_camerainfo, x

  if size(x,/n_dim) eq 0 then begin
    
    if size(x, /tname) eq 'STRING' then begin
      ;; Called with string
      camromnum = red_strreplace(x, 'cam', '') ; Remove "cam" if needed.
      camnum    = red_romannumber(camromnum)   ; The integer camera number.
    endif else begin
      ;; Called with number
      camromnum = red_romannumber(x) ; The roman camera number.
      camnum    = x
    endelse
    
    case camnum of
      4:  return, {romnum:'IV', $
                   defined:1B, $
                   model:'MegaPlus 1.6', $
                   xsize:1534L, $
                   ysize:1024L, $
                   pixelsize:9.0e-6, $
                   serialnumber:'22743JS', $
                   use:'', $
                   note:''}
      6:  return, {romnum:'VI', $
                   defined:1B, $
                   model:'MegaPlus 1.6', $
                   xsize:1534L, $
                   ysize:1024L, $
                   pixelsize:9.0e-6, $
                   serialnumber:'22099M4', $
                   use:'', $
                   note:''}
      8:  return, {romnum:'VIII', $
                   defined:1B, $
                   model:'MegaPlus 4.2i/10', $
                   xsize:2029L, $
                   ysize:2044L, $
                   pixelsize:9.0e-6, $
                   serialnumber:'62981G8CSY3B', $
                   use:'', $
                   note:'Broken internal shutter, needs external Uniblitz shutter.' $
                   + ' Stored in tower in Al case. "Lacking bits."'}
      9:  return, {romnum:'IX', $
                   defined:1B, $
                   model:'MegaPlus 4.2i/10', $
                   xsize:2029L, $
                   ysize:2044L, $
                   pixelsize:9.0e-6, $
                   serialnumber:'6298068CSY3E', $
                   use:'', $
                   note:'UV sensitive coating. Suspected to give slightly blurry images.'}
      10: return, {romnum:'X', $
                   defined:1B, $
                   model:'MegaPlus 1.6i', $
                   xsize:1534L, $
                   ysize:1024L, $
                   pixelsize:9.0e-6, $
                   serialnumber:'', $
                   note:'Spectroscopy/slit-jaw. UV Sensitive - blue plus chip'}
      11: return, {romnum:'XI', $
                   defined:1B, $
                   model:'MegaPlus 6.3i/10', $
                   xsize:3072L, $
                   ysize:2048L, $
                   pixelsize:9.0e-6, $
                   serialnumber:'64142S00FSYB', $
                   use:'', $
                   note:'UV Sensitive - blue plus chip'}
      12: return, {romnum:'XII', $
                   defined:1B, $
                   model:'MegaPlus 1.6i', $
                   xsize:1534L, $
                   ysize:1024L, $
                   pixelsize:9.0e-6, $
                   serialnumber:'61316 M6CRY1', $
                   use:'Spectroscopy/slit-jaw', $
                   note:''}
      13: return, {romnum:'XIII', $
                   defined:1B, $
                   model:'MegaPlus 1.6i', $
                   xsize:1534L, $
                   ysize:1024L, $
                   pixelsize:9.0e-6, $
                   serialnumber:'61317 M6CRY1', $
                   use:'Spectroscopy/slit-jaw', $
                   note:'Ugly orange peel pattern'}
      14: return, {romnum:'XIV', $
                   defined:1B, $
                   model:'MegaPlus 1.6i', $
                   xsize:1534L, $
                   ysize:1024L, $
                   pixelsize:9.0e-6, $
                   serialnumber:'64061EOOCSY3', $
                   use:'Spectroscopy/slit-jaw', $
                   note:'UV Sensitive - blue plus chip.'}
      15: return, {romnum:'XV', $
                   defined:1B, $
                   model:'MegaPlus II es1603', $
                   xsize:1536L, $
                   ysize:1024L, $
                   pixelsize:9.0e-6, $
                   serialnumber:'07000014 M', $
                   use:'Spectroscopy/slit-jaw', $
                   note:'Dark level problem. Cover edge, disable black subtraction.'}
      17: return, {romnum:'XVII', $
                   defined:1B, $
                   model:'MegaPlus II es1603', $
                   xsize:1536L, $
                   ysize:1024L, $
                   pixelsize:9.0e-6, $
                   serialnumber:'07000015 M', $
                   use:'Spectroscopy/slit-jaw', $
                   note:'Dark level problem. Cover edge, disable black subtraction'}
      18: return, {romnum:'XVIII', $
                   defined:1B, $
                   model:'Sarnoff CAM1M100', $
                   xsize:1024L, $
                   ysize:1024L, $
                   pixelsize:16.0e-6, $
                   serialnumber:'SAS 2', $
                   use:'CRISP', $
                   note:'Has new red AR coating. Transparent Sarnoff chip problem.'}
      19: return, {romnum:'XIX', $
                   defined:1B, $
                   model:'Sarnoff CAM1M100', $
                   xsize:1024L, $
                   ysize:1024L, $
                   pixelsize:16.0e-6, $
                   serialnumber:'SAS 3', $
                   use:'CRISP', $
                   note:'Has new red AR coating. Transparent Sarnoff chip problem.'}
      20: return, {romnum:'XX', $
                   defined:1B, $
                   model:'Sarnoff CAM1M100', $
                   xsize:1024L, $
                   ysize:1024L, $
                   pixelsize:16.0e-6, $
                   serialnumber:'SAS 1', $
                   use:'CRISP', $
                   note:'Has new red AR coating. Transparent Sarnoff chip problem.'}
      21: return, {romnum:'XXI', $
                   defined:1B, $
                   model:'MegaPlus II es4020', $
                   xsize:2048L, $
                   ysize:2048L, $
                   pixelsize:7.4e-6, $
                   serialnumber:'03000391M', $
                   use:'Blue beam', $
                   note:''}
      22: return, {romnum:'XXII', $
                   defined:1B, $
                   model:'MegaPlus II es4020', $
                   xsize:2048L, $
                   ysize:2048L, $
                   pixelsize:7.4e-6, $
                   serialnumber:'03000392M', $
                   use:'Blue beam', $
                   note:''}	 
      23: return, {romnum:'XXIII', $
                   defined:1B, $
                   model:'MegaPlus II es4020', $
                   xsize:2048L, $
                   ysize:2048L, $
                   pixelsize:7.4e-6, $
                   serialnumber:'03000404M', $
                   use:'Blue beam', $
                   note:''}
      24: return, {romnum:'XXIV', $
                   defined:1B, $
                   model:'MegaPlus II es4020', $
                   xsize:2048L, $
                   ysize:2048L, $
                   pixelsize:7.4e-6, $
                   serialnumber:'03000581M', $
                   use:'Blue beam', $
                   note:''}
      25: return, {romnum:'XXV', $
                   defined:1B, $
                   model:'Sarnoff CAM1M100', $
                   xsize:1024L, $
                   ysize:1024L, $
                   pixelsize:16.0e-6, $
                   serialnumber:'SAS 4', $
                   use:'CRISP', $
                   note:'New for 2008 season, has red AR coating. Transparent Sarnoff chip problem.'}
      26: return, {romnum:'XXVI', $
                   defined:1B, $
                   model:'MegaPlus II es4020', $
                   xsize:2048L, $
                   ysize:2048L, $
                   pixelsize:7.4e-6, $
                   serialnumber:'03000831M', $
                   use:'Blue beam', $
                   note:'Camera without window to improve fringing problems'}
      27: return, {romnum:'XXVII', $
                   defined:1B, $
                   model:'PointGrey GS3-U3-23S6M-C', $
                   xsize:1920L, $
                   ysize:1200L, $
                   pixelsize:5.86e-6, $
                   serialnumber:'15452652', $
                   use:'CHROMIS', $
                   note:''}
      28: return, {romnum:'XXVIII', $
                   defined:1B, $
                   model:'PointGrey GS3-U3-23S6M-C', $
                   xsize:1920L, $
                   ysize:1200L, $
                   pixelsize:5.86e-6, $
                   serialnumber:'14471636', $
                   use:'CHROMIS', $
                   note:''}
      29: return, {romnum:'XXIX', $
                   defined:1B, $
                   model:'PointGrey GS3-U3-23S6M-C', $
                   xsize:1920L, $
                   ysize:1200L, $
                   pixelsize:5.86e-6, $
                   serialnumber:'15452653', $
                   use:'CHROMIS', $
                   note:''}
      30: return, {romnum:'XXX', $
                   defined:1B, $
                   model:'PointGrey GS3-U3-23S6M-C', $
                   xsize:1920L, $
                   ysize:1200L, $
                   pixelsize:5.86e-6, $
                   serialnumber:'15452648', $
                   use:'CHROMIS', $
                   note:''}
      31: return, {romnum:'XXXI', $
                   defined:1B, $
                   model:'XIMEA MX262RG-GP-X8G3-MTP-LA', $
                   xsize:5120L, $
                   ysize:5120L, $
                   pixelsize:2.5e-6, $
                   serialnumber:'XWSMW2226000', $
                   use:'CRISP', $
                   note:''}
      32: return, {romnum:'XXXII', $
                   defined:1B, $
                   model:'XIMEA MX262RG-GP-X8G3-MTP-LA', $
                   xsize:5120L, $
                   ysize:5120L, $
                   pixelsize:2.5e-6, $
                   serialnumber:'XWSMW2226001', $
                   use:'CRISP', $
                   note:''}
      33: return, {romnum:'XXXIII', $
                   defined:1B, $
                   model:'XIMEA MX262RG-GP-X8G3-MTP-LA', $
                   xsize:5120L, $
                   ysize:5120L, $
                   pixelsize:2.5e-6, $
                   serialnumber:'XWSMW2226002', $
                   use:'CRISP', $
                   note:''}
      34: return, {romnum:'XXXIV', $
                   defined:1B, $
                   model:'XIMEA MX262RG-GP-X8G3-MTP-LA', $
                   xsize:5120L, $
                   ysize:5120L, $
                   pixelsize:2.5e-6, $
                   serialnumber:'XWSMW2326000', $
                   use:'CRISP', $
                   note:''}
      35: return, {romnum:'XXXV', $
                   defined:1B, $
                   model:'XIMEA MX262RG-GP-X8G3-MTP-LA', $
                   xsize:5120L, $
                   ysize:5120L, $
                   pixelsize:2.5e-6, $
                   serialnumber:'XWSMW2326001', $
                   use:'CRISP', $
                   note:''}
      36: return, {romnum:'XXXVI', $
                   defined:1B, $
                   model:'XIMEA MX203MG-SY-X4G3-FF', $
                   xsize:4512L, $
                   ysize:4512L, $
                   pixelsize:2.74e-6, $
                   serialnumber:'XINGF2448000', $
                   use:'CHROMIS', $
                   note:''}
      37: return, {romnum:'XXXVII', $
                   defined:1B, $
                   model:'XIMEA MX203MG-SY-X4G3-FF', $
                   xsize:4512L, $
                   ysize:4512L, $
                   pixelsize:2.74e-6, $
                   serialnumber:'XINGF2448002', $
                   use:'CHROMIS', $
                   note:''}
      38:  return, {romnum:'XXXVIII', $
                    defined:1B, $
                    model:'XIMEA MX203MG-SY-X4G3-FF', $
                    xsize:4512L, $
                    ysize:4512L, $
                    pixelsize:2.74e-6, $
                    serialnumber:'XINGF2448003', $
                    use:'CHROMIS', $
                    note:''}
      39:  return, {romnum:'XXXVIX', $
                    defined:1B, $
                    model:'XIMEA MX203MG-SY-X4G3-FF', $
                    xsize:4512L, $
                    ysize:4512L, $
                    pixelsize:2.74e-6, $
                    serialnumber:'XINGF2448004', $
                    use:'CHROMIS', $
                    note:''}
      40:  return, {romnum:'XL', $
                    defined:1B, $
                    model:'XIMEA MX203MG-SY-X4G3-FF', $
                    xsize:4512L, $
                    ysize:4512L, $
                    pixelsize:2.74e-6, $
                    serialnumber:'XINGF2448006', $
                    use:'CHROMIS', $
                    note:''}
      else: return, {romnum:'', $
                     defined:0B, $
                     model:'', $
                     xsize:0L, $
                     ysize:0L, $
                     pixelsize:0.0, $
                     serialnumber:'', $
                     use:'', $
                     note:''}
    endcase
  endif else begin

    Ncams = n_elements(x)
    
    caminfo = replicate({romnum:'', $
                         defined:0B, $
                         model:'', $
                         xsize:0L, $
                         ysize:0L, $
                         pixelsize:0.0, $
                         serialnumber:'', $
                         use:'', $
                         note:''}, Ncams)

    for i = 0, Ncams-1 do caminfo[i] = red_camerainfo(x[i])

    if size(x,/n_dim) gt 1 then begin
      ;; Not a one-dimensional array.
      caminfo = reform(caminfo, size(x, /dim))
    endif
    
    return, caminfo
    
  endelse

end
