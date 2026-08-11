; docformat = 'rst'

;+
; List of PRSTEP and associated keywords, as defined in the SOLARNET
; recommendations.
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
;    A list of defined PR*na keywords, sans the na part.
;
; :Keywords:
; 
;   count : out, optional, type=integer
;
;      The number of keywords in the list.
;   
; :History:
; 
;   2021-05-28 : MGL. First version.
; 
;-
function red_headerinfo_prkeys, count = count

  prkeys = ['PRSTEP', 'PRPROC', 'PRMODE', 'PRPARA', 'PRLIB', 'PRVER', 'PRREF', 'PRBRA', 'PRHSH', 'PRENV']
  count = n_elements(prkeys)

  return, prkeys
  
end
