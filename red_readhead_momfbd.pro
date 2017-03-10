; docformat = 'rst'

;+
; Return a fits-type header from a MOMFBD format file.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    fname : in, type=string
;
;       The name of the data file.
; 
; :History:
; 
;    2017-03-10 : MGL. Moved reading of headers from MOMFBD format
;                 files from red_readhead.pro.
; 
; 
; 
;-
function red_readhead_momfbd, fname

  compile_opt idl2

  mr = momfbd_read(fname, /names) ; Use /names to avoid reading the data parts
  mkhdr, header, 4, [mr.clip[0,1,1]-mr.clip[0,1,0]+1 $
                     , mr.clip[0,0,1]-mr.clip[0,0,0]+1]
  header = header[where(header ne replicate(' ',80))] ; Remove blank lines

  date_ave = mr.date + 'T' + mr.time
  sxaddpar, header, 'DATE-AVE', date_ave, ' ', before='COMMENT'

  return, header
        
end
