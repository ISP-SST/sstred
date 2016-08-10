; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    img : 
;   
;   
;   
;    x : 
;   
;   
;   
;    y : 
;   
;   
;   
; 
; :Keywords:
; 
;    aligned : in, optional, type=boolean
;
;       Get the location of the patch in coordinates after clipping and flipping the image.
;       These are the coordinates that, e.g., momfbd uses.
;
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2016-08-09 : Re-write to also handle channels with reversed align_clip
; 
; 
;-
function red_getclips, img, x, y, aligned=aligned

    
    channel=0     ; do ever need to concern ourselves with multiple channels here?
    
    ; patch coordinates are given relative to the cutout area (align-clip)
    ; of the reference channel
    clip = [ img.patch[x,y].xl, $
             img.patch[x,y].xh, $
             img.patch[x,y].yl, $
             img.patch[x,y].yh ]-1  ; momfbd coordinates are 1-based

    ; add the local shift (i.e. the integer part from the offset files)
    clip[0:1] += img.patch[x,y].dx[channel]
    clip[2:3] += img.patch[x,y].dy[channel]
    
    ; subtract the total shift induced by the momfbd code by shifting the subimages.
    clip[0:1] -= img.patch[x,y].offx
    clip[2:3] -= img.patch[x,y].offy

    if keyword_set(aligned) then return, clip

    if img.clip[0,0,0] gt img.clip[0,0,1] then begin    ;  x-reversed
        clip[0:1] = reverse(img.clip[0,0,0] - 1 - clip[0:1])
    endif else begin
        clip[0:1] += img.clip[0,0,0]-1
    endelse
    
    if img.clip[0,1,0] gt img.clip[0,1,1] then begin    ;  y-reversed
        clip[2:3] = reverse(img.clip[0,1,0] - 1 - clip[2:3])
    endif else begin
        clip[2:3] += img.clip[0,1,0]-1
    endelse


  return, clip
end
