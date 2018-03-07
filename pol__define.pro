PRO pol__define
nlc = 4
states = {pol, tfiles:strarr(nlc), rfiles:strarr(nlc),state:' ', $
            timg:ptrarr(nlc,/allocate_heap), rimg:ptrarr(nlc,/allocate_heap), $
            immt:ptr_new(), immr:ptr_new(), pref:' ', destretch:0B, wb:' ',$
            wbfiles:strarr(nlc), camt:' ', camr:' ', camwb:' ',$
            x0:0L, x1:0L, y0:0L, y1:0L, telog:' ', ftfiles:strarr(nlc), $
            frfiles:strarr(nlc), utflat:' ', urflat:' ', scan:' ', ftype:'', $
            tclip:lonarr(4), rclip:lonarr(4), xotfile:' ', yotfile:' ', xorfile:' ', yorfile:' '}
END
