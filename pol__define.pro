pro pol__define
  Nlc = 4
  states = { pol $
             , tfiles:strarr(Nlc) $
             , rfiles:strarr(Nlc) $
             , state:' ' $
             , timg:ptrarr(Nlc,/allocate_heap) $
             , rimg:ptrarr(Nlc,/allocate_heap) $
             , immt:ptr_new() $
             , immr:ptr_new() $
             , pref:' ' $
             , destretch:0B $
             , wb:' ' $
             , wbfiles:strarr(Nlc) $
             , camt:' ' $
             , camr:' ' $
             , camwb:' ' $
             , x0:0L, x1:0L, y0:0L, y1:0L $
             , telog:' ' $
             , ftfiles:strarr(Nlc) $
             , frfiles:strarr(Nlc) $
             , utflat:' ' $
             , urflat:' ' $
             , scan:' ' $
           }
end
