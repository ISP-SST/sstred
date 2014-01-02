PRO red::getproperty, item, value, HELP = hlp, PRINT = prnt, PTR = ptr

IF keyword_set(hlp) OR n_params() EQ 0 THEN BEGIN
    help, self, /obj, OUT = value
ENDIF ELSE BEGIN
    IF keyword_set(ptr) THEN $
      ok = execute('value = *self.'+item) $
    ELSE $
      ok = execute('value = self.'+item)
ENDELSE

IF keyword_set(prnt) OR n_params() LT 2 THEN $
  print, value, FORMAT = "((a))"

END
