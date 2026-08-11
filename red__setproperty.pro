PRO red::setproperty, item, value, PTR = ptr

cmd = 'self.'+item+' = value'
IF keyword_set(ptr) THEN cmd = '*'+cmd
ok = execute(cmd)

END
