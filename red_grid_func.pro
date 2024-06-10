FUNCTION red_grid_func, p, X = x, Y = y
mat = [[cos(p[1]), sin(p[1])], [-sin(p[1]), cos(p[1])]]
v = mat # (p[0]*x)
return, (y-v)[*]
END
