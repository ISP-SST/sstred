function minmax_filter, var
  n = n_elements(var)
  return, (sort(var))[1L:n-2]
end
