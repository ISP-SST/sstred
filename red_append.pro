pro red_append, array, data
    if ~n_elements(data) then return
    if n_elements(array) eq 0 then array = [ data ] else array = [ array, data ]
end
