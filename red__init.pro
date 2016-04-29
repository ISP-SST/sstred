function red::init, filename

    ; This function is called implicitly when an instance is created. 
        
    if n_elements(filename) eq 0 then filename = 'config.txt'
    
    if file_test(filename) then self->initialize, filename
    
    return,1
    
end
