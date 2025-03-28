; docformat = 'rst'

;+
; Split text into multiple lines, not exceeding a specified length. As
; a function, return the text. As a subroutine, print it.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
;   The reflowed text as a strarr.
; 
; :Params:
; 
;   txt : in, type="string or strarr"
; 
;     The text to be reflowed.
; 
; 
; :Keywords:
; 
;   indent :  in, optional, type=string
; 
;     This string, typically a few blanks, is added at the beginning
;     of the first line of text.
; 
;   width : in, optional, type=integer, default="min(60,Terminal width-1)"
;   
;     Make the lines no longer than this.
; 
; 
; :History:
; 
;   2021-04-15 : MGL. First version.
; 
;   2025-03-27 : MGL. New keyword indent.
; 
;-
function red_strflow, txt, width = width, indent = indent

  if n_elements(txt) eq 0 then return, ''
  if n_elements(width) eq 0 then width = ((TERMINAL_SIZE( ))[0]-1) <60
  if n_elements(indent) eq 0 then indent = '' 
  
  istr = strjoin(strtrim(strcompress(txt), 2), ' ')

  if strlen(istr) eq 0 then return, ''
  
  done = !false

  repeat begin

    pos = strpos(istr, ' ')
    if pos eq -1 then begin
      wrd = strtrim(istr, 2)
      done = !true
    endif else begin
      wrd = strmid(istr, 0, pos)
      istr = strmid(istr, pos+1)
    endelse

    case 1 of
        
      n_elements(line) eq 0 : begin
        ;; Start the first line
        line = indent + wrd
      end
      
      strlen(line) + strlen(wrd) + 1 gt width : begin
        ;; Need to make a line break
        red_append, ostr, line
        line = wrd
      end

      else : begin
        ;; Add to the line
        line += ' ' + wrd
      end
      
    endcase
    
  endrep until done

  if n_elements(ostr) eq 0 then return, line

  red_append, ostr, line
  return, ostr
  
end

pro red_strflow, txt, lun = lun, _ref_extra = extra  
  
  if n_elements(lun) eq 0 then begin
    print, red_strflow(txt, _strict_extra = extra), format = '(a0)'
  endif else begin
    printf, lun, red_strflow(txt, _strict_extra = extra), format = '(a0)'
  endelse
  
end

instring = 'The Duke of Cambridge and Duke of Sussex will walk behind their grandfather’s coffin at his funeral on Saturday. They will, however, be separated by the diplomatic presence of their cousin. The brothers, whose fractured relationship has not recovered since their last awkward encounter at Westminster Abbey a year ago, will be among the nine members of the royal family...'

instrings = ["(CNN) If you've been out driving on the eastern coast of Australia in" $
             , "the last few months, you might have seen Tom Drury. He would have been" $
             , "hard to miss, a 28-year-old with a droopy moustache and a backpack," $
             , "cruising along by the side of the Bruce Highway. Cruising, on a" $
             , "skateboard. What you probably wouldn't have known was that he was a" $
             , "long way from home, skating alone on an epic voyage of discovery that" $
             , "led him from Melbourne all the way north to Cairns, a 4,000-kilometer" $
             , "route on just four little wheels. "]

;openw, lun, 'tmp.txt', /get_lun
;red_strflow, instring, w = 40, lun = lun
;free_lun, lun
;
;stop

outstring1 = red_strflow(instring)
outstring2 = red_strflow(instrings, w = 40, indent = '  ')
outstring3 = red_strflow(instrings, w = 30)

print
hprint, outstring1
print
hprint, outstring2
print
hprint, outstring3
print

end
