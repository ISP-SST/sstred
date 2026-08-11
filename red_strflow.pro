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
;   before_line : in, optional, type=string
; 
;     Print this as a line before printing txt.
; 
;   after_line : in, optional, type=string
; 
;     Print this as a line after printing txt.
; 
;   indents : in, optional, type="string or number(s)"
; 
;     This string, is added at the beginning of the first line of
;     text. If array, first element indents first line, second element
;     indents all following lines. If number(s), indent by that many
;     blanks.
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
;   2025-03-27 : MGL. New keywords before_line, after_line. Keyword
;                indents can now be a 2-element array.
; 
;-
function red_strflow, txt $
                      , after_line = after_line $
                      , before_line = before_line $
                      , indents = indent $
                      , width = width 

  if n_elements(txt) eq 0 then return, ''
  if n_elements(width) eq 0 then width = ((TERMINAL_SIZE( ))[0]-1) <60
  if n_elements(indent) eq 0 then indent = ''

  case n_elements(indent) of

    0 : indents = ['', '']

    1 : begin
      case !true of
        isa(indent, /string) : indents = [indent, '']
        isa(indent, /number) : indents = [string(replicate(32b, round(indent))), '']
        else : stop
      endcase
    end

    2 : begin
      case !true of
        isa(indent, /string) : indents = indent
        isa(indent, /number) : indents = [string(replicate(32b, round(indent[0]))), $
                                          string(replicate(32b, round(indent[1])))  ]
        else : stop
      endcase
    end

    else : stop
    
  endcase
  
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
        line = indents[0] + wrd
      end
      
      strlen(line) + strlen(wrd) + 1 gt width : begin
        ;; Need to make a line break
        red_append, ostr, line
        line = indents[1] + wrd
      end

      else : begin
        ;; Add to the line
        line += ' ' + wrd
      end
      
    endcase
    
  endrep until done

  red_append, ostr, line

  if n_elements(before_line) gt 0 then ostr = [before_line, ostr]
  if n_elements(after_line) gt 0 then ostr = [ostr, after_line]
  
  
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

outstring1 = red_strflow("(CNN) If you've been out driving")
outstring2 = red_strflow(instrings, w = 40, indent = '  ')
outstring3 = red_strflow(instrings, indents = [7, 2])
outstring4 = red_strflow(instrings, indents = ['+ ', '-   '], before = '-----', after = ' + + + +')

print
hprint, outstring1
print
hprint, outstring2
print
hprint, outstring3
print
hprint, outstring4
print

end
