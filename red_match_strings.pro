; docformat = 'rst'

;+
; Examine an array of strings with respect to matching a regular
; expression. 
; 
; :Categories:
;
;    SST observations
; 
; 
; :Author:
; 
;    Mats Löfdahl, 2015-09-21
; 
; 
; :Params:
;
;    strings : in, type=strarr
;
;       The array of strings that you want to examine.
;
;    regex : in, type=string
;
;       The regular expression, the matching status of which you want
;       to examine the strings.
; 
; :Keywords:
;
;    indx_matching : out, optional, type=lonarr
;
;       The indices of the matching strings.
;
;    indx_nonmatching  : out, optional, type=lonarr
;
;       The indices of the non-matching strings.
;
;    matching_strings  : out, optional, type=lonarr
;
;       The matching strings.
;
;    nonmatching_strings  : out, optional, type=lonarr
;
;       The non-matching strings.
;
;    Nmatching  : out, optional, type=lonarr
;
;       The number of matching strings.
;
;    Nnonmatching  : out, optional, type=lonarr
;
;       The number of non-matching strings.
;
;    fold_case  : in, optional, type=boolean
;
;       The comparison is usually case sensitive. Setting the
;       FOLD_CASE keyword causes a case insensitive match to be done
;       instead. 
;
; :History:
;
;    2015-09-21 : MGL. Initial version.
;
;
;       
;
;
;    
;
;
;-
pro red_match_strings, strings $
                               , regex $
                               , indx_matching = indx_matching $
                               , indx_nonmatching = indx_nonmatching $
                               , matching_strings = matching_strings $
                               , nonmatching_strings = nonmatching_strings $
                               , Nmatching = Nmatching $
                               , Nnonmatching = Nnonmatching $
                               , fold_case = fold_case

  indx_matching = where(strmatch(strings, regex, fold_case = fold_case), Nmatching $
                        , complement = indx_nonmatching, Ncomplement = Nnonmatching)

  if arg_present(matching_strings) then matching_strings = strings[indx_matching]
  if arg_present(nonmatching_strings) then nonmatching_strings = strings[indx_nonmatching]

end
