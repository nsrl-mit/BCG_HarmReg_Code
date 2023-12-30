function tok = strntok(s,n,d)
%   take the nth token (substring) in input string s after being cut up by
%   delimiter(s) d.
sCut    = regexp(s,d,'split');
tok     = sCut{n};