function [cSorted,isi] = sortBySubstr(c,si,ds)
%   sorts a cell array of stereotyped strings of the form [ss1 d1 ss2 d2
%   ... dn ssn] where ss are substrings and d are delimiters (for
%   tokenization)
%   c   :: input cell array
%   si  :: sorting index (strings sorted (after conversion to num type) by
%       substring of index si)
%   ds  :: delimiters
cSplit      = regexp(c,ds,'split');
len         = length(cSplit);
i2s         = zeros(len,1);
for ci = 1:len
    i2s(ci) = str2double(cSplit{ci}{si});   %   "indeces2sort"
end
[~,isi]     = sort(i2s);
cSorted     = c(isi);