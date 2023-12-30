function choi = listInput(list)
%   give it a cell array of various choices. this function prints them all
%   out ~nicely and let's you pick from the list. choi is the list item
%   that you chose.
%   easy!
if ~iscell(list)
    error('ERR:     Input list must be a cell array.');
end
L       = length(list);
fprintf('\nplease choose from the following:\n')
for li = 1:L
    fprintf('\n(%d)\t%s',li,list{li})
end
chi     = input('\n\npick a number:  ');
choi    = list{chi};
%%  things to fix:
%   -   align all of the tabs
%   -   dealing more reasonably with heterogeneous data input types
%   -   ...?