%HASPAR Checks whether parallel toolbox is installed
%
%   Output:
%   TRUE if parallel toolbox is installed
%   haspar returns 1 if parallel toolbox is installed and 0 otherwise
%
%   See also configpar
%
%   Copyright 2011 Michael J. Prerau, Ph.D.
%   
%   Last modified 01/07/2011
%********************************************************************
function result = haspar()
%Check if the parallel toolbox has been installed
result=~isempty(ver('distcomp'));
end