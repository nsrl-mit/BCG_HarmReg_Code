function S = paesa(varargin)
%   PAESA (PreAllocate Empty Struct Array)
%       
%       PAESA(N) returns a struct array of length N with no defined fields
%       in each struct.
%       
%       PAESA(N,'field1','field2',...) returns a struct array of length N
%       with fields field1, field2,... initialized as empty values.
%       
%       PAESA('field1','field2',...) returns a structure with empty fields
%       field1, field 2,...
%       
%       PAESA returns an empty struct.
%
%
%   \////\
%   C:  1136.02.19.2015
%   U:  1136.02.19.2015

if nargin < 1                   %   no input arguments
    S = struct;
else
    if isnumeric(varargin{1});  %   if 1st argument is a number
        N = varargin{1};
        nStart = 2;             %   start looping at index 2
    else
        N = 1;
        nStart = 1;             %   1st argument is a field, length is 1
    end
    s = struct();
    for n = nStart:nargin
        s.(varargin{n}) = [];   %   tack each input to the struct
    end
    S = repmat(s,N,1);          %   replicate (way faster than a for-loop)
end