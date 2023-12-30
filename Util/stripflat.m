function fsf = stripflat(S,f)
%   STRIPFLAT takes a stuct array S with field f and returns an array fsf 
%   (field, stripped and flattened) s.t. fsf(k) = S(k).f. The size and type 
%   of f must be homogeneous throughout the struct S.
%
%
%   \////\
%   C: 1153.02.19.2015
%   U: 1153.02.19.2015

if ~isfield(S,f)
    throw(error('f is not a field in the struct array S'))
else
    N = length(S);
    dimF = size(S(1).f);
    for n = 1:N
        fsf(:,n) = S(n).(f);
    end
end
    