function [xn,indn] = findNearest(x,arr)
N = length(x);
for n = 1:N    
    dif = abs(arr-x(n));
    [~,indn(n)] = min(dif);
    xn(n) = arr(indn(n));
end