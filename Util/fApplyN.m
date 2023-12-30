function y = fApplyN(f,n,x)
if n > 1
    y = fApplyN(f,n-1,x);
else
    y = f(x);
end