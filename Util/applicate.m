function z = applicate(x,y,n)
w = fApplyN(@log,n,x)+fApplyN(@log,n,y);
z = fApplyN(@exp,n,w);