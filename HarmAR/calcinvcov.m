function [Winv, dvec, q] = calcinvcov(V,p,N,Linvt,Winv)   
% Calculate Inverse of AR Covariance Matrix via LDR
% V is NXps ensemble of residuals, Linv & Winv are initialized to I_N
% Source Code: Wasim Malik - 05/2011
% Modified to Speed up BottleNeck in Winv calculation
% Vectorized to enable speedy evaluation of > 1 AR coefficients
% Pavitra Krishnaswamy - 05/2012
%************************************************************************

%% Find AR Coefficients for All Cases
[b, E] = arburg_parvec(V,p);                                                %   Estimate p AR coefficient vectors at once
%b{r}(:) (size 1Xr) and E(r) (size 1X1) = outputs of arburg(V,r) 
    
%% Transpose of Linv matrix as per Levinson-Durbin Recursion
for r = 1:p                                                                 %   upper left corner block of Linv
	rind = r+1; cind = 1:r;                                                 %   row and column indices to write on
    coef = b{r}';                                                           %   coef is rX1
    Linvt(cind,rind) = coef;                                                %   Populate first p+1 rows of Linv [2, pp. 125]
end
for r = 1:N-p-1                                                             %   The rest of Linv
	rind = r+p+1; cind = r+1:r+p;                                           %   row and column indices to write on
    Linvt(cind,rind) = coef;                                                %   Write the AR coefficients for order ARord
end
    
%% Elements of Diagonal Matrix as per Levinson-Durbin Recursion
d1 = xcov(V,0,'unbiased');                                                  %   First entry r_xx[0]
dvec = [d1 E]; Dinv_d = 1./dvec;                                            %   1st row is d1, and 2nd to (p+1)th rows = E
q = E(p); c = 1/q;                                                          %   white noise error variance
        
%% Vectorized Inverse Covariance of AR(p) - Cholesky Form [2, pp. 121]
ind1 = 1:p+1; ind2 = p+2:N;                                                 %   to make block multiplication easy
A = Linvt(ind1,ind1);                                                       %   top left block of L_pr
B1 = Linvt(ind1,ind2);                                                      %   top right block of L_pr
B2 = sparse(Linvt(ind2,ind2));                                              %   bottom right block of L_pr
D = diag(Dinv_d);                                                           %   top left part of diag matrix
    
Winv(ind1,ind1) = A*D*A'+c*(B1*B1');                                        %   top left block of Winv
Winv(ind1,ind2) = c*(B1*B2');                                               %   top right block of Winv
Winv(ind2,ind1) = c*(B2*B1');                                               %   bottom left block of Winv
Winv(ind2,ind2) = c*(B2*B2');                                               %   bottom right block of Winv
end