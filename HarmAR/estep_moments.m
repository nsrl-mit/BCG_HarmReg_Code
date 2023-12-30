function [tr_E_bt_bt_N, tr_E_btm1_ATA_btm1_N, tr_E_bt_A_btm1_N] = ...
    estep_moments(A, R, N, Wb, beta_smooth, x_init, ErrCov)
% Use Estimates from Kalman Smoother to Compute Moments for E-Step
% Pavitra Krishnaswamy 12/2014
%************************************************************************

%% Initialize Matrices
% Read In Error Covariances
beta_init = x_init(1:R,1);
beta_smoothcov = ErrCov.Amps;
beta_crsmoothcov = ErrCov.Amps_cr;
% Initialize for Computations Here
init1 = zeros(R,R,N);	E_bt_bt_N = init1;      E_btm1_btm1_N = init1;          E_btm1_bt_N = init1;
init2 = zeros(1,N);     tr_E_bt_bt_N = init2;   tr_E_btm1_ATA_btm1_N = init2;   tr_E_bt_A_btm1_N = init2;

%% For First Time Point
E_bt_bt_N(:,:,1) = beta_smoothcov(:,:,1) + beta_smooth(:,1)*beta_smooth(:,1)';
E_btm1_btm1_N(:,:,1) = Wb + beta_init*beta_init';
E_btm1_bt_N(:,:,1) = beta_crsmoothcov(:,:,1) + beta_init*beta_smooth(:,1)';

%% Compute Expectations for Each Time Point
for i = 2:N
    E_bt_bt_N(:,:,i) = beta_smoothcov(:,:,i) + beta_smooth(:,i)*beta_smooth(:,i)';
    E_btm1_btm1_N(:,:,i) = beta_smoothcov(:,:,i-1) + beta_smooth(:,i-1)*beta_smooth(:,i-1)';
    E_btm1_bt_N(:,:,i) = beta_crsmoothcov(:,:,i) + beta_smooth(:,i-1)*beta_smooth(:,i)';
end

%% Compute Traces of Expectations for Each Time Point
for i = 1:N
    tr_E_bt_bt_N(i) = trace(E_bt_bt_N(:,:,i));
    tr_E_btm1_ATA_btm1_N(i) = trace(E_btm1_btm1_N(:,:,i)*(A'*A));
    tr_E_bt_A_btm1_N(i) = trace(E_btm1_bt_N(:,:,i)*A);
end

end