function [harmar_smooth, harm_smooth, ar_smooth, res_smooth, ...
          beta_smooth, beta_smoothcov, v_smooth, v_smoothcov, beta_crsmoothcov, v_crsmoothcov] ...
          = kalm_filt_smooth(y, Z, A, Psi, R, P, x_init, Wb, We, Irp, checkon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman filter and fixed interval smoother with random walk state model 
% augmented with AR process for estimating beta in colored noise
%
% MODEL
%=======================
% Model for each of t = 1,2,...,N timepoints
%-----------------------
% State Model:
%-----------------------
% beta(:,t) = A*beta(:,t-1) + wb(:,t)                       - R X 1
% v(:,t) = Psi*v(:,t-1) + epsilon(:,t)                      - P X 1
% x(:,t) = [beta(:,t); v(:,t)]                              - R+P X 1
% F = [A 0; 0 Psi]                                          - R+P X R+P
% x(:,t) = F*x(:,t-1) + w(:,t)                              - R+P X 1
% v and w are I.I.D Gaussian
%-----------------------
% Measurement Model:
%-----------------------
% y(t) = Z(t,:)*beta(:,t) + v(t) = [Z(t,:) Ip(t,:)]*x(:,t)  - 1 X 1
% H(t,:) = [Z(t,:) Ip(t,:)]                                 - 1 X R+P
% y(t) = H(t,:)*x(:,t)                                      - 1 X 1
% Measurement eqn noiseless - observation noise is folded into v
%
% INPUT:
%=======================
% y - data observations arranged as a column vector in time - N X 1
% Z - matrix of harmonic regressors                         - N X R
% A - drift of harmonic state                               - R X R
% Psi - state transition matrix for AR Akaike form          - P X P
% R is # harmonics                                          - 1 X 1
% P is # AR terms                                           - 1 X 1
% x_init - harmonic amplitude state and AR series at time 0 - R+P X 1
% Wb - state noise covariance for harmonics                 - R X R
% We - state noise covariance for AR series                 - P X P
% Irp - identity matrix of size R+P                         - R+P X R+P
%
% COMPUTE IN CODE:
%=======================
% F = [A 0; 0 Psi]                                          - R+P X R+P
% W - state noise covariance                                - R+P X R+P
% Ip - identity matrix for measurement equation             - N X P
% H(t,:) = [Z(t,:) Ip(t,:)]                                 - 1 X R+P
%
% OUTPUT:
%=======================
% beta_smooth - estimated harmonic amplitudes               - R X N
% beta_smoothcov - error covariance for harmamp estimates   - N X R X R 
% beta_smoothcrcov - error covariance 1 step lag cross terms- N X R X R
% v_smooth - estimated AR series                            - P X N
% v_smoothcov - error covariance for AR estimates           - N X P X P
% v_smoothcrcov - error covariance 1 step lag cross terms   - N X P X P
% harm_smooth - estimated harmonic fit                      - N X 1
% ar_smooth - estimated AR  signal                          - N X 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in Measurements and Initialize State Space Matrices
N = length(y); init = zeros(size(y));
F = [A zeros(R,P); zeros(P, R) Psi];                                        %   augmented state transition matrix
W = zeros(R+P,R+P);                                                         %   initialize augmented state noise covariance
W(1:R,1:R) = Wb; W(R+1:end, R+1:end) = We;                                  %   set augmented state noise covariance

%% Initialize State Estimate and Covariance
x_last = x_init;                                                            %   initialize augmented state
xcov_last = W;                                                              %   initialize augmented error covariance
x_predict = zeros(R+P,N); xcov_predict = zeros(R+P,R+P,N);                  %   initialize predictions
x_update = x_predict;  xcov_update =  xcov_predict;                         %   initialize updates

%% Initialize Output Variables
harmar_smooth = init; harm_smooth = init; ar_smooth = init; res_smooth = init;
beta_smooth = zeros(R,N); v_smooth = zeros(P,N); 
beta_smoothcov = zeros(R,R,N); v_smoothcov = zeros(P,P,N);
beta_crsmoothcov = zeros(R,R,N); v_crsmoothcov = zeros(P,P,N);

%% Kalman Filter - Forward Pass
for i = 1:N
    % Prediction Steps
    x_kalp  = F*x_last;                                                     %   predicted estimates    - R+P X 1
    xcov_kalp = F*xcov_last*F' + W;                                         %   predicted error covar	- R+P X R+P
    if rcond(xcov_kalp) < 1e-12
        %   we have a singular matrix here.
        condy = rcond(xcov_kalp);
        disp(condy)
    end
    x_predict(:,i) = x_kalp;                                                %   predicted state estimates - x_i|i-1 - R+P X 1
    xcov_predict(:,:,i) = xcov_kalp;                                        %   predicted state error covariance - Sigma_i|i-1 - R+P X R+P X 1
    
    % Kalman Innovations and Gain
    Hi = [Z(i,:) 1 zeros(1, P-1)];                                          %   meas process vector    - 1 X R+P
    innovi = y(i) - Hi*x_kalp;                                              %   innovation             - 1 X 1
    Ki = (xcov_kalp*Hi')/(Hi*xcov_kalp*Hi');                                %   Kalman gain            - R+P X 1
    
    % Update Steps
    x_last = x_kalp + Ki*innovi;                                            %   updated estimates      - R+P X 1 
    xcov_last = (Irp-Ki*Hi)*xcov_kalp;                                      %   updated error covar    - R+P X R+P    
    x_update(:,i) = x_last;                                                 %   updated state estimates - x_i|i - R+P X 1
    xcov_update(:,:,i) = xcov_last;                                         %   updated state error covariance - Sigma_i|i - R+P X R+P X 1
end
clear x_kalp xcov_kalp Hi Ki innovi x_last xcov_last

%% Kalman Smoother - Backward Pass
x_smooth(:,N) = x_update(:,N); xcov_smooth(:,:,N) = xcov_update(:,:,N); 
for j = N-1:-1:1
    xcov_last = squeeze(xcov_update(:,:,j));                                %   read in relevant updated cov
    xcov_kalp = squeeze(xcov_predict(:,:,j+1));                             %   read in relevant predicted cov
    xcov_sm = squeeze(xcov_smooth(:,:,j+1));                                %   read in relevant smoothed cov
    
    Jt = xcov_last*F'/xcov_kalp;                                            %   smoother "gain" - R+P X R+P
    x_smooth(:,j) = x_update(:,j) + Jt*(x_smooth(:,j+1)- x_predict(:,j+1)); %   smoothed state esimates - R+P X 1
    xcov_smooth(:,:,j) = xcov_last + Jt*(xcov_sm - xcov_kalp)*Jt';          %   smoothed error covar - R+P X R+P X 1
    xcov_crsmooth(:,:,j+1)= xcov_sm*Jt';                                    %   one-step cross covariance - R+P X R+P X 1
end
    
%% Store Final Results
for i = 1:N
    Hi = [Z(i,:) 1 zeros(1, P-1)];                                          %   meas process vector    - 1 X R+P
    harmar_smooth(i) = Hi*x_smooth(:,i);                                    %   harmonic + AR fit      - 1 X 1
    harm_smooth(i) = Hi(1,1:R)*x_smooth(1:R,i);                             %   harmonic fit           - 1 X 1
    ar_smooth(i) = Hi(1,R+1:end)*x_smooth(R+1:end,i);                       %   AR fit                 - 1 X 1
    res_smooth(i) = y(i)-harmar_smooth(i);                                  %   residual noise series  - 1 X 1
    
    beta_smooth(:,i) = x_smooth(1:R,i);                                     %   store smoothed harmamp estimates- R X 1
    beta_smoothcov(:,:,i) = xcov_smooth(1:R,1:R,i);                         %   store smoothed harmamp error covariance - R X R X 1
    v_smooth(:,i) = x_smooth(R+1:end,i);                                    %   store smoothed AR estimates - same as ar_hat - P X 1 
    v_smoothcov(:,:,i) = xcov_smooth(R+1:end,R+1:end,i);                    %   store smoothed AR error covariance - P X P X 1
    beta_crsmoothcov(:,:,i) = xcov_crsmooth(1:R,1:R,i);                     %   store smoothed harmamp error covariance cross terms - R X R X 1
    v_crsmoothcov(:,:,i) = xcov_crsmooth(R+1:end,R+1:end,i);                %   store smoothed AR error covariance cross terms - P X P X 1
end

% % Save Final Results
if checkon 
    save('test.mat','y',...
	'harmar_smooth','harm_smooth','ar_smooth','res_smooth',...
	'beta_smooth','beta_smoothcov','v_smooth','v_smoothcov',...
	'beta_crsmoothcov', 'v_crsmoothcov');
end
end