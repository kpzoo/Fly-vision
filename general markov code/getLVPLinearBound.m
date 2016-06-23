% Function to calculate the LVP bound in the linear encoding case that
% x2dot = lam(x1) = alpha*x1
function [relbnd rawbnd N1 N2] = getLVPLinearBound(alpha, x1mean, kdeath_x1, umean)

% Obtain parameters from inputs
tau1 = 1/kdeath_x1;
N1check = x1mean;
N1 = umean*tau1;
N2 = alpha*x1mean*tau1;

% Obtain linear encoding bound
relbnd = 2/(x1mean*(1 + sqrt(1 + 4*N2/N1)));
rawbnd = relbnd*(x1mean^2);
rawbnd_check = 2*(x1mean^2)/(x1mean*(1 + sqrt(1 + 4*N2/N1check)));

% Display results and check closeness of N1 calculations and effect on bound
disp(['Linear LVP bound: [rel abs] =  ' [num2str(relbnd) ' ' num2str(rawbnd)]]);
disp('********************************************************************');
disp(['The 2 estimates of N1 are: ' [num2str(N1) ' ' num2str(N1check)]]);
disp(['Corresponding abs bounds are: ' [num2str(rawbnd) ' ' num2str(rawbnd_check)]]);
disp('********************************************************************');

