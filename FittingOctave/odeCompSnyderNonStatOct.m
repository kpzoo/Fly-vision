% Modified to apply to the non-stationary case in which a z input is
% required to construct the lam matrix

% Function to solve ODEs on the Snyder filter between observed birth events
% with compensation for delayed observations (exponential) - note y is a
% column vector on input
function dy = odeCompSnyderNonStatOct(y, ts, Q, zt, dimQ, b, lam)
% OCTAVE: reordered inputs from odeSnyder(ts, y, ..) to odeSnyder(y, ts, ..)

% Boolean to check sizes
check_on = 0;

% Obtain lamhat estimate
lamcap = y'*lam*ones(dimQ, 1);
% lamhat = lamcap*eye(dimQ);
% OCTAVE: put sum command on lamhat
lamhat = sum(lamcap)*eye(dimQ);

% Solve differential equation set and check size
dy = y'*(Q + lamhat - lam);
dy = dy';

% Check array sizes if inner boolean set
if check_on
    disp(['Dim(y) = ' num2str(size(y))]);
    disp(['Dim(C) = ' num2str(size(C))]);
    disp(['Dim(dy) = ' num2str(size(dy'))]);
end


