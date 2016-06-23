% Note - use Kolmogorov forward equation as dP/dt = PQ as opposed to Neil's
% version and corrects for the sign change by Neil

% Function to solve ODEs on the Snyder filter between x2 birth events -
% uses Q matrix and the hybrid stochastic method but keep dq/dt vs dq/dI
function dy = odeSnyder(ts, y, B, coeff, S, birStr2)

% Boolean to check sizes
check_on = 0;

% Obtain C matrix from current estimate of lam
[lamhat x1hat] = calclamEst(birStr2, coeff, y', S);
C = B + lamhat*eye(size(S));

% Check array sizes if inner boolean set
if check_on
    disp(['Dim(y) = ' num2str(size(y))]);
    disp(['Dim(C) = ' num2str(size(C))]);
end

% Solve differential equation set and check size
dy = y'*C;
if check_on
    disp(['Dim(dy) = ' num2str(size(dy))]);
end
dy = dy';


