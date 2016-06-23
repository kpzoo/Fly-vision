% Note - use Kolmogorov forward equation as dP/dt = PQ as opposed to Neil's
% version and corrects for the sign change by Neil

% Function to solve ODEs on the Snyder filter between x2 birth events -
% uses Q matrix and the hybrid stochastic method but keep dq/dt vs dq/dI
function dy = odeSnyderOffline5(ts, y, B, coeff, S)

% Boolean to check sizes
check_on = 0;

% Check array sizes if inner boolean set
if check_on
    disp(['Dim(y) = ' num2str(size(y))]);
    disp(['Dim(C) = ' num2str(size(C))]);
end

% Solve differential equation set and check size
% y(end) = 1 - sum(y(1:end-1));
C = B + sum(coeff(1)*y'*S)*eye(length(diag(S)));
dy = y'*C;
% dy(end) = - sum(dy(1:end-1));
if check_on
    disp(['Dim(dy) = ' num2str(size(dy))]);
end
dy = dy';


