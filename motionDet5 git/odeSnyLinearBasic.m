% Function to solve linear ODEs on the Snyder filter between observed
% births with - note y is a column vector on input
function dy = odeSnyLinearBasic(ts, y, Q, lam)

% Boolean to check sizes
check_on = 0;

% Solve linear differential equation set and check size
dy = y'*(Q - lam);
dy = dy';

% Check array sizes if inner boolean set
if check_on
    disp(['Dim(y) = ' num2str(size(y))]);
    disp(['Dim(Q) = ' num2str(size(Q))]);
    disp(['Dim(dy) = ' num2str(size(dy'))]);
end