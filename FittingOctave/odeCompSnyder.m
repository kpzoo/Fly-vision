% Function to solve ODEs on the Snyder filter between observed birth events
% with compensation for delayed observations (exponential) - note y is a
% column vector on input
function dy = odeCompSnyder(ts, y, Q, lam, dimQ)

% Boolean to check sizes
check_on = 0;

% Obtain lamhat estimate
lamcap = y'*lam*ones(dimQ, 1);
lamhat = lamcap*eye(dimQ);

% Solve differential equation set and check size
dy = y'*(Q + lamhat - lam);
dy = dy';

% Check array sizes if inner boolean set
if check_on
    disp(['Dim(y) = ' num2str(size(y))]);
    disp(['Dim(C) = ' num2str(size(C))]);
    disp(['Dim(dy) = ' num2str(size(dy'))]);
end


