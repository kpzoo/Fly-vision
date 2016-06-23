% Assumes a certain form of SHM motion <--------------------------------
% Modified to include alpha for scaling the intensity

% Function to solve non-linear ODEs on the Snyder filter between observed
% births for a random variable modulating the observed Poisson process -
% note y is a column vector on input
function dy = odeSnyRVNonLin2(ts, y, rspace, w, alpha)

% Boolean to check sizes
check_on = 0;

% Calculate the time dependent lam and lamhat matrices
lam = alpha*diag((rspace + 1)*sin(w*ts) + max(rspace + 1));
lamhat = sum(y'*lam)*eye(length(y));

% Solve linear differential equation set and check size
dy = y'*(lamhat - lam);
dy = dy';

% Check array sizes if inner boolean set
if check_on
    disp(['Dim(y) = ' num2str(size(y))]);
    disp(['Dim(lam) = ' num2str(size(lam))]);
    disp(['Dim(dy) = ' num2str(size(dy'))]);
end