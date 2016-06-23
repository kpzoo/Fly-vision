% Simple function to calculate the lam matrix for a given z input
function lam = getLamZFn(z, dimQ, b)

% Calculation of lam for non-stationary solution
lamDiag = zeros(1, dimQ);
iV = 3;
for iL = 2:dimQ/2
    lamDiag(iV:iV+1) = b*((iL-1) - z);
    iV = iV + 2;
end
lam = diag(lamDiag);