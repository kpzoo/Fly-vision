% Function to calculate the lam estimate based on the functional form of
% lam(x1) which is specified via the x2 birth function
function [lamEst x1Est] = calclamEst(birStr, coeff, q, S)

% Birth rate options - must match other code <-----------------------------
birSet = {'cross', 'square', 'negexp'};
birType = strmatch(birStr, birSet, 'exact');
if isempty(birType)
    error('Input birth rate does not match any birth type');
end

% Calculate the estimate of x1
x1Est = sum(q*S);

% Using functional form of birth type of x2 calculate lamEst
switch(birType)
    case 1
        % Linear case - lam(x1) = a*x1 + b
        if length(coeff) ~= 2
            error('coeff of incorrect size for linear (cross) function');
        else            
            lamEst = coeff(1)*x1Est + coeff(2);
        end
    case 2
        % Square case - lam(x1) = a*x1^2 + b*x1 + c
        if length(coeff) ~= 3
            error('coeff of incorrect size for square function');
        else
            lamEst = coeff(1)*(x1Est^2) + coeff(2)*x1Est + coeff(3);
        end
    case 3
        % Negative exponential case - lam(x1) = ae^-bx1        
        if length(coeff) ~= 2 || coeff(2) <= 0
            error(['coeff of incorrect size or sign for negative exponential'... 
                '(cross) function']);
        else            
            lamEst = coeff(1)*exp(-coeff(2)*x1Est);
        end
end