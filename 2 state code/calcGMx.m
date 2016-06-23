% Note the correct calculation of G from A which is required for the
% non-linear cases is not developed yet <----------------------------------

% Function to calculate the G matrix which represents the perturbation
% caused by the birth events (jumps) of x2
function q = calcGMx(birStr, coeff, S, q, lamcap, iter)

% Birth rate options - must match other code <-----------------------------
birSet = {'cross', 'square', 'negexp'};
birType = strmatch(birStr, birSet, 'exact');
if isempty(birType)
    error('Input birth rate does not match any birth type');
end

% Code assumes diagonal S matrix
if ~isequal(diag(diag(S)), S)
    error('S matrix not diagonal');
end
I = eye(size(S));

% Determine G from birth type of x2 (functional form)
switch(birType)
    case 1
        % Linear case - lam(x1) = a*x1 + b
        if length(coeff) ~= 2
            error('coeff of incorrect size for linear (cross) function');
        else
            % Calculate the matrix form of lam matrix
            A = coeff(1)*S + coeff(2)*I;
            G = q*A;
        end
    case 2
        % Square case - lam(x1) = a*x1^2 + b*x1 + c
        if length(coeff) ~= 3
            error('coeff of incorrect size for square function');
        else
            % Calculate the matrix form of lam matrix
            A = coeff(1)*(S*S) + coeff(2)*S + coeff(3)*I;
            G = 'not developed yet';
        end
    case 3
        % Negative exponential case - lam(x1) = ae^-bx1
        if length(coeff) ~= 2 || coeff(2) <= 0
            error(['coeff of incorrect size or sign for negative exponential'...
                '(cross) function']);
        else
            % Calculate the matrix form of lam matrix
            A = coeff(1)*expm(-coeff(2)*S);
            G = 'not developed yet';
        end
end

% Obtain output from G and lamcap
q = G/lamcap;

% Check that q maintains the requirements of a probability distribution
if any(q < -10^-9)
    assignin('base', 'qerror', q);
    assignin('base', 'lamEst', lamcap);
    assignin('base', 'G', G);
    error(['q distribution has negative entries at i =' num2str(iter)]);
end
if max(abs(sum(q) - 1)) > 10^-9
    assignin('base', 'qerror', q);
    assignin('base', 'lamEst', lamcap);
    assignin('base', 'G', G);
    error(['q distribution does not sum to 1 at i = ' num2str(iter)]);
end
