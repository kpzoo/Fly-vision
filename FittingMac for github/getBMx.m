% Note: many of these algebraic operations are only valid under the
% assumption that S is diagonal <-----------------------------------------

% Function to calculate the functional form of lam = lam(x1) matrix which
% is then subtracted from Q to obtain B
function B = getBMx(Q, S, coeff, birStr)

% Check input dimensions, diagonality and calculate suitable identity matrix
if ~all(size(Q) == size(S))
    error('Dimension mismatch between Q and S');
end
if ~isequal(diag(diag(S)), S)
    error('S matrix not diagonal');
end
I = eye(size(S));

% Birth rate options - must match other code <-----------------------------
birSet = {'cross', 'square', 'negexp'};
birType = strmatch(birStr, birSet, 'exact');
if isempty(birType)
    error('Input birth rate does not match any birth type');
end

% Determine functional form from birth type of x2
switch(birType)
    case 1
        % Linear case - lam(x1) = a*x1 + b
        if length(coeff) ~= 2
            error('coeff of incorrect size for linear (cross) function');
        else
            % Calculate the matrix form of lam matrix
            A = coeff(1)*S + coeff(2)*I;
        end
    case 2
        % Square case - lam(x1) = a*x1^2 + b*x1 + c
        if length(coeff) ~= 3
            error('coeff of incorrect size for square function');
        else
            % Calculate the matrix form of lam matrix
            A = coeff(1)*(S*S) + coeff(2)*S + coeff(3)*I;
        end
    case 3
        % Negative exponential case - lam(x1) = ae^-bx1        
        if length(coeff) ~= 2 || coeff(2) <= 0
            error(['coeff of incorrect size or sign for negative exponential'... 
                '(cross) function']);
        else
            % Calculate the matrix form of lam matrix
            A = coeff(1)*expm(-coeff(2)*S);
        end
end

% Obtain B matrix as output
B = Q - A;