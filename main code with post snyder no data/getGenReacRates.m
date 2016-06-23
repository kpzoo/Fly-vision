% Function to obtain all the desired reaction rates in correspondence with
% the r_const protocol of [bir dea bir dea...]
function rdot = getGenReacRates(xhist, rprev, x, r_const, nReacs, reacType,...
    molecType, crossType, SlimSet, bulk)

% Basic functional forms of reaction rates
%   Case 0 - constant (includes zero rate)
%   Case 1 - linear self with or without restriced space
%   Case 2 - linear cross with or without restricted space


% Extract inputs
Smax = SlimSet.max;
Smin = SlimSet.min;

% Initialise output and separate even and odd r_const values
rdot = zeros(1, nReacs);

% Ensure the correct number of molecular identities are provided
if length(molecType) ~= nReacs
    assignin('base', 'molecType', molecType);
    error('The molecular indicator is of incorrect length');
end

% Calculate the reaction rates with molecType determining the
% dependent variables for the calculation
for j = 1:nReacs

    % Obtain molecular ids for the reaction to be evaluated
    mID = molecType(j);
    cID = crossType(j);
    
    % Assign rates based on type and functional definition
    switch(reacType(j))
        case 0
            % Constant constitutive rate
            rdot(j) = r_const(j);

        case 1
            % Linear in the molecule given by molecType (the reactant)
            if rem(j, 2) == 0
                % Death reaction
                rdot(j) = r_const(j)*(x(mID) - Smin(mID));
            else
                % Birth reaction - two versions based on Smax
                if ~isinf(Smax(mID))
                    rdot(j) = r_const(j)*(Smax(mID) - x(mID));
                else                    
                    rdot(j) = r_const(j)*x(mID);
                end
            end

        case 2
            % Linear in a non-reactant molecule - specified in crossType
            if rem(j, 2) == 0
                % Death reaction
                rdot(j) = r_const(j)*(x(cID) - Smin(mID));
            else
                % Birth reaction - assumes cross molecule well behaved so
                % no state restrictions are made here
                rdot(j) = r_const(j)*x(cID);
            end

        otherwise
            error('Unsupported birth type specified');
    end

    % Ensure that no positive rates are provided for reactions that could
    % result in defiance of the Smax and Smin boundaries
    if rem(j, 2) == 0
        % Death reaction correction
        if x(mID) - bulk(j) < Smin(mID)
            rdot(j) = 0;
        end
    else
        % Birth reaction correction
        if x(mID) + bulk(j) > Smax(mID)
            rdot(j) = 0;
        end
    end
end



% Check for negative rates or all rates at 0
if any(rdot < 0) || all(rdot == 0)
    assignin('base', 'rdot', rdot);
    error('The calculated rates are zero or negative');
end