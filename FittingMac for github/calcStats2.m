% Provides statistics: [max mean var min]

% Assumes zoh type behaviour between points

% Function to calculate statistics based on time averaging - applicable to
% cases with and without event times based on number of inputs
function [Amax Amean Avar Amin] = calcStats2(A, T)

% Determine from parameters the types of stats needed
if nargin == 1
    on_time = 0;
else
    % Check dimensions
    if length(A) ~= length(T)
        error('Inputs not of compatible dimensions');
    end
    on_time = 1;
end

% Time invariant statistics
Amax = max(A);
Amin = min(A);

% Calculate statistics without time
if ~on_time
    Amean = mean(A);
    Avar = var(A);
end

% Calculation for the statistics with respect to time
if on_time
    
    % Mean calculation
    dT = diff(T);
    integ_A = sum(dT.*A(1:end-1));
    integ_T = sum(dT);
    Amean = integ_A/integ_T;
    
    % Variance calculation
    A2 = A.*A;
    integ_A2 = sum(dT.*A2(1:end-1));
    Asec_mom = integ_A2/integ_T;
    Avar = Asec_mom - Amean^2;
    
end
