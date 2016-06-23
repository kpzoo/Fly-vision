% Simplified function to perform ODE solution in the simplest case of
% kbirth = kdeath and S = [0 1]
function dy = odeSimpleState(ts, y, k, a)


% Choose option based on input length
leny = length(y);

switch(leny)
    
    case 2
        % 2 state space ODEs with substitution
        dy(1) = k + y(1)*(a - 2*k) - a*(y(1)^2);
        dy(2) = -dy(1);
        
    case 3
        % 3 state space ODEs with substitution
        dy(1) = k*(y(2) - 2*y(1)) + a*y(1)*(2 - 2*y(1) - y(2));
        dy(2) = 2*k*(1 - 2*y(2)) + a*y(2)*(1 - 2*y(1) - y(2));
        dy(3) = -(dy(1) + dy(2));
        
        
    otherwise
        error(['No code to support state space of size = ' num2str(leny)]);
end

% Ensure column vector output
dy = dy';