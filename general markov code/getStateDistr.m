% Function to obtain the distribution P(x = a) based on the time spent in
% each state of value x = a <---------------- assumes ergodicity
function prob = getStateDistr(T, x, plotState)

% Ensure vector inputs of correct size
if ~all(size(T) == size(x))
    error('Inconsistent input dimensions');
end

% Remove transients from T vector - makes calculations simpler
T = T - T(1);

% Obtain state space from data and assign storage variables
xmax = max(x);
xmin = min(x);
state = xmin:xmax;
nStates = length(state);
tstate = zeros(1, nStates);

% Loop through possible states and obtain the time spent in each
for i = 1:nStates
    % Obtain indices of time instants of given state and remove case that
    % the last id is the last index of x
    id = find(x == state(i));
    id = id(id < length(T));
    
    % Time spent in state would be T(id+1) - T(id) = dT(id)
    dT = diff(T);
    dT = dT(id);
    tstate(i) = sum(dT); 
end

% Obtain probabilities via two methods and check validity
pr = tstate./max(T);
prcheck = tstate./sum(tstate);
if max(abs(pr - prcheck)) > 10^-9
    assignin('base', 'pr', pr);
    assignin('base', 'prcheck', prcheck);
    error('Probabilities seem inconsistent with check');
end

% Plot results if specified
if plotState
    figure;
    plot(state, pr, 'o-');
    xlabel('state');
    ylabel('probability');
    title('Probability distribution of holding times of states');
    grid;
end

% Assign outputs
prob.val = pr;
prob.state = state;