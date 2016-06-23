% Function that supplies rates based on the Q matrix of a Markov chain
% which represents x1's behaviour
function [rdotx1 bulkx1] = getRatesQMx(Qx1, x1, restrict)

% Essentially a lookup table on the Q matrix of x1 with the row
% corresponding to state x1 representing its reactions
nState = length(Qx1);
state = x1 + 1; % <----------------- x1 starts from 0 = 1st row
Qrow = Qx1(state, 1:nState);

% Obtain rates of birth and death reactions (removing negative diagonal)
rdotx1 = Qrow(Qrow >= 0);

% Obtain transitions of reactions;
bulk = (1:nState) - state;
bulkx1 = setdiff(bulk, 0);

% Depending on the restriction value decide how to treat outputs e.g. a
% rate of zero may be a reaction that simply went zero or a non-existent
% reaction that should not be included
if restrict ~= 0
    % Obtain restricted rates and bulk
    id = find(bulkx1 <= restrict & bulkx1 >= - restrict);
    rdotx1Old = rdotx1;
    bulkx1 = bulkx1(id);
    rdotx1 = rdotx1(id);
    
    % Check that the restriction does not include valid rates
    rdotLeft = setdiff(rdotx1Old, rdotx1);
    if any(rdotLeft > 0)
        assignin('base', 'rdotLeft', rdotLeft);
        error('It appears there are positive transitions outside the restricted set');
    end
end