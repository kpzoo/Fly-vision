% Function that supplies rates based on the Q matrix of a Markov chain
% which represents x1's behaviour
function [rdotx1 bulkx1] = getRatesQMx(Qx1, x1, restrict)

% Essentially a lookup table on the Q matrix of x1 with the row
% corresponding to state x1 representing its reactions
nState = length(Qx1);
state = x1 + 1; % <----------------- x1 starts from 0 = 1st row
if state < 1
    assignin('base', 'x1', x1);
    error('Incorrect state value');
end
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


%%%%%%%%%%%%%%%%%%%%%%%%%%CORRECTIONS TO MAKE%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Need to fix the rdot produced by this function so that it accounts for
% the cases x1 = 0 and nState-1 as at these 1 rate is missing and needs to
% be augmented into rdot as a 0 - e.g. for the +1 -1 case
if x1 == 0
    rdotx1 = [0 rdotx1];
end
if x1 == nState - 1
    rdotx1 = [rdotx1 0];
end
% Account for +1 -1 ordering vs ordering as seen in bulk variable here
rdotx1 = [rdotx1(2) rdotx1(1)];
% - Need to for higher transitions get proper Q matrices and format for the
% even death reactions and odd births