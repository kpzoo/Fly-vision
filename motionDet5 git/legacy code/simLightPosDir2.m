% Modified to account for more general light models with more generalised
% inputs that allow constuction of more complex schemes. Further uses
% updated function gillespieQGeneral which accounts for the intensity of
% the observed processes in a more general manner

% Function to perform simulation of the Markov light model in which each
% state codes for both position and direction
function [X T posStats Q] = simLightPosDir2(type, inpType, inten, Nev, nPos)

% Obtain the Q matrix for a M state (M/2 position) stimulus and the state
% distribution at equilibirum
posSpace = inpType.posSpace;
M = 2*nPos;
[Q Pi] = getLightModelQ(type, inpType);

% Obtain the number of reactions from the max jump in state found in Q -
% assumes a state structure of [0, 1, ... , m]
dimQ = length(Q);
jumpMax = 0;
for i = 1:dimQ
    % Check each row for highest difference between position of non-zero
    % rates and the current state
    Qrow = Q(i, 1:dimQ);
    idjump = find(Qrow > 0);
    currjumpMax = max(abs(idjump - i));
    jumpMax = max(currjumpMax, jumpMax);
end

% Assume 3 (based on pixels) poisson process modulated by reactions in Q and
% obtain the transit matrix based on the jumps assuming 1:jumpMax exists
nReacs = 2*jumpMax + nPos;
bulk = [1:jumpMax 1:jumpMax];
bulk = sort(bulk);
transitx1 = [bulk zeros(1, nPos)];
transitx1(2:2:end) = -transitx1(2:2:end);

% Calculate the further transit rows - one for each pixel
transitxn = [zeros(nPos, nReacs - nPos) eye(nPos)];
transit = [transitx1; transitxn];

% Set the input parameters to the Gillespie code noting that the number of
% molecules is x1 (information process) plus nPos observed processes
inpGill.N = Nev;
inpGill.Nstart = 1;
inpGill.len = 1 + nPos;
inpGill.Q = Q;
inpGill.transit = transit;
inpGill.nReacs = nReacs;
inpGill.alpha = inten;
inpGill.x0 = zeros(1, inpGill.len);
inpGill.nPos = nPos;
inpGill.type = type;

% Run the Gillespie algorithm for the Q matrix specified
[X, Xdot, T] = gillespieQGeneral(inpGill);
pos = X(:, 1);

% Obtain statistics of markov position and compare to theoretical
dT = diff(T);
posStats.posMean = sum(dT.*pos(1:end-1))/sum(dT);
pos2 = pos.*pos;
posStats.posVar = sum(dT.*pos2(1:end-1))/sum(dT);
posStats.posMeanTheo = sum(posSpace*Pi);
posSpace2 = posSpace.*posSpace;
posStats.posVarTheo = sum(posSpace2*Pi);