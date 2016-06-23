% Includes extra input to control the inclusion of a delay

% Modified to work with the Snyder functions that also produce positional
% estimates in addition to state and directional ones and renames posDirSpace
% to posDirSpace and recalculates the pos variable
% Also uses the getDirecfromStates2 function as more sensible

% Modified to work with the 2D MC that is the optimised model for usage
% with the reverse phi generated data (or optimally generated data instead)
% Also altered the function inputs

% Function to perform normal and deterministically delayed Snyder filtering
% and produce MSE for both the state and the optimal directional estimate
% assuming half the states code for 1 direction
function [mseDirecDel msePosDel mseDirec msePos lamInp] = snyderDelNoDel3(X, T, modType, delay, intenInp, nPos, nStates, Q, lamInp, modtype, wantDelay)

% Separate the information and vector observed processes <-----------------
x1 = X(:, 1);
pos = getPosfromStates(x1, nPos, nStates);
phot = X(:, 2:end);

% Extract state inputs and obtain intensity vector based on model type
% together with definition of the state space
posDirSpace = modType.posDirSpace;
intenSpace = modType.intenSpace;
switch(modtype)
    case 1
        % Constant dot intensity
        inten = intenInp.i0;
        Sx = diag(posDirSpace);
    case 2
        % Çontrast changes for 2D MC
        inten = [intenInp.i0 intenInp.ic];
        Sx = diag(0:nStates-1); % <------------------------------------------
end

% Obtain the times of the various photons from different pixels
Tphot = cell(1, nPos);
idphot = cell(1, nPos);
nEvents = 0;
for i = 1:nPos
    [Tphot{i}, perc, idphot{i}] = getEventTimes2(T, phot(:, i), 'birth');
    nEvents = nEvents + length(Tphot{i});
end
assignin('base', 'Tphot', Tphot);

% Obtain the full event set and ensure the times are increasing
Tevent = zeros(nEvents, 1);
offset = 0;
for i = 1:nPos
    if i > 1
        offset = length(idphot{i-1}) + offset;
    end
    range = 1:length(idphot{i});
    range = range + offset;
    Tevent(range) = T(idphot{i});
end
Tevent = sort(Tevent);

% Set the initial conditions for posteriors and states
% q0 = ones(1, nStates)/nStates;
q0 = rand(1, nStates);
q0 = q0/sum(q0);
x1cap0 = sum(q0*Sx, 2);


% Deterministically delay the observed events
TphotDel = cell(1, nPos);
for i = 1:nPos
    TphotDel{i} = Tphot{i} + delay;
end
TeventDel = Tevent + delay;

% Run Snyder on delayed photon data for vector observed processes and
% obtain the appropriate MSE values
calcBools = [1 1];
if wantDelay
    outDel = generalSnyFilterBasicDetDel3(q0, x1cap0, Q, lamInp, nPos, TphotDel, Sx, nEvents, TeventDel, delay, calcBools, posDirSpace, intenSpace, modtype);
else
    outDel = -1;
end
assignin('base', 'outDel', outDel);

% Run normal Snyder on undelayed data
outNorm = generalSnyFilterBasic3(q0, x1cap0, Q, lamInp, nPos, Tphot, Sx, nEvents, Tevent, calcBools, posDirSpace, intenSpace, modtype);
assignin('base', 'outNorm', outNorm);

% Obtain the actual directional estimate in the case that reverse phi is
% not used else obtain the direction expected from the positional data
% which involves the positional motion among nPos positions
plot_on = 1;
disp(' ');
if wantDelay
    disp(['With delay of ' num2str(delay) ' ms']);
    disp('*********************************************************************');
    [mseDirecDel msePosDel] = getMSEStateDirec(T, outDel, x1, pos, posDirSpace, inten, plot_on);
else
    mseDirecDel = -1; 
    msePosDel = -1;
end
disp('Without delay');
disp('*********************************************************************');
[mseDirec msePos] = getMSEStateDirec(T, outNorm, x1, pos, posDirSpace, inten, plot_on);


%% 

% Subfunction to obtain the position from the states
function pos = getPosfromStates(x1, nPos, nStates)

% Get the basic positional states and redundancy
basPos = 0:(nPos-1);
redRatio = nStates/nPos;

% Assign the redundant states to the position that they code
for i = basPos
    redSet = i:nPos:((redRatio-1)*nPos + i);
    x1(ismember(x1, redSet)) = i;
end
pos = x1;

% Check that the positions are correctly assigned
if ~all(ismember(unique(x1), basPos))
    assignin('base', 'x1PosErr', x1);
    error('The position assignment failed');
end

% Subfunction to calculate the state and directional MSE values
function [mseDirec msePos] = getMSEStateDirec(T, outFil, x1, pos, posDirSpace, inten, plot_on)

% Obtain 3 estimates of x1 statistics (with the appended x1 stream applied)
[x1Stats x1n x1capn Tn] = getAllStats(x1, outFil.x1capset, outFil.Tset, T, outFil.Qset, 0, -1);

% Estimate the dot position based on the MC states
[posStats posn poscapn Tn] = getAllStats(pos, outFil.poscapset, outFil.Tset, T, outFil.Qset, 0, -1);
if ~isnan(posStats.meth3(3))
    meanPos = posStats.meth3(1);
    msePos = posStats.meth3(3);
else
    disp('Method 3 failed for MSE calculation');
    meanPos = posStats.meth2(1);
    msePos = posStats.meth2(3);
end

% Obtain direction for the MC positions
[direc, ans] = getDirecfromStates2(x1, posDirSpace, '1D');
[direcReal, ans] = getDirecfromStates2(x1n, posDirSpace, '1D');

% Obtain the best estimate of the direction from the state estimate with
% direct application of the posterior
[dStats dn dcapn Tn] = getAllStats(direc, outFil.dcapset, outFil.Tset, T, outFil.Qset, 0, -1);
if ~isnan(dStats.meth3(3))
    meanDirec = dStats.meth3(1);
    mseDirec = dStats.meth3(3);
else
    disp('Method 3 failed for MSE calculation');
    meanDirec = dStats.meth2(1);
    mseDirec = dStats.meth2(3);
end

% Plot and save the trajectories of the solutions
plotend = ceil(0.25*length(Tn));
ran = 1:plotend;
if plot_on
    % Visualise the estimated and actual position trajectories
    figure;
    plot(Tn(ran), posn(ran), Tn(ran), poscapn(ran));
    legend('reference', 'estimate', 'location', 'best');
    xlabel('time');
    ylabel('position');
    title(['Optimal positional estimate for inten = ' num2str(1000*inten) 's^{-1}']);
    ylim([min(posn)-0.1 max(posn)+0.1]);
    
    % Visualise directional estimates (optimal)
    figure;
    stairs(Tn(ran), dn(ran), 'b');
    hold on
    plot(Tn(ran), dcapn(ran), 'r');
    hold off
    legend('reference', 'estimate', 'location', 'best');
    xlabel('time');
    ylabel('optimal direction');
    title(['Optimal direction estimate for inten = ' num2str(1000*inten) ' s^{-1}']);
    ylim([-1.1 1.1]);
end