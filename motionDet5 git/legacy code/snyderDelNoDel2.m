% Modified to work with the 2D MC that is the optimised model for usage
% with the reverse phi generated data (or optimally generated data instead)
% Also altered the function inputs

% Function to perform normal and deterministically delayed Snyder filtering
% and produce MSE for both the state and the optimal directional estimate
% assuming half the states code for 1 direction
function [mseDirecDel mseCalcDel mseDirec mseCalc lamInp] = snyderDelNoDel2(X, T, modType, delay, intenInp, nPos, nStates, Q, lamInp, modtype)

% Separate the information and vector observed processes <-----------------
pos = X(:, 1);
phot = X(:, 2:end);

% Extract state inputs and obtain intensity vector based on model type
% together with definition of the state space
posSpace = modType.posSpace;
intenSpace = modType.intenSpace;
switch(modtype)
    case 1
        % Constant dot intensity
        inten = intenInp.i0;
        Sx = diag(posSpace);
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
calcDirec = 1;
outFil = generalSnyFilterBasicDetDel2(q0, x1cap0, Q, lamInp, nPos, TphotDel, Sx, nEvents, TeventDel, delay, calcDirec, posSpace, intenSpace, modtype);
assignin('base', 'outDel', outFil);

% Run normal Snyder on undelayed data
outFil = generalSnyFilterBasic2(q0, x1cap0, Q, lamInp, nPos, Tphot, Sx, nEvents, Tevent, calcDirec, posSpace, intenSpace, modtype);
assignin('base', 'outNorm', outFil);

% Obtain the actual directional estimate in the case that reverse phi is
% not used else obtain the direction expected from the positional data
% which involves the positional chain in the assumed reverse phi model
plot_on = 1;
[mseDirecDel mseCalcDel] = getMSEStateDirec(T, outFil, pos, posSpace, inten, plot_on);
[mseDirec mseCalc] = getMSEStateDirec(T, outFil, pos, posSpace, inten, plot_on);

% Subfunction to calculate the state and directional MSE values
function [mseDirec mseCalc] = getMSEStateDirec(T, outFil, pos, posSpace, inten, plot_on)

% Obtain 3 estimates of x1 statistics (with the appended x1 stream applied)
[x1Stats x1n x1capn Tn] = getAllStats(pos, outFil.x1capset, outFil.Tset, T, outFil.Qset, 0, -1);
if ~isnan(x1Stats.meth3(3))
    meanCalc = x1Stats.meth3(1);
    mseCalc = x1Stats.meth3(3);
else
    disp('Method 3 failed for MSE calculation');
    meanCalc = x1Stats.meth2(1);
    mseCalc = x1Stats.meth2(3);
end

% Obtain direction for the MC positions
[direc, ans] = getDirecfromPos2(pos, posSpace, '1D');
[direcReal, ans] = getDirecfromPos2(x1n, posSpace, '1D');

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
if plot_on
%     % Visualise the estimated and actual state trajectories
%     figure;
%     plot(Tn, x1n, Tn, x1capn);
%     xlabel('time');
%     ylabel('state');
%     title(['Optimal state estimate for inten = ' num2str(1000*inten) 's^{-1}']);
    
    % Visualise directional estimates (optimal)
    figure;
    stairs(Tn, dn, 'b');
    hold on
    stairs(Tn, dcapn, 'r');
    hold off
    xlabel('time');
    ylabel('optimal direction');
    title(['Optimal direction estimate for inten = ' num2str(1000*inten) ' s^{-1}']);
end