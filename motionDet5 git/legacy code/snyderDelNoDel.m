% Function to perform normal and deterministically delayed Snyder filtering
% and produce MSE for both the state and the optimal directional estimate
% assuming half the states code for 1 direction
function [mseDirecDel mseCalcDel mseDirec mseCalc lamInp] = snyderDelNoDel(X, T, posSpace, delay, inten, nPos, Q)

% Separate the information and vector observed processes <-----------------
pos = X(:, 1);
phot = X(:, 2:end);

% Set various rate matrices based on assumed form of state MC <------------
lamDiagIndiv = zeros(nPos, length(posSpace));
for i = 1:nPos
    % Assume the coding of position and direction via a state then a
    % counting process represents 2 states e.g. 0 and 3 for space 0:5
    lamDiagIndiv(i, [i nPos+i]) = inten;
end
lamDiag = sum(lamDiagIndiv, 1);

% Set the input structures to the filtering function
lamInp.lamDiagIndiv = lamDiagIndiv;
lamInp.lamDiag = lamDiag;

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

% Obtain position state space matrix and set initial conditions
Sx = diag(posSpace);
% q0 = ones(1, length(posSpace))/length(posSpace);
q0 = rand(1, length(posSpace));
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
plot_on = 0;
outFil = generalSnyFilterBasicDetDel(q0, x1cap0, Q, lamInp, nPos, TphotDel, Sx, nEvents, TeventDel, delay, calcDirec, posSpace);
[mseDirecDel mseCalcDel] = getMSEStateDirec(T, outFil, pos, posSpace, inten, plot_on);
assignin('base', 'outDel', outFil);

% Run normal Snyder on undelayed data
outFil = generalSnyFilterBasic(q0, x1cap0, Q, lamInp, nPos, Tphot, Sx, nEvents, Tevent, calcDirec, posSpace);
[mseDirec mseCalc] = getMSEStateDirec(T, outFil, pos, posSpace, inten, plot_on);
assignin('base', 'outNorm', outFil);




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
[direc, ans] = getDirecfromPos(pos, posSpace);
[direcReal, ans] = getDirecfromPos(x1n, posSpace);

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
    % Visualise the estimated and actual state trajectories
    figure;
    plot(Tn, x1n, Tn, x1capn);
    xlabel('time');
    ylabel('state');
    title(['Optimal state estimate for inten = ' num2str(1000*inten) 's^{-1}']);
    
    % Visualise directional estimates (optimal)
    figure;
    stairs(Tn, dn, 'b');
    hold on
    stairs(Tn, dcapn, 'r');
    hold off
    xlabel('time');
    ylabel('optimal direction');
    title(['Optimal direction estimate for inten = ' num2str(1000*inten) 's^{-1}']);
end