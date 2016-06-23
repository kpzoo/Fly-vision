% Modified to account for full wave rectification and interleaving of
% normal and reversed phi as well as handle matching of models

% Modified to produce estimates of position, direction and states
% separately and renamed the old posSpace to posDirSpace

% Script to simulate the reverse phi problem in which the Snyder filter has
% a Q derived from a 2D MC but is fed with photon data derived from an
% illusory stimuli which keeps the Markov position but has a flashing dot
clear all
clc
close all

% Variables to control the execution of certain cells
cell1 = 1;
cell2 = 1;
savedata = 0;
testCon = 0;
fullWave = 0;
interleave = 0;
matching = 1;

%% Cell for simulating the reverse phi and filtering according to the model
if cell1
    % Set the initial simulation parameters which control the size of the MC
    % that is internal to the Snyder model
    posDirSpace = 0:5;
    intenSpace = 0:1;
    Nev = 10000;
    nPos = length(posDirSpace)/2;
    nStates = length(posDirSpace)*length(intenSpace);
    
    % Set the rates associated with the internal 2D MC with particular
    % attention to the relative speed of intensity and position changes
    a = 0.001;
    b = a/20;
    if ~matching
        h = 0.1*a;
    else
        % No flashing allowed so remove extra transitions
        h = 0;
    end
    
    % Set the delay and possible intensity values (ms and ms^-1)
    delay = 1*43.3;
    i0 = 10*0.001;
    if matching
        % As only white can have ic > i0 for better estimates
        ic = 1*i0;
    else
        % For flashing ic <= i0
        ic = i0;
    end
    
    % Set inputs to optimised filter light model
    modVal = 2;
    modType.a = a;
    modType.b = b;
    modType.c = a;
    modType.h = h;
    modType.posDirSpace = posDirSpace;
    modType.intenSpace = intenSpace;
    
    % Obtain the Q matrices for the Snyder filter
    if testCon
        testDir = [0 1];
        [Q Pi] = getLightModelQ3(modVal, modType, testDir);
    else
        [Q Pi] = getLightModelQ2(modVal, modType);
        testDir = [0 0];
    end
    
    % Obtain the Lam matrices based on intensities and matching boolean
    intenInp.i0 = i0;
    intenInp.ic = ic;
    lamInp = getLamManyStreams2(modVal, nPos, nStates, intenInp, matching);
    
    % Get the false data from the reverse phi stimulus which only matches the
    % real model in terms of the positional-directional parts
    datVal = 1;
    nStatesRevPhi = length(posDirSpace);
    test.testCon = testCon;
    test.testDir = testDir;
    [X T posStats Qphi] = simLightPosDir4(datVal, modVal, modType, intenInp, Nev, nPos, nStatesRevPhi, test, fullWave, matching);
    
    % Obtain the MSE via Snyder filtering
    [mseDirecDel msePosDel mseDirec msePos lamInp] = snyderDelNoDel3(X, T, modType, delay, intenInp, nPos, nStates, Q, lamInp, modVal);
end

% Save a subset of the simulation results with a unique name in a subfolder
% based on the position
if savedata
    clear outDel outNorm
    folder = ['revPos' num2str(nPos)];
    cd data
    try
        cd(folder);
        files = dir('*.mat');
        flen = length(files);
        id = flen + 1;
    catch excepFol
        mkdir(folder);
        cd(folder);
        id = 0;
    end
    save(['testPos' num2str(id)]);
    cd ..
    cd ..
else
    clear outDel outNorm
    save('revPhiMaxN');
end

% Quick check on model construction
if matching
    if all(Q == [Qphi zeros(length(Qphi)); zeros(length(Qphi)) Qphi])
        disp('The internal and external MC generators have a matching form');
    else
        disp('The infinitesimal generators do not match');
    end
end


%% Cell for plotting results
if cell2
    % View the intensity from the workspace Rset and T for a fixed interval as
    % well as the photons produced in that period
    if exist('Rset', 'var') && exist('T', 'var') && exist('pos', 'var') && exist('Tphot', 'var')
        % Obtain the light intensities, positions and intensity range
        Li = Rset(2:end, :);
        posit = pos(2:end);
        yrange = [min(min(Li)) max(max(Li))];
        if length(unique(yrange)) == 1
            noYLim = 1;
        else
            noYLim = 0;
        end
        
        % Plot across a fixed interval and intensity range
        plotend = 200;
        hFig = figure;
        set(hFig, 'Position', [100 100 800 800])
        hold on
        for j = 1:nPos
            % Plot the intensities as a column of subfigures
            subplot(nPos+1, 1, j);
            stairs(T(1:plotend), Li(1:plotend, j), 'linewidth', 3);
            if ~noYLim
                ylim(yrange);
            end
            ylabel(['I @ pos = ' num2str(j-1)]);
        end
        subplot(nPos+1, 1, nPos+1);
        h = stairs(T(1:plotend), posit(1:plotend), 'linewidth', 3);
        set(h, 'color', 'r');
        ylabel('pos switch');
        hold off
        xlabel('time');
        % Set a label in the centre
        set(gcf,'NextPlot','add');
        axes;
        h = title('Light intensities at different positions');
        set(gca,'Visible','off');
        set(h,'Visible','on');
        
        % Obtain the photon streams over the same interval as the above plots
        Tend = T(plotend);
        Tstart = Tend;
        for j = 1:nPos
            Tstart = min(Tstart, min(Tphot{j}));
        end
        % Plot across a fixed time range
        xrange = [Tstart Tend];
        hFig = figure;
        set(hFig, 'Position', [100 100 800 800])
        hold on
        for j = 1:nPos
            % Plot the streams as delta trains on a column
            subplot(nPos+1, 1, j);
            stem(Tphot{j}, ones(size(Tphot{j})), 'linewidth', 3);
            xlim(xrange);
            ylabel(['P(t) @ pos = ' num2str(j-1)]);
        end
        subplot(nPos+1, 1, nPos+1);
        h = stairs(T(1:plotend), posit(1:plotend), 'linewidth', 3);
        set(h, 'color', 'r');
        ylabel('pos switch');
        hold off
        xlabel('time');
        % Set a label in the centre
        set(gcf,'NextPlot','add');
        axes;
        h = title('Photon trains at different positions');
        set(gca,'Visible','off');
        set(h,'Visible','on');
    end
end