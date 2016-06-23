% Modified to account for more complex light models and to include a plot
% based on the equivalent percentage failures

% Script to generate Gillespie data and then filter it according both to a
% normal and deterministically delayed Snyder formulation
clear all
clc
close all

% Set the space for possible position and intensity values (note that
% the position space is twice the actual no. positions as the duplication
% of states is used to represent different directions
posSpace = 0:3;
intenSpace = 1;
Nev = 20000;
nPos = length(posSpace)/2;

% Set the intensity set of the dot implicitly and the explicit delay (ms)
% intenPm = 1:1:50;
intenPm = 10;
delay = 43.3;
nRuns = length(intenPm);

% Set the transition rates of the infinitesimal generator
a = 0.001;
b = a/5;
c = a;

% Set the intensity based on whether want fixed intensity or fixed alpha
fix_alpha = 1;
if fix_alpha
    inten = a*intenPm;
else
    inten = 0.001*intenPm; % removed dependence on a (intensity in ms^-1)
end

% Output variable declaration
mseDirec = zeros(size(intenPm));
mseDirecDel = zeros(size(intenPm));
mseCalc = zeros(size(intenPm));
mseCalcDel = zeros(size(intenPm));
lamInp = cell(size(intenPm));

% Set inputs to light model
type = 1;
inpType.a = a;
inpType.b = b;
inpType.c = c;
inpType.posSpace = posSpace;

% Set identifiers for data storage
if a < 1
    idnum = num2str(1/a);
    id = 'inv';
else
    idnum = num2str(a);
    id = 'norm';
end

% Set simulation file and folder names
if nRuns > 1
    folder = ['case' id idnum];
    figname = ['compDel' id idnum];
    savename = ['overall' id idnum];
    cd('data');
    mkdir(folder);
    cd ..
end

% Loop across the main functions for producing the Gillespie simulations
% and then running the Snyder filters
for i = 1:nRuns
    % Simulation of light model
    [X T posStats Q] = simLightPosDir2(type, inpType, inten(i), Nev, nPos);
    
    % Obtain the MSE via Snyder filtering
    [mseDirecDel(i) mseCalcDel(i) mseDirec(i) mseCalc(i) lamInp{i}] = snyderDelNoDel(X, T, posSpace, delay, inten(i), nPos, Q);
    
    % Save simulation data and display progress
    if nRuns == 1
        simname = 'simTest';
        save(simname);
    else
        cd('data');
        cd(folder);
        simname = ['sim' num2str(i)];
        save(simname);
        cd ..
        cd ..
    end
    disp(['Completed iteration ' num2str(i) ' of ' num2str(nRuns)])
end

% Calculate the percent of falures if one assumes that for time tp the
% error is max at 2 and 0 otherwise
percFail = 25*mseDirec;
percFailDel = 25*mseDirecDel;

% Plot and save the overall results
if nRuns > 1
    figure;
    plot(1000*inten, mseDirec, 1000*inten, mseDirecDel);
    legend('no delay', ['delay = ' num2str(delay)], 'location', 'best');
    xlabel('intensity (s^{-1})');
    ylabel('mse for directional estimates');
    title(['MSE for normal and determinsitically delayed observations, a = ' num2str(a)]);
    saveas(gcf, ['compDel' id idnum '.fig']);
    figure;
    plot(1000*inten, percFail, 1000*inten, percFailDel);
    legend('no delay', ['delay = ' num2str(delay)], 'location', 'best');
    xlabel('intensity (s^{-1})');
    ylabel('failure percent for directional estimates');
    title(['Failure percent for normal and determinsitically delayed observations, a = ' num2str(a)]);
    saveas(gcf, ['percDel' id idnum '.fig']);
    figure;
    plot(1000*inten, mseCalc, 1000*inten, mseCalcDel);
    legend('no delay', ['delay = ' num2str(delay)], 'location', 'best');
    xlabel('intensity (s^{-1})');
    ylabel('mse for state estimates');
    title(['State MSE for normal and determinsitically delayed observations, a = ' num2str(a)]);
    saveas(gcf, ['compStateDel' id idnum '.fig']);
    save(savename);
end