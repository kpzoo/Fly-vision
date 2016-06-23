% Script to generate Gillespie data and then filter it according both to a
% normal and deterministically delayed Snyder formulation
clear all
clc
close all

% Set the space for possible position and intensity values (note that
% the position space is twice the actual no. positions as the duplication
% of states is used to represent different directions
posSpace = 0:5;
intenSpace = 1;
Nev = 20000;
nPos = length(posSpace)/2;

% Set the intensity set of the dot implicitly and the explicit delay (ms)
% intenPm = 10*[2:2:50];
intenPm = 10;
delay = 43.3;
nRuns = length(intenPm);

% Set the transition rates of the infinitesimal generator
a = 0.01;
b = a/5;
c = a;
inten = 0.001*intenPm; % removed dependence on a (intensity in ms^-1)

% Output variable declaration
mseDirec = zeros(size(intenPm));
mseDirecDel = zeros(size(intenPm));
mseCalc = zeros(size(intenPm));
mseCalcDel = zeros(size(intenPm));
x1Stats = cell(size(intenPm));
x1StatsDel = cell(size(intenPm));
lamInp = cell(size(intenPm));

% Set simulation file and folder names
if nRuns > 1
    folder = ['case_2_' num2str(1/a)];
    figname = ['compDel_2_' num2str(1/a)];
    savename = ['overall_2_' num2str(1/a)];
    cd('data');
    mkdir(folder);
    cd ..
end

% Loop across the main functions for producing the Gillespie simulations
% and then running the Snyder filters
for i = 1:nRuns
    % Simulation of light model
    [X T posStats Q] = simLightPosDir(a, b, c, posSpace, inten(i), Nev, nPos);
    
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

% Plot and save the overall results
if nRuns > 1
    figure;
    plot(1000*inten, mseDirec, 1000*inten, mseDirecDel);
    legend('no delay', ['delay = ' num2str(delay)], 'location', 'best');
    xlabel('intensity (s^{-1})');
    ylabel('mse for directional estimates');
    title(['MSE for normal and determinsitically delayed observations, a = ' num2str(a)]);
    saveas(gcf, ['compDel_2_' num2str(1/a) '.fig']);
    figure;
    plot(1000*inten, mseCalc, 1000*inten, mseCalcDel);
    legend('no delay', ['delay = ' num2str(delay)], 'location', 'best');
    xlabel('intensity (s^{-1})');
    ylabel('mse for state estimates');
    title(['State MSE for normal and determinsitically delayed observations, a = ' num2str(a)]);
    saveas(gcf, ['compStateDel_2_' num2str(1/a) '.fig']);
    save(savename);
end