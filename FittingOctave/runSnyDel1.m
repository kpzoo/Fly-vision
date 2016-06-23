% Assumptions - only a 2 state symmetric x MC, only a single exponential
% delay, aim to compare normal Snyder, distorted input Snyder and 
% compensated Snyder

% Main code to run a Snyder filter modified to account for delayed photon
% streams via state appendages and QBD theory - code derived from
% runDistort2 from the Snyder Test folders
clear all
clc
close all

% profile clear
% profile on

% Specify location inputs
locfolder1 = 'doubleRange';
locfolder2 = 'test1';
savename = 'linSny1';

% Obtain source files - Gillespie simulation data
cd(locfolder1);
cd(locfolder2);
files = dir('*.mat');
cd ..
cd ..
flen = length(files);

% Set options on distortion structure
lenTh = 1;
noiseMeth = 3*ones(1, lenTh);
paramDistr = 1000;
probDel = 0; % set to 0.34 for QE
probIns = 0;
delayDistr = ones(1, lenTh);
noiseTraits = cell(lenTh, 1);
for i = 1:lenTh
    noiseStruc.meth = noiseMeth(i);
    noiseStruc.paramDistr = paramDistr(:, i);
    noiseStruc.delayDistr = delayDistr(i);
    noiseStruc.probDel = probDel(i);
    noiseStruc.probIns = probIns(i);
    noiseTraits{i} = noiseStruc;
end

% Set control parameters including a cap on the photons allowed
adjust = 1;
limPhoton = 100;

% Declare loop variables
MSEreal1 = zeros(1, flen);
MSEreal2 = zeros(1, flen);
MSEest = zeros(lenTh, flen);
MSEcomp = zeros(1, flen);
lenTreal1 = zeros(1, flen);
lenTreal2 = zeros(1, flen);
lenTest1 = zeros(lenTh, flen);
lenTest2 = zeros(lenTh, flen);
Treal1 = cell(1, flen);
Treal2 = cell(1, flen);
Test1 = cell(1, flen);
Test2 = cell(1, flen);
beta = zeros(1, flen);
gamma = zeros(1, flen);
x1Ests = cell(1, flen);
Qs = cell(1, flen);
Lam = cell(1, flen);
Sy = cell(1, flen);
Sx = cell(1, flen);


% Main loop loads every file in turn and runs the Snyder filter on it with
% the appropriate distortion setting and then runs the compensated Snyder
% for comparison and stores the results
for i = 1:flen
    % Distort appropriate Gillespie data and store data in output - normal
    % Snyder with and without distorted input
    simname = files(i).name;
%     outDist = cleanDistort2(locfolder1, locfolder2, adjust, noiseTraits, limPhoton, simname);
%     
%     % Assign suitable outputs for normal Snyder filters
%     MSEreal1(i) = outDist.mseAct;
%     MSEest(:, i) = outDist.mseEst;
%     Test1{i} = outDist.Test;
%     Treal1{i} = outDist.Treal;
%     lenTreal1(i) = outDist.lenTreal;
%     lenTest1(:, i) = outDist.lenTest;
%     beta(i) = outDist.beta;
%     gamma(i) = outDist.gamma;
    
    % Apply compensated Snyder with the appropriately delayed photon stream
    if ~all(noiseMeth == 3)
        error('Delayed Snyder can only handle exponential delays at the moment');
    else
%        outComp = delSny(locfolder1, locfolder2, adjust, noiseTraits, limPhoton, simname);
       outComp = delSnyNonStat(locfolder1, locfolder2, adjust, noiseTraits, limPhoton, simname);
    end
    
    % Assign further outputs for compensated Snyder filters
    MSEreal2(i) = outComp.mseAct;
    MSEcomp(:, i) = outComp.mseEst;
    Test2{i} = outComp.Test;
    Treal2{i} = outComp.Treal;
    lenTreal2(i) = outComp.lenTreal;
    lenTest2(:, i) = outComp.lenTest;
    x1Ests{i} = outComp.x1StatsEst;
    Qs{i} = outComp.Qs;
    Lam{i} = outComp.Lam;
    Sx{i} = outComp.Sx;
    Sy{i} = outComp.Sy;
    
    % Display progress
    disp('***************************************************************');
    disp(['Completed iteration ' num2str(i) ' on ' simname]);
    disp('***************************************************************');
end

% Save complete data
cd(locfolder1);
cd(locfolder2);
cd('complete');
save(savename);
cd ..
cd ..
cd ..

% Sort data in ascending beta with fixed gamma
if ~all(gamma == gamma(1))
    error('Simulations done over a non-constant gamma');
else
    gam = gamma(1);
    [beta isort] = sort(beta);
    MSEreal1 = MSEreal1(isort);
    MSEreal2 = MSEreal2(isort);
    MSEest = MSEest(:, isort);
    MSEcomp = MSEcomp(:, isort);
    Test1 = Test1{isort};
    Test2 = Test2{isort};
    Treal1 = Treal1{isort};
    Treal2 = Treal2{isort};
    lenTest1 = lenTest1(:, isort);
    lenTest2 = lenTest2(:, isort);
    lenTreal1 = lenTreal1(isort);
    lenTreal2 = lenTreal2(isort);
end

% Save complete and sorted data
cd(locfolder1);
cd(locfolder2);
cd('complete');
save(savename);
cd ..
cd ..
cd ..

% % Plot the distorted, compensated and real MSE
% figure;
% plot(beta, MSEreal1, 'k--');
% hold on
% plot(beta, MSEest, 'b--');
% plot(beta, MSEcomp, 'r--');
% hold off
% xlabel('beta');
% ylabel('mse estimates');
% legend('normal no distortion', 'normal with distortion', 'compensated', 'location', 'best');
% title('Comparison of Snyder estimates with compensation for delayed photons');

% Plot the estimate comparison for delay when only 1 file
if flen == 1 && all(noiseMeth == 3)
    figure;
    plot(TnAct, x1capnAct, Tnn, x1capn);
    xlabel('time');
    ylabel('x1 estimates');
    title(['Comparison of pure Snyder and delayed Snyder when the eta = ' num2str(paramDistr)]);
    legend('pure Snyder', 'delayed Snyder');
end
    
% profile off 
% profile report
