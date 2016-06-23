% Main script to run strip and dot reverse phi experiments with aims of
% matching Tuthill and other stimuli
clear all
clc
close all

% Set main booleans to control code
stimType = 1;       % controls dots or strips
fullWave = 0;       % controls full wave rectification
wantDelay = 0;      % controls running of delayted Snyder
plotPos = 0;        % plots the positional intensities and photons 
plotChess = 0;      % plots the chess board type mode
testCon = 0;        % removes some connections on the MC

% Set the spaces of the internal model and calculate states
posDirSpace = 0:5;
intenSpace = 0:1;
nIntens = length(intenSpace);
nPos = length(posDirSpace)/2;
nStates = length(posDirSpace)*length(intenSpace);
disp(['Internal model has ' num2str(nStates) ' states and ' num2str(nPos) ' positions']);
disp(['Internal model expects ' num2str(nIntens) ' intensities']);

% Set the no. events and the main internal rates (ms and ms^-1)
a = 0.001;
b = a/20;
c = a;
h = a;
Nev = 20000;

% Set the intensity values based on extremes with a linear gradient and
% convert to contrasts based on median intensity
imin = 0;
imax = 20*a;
iset = linspace(imin, imax, nIntens);
ibgnd = median(iset);  
cset = iset - ibgnd;

% Fullwave rectify the contrasts if required
if fullWave
    cset = abs(cset);
    iset = cset + ibgnd;
    disp('Contrasts have been full wave rectified');
end

% Collect model and MC state parameters
intMod.a = a;
intMod.b = b;
intMod.c = c;
intMod.h = h;
intMod.posDirSpace = posDirSpace;
intMod.intenSpace = intenSpace;
intMod.nStates = nStates;
intMod.nPos = nPos;
intMod.nIntens = nIntens;
intMod.model = 2;

% Collect intensity inputs
intenInp.imin = imin;
intenInp.imax = imax;
intenInp.iset = iset;
intenInp.ibgnd = ibgnd;
intenInp.cset = cset;

% Obtain the Q matrix for the internal model with different connections
testDir = [0 0];
[Q Pi] = getQLightExp(intMod, testDir, nIntens);

% Obtain the positional rate matrices (lam) based on if using bars or a dot
lamInp = getLamSetExp(stimType, intMod, intenInp);
