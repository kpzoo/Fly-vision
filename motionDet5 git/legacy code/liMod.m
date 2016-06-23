% Simple code to simulation motion of a dot across pixel positions with the
% aim of maintaining a Markov format
clear all
clc
close all

% Set the spaces of the modulating variables
intenSpace = 1;
posSpace = 0:1;
velSpace = -1:1;

% The position controls the location of the point intensities which then
% leads to the overall intensity I. This is indexed based on the no.
% velocity changes used
Nev = 500000;
I = zeros(Nev, length(intenSpace), length(posSpace));

% Simulate a symmetrical velocity markov chain via Gillespie algorithm with
% the velocity having a fixed state space
inpGill.N = Nev;
inpGill.Nstart = 1;
inpGill.len = 1;
inpGill.x0 = velSpace(1);
inpGill.stateMin = min(velSpace);
inpGill.stateMax = max(velSpace);
inpGill.r_const = [0.1 0.1];
inpGill.transit = [1 -1];
[vel veldot T] = simpleGill(inpGill);

% Obtain statistics of markov velocity
dT = diff(T);
velMean = sum(dT.*vel(1:end-1))/sum(dT);
vel2 = vel.*vel;
velVar = sum(dT.*vel2(1:end-1))/sum(dT);

% Obtain the positions via integration of the velocities across the
% continuous components between transitions
pos1 = cumtrapz(T, vel);
pos2 = zeros(size(T));
pos2(2:end) = dT.*vel(1:end-1);

% Obtain a position estimate accounting for the zero velocity terms
pos3 = zeros(size(T));
for i = 1:length(dT)
    if vel(i) ~= 0
        pos3(i+1) = dT(i)*vel(i);
    else
        pos3(i+1) = pos3(i);
    end
%     pos3(i+1) = dT(i)*vel(i) + pos3(i);
end

% Convert the raw position into discrete 0 and 1 values
% pos3Raw = pos3;
% pos3 = sign(pos3);
% pos3(pos3 == -1) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%% ISSUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It is likely that the problem results from failing to sample the velocity
% wave and then using these interpolated points to get position - if this
% is done a sensible positional estimate should be obtained

% Secondly by using a little E(ds/dt) = d(E(s))/dt and similarly for
% integration one can see that a stationary velocity process gives rise to
% a stationary positional process