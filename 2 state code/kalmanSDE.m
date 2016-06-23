% Code assumes kb = kd = k <---------------------------------------------
% Some kalman filtering code to compare to the Snyder filter
clear all
clc
close all

% Set parameters
k = 10;
alpha = 10;

% Define set of SDEs which become the system and observation equations
