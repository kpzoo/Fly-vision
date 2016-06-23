% Diagnostic function to calculate the time between all types of events in
% order to evaluate relative speed of x1 and x2 events
function eventStats = getRelMeanEventTimes(T, x1, x2)

% Obtain birth event times of both species
type = 'birth';
[Tbir1 percbir1] = getEventTimes(T, x1, type);
[Tbir2 percbir2] = getEventTimes(T, x2, type);

% Obtain death event times of x1
type = 'death';
[Tdea1 percdea1] = getEventTimes(T, x1, type);
[Tdea2 percdea2] = getEventTimes(T, x2, type);

% Obtain both event times of both species
type = 'both';
[Tboth1 percboth1] = getEventTimes(T, x1, type);
[Tboth2 percboth2] = getEventTimes(T, x2, type);

% Calculate time differences between event types for x1 and x2
dTbir1 = diff(Tbir1);
dTbir2 = diff(Tbir2);
dTdea1 = diff(Tdea1);
dTdea2 = diff(Tdea2);
dTboth1 = diff(Tboth1);
dTboth2 = diff(Tboth2);

% Obtain the mean time difference between events
mTbir1 = mean(dTbir1);
mTbir2 = mean(dTbir2);
mTdea1 = mean(dTdea1);
mTdea2 = mean(dTdea2);
mTboth1 = mean(dTboth1);
mTboth2 = mean(dTboth2);

% Obtain speed ratios for x2/x1 total events, x1bir/x1dea and x2bir/x1bir -
% as speed comparisons the mean time ratio is inverted
rx2x1 = (mTboth2/mTboth1)^(-1);
rx1bx1d = (mTbir1/mTdea1)^(-1);
rx2bx1b = (mTbir2/mTbir1)^(-1);


% Assign outputs
eventStats.rx2x1 = rx2x1;
eventStats.rx1bx1d = rx1bx1d;
eventStats.rx2bx1b = rx2bx1b;
eventStats.percbir = [percbir1 percbir2];
eventStats.percdea = [percdea1 percdea2];
eventStats.percboth = [percboth1 percboth2];
eventStats.meanbir = [mTbir1 mTbir2];
eventStats.meandea = [mTdea1 mTdea2];
eventStats.meanboth = [mTboth1 mTboth2];
