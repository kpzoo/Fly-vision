% Function to collect the 3 estimates of x1 statistics
function [x1Stats x1n x1capn Tn] = getAllStats(x1, x1Set, Tset, T, Qset, rem_trans, rawbnd)

% First method
[x1n x1capn Tn x1stats1] = getIntensityStatsFull2(x1, x1Set, Tset, T, Qset, rem_trans);

% Second method
ex1n = x1n - x1capn;
ex1Sqn = ex1n.*ex1n;
x1stats2.meanErr = trapz(Tn, ex1n)/(max(Tn) - min(Tn));
x1stats2.mseErr = trapz(Tn, ex1Sqn)/(max(Tn) - min(Tn));
x1stats2.varErr = x1stats2.mseErr - x1stats2.meanErr^2;  

% Third method (modified to use faster getStatsAltMeths2)
x1stats3 = getStatsAltMeths2(Tn, x1n, x1capn, x1, T);

% Collect all data into a single structure
x1Stats.order = x1stats3.order;
x1Stats.meth1 = [x1stats1.meanErr x1stats1.varErr x1stats1.mseErr];
x1Stats.meth2 = [x1stats2.meanErr x1stats2.varErr x1stats2.mseErr];
x1Stats.meth3 = x1stats3.interpMeth;

% Add another field if the bound is provided - input -1 if no bound
if rawbnd ~= -1
    x1Stats.psi = [x1Stats.meth1(3) x1Stats.meth2(3) x1Stats.meth3(3)]/rawbnd;
end