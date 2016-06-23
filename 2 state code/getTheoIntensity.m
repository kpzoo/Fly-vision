% Function to get theoretical intensity statistics as opposed to using the
% actual simulations of lam - this code assumes that the intensity relation
% is linear: intensity = kgain*x1 <----------------------------------------
function [inmean invar] = getTheoIntensity(Pi, state, kgain)

% Obtain simple MC statistics from discrtet state distribution
xmean = sum(Pi.*state);
xms = sum(Pi.*(state.^2));
xvar = xms - xmean^2;

% Obtain intensity statistics from linear relation
inmean = kgain*xmean;
invar = (kgain^2)*xvar;