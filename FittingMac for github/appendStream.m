% NOTE - need to account for repeated values between T and TphEst  as want 
% to remove repeated times yet cannot exactly due to multiple reactions <------------------

% This function does not impose any restrictions on the dimensionality of X
% and thus does not require X = [x1 x2] only but can include additional
% molecules to allow multi-stage delays

% Function that creates a coherent stream that accounts for the new z
% events relayed as delayed photon times
function [Tnew Xnew] = appendStream(T, X, TphEst)

% Obtain observed counting process z
TphEst = sort(TphEst);
z = 1:length(TphEst);
z = z';

% Append the delayed photon times into the time stream and sort
Tnew = [T; TphEst];
[Tnew idsort] = sort(Tnew);
applen = length(TphEst);

% Apply sorting indices to X values and reformat suitably
[nr nc] = size(X);
X = [X; zeros(applen, nc)];
X = X(idsort, :);

% Replace appended zeros with previous molecular numbers for all X species
idcurr = zeros(size(TphEst));
idcurr(1) = find(Tnew == TphEst(1));
for i = 2:length(TphEst)+1
    % Obtain current photon event id
    if i ~= length(TphEst)+1
        idcurr(i) = find(Tnew == TphEst(i));
    end
    
    % Indices at the sort id must have previous value as no events
    % occurred at this point (at i = 1 its zero)
    if idcurr(i-1) ~= 1
        X(idcurr(i-1), 1:nc) = X(idcurr(i-1)-1, 1:nc);
    end
end

% Make z stream coherent with the other molecular species event times
z = [z; zeros(length(T), 1)];
Ttemp = [TphEst; T];
[Ttemp idsortz] = sort(Ttemp);
z = z(idsortz);

% Ensure z is consistent by matching previous values (counting process)
idcurrz = zeros(size(T));
idcurrz(1) = find(Ttemp == T(1));
for i = 2:length(T)+1
    % Obtain current photon event id
    if i ~= length(T)+1
        idcurrz(i) = find(Ttemp == T(i));
    end
    
    % Indices at the sort id must have previous value as no events
    % occurred at this point (at i = 1 its zero)
    if idcurrz(i-1) ~= 1
        z(idcurrz(i-1)) = z(idcurrz(i-1)-1);
    end
end 

% Append X with z vector and decrement the undelayed counting process at
% the times of z jumps due to stoichiometric coupling
X(:, nc) = X(:, nc) - z;
Xnew = [X(:, 1:nc) z];
    
