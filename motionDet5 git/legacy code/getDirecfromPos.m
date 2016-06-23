% Warning - this is not an optimal way of obtaining a directional estimate
% from a state estimate - use the posterior distribution in that case

% Simple function to produce directional estimates assuming that the first
% half of states code left and the second half right
function [direc discrimin] = getDirecfromPos(pos, posSpace)

% Check that the position is within the posSpace
maxPos = max(pos);
minPos = min(pos);
if maxPos > posSpace(end) || minPos < posSpace(1)
    assignin('base', 'pos', pos);
    assignin('base', 'posSpace', posSpace);
    error('The position has values outside the limits of its space');
end

% Obtain state which gives the threshold between the directions and ensure
% that the number of states is even
if rem(length(posSpace), 2) == 1
    assignin('base', 'posSpace', posSpace);
    error('Code assumes an even no. states');
end
discrimin = posSpace(length(posSpace)/2);

% Convert position to direction under assumption of states representing
% both position and direction
direc = pos;
direc(pos > discrimin) = 1;
direc(pos <= discrimin) = -1;