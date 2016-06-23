% Function converts the Q matrix to an adjacency matrix A and plots the
% corresponding Markov chain graph
function [A coords] = plotMCGraph(Q, states)

% Default state labels if none are provided
nStates = length(Q);
if nargin == 1
    states = 0:(nStates-1);
end

% Obtain adjaceny matrix from Q without self loops
[rid cid] = find(Q > 0);
A = zeros(size(Q));
for i = 1:length(rid)
    A(rid(i), cid(i)) = 1;
end

% Obtain number of nodes and connections
nNodes = nStates/2;
nCons = length(rid);

% Set positions for the nodes (states) for a layout in which half the
% states are above the other half
horizSpace = 10;
initx = 10;
finx = initx + (nNodes-1)*horizSpace;
xpos = initx:horizSpace:finx;
xpos = [xpos xpos];
vertSpace = 50;
inity = 50;
finy = inity + vertSpace;
ypos = [inity*ones(1, nNodes) finy*ones(1, nNodes)];

% Check the size of the position vectors
if length(xpos) ~= nStates || length(ypos) ~= nStates
    assignin('base', 'xpos', xpos);
    assignin('base', 'ypos', ypos);
    error('Incorrect number of node coordinates');
end

% Define labels from states
stateLabel = cell(1, nStates);
for i = 1:nStates
    stateLabel{i} = num2str(states(i));
end

% Plot the connections with scaled arrows in turn
scale = 1;
figure;
hold on
for i = 1:nCons
    % Obtain the start and end coordinates of each connection
    startNode = rid(i);
    endNode = cid(i);
    startx = xpos(startNode);
    starty = ypos(startNode);
    endx = xpos(endNode);
    endy = ypos(endNode);
    % plot([startx endx], [starty endy], 'r', 'LineWidth', 5);
    
    % Plot vectors with the appropriate directions
    quiver(startx, starty, scale*(endx-startx), scale*(endy-starty), 'r', 'filled', 'LineWidth', 2);
end

% Plot the nodes at suitable positions
for i = 1:nStates
    plot(xpos(i), ypos(i), 'bo', 'MarkerSize', 25, 'MarkerFaceColor', [.49 1 .63]);
    text(xpos(i) - 0.2, ypos(i), stateLabel{i}); 
end
coords = [xpos' ypos'];

% Polish the plots
xlim([min(xpos) - horizSpace, max(xpos) + horizSpace]);
ylim([min(ypos) - horizSpace, max(ypos) + horizSpace]);

% Subfunction that plots different types of curves based on the nodes
function plotNodeConnections(N1, N2, x1, y1, x2, y2, dh)

% If the nodes are adjacent states then a straight line is sufficient
% provided the nodes also code the same direction (same horizontal)
if abs(N2 - N1) == 1 && (x2 - x1 == dh)
    quiver(x1, y1, x2-x1, y2-y1, 'r', 'filled', 'LineWidth', 2);
end

% If the nodes are not adjacent then use a parabolic arc constructed using
% a parametric variable t
if abs(N2 - N1) ~= 1
    t = (0:0.01:1)';
    % Straight line bwtween the points
    L1 = (1-t)*N1 + t*N2;
    % Add a parabolic offset as a function of C (C = 0 is a straight line)
    C = 0.5;
    N = (N1-N2)*[0 -1;1 0];
    L2 = L1 + C*(t.*(1-t))*N;
end 