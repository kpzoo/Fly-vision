% Removed the kb and kd parameters as they are not useful for general MC
% formulations

% Major correction to work with the full ODE solutions and the
% corresponding times which are obtained from filterHybrid3 

% Modification to include plot with lamcap referenced to x2 events 
% Modification to include an extra plot for the x and xcap estimate
% Modification to include relative speed in title and to plot the error
% curve and to save figures if boolean demands it

% Function to produce plots of modulating intensity estimates and a view of
% how the posterior distribution evolves with KL divergence calculation
function plotPosteriorMC(Tn, q, lamn, lamcapn, x1n, x1capn, T, x2, S, kg, reldyn, save_sim, plot_on)

% No plots desired
if ~plot_on
    return;
end

% Check input dimensions
[nr nc] = size(q);
if length(lamn) ~= length(lamcapn) || length(lamcapn) ~= length(Tn)
    error('Input dimensions mismatched');
end

% Determine title boolean and naming of plots boolean
if nargin == 5
    tit_id = 0;
elseif nargin >= 8
    tit_id = 1;
else
    error(['Inconsistent number of arguments, nargin = ' num2str(nargin)]);
end
if kg < 1
    plotname = 0;
else
    plotname = 1;
end

% Plot the actual and estimated intensity
figure;
plot(Tn, [lamn lamcapn]);
xlim([min(Tn) max(Tn)])
legend('actual intensity', 'estimate', 'location', 'best');
xlabel('time');
ylabel('modulating intensity');
if tit_id
    title(['Intensity comparison at [kg relspd] = ' num2str(kg) ' ' num2str(reldyn)]);
else
    title('Comparison of actual and estimated intensity');
end
if save_sim
    if plotname
        saveas(gcf, ['snyderA_' num2str(kg)]);
    else
        saveas(gcf, 'snyderAinv.fig');
    end
end

% Plot the error curve between actual and estimated intensity
figure;
plot(Tn, lamn - lamcapn);
xlim([min(Tn) max(Tn)])
xlabel('time');
ylabel('estimation error');
if tit_id
    title(['Intensity error curve at [kg relspd] = ' num2str(kg) ' ' num2str(reldyn)]);
else
    title('Error between actual and estimated intensity');
end
if save_sim
    if plotname
        saveas(gcf, ['snyderD_' num2str(kg)]);
    else
        saveas(gcf, 'snyderDinv.fig');
    end
end

% Plot the actual and estimated population of x1
figure;
plot(Tn, [x1n x1capn]);
xlim([min(Tn) max(Tn)])
legend('actual x1', 'estimate', 'location', 'best');
xlabel('time');
ylabel('x1 population number');
if tit_id
    title(['Population comparison at [kg relspd] = ' num2str(kg) ' ' num2str(reldyn)]);
else
    title('Comparison of actual and estimated x1 population');
end
if save_sim
    if plotname
        saveas(gcf, ['snyderE_' num2str(kg)]);
    else
        saveas(gcf, 'snyderEinv.fig');
    end
end

% Characterise the closeness of the posterior densities via KL divergence
qKL = zeros(1, nr-1);
for i = 2:nr
    for j = 1:nc
        % Account for 0log0 and 0log(1/0) limits
        if q(i, j) > 10^-9 && q(i-1, j) > 10^-9
            qKL(i-1) = q(i, j)*log(q(i, j)/q(i-1, j)) + qKL(i-1);
        end
    end
end
assignin('base', 'qKL', qKL);

% Plot the successive KL divergence
figure;
plot(1:length(qKL), qKL);
xlim([1 length(qKL)]);
xlabel('starting iteration of calculation');
ylabel('KL divergence');
if tit_id
    title(['Successive KL divergence at [kg relspd] = ' num2str(kg) ' ' num2str(reldyn)]);
else
    title('Successsive KL divergence across iterations');
end
if save_sim
    if plotname
        saveas(gcf, ['snyderB_' num2str(kg)]);
    else
        saveas(gcf, 'snyderBinv.fig');
    end
end

% Plot a sample of at most 10 posterior distributions
state = S(S>0);
state = [0 state'];
if nr > 10
    inter = round(nr/10);
    seq = 1:inter:nr;
else
    seq = 1:nr;
end
% assignin('base', 'q', q);
% assignin('base', 'seq', seq);

figure;
plot(state, q(seq(1), :), 'ro-');
hold on
for i = 2:length(seq)-1
    plot(state, q(seq(i), :));
end
plot(state, q(seq(length(seq)), :), 'ko-');
hold off
xlabel('discrete state space');
ylabel('posterior density');
if tit_id
    title(['Posterior density at [kg relspd] = ' num2str(kg) ' ' num2str(reldyn)]);
else
    title('Comparison of posterior density across iterations');
end
if save_sim
    if plotname
        saveas(gcf, ['snyderC_' num2str(kg)]);
    else
        saveas(gcf, 'snyderCinv.fig');
    end
end

% Obtain and plot lamcap with ODE points with event times of x2
[Tevent perc] = getEventTimes(T, x2, 'birth');
% figure;
% stairs(Tn, lamcapn);
% set(gca, 'XTick', Tevent);
% set(gca, 'XGrid', 'on');
% set(gca, 'XTickLabel', []);
% xlabel('x2 birth event times');
% ylabel('lamcap');
% title('Behaviour of lamcap with reference to x2 birth times');