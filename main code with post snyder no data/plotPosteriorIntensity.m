% Modification to include relative speed in title and to plot the error
% curve and to save figures if boolean demands it

% Function to produce plots of modulating intensity estimates and a view of
% how the posterior distribution evolves with KL divergence calculation
function qKL = plotPosteriorIntensity(T, q, lam, lamcap, S, kg, kb, kd, reldyn, save_sim, plot_on)

% Terminate if plots not desired
if ~plot_on
    qKL = 0;
    return;
end

% Check input dimensions
[nr nc] = size(q);
if length(T) ~= nr || length(lam) ~= length(lamcap)...
        || length(lamcap) ~= length(T)
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
if any([kg kb kd] < 1)
    plotname = 0;
else
    plotname = 1;
end

% Plot the actual and estimated intensity
figure;
plot(T, [lam lamcap]);
xlim([min(T) max(T)])
legend('actual intensity', 'estimate', 'location', 'best');
xlabel('time');
ylabel('modulating intensity');
if tit_id
    title(['Intensity comparison at [kg kb kd relspd] = ' num2str(kg) ' '...
        num2str(kb) ' ' num2str(kd) ' ' num2str(reldyn)]);
else
    title('Comparison of actual and estimated intensity');
end
if save_sim
    if plotname
        saveas(gcf, ['snyderA_' num2str(kg) '_' num2str(kb) '_' num2str(kd)]);
    else
        saveas(gcf, 'snyderAinv.fig');
    end
end

% Plot the error curve between actual and estimated intensity
figure;
plot(T, lam - lamcap);
xlim([min(T) max(T)])
xlabel('time');
ylabel('estimation error');
if tit_id
    title(['Intensity error curve at [kg kb kd relspd] = ' num2str(kg) ' '...
        num2str(kb) ' ' num2str(kd) ' ' num2str(reldyn)]);
else
    title('Error between actual and estimated intensity');
end
if save_sim
    if plotname
        saveas(gcf, ['snyderD_' num2str(kg) '_' num2str(kb) '_' num2str(kd)]);
    else
        saveas(gcf, 'snyderDinv.fig');
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

% Plot the successive KL divergence
figure;
plot(qKL);
xlabel('starting iteration of calculation');
ylabel('KL divergence');
if tit_id
    title(['Successive KL divergence at [kg kb kd relspd] = ' num2str(kg) ' '...
        num2str(kb) ' ' num2str(kd) ' ' num2str(reldyn)]);
else
    title('Successsive KL divergence across iterations');
end
if save_sim
    if plotname
        saveas(gcf, ['snyderB_' num2str(kg) '_' num2str(kb) '_' num2str(kd)]);
    else
        saveas(gcf, 'snyderBinv.fig');
    end
end

% Plot a sample of at most 10 posterior distributions
state = S(S>0);
state = [0 state'];
if nr > 10
    inter = nr/10;
    seq = 1:inter:nr;
else
    seq = 1:nr;
end
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
    title(['Posterior density at [kg kb kd relspd] = ' num2str(kg) ' '...
        num2str(kb) ' ' num2str(kd) ' ' num2str(reldyn)]);
else
    title('Comparison of posterior density across iterations');
end
if save_sim
    if plotname
        saveas(gcf, ['snyderC_' num2str(kg) '_' num2str(kb) '_' num2str(kd)]);
    else
        saveas(gcf, 'snyderCinv.fig');
    end
end