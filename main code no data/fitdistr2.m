% Modified to account for data like inter-departure times which cannot be
% put into bins with min(data):max(data)

% Function to fit a distribution to data - assumes unit values between min
% and max of data
function [param m v stat] = fitdistr2(data, plot_on, distr, nbins)

% Obtain the bins from the data and the expected counts
if nbins == 0
    g = min(data):max(data);
    [obsCounts bins] = hist(data, g);
else
    [obsCounts bins] = hist(data, nbins);
end
cum_obsCounts = cumsum(obsCounts);
n = sum(obsCounts);

% Fit data to distrbution of choice and obtain counts and statistics
type = {'poiss', 'exp', 'gam', 'nbin', 'gauss'};
id = find(strcmpi(type, distr));
if isempty(id)
    error('Not a supported distribution entered');
end
switch(id)
    case 1
        % Poisson distribution
        param = poissfit(data);
        expCounts = n*poisspdf(bins, param);
        cum_expCounts = n*poisscdf(bins, param);
        [m v] = poisstat(param);
    case 2
        % Exponential distribution
        param = expfit(data);
        expCounts = n*exppdf(bins, param);
        cum_expCounts = n*expcdf(bins, param);
        [m v] = expstat(param);
    case 3
        % Gamma distribution
        param = gamfit(data);
        a = param(1);
        b = param(2);
        expCounts = n*gampdf(bins, a, b);
        cum_expCounts = n*gamcdf(bins, a, b);
        [m v] = gamstat(a, b);
    case 4
        % Negative binomial
        param = nbinfit(data);
        r = param(1);
        p = param(2);
        expCounts = n*nbinpdf(bins, r, p);
        cum_expCounts = n*nbincdf(bins, r, p);
        [m v] = nbinstat(r, p);
        
    case 5
        % Gaussian distribution
        [mu sigma] = normfit(data);
        param(1) = mu;
        param(2) = sigma;
        expCounts = n*normpdf(bins, mu, sigma);
        cum_expCounts = n*normcdf(bins, mu, sigma);
        [m v] = normstat(mu, sigma);
end

% Chi squared goodness of fit
[h1 p1 stat1] = chi2gof(bins, 'ctrs', bins, 'frequency', 1000*obsCounts/(sum(obsCounts)), ...
    'expected', 1000*expCounts/(sum(expCounts)), 'nparams', 1);
[h2 p2 stat2] = chi2gof(bins, 'ctrs', bins, 'frequency', 1000*cum_obsCounts/(sum(cum_obsCounts)), ...
    'expected', 1000*cum_expCounts/(sum(cum_expCounts)), 'nparams', 1);
stat.h = [h1 h2];
stat.p = [p1 p2];
stat.stats{1} = stat1;
stat.stats{2} = stat2;

% Plot the fit
if plot_on
    figure;
    plot(bins, obsCounts, bins, expCounts);
    hold on
    legend('observed', 'expected', 'Location', 'Best');
    title(['Comparison to ' distr ' pdf for param = ' num2str(param)]);
    xlabel('data');
    ylabel('frequency');
    grid;
    hold off
    figure;
    plot(bins, cum_obsCounts, bins, cum_expCounts);
    hold on
    legend('observed', 'expected', 'Location', 'Best');
    title(['Comparison to ' distr ' cdf for param = ' num2str(param)]);
    xlabel('data');
    ylabel('cumulative frequency');
    grid;
    hold off
end