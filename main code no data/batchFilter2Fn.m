% Functional form for doubleBatch Gillespie simulations to run a batch of 
% filter simulations across various parameters and log the compound results

% Modified to include skip_Snyder and extra outputs to store the speed of
% the x1 and x2 reactions i.e. rates

% Modification to account for 3 mse calculation methods
% Modification to work with the correct filter code of runDSPP2
% Modification to compare filter mse results agains the linear LVP bound
% and plot the corresponding curves
function batchFilter2Fn(kbirth, simroot, beta_ratio, ind)

% Set batch numbers to determine simulations
batchNo = [1 0];

% Set relevant booleans and inputs to control code
inpB.save_sim = [1 0];
inpB.markov = 1;
inpB.skip_SSA = 0;
inpB.skip_Snyder = 0;
inpB.profile = 0;
inpB.plot_on = [0 0 0];
% inpB.simname = 'simF1';

% Set some initial simulation parameters
inpB.kbirth = kbirth;
inpB.kdeath = inpB.kbirth;
inpB.Nset = [20000 100000];
inpB.Slim = [0 1];

% Set initial birth and death rates and types
inpB.coeff = [0 0];
inpB.avgR = [inpB.coeff(1) inpB.coeff(1)];
inpB.birType = [7 3];
inpB.deaType = [2 0];

%%
% Cell to batch simulate a set encoder for various gain values and
% then at various state space sizes

% Set encoding function and birth-death forms for specific case as well as
% loop controlling variable
inpB.birType = [7 3];
inpB.deaType = [2 0];
% kgainSet = [0.0005:0.0005:0.01];
kgainSet = beta_ratio*inpB.kbirth;
lenB1 = length(kgainSet);

% Declare loop variables
msex1a1 = zeros(1, lenB1);
msex1a2 = zeros(1, lenB1);
msex1a3 = zeros(1, lenB1);
bnd1a = zeros(1, lenB1);
lamMeanSet1 = zeros(1, lenB1);
rateSet1 = cell(1, lenB1);
molecSet1 = cell(1, lenB1);

% Loop across the kgain = coeff(1) parameter defined as the set below and
% alter N value based on gain
if any(batchNo == 1)
    disp(['Starting batch simulation: ' num2str(1)]);
    for i = 1:lenB1

        % Alter coeff based on gain setting and run filter
        inpB.coeff(1) = kgainSet(i);
        inpB.simname = [simroot 'No' num2str(i)];
        outB = runDSPP2Fn(inpB);

        % Assign outputs and disp progress
        if ~inpB.skip_Snyder
            msex1a1(i) = outB.x1Stats.meth1(3);
            msex1a2(i) = outB.x1Stats.meth2(3);
            msex1a3(i) = outB.x1Stats.meth3(3);
            bnd1a(i) = outB.rawbnd;
        end
        lamMeanSet1(i) = outB.x2rateMean;
        rateSet1{i} = outB.rate;
        molecSet1{i} = outB.molec;
        disp(['Finished iteration: ' num2str(i) ' of ' num2str(lenB1)]);

    end
    disp(['Completed batch simulation: ' num2str(1)]);

    if ~inpB.skip_Snyder
        % Plot results for mse and bound
        figure;
        plot(kgainSet, msex1a1, 'bo-', kgainSet, bnd1a, 'ro-');
        xlabel('kgain');
        ylabel('mse of x1 or normalised lam mse and bound');
        legend('mse', 'bound', 'location', 'best');
        title(['Evolution of mse and bound with kgain for ' outB.birth{2} ' at Smax = ' num2str(inpB.Slim(2))]);
        saveas(gcf, 'batch1a.fig');

        % Plot results for mse and bound on log scale
        figure;
        semilogx(kgainSet, msex1a1, 'bo-', kgainSet, bnd1a, 'ro-');
        xlabel('kgain');
        ylabel('mse of x1 or normalised lam mse and bound');
        legend('mse', 'bound', 'location', 'best');
        title(['Evolution of mse and bound with kgain for ' outB.birth{2} ' at Smax = ' num2str(inpB.Slim(2))]);
        saveas(gcf, 'batch1aS.fig');

        % Plot ratio of mse to bound
        figure;
        plot(kgainSet, msex1a1./bnd1a, 'bo-');
        xlabel('kgain');
        ylabel('mse/bound for x1');
        title(['Relative performance mse/bound with kgain for ' outB.birth{2} ' at Smax = ' num2str(inpB.Slim(2))]);
        saveas(gcf, 'batch1b.fig');

        % Plot ratio of mse to bound on log scale
        figure;
        semilogx(kgainSet, msex1a1./bnd1a, 'bo-');
        xlabel('kgain');
        ylabel('mse/bound for x1');
        title(['Relative performance mse/bound with kgain for ' outB.birth{2} ' at Smax = ' num2str(inpB.Slim(2))]);
        saveas(gcf, 'batch1bS.fig');

        % Some further plots comparing the mse estimates
        figure;
        plot(kgainSet, [msex1a1' msex1a2' msex1a3']);
        xlabel('kgain');
        ylabel('mse estimates');
        title(['Comparison of mse estimates at Smax = ' num2str(inpB.Slim(2))]);
        saveas(gcf, 'mseComp1.fig');
    end
end

% Set changed booleans and control inputs for batch across state space
% SmaxSet = [10 100];
SmaxSet = 5:5:80;
lenB2 = length(SmaxSet);
inpB.Nset = [30000 31000];
inpB.coeff(1) = 1000;

% Declare loop variables
msex1b1 = zeros(1, lenB2);
msex1b2 = zeros(1, lenB2);
msex1b3 = zeros(1, lenB2);
bnd1b = zeros(1, lenB2);
lamMeanSet2 = zeros(1, lenB2);
rateSet2 = cell(1, lenB2);
molecSet2 = cell(1, lenB2);

% Loop across the kgain = coeff(1) parameter defined as the set below and
% alter N value based on gain
if any(batchNo == 2)
    disp(['Starting batch simulation: ' num2str(2)]);
    for i = 1:lenB2
        % Alter Smax based on space setting and run filter
        inpB.Slim(2) = SmaxSet(i);
        outB = runDSPP2Fn(inpB);

        % Assign outputs and disp progress
        if ~inpB.skip_Snyder
            msex1b1(i) = outB.x1Stats.meth1(3);
            msex1b2(i) = outB.x1Stats.meth2(3);
            msex1b3(i) = outB.x1Stats.meth3(3);
            bnd1b(i) = outB.rawbnd;
        end
        lamMeanSet1(i) = outB.x2rateMean;
        rateSet1{i} = outB.rate;
        molecSet1{i} = outB.molec;
        disp(['Finished iteration: ' num2str(i) ' of ' num2str(lenB2)]);

    end
    disp(['Completed batch simulation: ' num2str(2)]);

    if ~inpB.skip_Snyder
        % Plot results for mse and bound
        figure;
        plot(SmaxSet, msex1b1, 'bo-', SmaxSet, bnd1b, 'ro-');
        xlabel('Smax');
        ylabel('mse of x1 or normalised lam mse and bound');
        legend('mse', 'bound', 'location', 'best');
        title(['Evolution of mse and bound with Smax for ' outB.birth{2} ' at kgain = ' num2str(inpB.coeff(1))]);
        saveas(gcf, 'batch2a.fig');

        % Plot results for mse and bound on log scale
        figure;
        semilogx(SmaxSet, msex1b1, 'bo-', SmaxSet, bnd1b, 'ro-');
        xlabel('Smax');
        ylabel('mse of x1 or normalised lam mse and bound');
        legend('mse', 'bound', 'location', 'best');
        title(['Evolution of mse and bound with Smax for ' outB.birth{2} ' at kgain = ' num2str(inpB.coeff(1))]);
        saveas(gcf, 'batch2aS.fig');

        % Plot ratio mse/bound
        figure;
        plot(SmaxSet, msex1b1./bnd1b, 'bo-');
        xlabel('Smax');
        ylabel('mse/bound of x1');
        title(['Relative performance mse/bound with Smax for ' outB.birth{2} ' at kgain = ' num2str(inpB.coeff(1))]);
        saveas(gcf, 'batch2b.fig');

        % Plot ratio mse/bound on log scale
        figure;
        semilogx(SmaxSet, msex1b1./bnd1b, 'bo-');
        xlabel('Smax');
        ylabel('mse/bound of x1');
        title(['Relative performance mse/bound with Smax for ' outB.birth{2} ' at kgain = ' num2str(inpB.coeff(1))]);
        saveas(gcf, 'batch2bS.fig');

        % Some further plots comparing the mse estimates
        figure;
        plot(SmaxSet, [msex1b1' msex1b2' msex1b3']);
        xlabel('Smax');
        ylabel('mse estimates');
        title(['Comparison of mse estimates at kgain = ' num2str(inpB.coeff(1))]);
        saveas(gcf, 'mseComp2.fig');
    end
end


% Save results after removing unnecessary variables
clear Tset Yset outB
save(['batchDBndFil' num2str(ind)]);