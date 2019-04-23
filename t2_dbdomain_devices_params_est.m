%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shalin Shah
% Affiliation: Dept. of Electrical & Computer Engineering, Duke University
% Email: shalin.shah@duke.edu
% Last modified: 01/08/2019
% Matlab version used: R2017a
%
% Description: This code tunes number of domains of DNA devices, generates 
% barcodes, analyze signals and generates scatter plot with estimated parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; clf;

DEBUGGING = 0;
NO_OF_SAMPLES = 100;
NO_OF_DEVICES = 16;
BASE_DEVICE_NO = 4;
COLLECTION_TIME = 1800;
legend_str = [];

% device avg times for 8 nt, 9 nt, 10 nt
dev_on_times = [0.005, 0.05, 0.55, 9.0];
% 16 distinct colors
C = [0,0,1;1,0,0;0,1,0;0,0,0.172413793103448;1,0.103448275862069,0.724137931034483;...
    1,0.827586206896552,0;0,0.344827586206897,0;...
    0.517241379310345,0.517241379310345,1;0.620689655172414,0.310344827586207,0.275862068965517;...
    0,1,0.758620689655172;0.517241379310345,0.58  6206896551724;...
    0,0,0.482758620689655;0.586206896551724,0.827586206896552,0.310344827586207;...
    0.965517241379310,0.620689655172414,0.862068965517241;...
    0.827586206896552,0.0689655172413793,1;...
    0.482758620689655,0.103448275862069,0.413793103448276];

single_blink_data = zeros(NO_OF_SAMPLES, NO_OF_DEVICES);
double_blink_data = zeros(NO_OF_SAMPLES, NO_OF_DEVICES);
group = zeros(NO_OF_SAMPLES, NO_OF_DEVICES);
for i=1:NO_OF_DEVICES
    group(:, i) = i;
end

for k = 1:NO_OF_SAMPLES
    single_blink = []; double_blink = [];
    device_count = 1;
    for i = 3:BASE_DEVICE_NO
        for j = i:BASE_DEVICE_NO
            fprintf('Generating device %d, %d\n', i+6, j+6);
             if k == 1 
                legend_str = [legend_str; sprintf('Generating device %d, %d\n', i, j)];
             end
            
            [single_blink_est, double_blink_est] = ...
                simulate(dev_on_times(i), dev_on_times(j), COLLECTION_TIME, DEBUGGING);
            single_blink_data(k, device_count) = (single_blink_est.mean);
            double_blink_data(k, device_count) = (double_blink_est.mean);
            figure(1); hold on;
            p = plot(single_blink_data(group==device_count), ...
                double_blink_data(group==device_count),'.','markersize', 10,...
                'color', C(device_count, :)); 
            
            device_count = device_count + 1;
        end
    end
end


box on; set(gca,'fontsize',20); set(gca,'fontweight','bold');
set(gca,'linew',3.0); set(gca,'yscale','log'); set(gca,'xscale','log');
%legend(legend_str);
%% Andor iXon 897 full-frame rate 56 Hz
p = plot([1e-3, 1e-3], [1e-4, 1e2], ':'); 
set(p,'Color', [0.01 0.01 0.01]); set(p,'linew',3.0);
p = plot([1e-3, 1e2], [1e-3, 1e-3], ':');
set(p,'Color', [0.01 0.01 0.01]); set(p,'linew',3.0);


function [single_blink_est, double_blink_est] = simulate(avg_on_timeA, avg_on_timeB, TOTAL_TIME, DEBUGGING)
    
    %% Rate constants for imager strand from Jungmann et. al (2010)
    % 8nt - 20 sec^-1, 9nt - 2 sec^-1, 10nt - 0.11 sec^-1 off-rate
    imagerC = 5;                                   % [nM]
    kon = 1e6;                                      % molar^-1 sec^-1
    off_rateA = 1/avg_on_timeA;                     % sec^-1
    off_rateB = 1/avg_on_timeB;                     % sec^-1
    on_rate = kon * imagerC * 1e-9;                 % sec^-1 

    [time, signal] = simulate_reaction(off_rateA, off_rateB, on_rate, TOTAL_TIME, DEBUGGING);

    %% Generate statistics for different length of time-signal
    figure(2);
%     subplot(1, 2, 1);
    [single_blink_time, double_blink_time, ~] = genOnOffStats(time, signal, time(end));
    if ~isempty(single_blink_time)
%         histfit(single_blink_time, 3, 'exponential');
        single_blink_est = fitdist(single_blink_time, 'exponential');
    %     xlabel('single blink time'); ylabel('counts'); 
    %     legend(sprintf('Capture time: %d\nImager Conc.: %d nM\n', TOTAL_TIME, imagerC));
    else
        single_blink_est = makedist('exponential', 'mu', 1e-3);
    end
    
    if ~isempty(double_blink_time)
        double_blink_est = fitdist(double_blink_time, 'exponential');
    else
        double_blink_est = makedist('exponential', 'mu', 1e-3);
    end

    if DEBUGGING
        figure(3);
        subplot(2, 2, 1); histogram(single_blink_time, 4); title('single-blink plot');
        subplot(2, 2, 2); histogram(double_blink_time, 4); title('double-blink plot');
        subplot(2, 2, [3, 4]); stairs(time, signal, 'LineWidth', 2.0); xlabel('time (s)'); ylabel('Intensity'); title('intensity plot'); 
        fprintf('Mean single-blink time: %f s\nMean double-blink time: %f s\n', mean(single_blink_time), mean(double_blink_time));

        n_signal = interpolate(10, signal, time, TOTAL_TIME);
    end
end

function [t_ssa, signal] = simulate_reaction(off_rateA, off_rateB, on_rate, TOTAL_TIME, DEBUGGING)
    %% Reaction network:
    %  1. imagerA + docker {off}<-->{on} imagerA-docker
    %  2. imagerB + docker {off}<-->{on} imagerB-docker
    %  3. imagerA-docker + imagerB {off}<-->{on} imagerA-docker-imagerB
    %  4. imagerB-docker + imagerA {off}<-->{on} imagerA-docker-imagerB

    %% Create Decaying-Dimerizing Model
    model = sbiomodel('Tuning the number of domains');

    %% Enter Reactions
    r1 = addreaction(model, 'I + D <-> I-Da');
    r2 = addreaction(model, 'I + D <-> I-Db');
    r3 = addreaction(model, 'I-Da + I <-> I-Dab');
    r4 = addreaction(model, 'I-Db + I <-> I-Dab');

    %% Set Reactions to be MassAction
    kl1 = addkineticlaw(r1, 'MassAction');
    kl2 = addkineticlaw(r2, 'MassAction');
    kl3 = addkineticlaw(r3, 'MassAction');
    kl4 = addkineticlaw(r4, 'MassAction');

    %% Add Rate Constants for Each Reaction
    addparameter(kl1, 'k_on', 'Value', on_rate);
    addparameter(kl1, 'k_off_a', 'Value', off_rateA);
    addparameter(kl2, 'k_on', 'Value', on_rate);
    addparameter(kl2, 'k_off_b', 'Value', off_rateB);
    addparameter(kl3, 'k_on', 'Value', on_rate);
    addparameter(kl3, 'k_off_b', 'Value', off_rateB);
    addparameter(kl4, 'k_on', 'Value', on_rate);
    addparameter(kl4, 'k_off_a', 'Value', off_rateA);


    %% Set the Kinetic Law Constants for Each Kinetic Law.
    kl1.ParameterVariableNames = {'k_on', 'k_off_a'};
    kl2.ParameterVariableNames = {'k_on', 'k_off_b'};
    kl3.ParameterVariableNames = {'k_on', 'k_off_b'};
    kl4.ParameterVariableNames = {'k_on', 'k_off_a'};

    %% Specify Initial Amounts of Each Species
    model.species(1).InitialAmount = 10;         % I
    model.species(2).InitialAmount = 1;          % D


    %% Display the Completed Model Objects
    model;

    %% Display the Reaction Objects
    model.Reactions;


    %% Display the Species Objects
    model.Species;

    %% Get the Active Configuration Set for the Model.
    cs = getconfigset(model,'active');

    %% Simulate Model Using SSA Stochastic Solver and Plot
    % tfinal = 30, logging every 10th datapoint.
    cs.SolverType = 'ssa';
    cs.StopTime = TOTAL_TIME;
    solver = cs.SolverOptions;
    solver.LogDecimation = 1;
    cs.CompileOptions.DimensionalAnalysis = false;
    [t_ssa, x_ssa] = sbiosimulate(model);

    if DEBUGGING
        subplot(4, 1, 1); stairs(t_ssa, x_ssa(:, 1:2),'LineWidth', 2);
        legend('imager', 'docker')
        axis([0 cs.StopTime 0 model.Species(1).InitialAmount])
        subplot(4, 1, 2); stairs(t_ssa, x_ssa(:, 3:5),'LineWidth', 2);
        legend('I-Da', 'I-Db', 'I-Dab')
        axis([0 cs.StopTime 0 model.Species(1).InitialAmount])
        subplot(4, 1, 3); stairs(t_ssa, x_ssa(:, 3) + x_ssa(:, 4) + 2*x_ssa(:, 5),'LineWidth', 2);
        axis([0 cs.StopTime 0 model.Species(1).InitialAmount])
        subplot(4, 1, 4); stairs(t_ssa, x_ssa(:, 2) + x_ssa(:, 3) + x_ssa(:, 4) + x_ssa(:, 5),'LineWidth', 2);
        axis([0 cs.StopTime 0 model.Species(1).InitialAmount])
    end
    signal = x_ssa(:, 3) + x_ssa(:, 4) + 2*x_ssa(:, 5);
    iSignal = interpolate(
end

function [single_blink, double_blink, off_times] = genOnOffStats(time, signal, trun_time)
    
    CAPTURE_TIME = 0.001;
    %% analyse the on off-time
    off_times = [];
    %% This will be total on-time
    single_blink = [];
    double_blink = [];

    last_off_step = 1;

    %% only the state change, time is recorded except last state
    for step = 2 : length(signal)-1
        if time(step) < trun_time
            cur_step_len = time(step) - time(step-1);
            if signal(step) == 1 && signal(step-1) == 0 && cur_step_len >= CAPTURE_TIME
                off_times = [off_times; cur_step_len];
                last_off_step = step;
            elseif signal(step) == 1 && signal(step-1) == 2 && cur_step_len >= CAPTURE_TIME
                double_blink = [double_blink; cur_step_len];
            elseif signal(step) == 0 && signal(step-1) == 1 && cur_step_len >= CAPTURE_TIME
                    single_blink = [single_blink; time(step)-time(last_off_step)];
            end
        else
            break
        end
    end
end

function nSignal = interpolate(freq, signal, time)
    % This function takes the original ODE signal containing 0 and 1 and  
    % interpolates it with given frequency
    %
    % freq is in hertz (eg. 10 Hz)
    %
    % signal is the original CTMC signal with only 0 and 1
    %
    % time is the time vector of the original signal
    % 
     
    nSignal = zeros(1, time(end)*freq);
    
    for step = 2:length(time)-1
        if signal(step) == 0 && signal(step-1) == 1
            nSignal(int32(time(step-1))*100:int32(time(step))*100) = 1.0;
        elseif signal(step) == 1 && signal(step-1) == 2
            nSignal(int32(time(step-1))*100:int32(time(step))*100) = 2.0;
        elseif signal(step) == 2 && signal(step-1) == 1
            nSignal(int32(time(step-1))*100:int32(time(step))*100) = 1.0;
        end
    end
    figure;
    subplot(2, 1, 1); stairs(time, signal);
    time = 0:TOTAL_TIME/length(nSignal):TOTAL_TIME;
    subplot(2, 1, 2); plot(nSignal+ randn(size(nSignal))/10);
    
%     figure(3); plot(time(2:end), n_signal+ randn(size(n_signal))/10);
%     xlabel('time (seconds)'); ylabel('Intensity (a.u.)');
end