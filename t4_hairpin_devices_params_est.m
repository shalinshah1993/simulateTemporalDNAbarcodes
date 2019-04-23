%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shalin Shah
% Affiliation: Dept. of Electrical & Computer Engineering, Duke University
% Email: shalin.shah@duke.edu
% Last modified: 01/08/2019
% Matlab version used: R2017a
%
% Description: This code tunes length of haripin DNA devices, generates 
% barcodes, analyze signals and generates scatter plot with estimated dark 
% time parameter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; clf;

DEBUGGING = 0;
NO_OF_SAMPLES = 100;
NO_OF_DEVICES = 5;
COLLECTION_TIME = 3600;

%% device avg times
hp_on_time = [0.001, 0.01, 0.06, 0.6, 9.0];
reporter_on_time = 0.11;

%% 16 distinct colors
C = [0,0,1;1,0,0;0,1,0;0,0,0.172413793103448;1,0.103448275862069,0.724137931034483;...
    1,0.827586206896552,0;0,0.344827586206897,0;...
    0.517241379310345,0.517241379310345,1;0.620689655172414,0.310344827586207,0.275862068965517;...
    0,1,0.758620689655172;0,0.517241379310345,0.586206896551724;...
    0,0,0.482758620689655;0.586206896551724,0.827586206896552,0.310344827586207;...
    0.965517241379310,0.620689655172414,0.862068965517241;...
    0.827586206896552,0.0689655172413793,1;...
    0.482758620689655,0.103448275862069,0.413793103448276];

on_time_data = zeros(NO_OF_SAMPLES, NO_OF_DEVICES);
off_time_data = zeros(NO_OF_SAMPLES, NO_OF_DEVICES);
group = zeros(NO_OF_SAMPLES, NO_OF_DEVICES);

for i=1:NO_OF_DEVICES
    group(:, i) = i;
end

legend_str = [];
for k = 1:NO_OF_SAMPLES
    device_count = 1;
    fprintf('Generating sample %d\n', k);
    for i = 1:NO_OF_DEVICES
            if k == 1 
                legend_str = [legend_str; sprintf('Generating device %d, %d\n', i, j)];
            end
            
            [on_time_est, off_time_est] = ...
                simulate(reporter_on_time, hp_on_time(i), COLLECTION_TIME, DEBUGGING);
            on_time_data(k, device_count) = (on_time_est.mean);
            off_time_data(k, device_count) = (off_time_est.mean);
                       
            figure(1); 
            hold on;
            plot(off_time_data(group==device_count), ...
                '.','markersize',10, 'color', C(device_count, :)); 
            
            device_count = device_count + 1;
    end
end

box on; set(gca,'fontsize',20); set(gca,'fontweight','bold');
set(gca,'linew',3.0); 
set(gca,'yscale','log');

function [on_time_est, off_time_est] = simulate(avg_on_time, h_avg_on_time, TOTAL_TIME, DEBUGGING)
    
    %% Rate constants for imager from Jungmann et. al (2010)
    %% Rate constants for hairpins from Tsukanov et. al (2013)
    % 8nt - 20 sec^-1, 9nt - 2 sec^-1, 10nt - 0.11 sec^-1 off-rate
    imagerC = 10;                                 % [nM]
    kon = 1e6;                                    % molar^-1 sec^-1
    off_rate = 1/avg_on_time;                     % sec^-1
    on_rate = kon * imagerC * 1e-9;               % sec^-1 
    h_off_rate = 1/h_avg_on_time;                 % sec^-1 
    h_on_rate = 100;                              % sec^-1 

    [time, signal] = simulate_reaction(off_rate, on_rate, h_off_rate, h_on_rate, TOTAL_TIME, DEBUGGING);

    %% Generate statistics for different length of time-signal
    [on_time, off_time] = genOnOffStats(time, signal);
    if ~isempty(on_time)
        on_time_est = fitdist(on_time, 'exponential');
    else
        on_time_est = makedist('exponential', 'mu', 1e-3);
    end
    
    if ~isempty(off_time)
        off_time_est = fitdist(off_time, 'exponential');
    else
        off_time_est = makedist('exponential', 'mu', 1e-3);
    end

%      figure(2);
%      subplot(2, 2, 1); histogram(on_time, 4); title('on-time plot');
%      subplot(2, 2, 2); histogram(off_time, 4); title('off-time plot');
%      fprintf('Mean on-time: %f s\nMean off-time: %f s\n', mean(on_time_est), mean(off_time_est));

    %n_signal = interpolate(10, signal, time, TOTAL_TIME);
end

function [t_ssa, signal] = simulate_reaction(off_rate, on_rate, h_off_rate, h_on_rate, TOTAL_TIME, DEBUGGING)
    %% Reaction network:
    %  1. imager + hairpinO {off}<-->{on} imager-hairpinO
    %  2. hairpinO {h_off}<-->{h_on} hairpinC

    %% Create Decaying-Dimerizing Model
    model = sbiomodel('Tuning the number of domains');

    %% Enter Reactions
    r1 = addreaction(model, 'I + Ho <-> I-Ho');
    r2 = addreaction(model, 'Ho <-> Hc');
    
    %% Set Reactions to be MassAction
    kl1 = addkineticlaw(r1, 'MassAction');
    kl2 = addkineticlaw(r2, 'MassAction');
    
    %% Add Rate Constants for Each Reaction
    addparameter(kl1, 'k_on', 'Value', on_rate);
    addparameter(kl1, 'k_off', 'Value', off_rate);
    addparameter(kl2, 'h_on', 'Value', h_on_rate);
    addparameter(kl2, 'h_off', 'Value', h_off_rate);


    %% Set the Kinetic Law Constants for Each Kinetic Law.
    kl1.ParameterVariableNames = {'k_on', 'k_off'};
    kl2.ParameterVariableNames = {'h_on', 'h_off'};
    
    %% Specify Initial Amounts of Each Species
    model.species(1).InitialAmount = 5;         % I
    model.species(4).InitialAmount = 1;          % D

    %% Display the Completed Model Objects
    model;

    %% Display the Reaction Objects
    model.Reactions;

    %% Display the Species Objects
    model.Species;

    %% Get the Active Configuration Set for the Model.
    cs = getconfigset(model,'active');

    %% Simulate Model Using SSA Stochastic Solver and Plot
    cs.SolverType = 'ssa';
    cs.StopTime = TOTAL_TIME;
    solver = cs.SolverOptions;
    solver.LogDecimation = 1;
    cs.CompileOptions.DimensionalAnalysis = false;
    [t_ssa, x_ssa] = sbiosimulate(model);

    if DEBUGGING
        subplot(2, 1, 1); stairs(t_ssa, x_ssa(:, [1, 3]),'LineWidth', 2);
        legend('imager', 'output')
        axis([0 cs.StopTime 0 model.Species(1).InitialAmount])
        subplot(2, 1, 2); stairs(t_ssa, x_ssa(:, [2, 4]),'LineWidth', 2);
        legend('hairpin open', 'hairpin closed')
        axis([0 cs.StopTime 0 model.Species(1).InitialAmount])
    end
    signal = x_ssa(:, 3);
end

function [on_times, off_times] = genOnOffStats(time, signal)
    DETECT_LIMIT = 0.001;
    
    %% analyse the on off-time
    off_times = [];
    on_times = [];

    last_off_step = 1;
    last_on_step = 1;

    %% only the state change, time is recorded except last state
    for step = 2 : length(signal)
        if step == length(signal)
            if signal(step) == 1
                on_times = [on_times; time(step)-time(last_off_step)];
            else
                off_times = [off_times; time(step)-time(last_on_step)];
            end
        end

        if signal(step) == 1 && signal(step-1) == 0
            off_times = [off_times; time(step)-time(last_on_step)];
            last_off_step = step;    
       elseif signal(step) == 0 && signal(step-1) == 1
            on_times = [on_times; time(step)-time(last_off_step)];
            last_on_step = step;
        end
    end
end

function n_signal = interpolate(freq, signal, time, TOTAL_TIME)
     %% Fill the values 
    n_signal = zeros(1, time(end)*freq);
    n_time = 1:time(end)*freq;
    last_state = 0;
    for step = 2:length(time)-1
        if signal(step) == 0 && signal(step-1) == 1
            n_signal(int32(time(step-1))*100:int32(time(step))*100) = 1.0;
        elseif signal(step) == 1 && signal(step-1) == 2
            n_signal(int32(time(step-1))*100:int32(time(step))*100) = 2.0;
        elseif signal(step) == 2 && signal(step-1) == 1
            n_signal(int32(time(step-1))*100:int32(time(step))*100) = 1.0;
        end
    end
    figure(4);
    subplot(2, 1, 1); stairs(time, signal);
    time = 0:TOTAL_TIME/length(n_signal):TOTAL_TIME;
    subplot(2, 1, 2); plot(n_signal+ randn(size(n_signal))/10);
    
    figure(3); plot(time(2:end), n_signal+ randn(size(n_signal))/10);
    xlabel('time (seconds)'); ylabel('Intensity (a.u.)');
end