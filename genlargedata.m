%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shalin Shah
% Affiliation: Dept. of Electrical & Computer Engineering, Duke University
% Email: shalin.shah@duke.edu
% Last modified: 01/08/2019
% Matlab version used: R2017a
%
% Description: This code generate nSample for nDevice. Each simulation
% experiment will be exprmntDrtn long and will assume a DNA nanostructure
% with maxStapleCount staples. The unbinding rates are adopted from
% Jungmann et al. (2010). Output feature vector for each device is stored
% in a csv file each row is a data vector and first maxStapleCount columns
% are feature vector. Rest of the columns indicate device label.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear workspace to avoid variable conflicts
clear; clc;
fprintf('RUNNING CODE TO SIMULATE TEMPORAL DNA BARCODES...\n');
% no of repeats for the simulation experiments
nSample = 100;
% total device combinations for the given number staples extension
nDevice = 56;
% length of temporal barcode collected using simulation experiment
exprmntDrtn = 5*3600;
% nanostructure has these many extended staples acting as device
maxStapleCount = 5;
% device avg times for 8 nt, 9 nt, 10 nt
deviceOnTime = [0.01, 0.06, 0.55, 9.0];
% training data X and training data output t
X = [];

% each number indicates (#|7|, #|8|, #|9|, #|10|)
% This is sample solution for low reporter regime
% sampleSoln = [1 0 0 0;
%               2 1 0 0;
%               1 1 0 0;
%               1 2 0 0;
%               0 4 0 0;
%               1 2 1 0;
%               0 1 1 0;
%               0 0 4 0;
%               3 0 0 1;
%               0 0 2 1;
%               0 0 2 2;
%               0 0 0 1];    
% This is sample solution for high reporter regime
% sampleSoln = [1 0 0 0;
%               2 1 0 0;
%               1 1 0 0;
%               1 2 0 0;
%               0 4 0 0;
%               1 2 1 0;
%               0 1 1 0;
%               0 0 4 0;
%               3 0 0 1;
%               0 0 2 1;
%               0 0 2 2;
%               0 0 0 1];    
% SET MAXSTAPLECOUNT TO 4

% simulate temporal barcodes and assign ground truth labels
for iSample = 1:nSample
    fprintf('Generating sample %d\n', iSample);
    iDevice = 1;
    for iSeven = 0 : maxStapleCount
        for iEight = 0 : maxStapleCount
            for iNine = 0 : maxStapleCount
                for iTen = 0 : maxStapleCount
                    if iEight + iNine + iTen ~= maxStapleCount
                        continue
                    end
                    
                    time = [repmat(deviceOnTime(1), [1 iSeven]) ...
                            repmat(deviceOnTime(2), [1 iEight]) ...
                            repmat(deviceOnTime(3), [1 iNine]) ...
                            repmat(deviceOnTime(4), [1 iTen])];

                    % generate a nDevice bit vector set device # 1
                    n = zeros(1, nDevice); n(iDevice) = 1;
                    y = [simulate(time, exprmntDrtn)' n];
                    % store all the training data in X
                    X = [X; y];
                    iDevice = iDevice + 1;
                end
            end
        end
    end
end

fprintf('Finished simulation experiment, saving data...\n');
csvwrite('traindata500.csv', X)

function onTimeEst = simulate(avgOnTime, TOTAL_TIME)

    % Rate constants for imager from Jungmann et. al (2010)
    % Rate constants for hairpins from Tsukanov et. al (2013)
    % 8nt - 20 sec^-1, 9nt - 2 sec^-1, 10nt - 0.11 sec^-1 off-rate
    imagerC = 30;                                 % [nM]
    kon = 1e6;                                    % molar^-1 sec^-1
    on_rate = kon * imagerC * 1e-9;               % sec^-1
    nDevice = length(avgOnTime);
    
    % interpolate time signal to capture rate of 1 ms
    outputTime = (1:TOTAL_TIME*1000); 
    % There are upto K devices attached to the nanostructure
    outputSignal = 0;
    parfor iDevice = 1 : nDevice
        off_rate = 1/avgOnTime(iDevice);                     
        if ~isnan(off_rate)
            [time, signal] = simulate_reaction(off_rate, on_rate,  TOTAL_TIME);
            outputSignal = outputSignal + interpolate(time, signal);
        end
    end

%     UNCOMMENT THIS CODE TO SEE A SAMPLE PLOT OF TEMPORAL BARCODE    
%     plot(oTime, output);
%     figure1 = figure('Color', [1 1 1], 'PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
%     axes1 = axes('Parent', figure1, 'FontSize', 16);
%     xlabel(axes1, 'Time (seconds)');
%     ylabel(axes1, 'Normalized Intensity (a.u.)');
%     plot(oTime, output + randn(size(output))/10, 'Parent', axes1, 'LineWidth', 2, 'Color', [0 0 0]);
%     ADD A BREAKPOINT HERE TO SEE THE PLOT

    % Generate statistics for different length of time-signal
    onTime = genOnOffStats(outputTime, outputSignal, nDevice);
    onTimeEst = zeros(nDevice, 1);
    for iDevice = 1 : nDevice
        if ~isempty(onTime{iDevice})
            estDist = fitdist(onTime{iDevice}, 'exponential');
        else
            estDist = makedist('exponential', 'mu', 1e-3);
        end
        onTimeEst(iDevice) = estDist.mean;
    end
end

function [t_ssa, signal] = simulate_reaction(off_rate, on_rate, TOTAL_TIME)
    %% Reaction network:
    %  1. imager + device {off}<-->{on} imager-device

    %% Create Decaying-Dimerizing Model
    model = sbiomodel('Tuning the nanostructure');

    %% Enter Reactions
    r1 = addreaction(model, 'I + D <-> I-D');

    %% Set Reactions to be MassAction
    k1 = addkineticlaw(r1, 'MassAction');

    %% Add Rate Constants for Each Reaction
    addparameter(k1, 'k_on', 'Value', on_rate);
    addparameter(k1, 'k_off', 'Value', off_rate);


    %% Set the Kinetic Law Constants for Each Kinetic Law.
    k1.ParameterVariableNames = {'k_on', 'k_off'};

    %% Specify Initial Amounts of Each Species
    model.species(1).InitialAmount = 1;          % I
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
    cs.SolverType = 'ssa';
    cs.StopTime = TOTAL_TIME;
    solver = cs.SolverOptions;
    solver.LogDecimation = 1;
    cs.CompileOptions.DimensionalAnalysis = false;
    [t_ssa, x_ssa] = sbiosimulate(model);
    
    signal = x_ssa(:, 3);
end

function onTimes = genOnOffStats(time, signal, K)
    % find the on-time for all the states
    onTimes = cell(1, K);
    for ik = 1 : K
        on_times = [];
        last_off_step = 1;
        iSignal = signal;
        
        % set all the lower values to zero and higher values to 1
        iSignal(iSignal < ik) = 0;
        iSignal(iSignal >= ik) = 1;

        % only the state change, time is recorded except last state
        for step = 2 : length(iSignal)
            if step == length(iSignal)
                if iSignal(step) > 0
                    on_times = [on_times; time(step)-time(last_off_step)];
                end
            end

            if iSignal(step) > 0 && iSignal(step-1) == 0
                last_off_step = step;
           elseif iSignal(step) == 0 && iSignal(step-1) > 0
                on_times = [on_times; time(step)-time(last_off_step)];
            end
        end

        onTimes{ik} = on_times.*0.001;
    end
end

function output = interpolate(time, signal)
    % This function interpolates signal to 1 ms capture rate. The output is
    % set to end at the last time stamp.
    %
    % time is the time stamps for state chain for reaction
    %
    % signal is 0 and 1 in alternate fashion

    capRate = 0.001;
    output = [];
    for t = 2:length(time)
        cur = interp1([time(t-1) time(t)], [signal(t-1) signal(t-1)], ...
                                                    time(t-1):capRate:time(t));
        output = cat(2, output, cur);
    end

    output = output(1:time(end)*(1/capRate));
end
