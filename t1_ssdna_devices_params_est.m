%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shalin Shah
% Affiliation: Dept. of Electrical & Computer Engineering, Duke University
% Email: shalin.shah@duke.edu
% Last modified: 01/08/2019
% Matlab version used: R2017a
%
% Description: This code tunes length of DNA devices, generates barcodes,
% analyze signals and generates scatter plot with estimated parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; 
addpath(genpath(strcat(userpath,'\distinguishable_colors')));

dev_on_times = [0.005, 0.05, 0.55, 9.0];
NO_OF_SAMPLES = 100;
NO_OF_DEVICES = 4;
CAPTURE_TIME = 30*60;

C = distinguishable_colors(NO_OF_DEVICES, 'w');
single_blink_data = zeros(NO_OF_SAMPLES, NO_OF_DEVICES);
group = zeros(NO_OF_SAMPLES, NO_OF_DEVICES);
for i=1:NO_OF_DEVICES
    group(:, i) = i;
end


for j = 1:1
    full_set = []; half_set = []; less_set = []; 
    device_count = 1;
    for i = 4:NO_OF_DEVICES
        [full_est, half_est, less_est] = simulate(dev_on_times(i), CAPTURE_TIME);
        full_set = [full_set, full_est];
        half_set = [half_set, half_est];
        less_set = [less_set, less_est];
        
        single_blink_data(j, device_count) = (full_est.mean);
        figure(2); hold on;
        p = plot(single_blink_data(group==device_count), ...
            '.','markersize', 20,...
            'color', C(device_count, :)); 
        device_count = device_count + 1;
    end
    %plotMultiDeviceParam([less_set; half_set; full_set]);
end
grid on; box on; set(gca,'fontsize',20); set(gca,'fontweight','bold');
set(gca,'linew',3.0); set(gca,'yscale','log'); %axis([0 2 -0.25 0.25]);


function [full_est, half_est, less_est] = simulate(avg_on_time, TOTAL_TIME)
    import Gillespie.*

    %% Reaction network:
    %  1. imager + docker ->{on_rate} imager-docker
    %  2. imager-docker ->{off_rate} imager + docker

    % Rate constants for 10 nt imager strand from Jungmann et. al (2010)
    % 8nt - 20 sec^-1, 9nt - 2 sec^-1, 10nt - 0.11 sec^-1 off-rate
    imagerC = 30;                                 % [nM]
    kon = 1e6;                                    % molar^-1 sec^-1
    p.off_rate = 1/avg_on_time;                   % sec^-1
    p.on_rate = kon * imagerC * 1e-9;             % sec^-1 

    % Initial state
    tspan = [0, TOTAL_TIME];           % Observe for 10 minutes
    x0    = [1, 1, 0];           % docker, imager, imager-docker

    %% Specify reaction network
    pfun = @propensities_2state;
    stoich_matrix = [ -1  -1   1   % R1. binding              
                      1    1  -1]; % R2. unbinding

    %% Run simulation
    [time, signal] = directMethod(stoich_matrix, pfun, tspan, x0, p);
    signal = signal(:,3);
    [nTime, nSignal] = interpolate(1000, signal, time, 0);

%     binSize = 2;
%     x = reshape(nSignal, [binSize length(nSignal)/binSize]);
%     binnedSignal = sum(x, 1);
%     binnedSignal = binnedSignal + randn(size(binnedSignal))/2;
    
%     figure;
%     subplot(3, 1, 1); 
%     plot(nSignal);
%     subplot(3, 1, 2); 
%     plot(binnedSignal)
%     subplot(3, 1, 3);
%     histogram(binnedSignal, 50); set(gca, 'YScale', 'log')
    
    %% Generate statistics for different length of time-signal
    [full_on_time, ~] = genOnOffStats(time, signal, time(end));
    if ~isempty(full_on_time)
     %histfit(full_on_time, 4, 'exponential');
     full_est = fitdist(full_on_time, 'exponential');
    else
     full_est = makedist('exponential', 'mu', 1e-4);
    end
    
    [half_on_time, ~] = genOnOffStats(time, signal, time(end)/2);
    if ~isempty(half_on_time)
    %histfit(half_on_time, 4, 'exponential');
    half_est = fitdist(half_on_time, 'exponential');
    else
     half_est = makedist('exponential', 'mu', 1e-4);
    end
 
    [lim_on_time, ~] = genOnOffStats(time, signal, time(end)/6);
    if ~isempty(lim_on_time)
        less_est = fitdist(lim_on_time, 'exponential');
    else
     less_est = makedist('exponential', 'mu', 1e-4);
    end
end

function a = propensities_2state(x, p)
    %% Return reaction propensities given current state x
    imager = x(2);
    imager_docker = x(3);

    % Order of the elements in a should match the order of reactions in stoch
    a = [p.on_rate * imager;                     % R1. binding
         p.off_rate * imager_docker];            % R2. unbinding
end

function [on_times, off_times] = genOnOffStats(time, signal, trun_time)
    %% analyse the on off-time
    off_times = [];
    on_times = [];
    %% only the state change, time is recorded except last state
    last_step = signal(1);
    
    for step = 2 : length(signal)-1
        if time(step) < trun_time
            cur_step_len = time(step) - time(step-1);
            if cur_step_len > 0.3
                if signal(step) == 1
                    off_times = [off_times; cur_step_len];
                else
                    on_times = [on_times; cur_step_len];
                end
            else 
                break
            end
        else
            break
        end
    end
end

function plotEstParms(full_est, half_est, less_est)
    est_conf_less = paramci(less_est);
    est_conf_half = paramci(half_est);
    est_conf_full = paramci(full_est);
    
    est_mean = [less_est.mean, half_est.mean, full_est.mean];
    est_confI = [est_conf_less,est_conf_half, est_conf_full];
    
    hb = bar([1 2 3], est_mean', 0.5);
    set(hb(1), 'FaceColor',[0.85 0 0],  'LineWidth', 1.0);

    set(gca, 'FontSize',12,'XTick',[1 2 3],'XTickLabel',{'10 minutes','30 minutes','1 hour' });
    ylabel('Estimated Parameter', 'FontSize', 16)
    xlabel('Signal Length', 'FontSize', 16)
    
    hold on;
    errorbar([1, 2, 3], est_mean , ...
        [est_mean(1)-est_conf_less(1), est_mean(2)-est_conf_half(1), est_mean(3)-est_conf_full(1)],...
        [-est_mean(1)+est_conf_less(2), -est_mean(2)+est_conf_half(2), -est_mean(3)+est_conf_full(2)],...
        'ok', 'LineWidth', 2.0, 'CapSize', 8);
end

function plotMultiDeviceParam(esti_params)
    % This method generates a bar graph for estimated parameters with error
    % bars showing 95% confidence intervals
    figure(1); clf; hold on;
    set(gca, 'FontSize',12,'XTick',[1 2 3],'XTickLabel',{'10 minutes','30 minutes','1 hour' });
    ylabel('Estimated Parameter', 'FontSize', 16)
    xlabel('Signal Length', 'FontSize', 16)
    est_mean = [];
    est_conf = [];
    for k1 = 1:3
        est_mean = ([est_mean; esti_params(1, k1).mean, esti_params(2, k1).mean, esti_params(3, k1).mean]);
        est_conf = ([est_conf; paramci(esti_params(1, k1)), paramci(esti_params(2, k1)), paramci(esti_params(3, k1))]);

    end     
    hb = bar([1 2 3], est_mean', 1.0);        
    set(hb(1), 'FaceColor',[0.956, 0.760, 0.760],  'LineWidth', 2.0);
    set(hb(2), 'FaceColor', [0.698, 0.133, 0.133],  'LineWidth', 2.0);
    set(hb(3), 'FaceColor', [0.501, 0, 0],  'LineWidth', 2.0);
    

   for k1 = 1:3
       errorbar([1, 2, 3]+.22*(k1-2.0), est_mean(k1,:) , ...
          [est_mean(k1,1)-est_conf(k1*2-1,1), est_mean(k1,2)-est_conf(k1*2-1,2), est_mean(k1,3)-est_conf(k1*2-1,3)],...
          [-est_mean(k1,1)+est_conf(k1*2,1), -est_mean(k1,2)+est_conf(k1*2,2), -est_mean(k1,3)+est_conf(k1*2,3)],...
          'ok', 'LineWidth', 2.0, 'CapSize', 8);
   end
   legend('Device length: 08 nt', 'Device length: 09 nt', 'Device length: 10 nt');
   set(gca,'yscale','log'); grid on; box on;
end

function [nTime, nSignal] = interpolate(freq, signal, time, toDisp)
    % This function takes the original ODE signal containing 0 and 1 and  
    % interpolates it with given frequency
    %
    % freq is in hertz (eg. 10 Hz)
    %
    % signal is the original CTMC signal with only 0 and 1
    %
    % time is the time vector of the original signal
    % 
    
    % new signal will be more length determined by frequency
    nSignal = zeros(1, time(end)*freq);
    nTime = 1/freq:1/freq:time(end);
    
    % add bleaching effect
    avgBleachTime = 0.5;
    stdBleachTime = 0.1;
    
    for step = 2:length(time)
        if signal(step) == 0 && signal(step-1) == 1
            
            % if bleaching is faster than on-time total on-time reduces
            bleachTime = normrnd(avgBleachTime, stdBleachTime);
            if bleachTime < time(step) - time(step-1)
                time(step) = bleachTime;
            end
            nSignal(int32(time(step-1)*freq):int32(time(step)*freq)) = 1;
        end
    end
    % edge case of signal ending with 1 never going to 0
    if signal(end) == 1
         nSignal(int32(time(end-1)*freq):int32(time(end)*freq)) = 1;
    end
    
    if toDisp
        figure;
        subplot(3, 2, 1); 
        stairs(time, signal);
        subplot(3, 2, 2); 
        histogram(signal);

        subplot(3, 2, 3); 
        plot(nTime, nSignal);
        subplot(3, 2, 4); 
        histogram(nSignal);

        subplot(3, 2, 5); 
        plot(nTime, nSignal+ randn(size(nSignal))/10);
        subplot(3, 2, 6); 
        histogram(nSignal+ randn(size(nSignal))/10)
    end
end