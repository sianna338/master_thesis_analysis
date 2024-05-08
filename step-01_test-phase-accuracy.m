%%%% this tests the accurary of phase dependent stimulation for single and paired pulses
%% preliminaires
clc 
clear all; 
close all; 

% fieldtrip
path_ft   = 'C:\Users\siann\Downloads\fieldtrip-20231113\fieldtrip-20231113';
addpath(path_ft);
ft_defaults;

sf = 5000; 

datapath = 'E:\spindle_ppTMS\EEG'

subjects = {'sub-09'}
session = {'ses-adapt'}
% define conditions
condition_peak = 'S155'
condition_trough = 'S156'
condition_rising = 'S157'
condition_falling = 'S158'
condition_sp_free = 'S159'

% individual spindle frequency
sp_freq = 11.84
%% main analysis
% look at markers in dataset and segment based on markers
for isub=1:length(subjects)
    for ises=1:length(session)
        cfg = [];
        cfg.datafile = [datapath, filesep, subjects{isub}, filesep, session{ises}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{ises}, '.eeg']
        cfg.headerfile = [datapath, filesep, subjects{isub}, filesep, session{ises}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{ises}, '.vhdr']
        cfg.continous = 'yes';
        cfg.trialdef.prestim = 1
        cfg.trialdef.poststim = -0.004 % end segment before TMS pulse to avoid confounding by artifact
        cfg.trialdef.eventtype = 'Stimulus';
        cfg.trialdef.eventvalue = {condition_peak, condition_trough, condition_falling,...
            condition_rising, condition_sp_free}
        cfg = ft_definetrial(cfg);
        trial_matrix = cfg.trl
   

        % read-in data 
        % cfg = []; 
        cfg.channel = {'C4', 'TP9'};
        cfg.reref = 'yes'
        cfg.refchannel = 'TP9'
        % cfg.bpfilter = 'yes'
        % cfg.bpfiltord = 2
        % cfg.bpfreq = [10.9 14.9]
        % cfg.bpfiltdir = 'twopass'
        data_raw_c4= ft_preprocessing(cfg);
        save([datapath, filesep, subjects{isub}, filesep, session{ises}, filesep, 'data_', subjects{isub}, '_', session{ises}, '_C4'], "data_raw_c4", '-v7.3')
    end
end 

% filter the data in the spindle frequency range
for itrial=1:length(data_raw_c4.trial)
    data_raw_c4.trial_filtered{itrial} = bandpass(data_raw_c4.trial{itrial}(1,:),[sp_freq-2 sp_freq+2],sf)
end 

% plot all the trials
for itrial=1:length(data_raw_c4.trial)
    figure;
    plot(data_raw_c4.trial_filtered{itrial})
end 

%% extract the phase for each trial
xh = zeros(length(data_raw_c4.trial),length(data_raw_c4.trial{1}(1,:)))
xphase = zeros(length(data_raw_c4.trial),length(data_raw_c4.trial{1}(1,:)))
xphase_unwrap = zeros(length(data_raw_c4.trial),length(data_raw_c4.trial{1}(1,:)))
for itrial = 1:length(data_raw_c4.trial)
    xh(itrial, :) = hilbert(data_raw_c4.trial_filtered{itrial}(1,:));
    xphase_unwrap(itrial, :) = (unwrap(angle(xh(itrial,:))));
    xphase(itrial, :) = angle(xh(itrial,:))
end 

% get phase for each condition 
conditions = {condition_peak, condition_trough, condition_rising,...
            condition_falling, condition_sp_free}
trials = cell(1, numel(conditions));
trials_unwrap = cell(1, numel(conditions));
for i = 1:numel(conditions)
    trials{i} = xphase(data_raw_c4.trialinfo == str2double(conditions{i}(2:end)),:)
    trials_unwrap{i} = xphase_unwrap(data_raw_c4.trialinfo == str2double(conditions{i}(2:end)),:)
end

trials_peak = trials{1} 
trials_trough = trials{2}
trials_rising = trials{3} 
trials_falling = trials{4} 
trials_post = trials{5}
trials_peak_unwrap = trials_unwrap{1} 
trials_trough_unwrap = trials_unwrap{2}
trials_falling_unwrap = trials_unwrap{3} 
trials_rising_unwrap = trials_unwrap{4}
trials_post_unwrap = trials_unwrap{5}
%% plot 
% for not unwrapped phase 
plot_positions = [0.5, 1, 1.5, 2, 2.5]
colors = {'.g', '.b', '.r', '.k', '.y'}
figure;
for i=1:numel(conditions)
    for itrial = 1:size(trials{i},1)
        h(i) = plot(plot_positions(i), trials{i}(itrial, end), colors{i}, 'MarkerSize', 10) % phase 0.004s before TMS pulse to avoid artifact
        hold on; 
    end 
end 

xlim([0.4 2.6])
ylim([-pi pi]); 
xlabel('targeted phase', 'FontSize', 12); 
ylabel('estimated phase (radians)', 'FontSize', 12); 
title('Estimated Phase vs. Targeted Phase', 'FontSize', 14); 
h = legend(h, {'Peak', 'Trough', 'Rising', 'Falling', 'Post'}, ...
           'FontSize', 10, 'Location', 'best'); 
grid on; 


% for unwrapped phase 
% plot_positions = [0.5, 1, 1.5, 2]
% colors = {'.g', '.b', '.r', '.k'}
trials_list = {trials_peak_unwrap, trials_falling_unwrap, trials_trough_unwrap, trials_rising_unwrap}
figure;
for i=1:numel(trials_list)
    for itrial = 1:size(trials_list{i},1)
        h(i) = plot(plot_positions(i), trials_list{i}(itrial, end), colors{i}, 'MarkerSize', 10)
        hold on; 
    end 
end  
hold on; 
t = linspace(0, 2*pi, 1000);
x = (mean(trials_list{1}(:, end))-mean(trials_list{3}(:, end)))/2 * cos(t) + (mean(trials_list{1}(:, end))+mean(trials_list{3}(:, end)))/2; % create cosine wave for reference 
plot(t/pi + 0.5, x, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]); 
xlim([0.4 2.6])
xlabel('targeted phase', 'FontSize', 12); 
ylabel('estimated phase (radians)', 'FontSize', 12); 
title('Comparison of Estimated Phase (unwrapped) with Cosine Wave', 'FontSize', 14); 
h = legend(h, {'Peak', 'Falling', 'Trough', 'Rising', 'Post'}, ...
           'FontSize', 10, 'Location', 'best')
grid on; 

%% circular plots
figure;
titles = {'peak stimulation', 'trough stimulation', 'rising flank stimulation', ...
    'falling flank stimulation'}
for i=1:numel(conditions)-1
    subplot(2,2,i)
    std_peak = circ_std(trials{i}(:, end), [], [], 1)
    mean_peak = circ_mean(trials{i}(:, end))
    circ_plot(trials{i}(:, end), 'pretty', 'bo', true,'linewidth',2,'color','b')
    hold on; 
    text(1.3, 0.8, sprintf('SD: %.2f\nMean: %.2f', std_peak, mean_peak), 'FontSize', 10);
    hold off; 
    title(titles{i}, 'Position',[0.2 1.27])
    % text(3, 0.2, ['RMS = 12.6'])
end 